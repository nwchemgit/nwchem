module wfkvars
! Stuff from Abinit WFK file
  character*6 :: codvsn
  integer :: headform,fform
  character*132 :: title
  double precision :: znuclpsp,zionpsp
  integer :: npsp,pspso,pspdat,pspcod,pspxc,lmax,lloc
  integer,allocatable :: kg(:,:)
  double precision,allocatable :: eigen(:)
  double complex, allocatable :: cg(:)
  double precision, allocatable :: kpt(:,:),occ(:)
  integer, allocatable :: indxkbnd(:),indxkpw(:),indxkcg(:)
  integer :: bantot,date,intxc,ixc,natom,ngfft(3),nkpt,&
& nspden,nspinor,nsppol,nsym,ntypat,occopt,pertcase,usepaw
  double precision :: acell(3),ecut,ecutdg,ecutsm,ecut_eff,qptn(3), &
&                     stmbias,tphysel,tsmear
  integer,allocatable :: istwfk(:),nband(:),npwarr(:),so_typat(:),&
& symafm(:),symrel(:,:,:),typat(:)
  double precision, allocatable :: tnons(:,:),znucltypat(:)
  double precision :: residm,etotal,fermie
  double precision,allocatable :: xred(:,:)
  integer :: npw,nbnd,npwt,ncg
  integer :: ngkpt(3)
  integer :: kptopt
  double precision :: shiftk(3)
end module wfkvars

!**************************************************************************

module kssvars
! Stuff from Abinit KSS file
  integer :: nsym2,nbandkss,npwkss,ishm,mpsang
  integer,allocatable :: symrel2(:,:,:),gbig(:,:)
  real,allocatable :: tnons2(:,:)
  integer,allocatable :: shlim(:)
  double precision,allocatable :: vkbsign(:,:),vkb(:,:,:,:),vkbd(:,:,:,:)
  double precision,allocatable :: en(:,:)
  double complex, allocatable :: wfg(:,:,:)
end module kssvars

!**************************************************************************

module pawvars
! Most of this is read from PSW files used by Abinit.
  integer :: ntypepaw
  integer :: nmeshmax=1,nbasismax=1,nlmnmax=1,ntlmnmax=1,llmaxmax=1,llocmax=1
  integer :: gridsizemax=1,ntlmn
  character(80), allocatable :: fnamepaw(:)
  integer, allocatable :: nmesh(:),nbasis(:),nlmn(:),llmax(:),lloc(:)
  integer, allocatable :: iphigrid(:),iprojgrid(:),icoregrid(:),ivgrid(:)
  integer, allocatable :: lorbital(:,:),llorb(:,:),mmorb(:,:),iorbno(:,:)
  integer, allocatable :: gridtype(:,:),gridsize(:,:)
  double precision, allocatable :: rstep(:,:),lstep(:,:)
  double precision, allocatable :: rphigrid(:,:),rprojgrid(:,:), &
&                                  rcoregrid(:,:),rvgrid(:,:)
  double precision, allocatable :: phi(:,:,:),tphi(:,:,:),tproj(:,:,:), &
&                                  coredens(:,:),tcoredens(:,:),vhntzc(:,:)
  double precision, allocatable :: dij0(:,:),rhoij0(:,:)
  double complex, allocatable :: projwf(:,:,:,:)
! projwf: inner product of projector with plane wave orbital, calculated
end module pawvars

!**************************************************************************

module pwx
  integer :: npwx
  integer, allocatable :: ipwx(:,:)
  double precision, allocatable :: engpw(:)
end module pwx


!**************************************************************************

module qpoints
  integer :: nqpt
  double precision, allocatable :: qpt(:,:)
end module qpoints

!**************************************************************************

module geometry
  double precision :: rprimd(3,3),rmet(3,3),bmet(3,3),blat(3,3),vol,vbz,bmetinv(3,3)
! rprimd: primitive direct lattice vectors
! blat: reciprocal lattice vectors
! vol: volume of primitive unit cell
! vbz: volume of 1st Brillouin zone
! bmet: metric tensor used for expanding dot products of quantities expressed in basis of reciprocal lattice vectors
! rmet: metric tensor used for expanding dot products of quantities expressed in basis of direct lattice vectors
end module geometry

!**************************************************************************

program mkprm
use wfkvars
use kssvars
use pwx
use qpoints
use pawvars
use geometry
use par
implicit none
character(80) :: fnamewf1,fnameinp1,fnameq,fname5,fnamevxc,fnamegwvxc,fnamecore,bindir,blank
character(80) :: fnamekss1,fnamekss2
character(2) :: chq1,chq2,chq3,chpw1,chpw2,chat,chbnd1,chbnd2
character(3) :: chpwx1,chpwy1,chpwz1,chpwx2,chpwy2,chpwz2,chpwl,chpwxt
character(11) :: strpw1,strpw2
character(200) :: fileout1,fileout2,fileout3,dirout,fileout4,fileout5,fileoutr
 integer :: mode
 integer, allocatable :: iptind(:,:,:),iptindfc(:,:,:,:),iptindbc(:,:,:), &
& ipttr(:,:,:),iptsym(:,:,:),idg(:,:,:,:),syminv(:,:,:)
 integer :: ncband,nbcore,nbocc,nbsev,nbsec,nwpt,ngcut,nvalel,ios,ios2,iinv,nsswpt
 integer :: nxtwpt,nxtqpt
 double precision, allocatable :: xtqpt(:,:)
 integer, allocatable :: xtqpt2qpt(:),qpt2xtqpt(:),qpt2xtg(:,:)
! xtqpt2qpt - input index of extended q-point grid, get index of standard q-point grid
! qpt2xtqpt - input index of standard q-point grid, get index of extended q-point grid at lower corner of cube containing standard q-point
! qpt2xtg - input index of standard q-point grid and direction, get shift in plane wave to go from qpt2xtqpt value to true corner
 integer :: test_bands_se(2), test_bands_pol(4)
 double precision :: wmax,xtwmax,qq(3),qs(3),qp(3),qps(3),vq,qq2,qp2,xck(3),xckq(3),qqmag
 double precision :: kv(3),kr(3)
 integer :: ix,iy,iz,ib,ik,ig,isum,igmx(3),igmn(3),ikk(3),ikx(3),iks(3),ngrid(3)
 integer :: igkmx(3),igkmn(3)
 integer, allocatable :: nsymk(:),symk(:,:),nqsym(:),qsym(:,:)
 integer, allocatable :: ihlf(:)
 integer, allocatable :: iwfg(:,:,:,:,:)
 integer, allocatable :: igndx(:,:,:,:)
 integer, allocatable :: igkndx(:,:,:)
 integer, allocatable :: ikndx(:,:,:),iqndx(:,:,:),ixtqndx(:,:,:)
 integer, allocatable :: isymndx(:,:,:),iqsymndx(:,:,:)
 integer, allocatable :: isymg(:,:,:),iqsymg(:,:,:)
 integer, allocatable :: lvtrans(:,:,:,:)
 integer :: npwlf,npwc,npwct,npwlft
 double complex, allocatable :: wfb(:,:,:,:),wfx(:)
 double complex :: wf
 double precision :: pi,rr(3),rrmag
 double complex, allocatable :: pola(:,:,:),polb(:,:,:),lossfn(:,:,:), &
& dielf(:),dielf2(:,:),lossfn2(:,:),coredielf(:),lossfn_loc(:,:,:)
! pola - hermitian part of vq * polarization 
! polb - antihermitian part of vq * polarization 
 double complex, allocatable :: xtpola(:,:,:),xtpolb(:,:,:),xtlossfn(:,:,:)
! extended dielectric response
 double complex :: rpola(3,3),rpolb(3,3),cpola(3,3),cpolb(3,3),rpoldum
 double complex :: bmet_t(3,3)
 double precision, allocatable :: secspec(:)
 double precision :: ww,dw,dxtw,hart,bohr,atomictime,mp,alpha
 double precision, allocatable :: enrgy(:),egrad(:,:),vqgrid(:)
 double precision, allocatable :: rhor(:),vxc(:),exc(:,:)
 integer :: ii,jj,kk,ll,iw,ikpt,lkpt,iqpt,iqpt0,ibnd,iband,jband,iiq,iq,jq
! ikpt - index of k point grid
! iqpt - index of q on q-point grid
! iqpt0 - index of q=(/0,0,0/) on q-point grid
 integer :: ipw,iipw,ipw1,ipw2,isym,ie,ikpt2,ipwx1,ipwy1,ipwz1,ipwx2,ipwy2,ipwz2
 integer :: ixx(3),ixp2(3),ixp1(3),ixm1(3),ixm2(3),itqp2,itqp1,itqm1,itqm2
 double precision :: poldum
 integer :: iatom,ilmn,jlmn,jkpt,jcg
 integer :: iocc,iunocc,igg(3),iqq(3),iqp(3),jka(3),jkk(3),iqs(3)
 integer :: igglf(3),igg0(3),igb(3)
 integer :: iisym,isign
 double precision :: omegap2,edens,gsum,speed,rate,atomdens,mfp,xsect,xkf,Ef,qtf,Eq,qkf,broadening,Eg
 double precision, allocatable :: fsum(:,:)
 double precision :: epwlf,epwc,epwx,epwg,xk(3),ek
 integer :: nln,ioff,ioffmx,jloop,ilsym,icg
 integer :: nprefsym,iprefsym(20)
 integer, allocatable :: ipsymndx(:)
 double precision :: exch
 integer, allocatable :: isecalculated(:,:),ipwsecalculated(:,:)
 integer :: tcalc,tncband,tnpwx,tnqpt,t_test_bands_se(2)
 double precision, allocatable :: xse(:,:),se(:,:)
 double complex, allocatable :: cse(:,:),zz(:,:)
 double precision, allocatable :: pwxse(:,:),pwse(:,:)
 double complex, allocatable :: pwcse(:,:),pwzz(:,:)
 double precision, allocatable :: ploss(:,:)
 double precision :: pwse_x, pwxse_x
 double complex :: pwcse_x, pwzz_x, pwdcse_x, pwdcse
 double precision :: xtpwse_x, xtpwxse_x
 double complex :: xtpwcse_x, xtpwzz_x, xtpwdcse_x
 double complex :: pse
 double complex :: wfval
 double precision :: strate,stpow,mass,vmax,cmin,bg
 integer, allocatable :: ipwndx(:,:),invpwndx(:,:,:),pwsymndx(:,:),invpw2ndx(:,:)
! ipwndx - given index of dual plane wave iipw, ipwndx(1,iipw) gives 1st plane wave, ipwndx(2,iipw) gives second plane wave
! invpwndx - given plane wave coordinates in integer multiples of reciprocal lattice vectors, finds plane wave index
! invpw2ndx - given pair of plane waves, finds index of dual plane wave
 double complex, allocatable :: pwmatel(:,:,:,:,:,:,:),tpwmatel(:,:,:,:,:,:,:)
! pwmatel - density matrix elements between PAW atomic orbitals
! tpwmatel - density matrix elements between PAW auxilliary atomic orbitals
 double complex, allocatable :: pwjmatel(:,:,:,:),tpwjmatel(:,:,:,:)
! pwjmatel - current matrix elements between PAW atomic orbitals at zero momentum transfer
! tpwjmatel - current matrix elements between PAW auxilliary atomic orbitals at zero momentum transfer
 double complex, allocatable :: aoovlp(:,:,:),taoovlp(:,:,:)
 double complex :: pwovlp,pawovlp,ovlp
 integer :: npwndx,ntpwndx,napwndx,pwcmax(3),pwcmin(3)
 double precision :: rat(3),xdummy,vdummy(100),mtxdummy(3,3),mtxdummy2(3,3),mtxdummy3(3,3),vtrdummy(3),vtrdummy2(3)
 double complex :: cvdummy(100),cmtxdummy(3,3)
 integer :: iat,inum,imult, iat_start, iat_end, iqpt_start, iqpt_end, iproc, &
 unit_num
 integer :: nfft
 integer :: iscissors,igw,sebandlo,sebandhi,sekptlo,sekpthi,iprtsym, &
& iprtdielf(5),iprtlosfn(5),iprtxtdielf(5),iprtxtlosfn(5),irdloss(2),ivxc,igwvxc,icore,itetrahedron
 integer :: pwlo, pwhi
 integer :: pwse_pwlo, pwse_pwhi
 integer :: pwse_nmin, pwse_nmax
 double precision :: pwse_mass
 double precision :: pwse_kstp(3)
 integer :: ipaw,ipolcalc,ipolpw(2),ipolqpt(2),iprtpolpw
 double precision :: scissors
 character*3 :: chkpt
 integer :: info,illpk,iulpk,neig
 double precision :: vllpk,vulpk,abstol
 double precision :: ssnrg, ssdnrg
 integer :: sskpt, sspw
 double precision :: sum
 integer :: gwbandmax
 integer :: qtensor
 double precision :: ktemp(3),kmag
 integer :: nxtgkpt(3)
 integer :: calc_dielf, calc_se, calc_ploss, calc_xtdielf, calc_pwse, calc_pwploss, calc_scspec
 integer :: OMP_GET_THREAD_NUM
 common /feg/ qqmag,Eq,Ef,qtf,xkf,qkf,edens,broadening,Eg,omegap2,pi

 call par_begin

 nxtgkpt = (/2,2,2/)
 open(unit=15, file="material_info", status='unknown')
! define some constants
 pi=acos(-1.d0)
 hart=27.2113845d0            ! atomic unit of energy = 27.2 eV
 bohr=0.5291772108d-10        ! atomic unit of distance = 0.529E-10 meters
 atomictime=2.4188843265d-17  ! atomic unit of time = 2.419E-17 seconds
 mp=1836.15267247             ! proton mass in atomic units
 alpha=1.d0/137.035999679d0   ! fine structure constant
 igg0=(/0,0,0/)
 iqpt0=-1  ! if the q-point grid does not contain q=(/0,0,0/), iqpt0 remains negative throughout calculation.

! read the input file "input.dat"
 call rdinput(calc_dielf, calc_se, calc_ploss, calc_xtdielf, calc_pwse, calc_pwploss, calc_scspec, &
& nbcore,nbocc,nbsev,nbsec,nwpt,wmax,nxtwpt,xtwmax,ncband,epwx,epwc,epwlf,epwg, &
& test_bands_se,test_bands_pol, &
& ngrid,nvalel,nprefsym,iprefsym, &
& fnamewf1,fnameinp1,fnamekss1,fnamekss2,fnameq,fnamevxc,fnamegwvxc,fnamecore,bindir, &
& iscissors,scissors,igw,sebandlo,sebandhi,sekptlo,sekpthi,pwse_pwlo,pwse_pwhi, &
& pwse_mass, pwse_nmin, pwse_nmax, pwse_kstp, &
& iprtsym,iprtdielf,iprtlosfn,iprtxtdielf,iprtxtlosfn,ipaw,ipolcalc,ipolpw,ipolqpt, &
& iprtpolpw,irdloss,ivxc,igwvxc,icore,itetrahedron, &
& sskpt, sspw, ssnrg, ssdnrg)
 wmax=wmax/hart
 xtwmax=xtwmax/hart
 dw=wmax/dble(nwpt)
 dxtw=xtwmax/dble(nxtwpt)

 iks=(/1,0,0/)
! Read the Abinit wave function
 write(*,*) 'Reading wave functions...'
 call rdwfk(fnamewf1,fnameinp1,ipaw,ncband)
! read q-points that will be used.
 call rdqpts(fnameq)
 write(6,'(a)') 'wave functions read'
 if (iprtdielf(5).eq.0) iprtdielf(5)=nqpt
 if (iprtlosfn(5).eq.0) iprtlosfn(5)=nqpt
 if (sekpthi.eq.0) sekpthi=nqpt

! compute geometrical parameters
 call mkgeom(pi)

 write(6,'(a,f10.6,a)') "Bare plasma frequency ",sqrt(4.d0*pi*(nbocc-nbcore)/vol)*hart," eV"
 if(master) write(15,'(a,f10.6,a)') "Bare plasma frequency ",sqrt(4.d0*pi*(nbocc-nbcore)/vol)*hart," eV"

 nfft=ngfft(1)*ngfft(2)*ngfft(3)
! If we want to compute electron density or exchange-correlation from abinit
 allocate(rhor(nfft),vxc(nfft))
! call rdden(fname5,ipaw,ngfft,npsp,nspden,rhor)
 if (ivxc.eq.1.or.igwvxc.eq.1) then
   allocate(exc(ncband,nqpt))
 endif
 if (ivxc.eq.1) then
   call rdvxc(fnamevxc,ipaw,ngfft,npsp,nspden,vxc)
 endif
 if (igwvxc.eq.1) then
   call rdgwvxc(fnamegwvxc,ncband,nqpt,exc,gwbandmax)
 endif

! read PAW files
 if (ipaw.ne.0) then
   call readpaw()
   write(6,*) "PAW files read"
 else
   nlmnmax=1
 endif
 
! Some free-electron electronic properties
 edens=dble(nvalel)/vol
 atomdens=dble(natom)/vol
 omegap2=4.d0*pi*edens
 xkf=(edens*(3.d0*pi*pi))**(1.d0/3.d0)
 Ef = 0.5d0*xkf*xkf
 qtf = 6*pi*edens/Ef
!write(*,*) xkf,Ef*hart,omegap2*hart

!do ii=1,natom
!  write(6,'(2i4,3x,a)') ii,typat(ii),trim(fnamepaw(typat(ii)))
!enddo
!

! Symmetry indices
 write(*,*) 'Handling symmetry...'
 allocate(ipsymndx(nsym))
 call mkipsymndx(nprefsym,iprefsym,nsym,ipsymndx)

 allocate(nqsym(nqpt),qsym(nqpt,nsym*2))
 allocate(iqndx(ngkpt(1),ngkpt(2),ngkpt(3)))
 allocate(ixtqndx(nxtgkpt(1),nxtgkpt(2),nxtgkpt(3)))
 allocate(iqsymndx(ngkpt(1),ngkpt(2),ngkpt(3)))
 allocate(iqsymg(3,nqpt,nsym))
 call fndsymk(symrel,ipsymndx,qpt,nqpt,nsym,ngkpt,(/0.d0,0.d0,0.d0/),1,nqsym,qsym,iqndx,iqsymndx,iqsymg)
 do iqpt=1,nqpt
   if(qpt(1,iqpt).eq.0.d0.and.qpt(2,iqpt).eq.0.d0.and.qpt(3,iqpt).eq.0.d0) then
     iqpt0=iqpt
   endif
 enddo

 allocate(nsymk(nkpt),symk(nkpt,nsym*2))
 allocate(ikndx(ngkpt(1),ngkpt(2),ngkpt(3)))
 allocate(isymndx(ngkpt(1),ngkpt(2),ngkpt(3)))
 allocate(isymg(3,nkpt,nsym))
 call fndsymk(symrel,ipsymndx,kpt,nkpt,nsym,ngkpt,shiftk,kptopt,nsymk,symk,ikndx,isymndx,isymg)
!stop
 allocate(syminv(3,3,nsym))
 call invsym(nsym,symrel,syminv)
 allocate(lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3)))
 call flvtrans(symrel,shiftk,kpt,ngkpt,nsym,nkpt,nsymk,symk,ikndx,isymndx,lvtrans)

! isum=0
! do ik=1,nkpt
!   write(6,'(2i3,2x,48i3)') ik,nsymk(ik),symk(ik,1:nsymk(ik))
!   isum=isum+nsymk(ik)
! enddo
! write(6,*) isum,ngkpt(1)*ngkpt(2)*ngkpt(3)
!stop

 if (iprtsym.ne.0) call prtsym(nsym,symrel,syminv,tnons,ngkpt,isymndx,ikndx)

 write(6,'(a)') 'k point and symmetry indices made'

! indices of transitions between reciprocal cells
 call figmxn(kg,npwt,igmx,igmn)
! write(6,'(3i3)') igmx(1:3)
! write(6,'(3i3)') igmn(1:3)
! call fprojwf(igmn,igmx,ncband)
 allocate(igndx (nkpt ,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3)))
 call mkigndx(kg ,indxkpw ,npwarr ,igmx,igmn,nkpt ,npwt ,igndx )

! ikpt=ikndx(5,5,10)
! do ii=igmn(1),igmx(1)
! write(6,'(i3,a)') ii,'------------------------------------------------------'
! do jj=igmn(2),igmx(2)
!   write(6,'(11i5)') igndx(ikpt,ii,jj,igmn(3):igmx(3))
! enddo
! enddo
! stop

 allocate(ihlf(nkpt))
 call fhlf(kpt,nkpt,ihlf)

 write(6,'(a)') 'plane wave indices made'

! ipwx(:,1:npwx) plane waves used for exchange
! ipwx(:,1:npwc) plane waves used for correlation
! ipwx(:,1:npwlf) plane waves used in local field block
 call mkpwx(epwx,bmet)
 if (pwse_pwhi.eq.0) pwse_pwhi=npwx
 npwc=0
 npwlf=0
 do ipw=1,npwx
   if (engpw(ipw).le.epwc) npwc=ipw
   if (engpw(ipw).le.epwlf) npwlf=ipw
 enddo
! write(6,*) npwx
! do ipw=1,npwx
!   write(6,'(i4,3x,3i3,3x,f14.6)') ipw,ipwx(1:3,ipw),engpw(ipw)
! enddo
! write(6,*) npwc
! do ipw=1,npwc
!   write(6,'(i4,3x,3i3,3x,f14.6)') ipw,ipwx(1:3,ipw),engpw(ipw)
! enddo
! write(6,*) npwlf  !,nqpt,nwpt,nwpt*nqpt*npwc,nwpt*nqpt*npwc*npwlf
! do ipw=1,npwlf
!   write(6,'(i4,3x,3i3,3x,f14.6)') ipw,ipwx(1:3,ipw),engpw(ipw)
! enddo
!stop
 pwcmax(1)=maxval(ipwx(1,1:npwc))
 pwcmax(2)=maxval(ipwx(2,1:npwc))
 pwcmax(3)=maxval(ipwx(3,1:npwc))
 pwcmin(1)=minval(ipwx(1,1:npwc))
 pwcmin(2)=minval(ipwx(2,1:npwc))
 pwcmin(3)=minval(ipwx(3,1:npwc))
 allocate(invpwndx(pwcmin(1):pwcmax(1),pwcmin(2):pwcmax(2),pwcmin(3):pwcmax(3)),pwsymndx(npwc,2*nsym))
 call mkinvpwndx(ipwx(:,1:npwc),npwc,invpwndx,pwcmax,pwcmin)
 call mkpwsymndx(ipwx(:,1:npwc),npwc,syminv,nsym,invpwndx,pwcmin,pwcmax,pwsymndx)

 npwndx=(npwlf*(npwlf+1))/2+npwc-npwlf
 ntpwndx=npwlf*npwlf+npwc-npwlf
 napwndx=ntpwndx+npwx-npwc
 if (iprtdielf(3).eq.0) iprtdielf(3)=npwndx
 if (iprtlosfn(3).eq.0) iprtlosfn(3)=npwndx
 if (iprtlosfn(3).eq.-1) iprtlosfn(3)=ntpwndx
 if (iprtxtdielf(3).eq.0) iprtxtdielf(3)=npwx
 if (iprtxtlosfn(3).eq.0) iprtxtlosfn(3)=npwx
 if (iprtpolpw.gt.0) then
   write(6,'(3i6)') npwlf,npwc,npwx
   write(6,'(3i6)') npwndx,ntpwndx,napwndx
   call par_end
   stop
 endif
 allocate(ipwndx(2,napwndx))
 call mkpwndx(npwx,npwc,npwlf,ipwndx,npwndx,ntpwndx,napwndx)
! do iipw=1,ntpwndx
!!   write(6,'(i3,5x,3(3i3,5x))') ipw,ipwx(:,ipwndx(1,ipw)),ipwx(:,ipwndx(2,ipw)), &
!!&        ipwx(:,ipwndx(1,ipw))-ipwx(:,ipwndx(2,ipw))
!   write(6,'(3i5,2(3x,3i5))') iipw,ipwndx(1,iipw),ipwndx(2,iipw),ipwx(:,ipwndx(1,iipw)),ipwx(:,ipwndx(2,iipw))
! enddo
!! write(6,*) npwndx,npwc*npwlf,npwc,npwlf
!stop
 allocate(invpw2ndx(npwx,npwx))
 call mkinvpw2ndx(ipwndx,napwndx,npwx,invpw2ndx)

! Write out file of momenta used - q-points and reciprocal lattice plane waves
 if (master) then
   open(unit=44,file='momenta.dat',status='unknown')
   write (44,*) "Reciprocal Lattice Vectors"
   write (44,*) 
   do ii=1,3
     write(44,'(a,i1,a,3(f9.6,a))') "K",ii," = (",blat(ii,1),",",blat(ii,2),",",blat(ii,3),")"
   enddo
   write (44,*)
   write (44,*) "Reciprocal Lattice Coordinates"
   write (44,*)
   write (44,*) "q-points"
   do iqpt=1,nqpt
     write(44,'(a,i3,a,3(f9.6,a))') "qpt: ",iqpt,", q  = (",qpt(1,iqpt),",",qpt(2,iqpt),",",qpt(3,iqpt),")"
   enddo
   write (44,*)
   write (44,*) "plane waves"
   do ipw=1,npwx
     write(44,'(a,i3,a,3(i6,a)f9.6)') "pw: ",ipw, &
  & ", q  = (",ipwx(1,ipw),",",ipwx(2,ipw),",",ipwx(3,ipw),"), q^2 = ",engpw(ipw)
   enddo
   write (44,*)
   write (44,*) "Cartesian Coordinates"
   write (44,*)
   write (44,*) "q-points"
   do iqpt=1,nqpt
     ktemp = qpt(1,iqpt)*blat(1,:) + qpt(2,iqpt)*blat(2,:) + qpt(3,iqpt)*blat(3,:)
     kmag = ktemp(1)*ktemp(1)+ktemp(2)*ktemp(2)+ktemp(3)*ktemp(3)
     write(44,'(a,i3,a,3(f9.6,a),f9.6,a,f9.6)') "qpt: ",iqpt,", q  = (",ktemp(1),",",ktemp(2),",",ktemp(3),"), q^2 = ",kmag, &
  & " |q| = ",sqrt(kmag)
   enddo
   write (44,*)
   write (44,*) "plane waves"
   do ipw=1,npwx
     ktemp = ipwx(1,ipw)*blat(1,:) + ipwx(2,ipw)*blat(2,:) + ipwx(3,ipw)*blat(3,:)
     kmag = ktemp(1)*ktemp(1)+ktemp(2)*ktemp(2)+ktemp(3)*ktemp(3)
     write(44,'(a,i3,a,3(f9.6,a),f9.6,a,f9.6)') "pw: ",ipw, &
  & ", q  = (",ktemp(1),",",ktemp(2),",",ktemp(3),"), q^2 = ",kmag, &
  & " |q| = ",sqrt(kmag)
   enddo
   write (44,*)
   write (44,*) "dual plane wave basis"
   write (44,*) "dual coordinate     plane wave set"
   do iipw=1,npwndx
     write(44,'(i8,11x,i5,a,i5)') iipw,ipwndx(1,iipw),",",ipwndx(2,iipw)
   enddo
   close(44)
 end if

! Adjust DFT eigen-energies
 write(6,'(a)') 'applying band energy corrections'

 allocate(enrgy(bantot))
 if (igw.ne.0) then
   call fndgwc(bindir,nbcore,ncband,eigen,enrgy,bantot, &
&   kpt,nkpt,indxkbnd)
   if (iscissors.ne.0) then
     do iband=nbocc+1,ncband
       do ikpt=1,nkpt
         ib=indxkbnd(ikpt)+iband
         enrgy(ib)=enrgy(ib)+scissors
       enddo
     enddo
   endif
   if(master) open(unit=44,file='energies.dat',status='unknown')
   do iband=1,ncband
     do iqpt=1,nqpt
       iqq=nint(qpt(:,iqpt)*ngkpt)+ngkpt/2
       ikpt=ikndx(iqq(1),iqq(2),iqq(3))
       if(master) write(44,'(2i5,3x,3i3,3f12.4)') iband,ikpt,nint(kpt(:,ikpt)*ngkpt), &
& eigen(indxkbnd(ikpt)+iband)*hart,enrgy(indxkbnd(ikpt)+iband)*hart, &
& (enrgy(indxkbnd(ikpt)+iband)-eigen(indxkbnd(ikpt)+iband))*hart
     enddo
   enddo
   close(44)
   write(6,'(a)') 'plasmon pole GW corrections read and applied'
 elseif (iscissors.ne.0) then
   if(master) open(unit=44,file='energies.dat',status='unknown')
   do iband=1,ncband
     do ikpt=1,nkpt
       ib=indxkbnd(ikpt)+iband
       if (iband.le.nbocc) then
         enrgy(ib)=eigen(ib)
       else
         enrgy(ib)=eigen(ib)+scissors
       endif
     enddo
     do iqpt=1,nqpt
       iqq=nint(qpt(:,iqpt)*ngkpt)+ngkpt/2
       ikpt=ikndx(iqq(1),iqq(2),iqq(3))
       if(master) write(44,'(2i5,3x,3i3,3f12.4)') iband,ikpt,nint(kpt(:,ikpt)*ngkpt), &
& eigen(indxkbnd(ikpt)+iband)*hart,enrgy(indxkbnd(ikpt)+iband)*hart, &
& (enrgy(indxkbnd(ikpt)+iband)-eigen(indxkbnd(ikpt)+iband))*hart
     enddo
   enddo
   if(master) close(44)
   write(6,'(a)') 'scissors operator energy corrections applied'
 else
   do ib=1,bantot
     enrgy(ib)=eigen(ib)
   enddo
   write(6,'(a)') 'no band energy corrections to apply'
 endif
! write(6,'(a)') 'Energy analysis'
 vmax=-1.d99
 cmin=1.d99
 do ikpt=1,nkpt
   vmax=max(vmax,enrgy(indxkbnd(ikpt)+nbocc))
   cmin=min(cmin,enrgy(indxkbnd(ikpt)+nbocc+1))
 enddo
 Eg = cmin-vmax
 write(6,'(a,f10.6,a)') 'minimum bandgap ',(cmin-vmax)*hart,' eV'
 write(15,'(a,f10.6,a)') 'minimum bandgap ',(cmin-vmax)*hart,' eV'
! ikpt=ikndx(ngkpt(1)/2,ngkpt(2)/2,ngkpt(3)/2)
!write(6,*) ikpt
! write(6,'(a,f10.6)') 'gamma bandgap ',(enrgy(indxkbnd(ikpt)+nbocc+1)-enrgy(indxkbnd(ikpt)+nbocc))*hart
! do iband=1,ncband
!   write(6,*) iband,eigen(indxkbnd(ikpt)+iband)
! enddo
!stop

! gradient of band energies in momentum space
 allocate(egrad(3,bantot))
 call fegrad(enrgy,bantot,rprimd,ngkpt,kpt,nkpt,shiftk,ikndx,indxkbnd,ncband,egrad)
 if (master) then
   open(unit=44,file='dispersion.dat',status='unknown')
     write(44,*) "band   k-point  gradient"
     do iband=1,ncband
       do ikpt=1,nkpt
         write(44,'(i5,2x,i5,2x,3es12.5)') iband,ikpt,egrad(:,indxkbnd(ikpt)+iband)
       enddo
     enddo
   close(44)
 end if

! Compute arrays for using PAW
 allocate(projwf(natom,nlmnmax,ngkpt(1)*ngkpt(2)*ngkpt(3),nbcore+1:ncband))
 allocate(pwmatel(ntypepaw,nlmnmax,nlmnmax,npwx,ngkpt(1),ngkpt(2),ngkpt(3)), &
&        tpwmatel(ntypepaw,nlmnmax,nlmnmax,npwx,ngkpt(1),ngkpt(2),ngkpt(3)))
 allocate(pwjmatel(3,ntypepaw,nlmnmax,nlmnmax), &
&        tpwjmatel(3,ntypepaw,nlmnmax,nlmnmax))
 allocate(aoovlp(ntypepaw,nlmnmax,nlmnmax),taoovlp(ntypepaw,nlmnmax,nlmnmax))
 if (ipaw.ne.0) then
! if (.false.) then
   !===============================
   !i think we want to do this loop with mpi .....
   call MPE_DECOMP1D(natom,numprocs,my_rank,iat_start,iat_end)
   !do iat=1,natom
   do iat=iat_start,iat_end
     write(chat,'(i2)') iat
     if (chat(1:1).eq.' ') chat(1:1)='0'
     fileout4=trim(bindir)//"projwf"//chat//".bin"
     open(unit=30,file=fileout4,status="old",form='unformatted',iostat=ios)
     if (ios.ne.0) then
!     if (.true.) then
       write(6,'(a)') 'Computing atomic orbital expansion coefficients, atom '//chat
       call ftprojwf(iat,ncband,nbcore,ihlf,igndx,igmn,igmx,ikndx,isymndx,syminv,isymg)
       open(unit=30,file=fileout4,status="unknown",form='unformatted')
       write(30) projwf(iat,:,:,:)
     !else
     !  write(6,'(a)') 'Reading atomic orbital expansion coefficients, atom '//chat
     !  read(30) projwf(iat,:,:,:)
     endif
     close (30)
   enddo
   call par_barrier
   do iat=1,natom
     write(chat,'(i2)') iat
     if (chat(1:1).eq.' ') chat(1:1)='0'
     fileout4=trim(bindir)//"projwf"//chat//".bin"
     open(unit=30,file=fileout4,status="old",form='unformatted')
     write(6,'(a)') 'Reading atomic orbital expansion coefficients, atom '//chat
     read(30) projwf(iat,:,:,:)
     close (30)
   enddo

   !===============================
!   iband=1
!   ikpt=2
!   iat=2
!   rr=(/0.1,0.1,0.1/)
!   call pwrwf(rr,iband,ikpt,wfval,ihlf,igndx,igmx,igmn)
!   write(6,*) wfval
!   call prjrtwf(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   call prjrpwf(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   call prjrwf(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   call pwrwfq(rr,iband,ikpt,wfval)
!   write(6,*) wfval
!   call prjrwfq(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   write(6,*)
!   iat=1
!   rr=(/-0.1,-0.1,-0.1/)
!   call pwrwf(rr,iband,ikpt,wfval,ihlf,igndx,igmx,igmn)
!   write(6,*) wfval
!   call prjrtwf(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   call prjrpwf(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   call prjrwf(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   call pwrwfq(rr,iband,ikpt,wfval)
!   write(6,*) wfval
!   call prjrwfq(rr,iat,iband,ikpt,wfval)
!   write(6,*) wfval
!   do ii=1,nlmn
!     write(6,'(3i3,2f20.10)') iorbno(ii),llorb(ii),mmorb(ii),projwf(iat,ii,ikpt,iband)
!   enddo
!stop

   qq=(/0.d0,0.d0,0.d0/)
   call fpwmatel(qq,igg0,aoovlp,taoovlp)
   write(6,'(a)') 'Computing atomic orbital current matrix elements'
   call fpwjmatel(pwjmatel,tpwjmatel)
   call MPE_DECOMP1D(npwx,numprocs,my_rank,iat_start,iat_end)
   do ipw=iat_start,iat_end
     if (ipw.lt.100) then
       write(chpw1,'(i2)') ipw
       if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
       fileout4=trim(bindir)//'pwmatel'//chpw1//'.bin'
     else
       write(chpwl,'(i3)') ipw
       fileout4=trim(bindir)//'pwmatel'//chpwl//'.bin'
     endif
     open(unit=30,file=fileout4,status="old",form='unformatted',iostat=ios)
     if (ios.ne.0) then
       write(6,'(a,i4,a,i4)') 'Computing atomic orbital density matrix elements, plane wave: ',ipw,' of ',npwx
       !============================
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz,iy,ix,iqq,iks,qq,qs)
       !!$OMP DO SCHEDULE(GUIDED)
       !$OMP DO SCHEDULE(dynamic)
       do ix=1,ngkpt(1)
       do iy=1,ngkpt(2)
       do iz=1,ngkpt(3)
         iqq=(/ix,iy,iz/)
         iks=iqq-ngkpt/2
         call fqq(shiftk,shiftk,iks,ngkpt,qq,qs)
         call fpwmatel(qq,ipwx(:,ipw),pwmatel(:,:,:,ipw,ix,iy,iz),tpwmatel(:,:,:,ipw,ix,iy,iz))
!         write(6,'(3i3,2x,3i3,2x,3f6.2)') iqq,ipwx(:,ipw),qq
!         do ii=1,nlmn
!           write(6,'(14f6.2)') dble(pwmatel(ii,:,ipw,ix,iy,iz))
!         enddo
!         write(6,*)
!         do ii=1,nlmn
!           write(6,'(14f6.2)') dble(tpwmatel(ii,:,ipw,ix,iy,iz))
!         enddo
!         write(6,*)
       enddo
       enddo
       enddo
       !$OMP END DO 
       !$OMP END PARALLEL
       !============================
       open(unit=30,file=fileout4,status="unknown",form='unformatted')
       write(30) pwmatel(:,:,:,ipw,:,:,:)
       write(30) tpwmatel(:,:,:,ipw,:,:,:)
     endif
     close(30)
   enddo
   !wait for everyone to finish their matrix elements before reading...
   call par_barrier()
   do ipw=iat_start,iat_end
     if (ipw.lt.100) then
       write(chpw1,'(i2)') ipw
       if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
       fileout4=trim(bindir)//'pwmatel'//chpw1//'.bin'
     else
       write(chpwl,'(i3)') ipw
       fileout4=trim(bindir)//'pwmatel'//chpwl//'.bin'
     endif
     open(unit=30,file=fileout4,status="old",form='unformatted')
     !write(6,'(a,i4)') 'Reading atomic orbital matrix elements, plane wave: ',ipw
     read(30) pwmatel(:,:,:,ipw,:,:,:)
     read(30) tpwmatel(:,:,:,ipw,:,:,:)
     close (30)
   end do
!   ipw=invpwndx(0,0,1)
!   write(6,*) ipw
!   do ii=1,nlmn
!     write(6,'(14f6.2)') dble(pwmatel(ii,:,ipw,5,5,1))
!   enddo
!   write(6,*)
!   do ii=1,nlmn
!     write(6,'(14f6.2)') dble(tpwmatel(ii,:,ipw,5,5,1))
!   enddo
!   ikpt=100
!   jkpt=100
!   iband=16
!   jband=16
!   call fovlp(ikpt,iband,jkpt,jband,ipaw,nkpt,ncg,nlmn,natom,ncband,cg,indxkcg,npwarr,projwf,aoovlp,taoovlp,pwovlp,pawovlp,ovlp)
!stop
!   if (iprtsym.ne.0) call chkwfsym(nsym,syminv,symrel,natom,ihlf,nkpt,igndx,igmx,igmn)
!   ipw=1
!!   qq=(/-0.1d0,-0.1d0,0.2d0/)
!   qq=(/0.0d0,0.0d0,0.0d0/)
!   call fpwmatel(qq,ipwx(:,ipw),aoovlp,taoovlp)
 endif
!stop
call par_barrier()

 if (ivxc.eq.1) then
!   write(6,*)
!   write(6,*) "band  k point    VXC"
   ikpt=1
   do iband=1,ncband
     call mkxclda(iband,ikpt,ngfft,vxc,vol,ihlf,nkpt,igndx,igmx,igmn,natom, &
&       exc(iband,ikpt))
!     write(6,'(i4,2x,i4,2x,f12.6,2x,f12.6)') iband,ikpt,exc(iband,ikpt),exc(iband,ikpt)*hart
   enddo
!   write(6,*)
 endif

!!!!!! DEBUG
! igglf=(/0,0,0/)
! iks=(/4,0,0/)
! ikk=(/1,5,7/)
! iocc=4
! iunocc=5
! iqs=iks+ngkpt/2
! ipw=invpwndx(igglf(1),igglf(2),igglf(3))
! !do ix=1,ngkpt(1)
! !do iy=1,ngkpt(2)
! !do iz=1,ngkpt(3)
! do ix=1,1
! do iy=5,5
! do iz=7,7
! ikk=(/ix,iy,iz/)
! call testME(igglf,iks,ikk,iocc,iunocc, &
!& pi,nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt,ncg,bantot,npwt,nlmn, &
!& nwpt,nbcore,nbocc,ncband,ngkpt,natom,projwf,kg,enrgy,cg, &
!& indxkpw,indxkbnd,indxkcg,npwarr,kpt,nsymk,symk,symrel,syminv, &
!& ihlf,lvtrans,ntypepaw, &
!& pwjmatel,tpwjmatel, &
!& pwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)), & 
!& tpwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)), &
!& igndx,ikndx,isymndx,isymg)
! enddo
! enddo
! enddo
! stop
!!!!!! DEBUG

!!!!!! Calculate secondary particle spectrum
 if (calc_scspec.ne.0.and.sspw.gt.0) then
   write(6,*) "Calculating secondary particle spectrum at "
   write(6,'(a,f10.6,a)') "energy: ",(ssnrg*hart)," eV"
   nsswpt=wmax/ssdnrg
   allocate(secspec(nsswpt))
   do ii=1,3
     iks(ii)=nint(qpt(ii,sskpt)*ngkpt(ii))
     iqq(ii)=iks(ii)+ipwx(ii,sspw)*ngkpt(ii)
     iqp(ii)=iks(ii)+ipwx(ii,sspw)*ngkpt(ii)
   enddo
   iqs=iks+ngkpt/2
   call fqq(shiftk,shiftk,iqq,ngkpt,qq,qs)
   call fqq(shiftk,shiftk,iqp,ngkpt,qp,qps)
   if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
     write(6,'(2(a,3(f9.6,a)))') "q  = (",qq(1),",",qq(2),",",qq(3),"), ","q' = (",qp(1),",",qp(2),",",qp(3),")"
     call mksecspecj(qq,ipwx(:,sspw),ipwx(:,sspw),iks,vol,pi,ssnrg,ssdnrg, &
&     bmet,test_bands_pol, &
&     nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt, &
&     ncg,bantot,npwt,nlmnmax, &
&     nsswpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
&     xred, &
&     projwf, &
&     kg, &
&     enrgy, &
&     cg, &
&     indxkpw, &
&     indxkbnd, &
&     indxkcg, &
&     npwarr, &
&     kpt, &
&     nsymk,symk,symrel,syminv, &
&     ihlf,lvtrans, &
&     ntypepaw, &
&     pwjmatel,tpwjmatel, &
&     pwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)), &
&     pwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)), &
&     igndx,ikndx,isymndx,isymg, &
&     nband,secspec)
   else 
!     call fqq(shiftk,shiftk,iqq,ngkpt,qq,qs)
!     call fqq(shiftk,shiftk,iqp,ngkpt,qp,qps)
     write(6,'(2(a,3(f9.6,a)))') "q  = (",qq(1),",",qq(2),",",qq(3),"), ","q' = (",qp(1),",",qp(2),",",qp(3),")"
     call mksecspec(qq,ipwx(:,sspw),ipwx(:,sspw),iks,vol,pi,ssnrg,ssdnrg, &
&     test_bands_pol, &
&     nsppol,shiftk,npw,igmx,igmn,ipaw,nsym,nkpt, &
&     ncg,bantot,npwt,nlmnmax, &
&     nsswpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
&     xred, &
&     projwf, &
&     kg, &
&     enrgy, &
&     cg, &
&     indxkpw, &
&     indxkbnd, &
&     indxkcg, &
&     npwarr, &
&     kpt, &
&     nsymk,symk,symrel,syminv, &
&     ihlf,lvtrans, &
&     ntypepaw, &
&     pwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)), &
&     pwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,sspw,iqs(1),iqs(2),iqs(3)), &
&     igndx,ikndx,isymndx,isymg, &
&     nband,secspec)
   endif
   if(master) then
     open(unit=26,file='secspec.dat',status='unknown',iostat=ios)
     write(26,*) "#    energy                 probability         cumulative"
     write(26,*) "#electron  hole             distribution    distribution function"
     sum=0.d0;
     do iw=1,nsswpt
       ww=iw*ssdnrg
       sum=sum+secspec(iw)*ssdnrg
       write(26,'(2(f6.2,3x),2f20.6)') (ww-vmax)*hart,(ssnrg-ww+vmax)*hart,secspec(iw)/hart,sum
     enddo
     close(26)
   end if
   write(6,*) "Secondary particle spectrum completed"
 endif

!!!!!! Calculate polarizability
 if (calc_dielf.gt.0.or.calc_se.gt.0.or.calc_ploss.gt.0.or.calc_pwse.gt.0.or.calc_pwploss.gt.0) then
   if (ipolpw(2).le.0) ipolpw(2)=npwndx
   if (ipolqpt(2).le.0) ipolqpt(2)=nqpt
   inum=4
   imult=12
   allocate(vqgrid(nqpt),pola(nwpt,nqpt+9,npwndx),polb(nwpt,nqpt+9,npwndx))
!   q-points greater than nqpt describe the 
!   q->(0,0,0) part of the polarization tensor
   qtensor = 0
   do iqpt=ipolqpt(1),ipolqpt(2)
     do ii=1,3
       iks(ii)=nint(qpt(ii,iqpt)*ngkpt(ii))
     enddo
     if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
       qtensor = iqpt
     endif
   enddo
 !=======================================
 !mpi for this loop?
 !why yes, i think so.
 call MPE_DECOMP1D(ipolpw(2)-ipolpw(1)+1,numprocs,my_rank,iat_start,iat_end)
 iat_start=iat_start+ipolpw(1)-1
 iat_end=iat_end+ipolpw(1)-1
 !write(*,*) 'proc ',my_rank,' will do plane waves from ',iat_start&
 !,' to ',iat_end,'.'
 !write(*,*) 'Global plane wave bounds: ',ipolpw(1),ipolpw(2)
 !stop
 !do iipw=ipolpw(1),ipolpw(2)
 do iipw=iat_start,iat_end
   ipw1=ipwndx(1,iipw)
   ipw2=ipwndx(2,iipw)
   write(chpw1,'(i2)') ipw1
   if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
   write(chpw2,'(i2)') ipw2
   if (chpw2(1:1).eq.' ') chpw2(1:1)='0'
   ipwx1=ipwx(1,ipwndx(1,iipw))
   ipwy1=ipwx(2,ipwndx(1,iipw))
   ipwz1=ipwx(3,ipwndx(1,iipw))
   ipwx2=ipwx(1,ipwndx(2,iipw))
   ipwy2=ipwx(2,ipwndx(2,iipw))
   ipwz2=ipwx(3,ipwndx(2,iipw))
   !write(*,*) 'Top of loop 1'
   !write(*,*) ipwx1,ipwy1,ipwz1,ipwx2,ipwy2,ipwz2
   !call system('date')
   !write(*,*) 'my_rank: ',my_rank
   write(chpwx1,'(i3)') abs(ipwx1)
   if (chpwx1(2:2).eq.' ') chpwx1(2:2)='0'
   if (ipwx1.ge.0) then
     chpwx1(1:1)='+'
   else
     chpwx1(1:1)='-'
   endif
   write(chpwy1,'(i3)') abs(ipwy1)
   if (chpwy1(2:2).eq.' ') chpwy1(2:2)='0'
   if (ipwy1.ge.0) then
     chpwy1(1:1)='+'
   else
     chpwy1(1:1)='-'
   endif
   write(chpwz1,'(i3)') abs(ipwz1)
   if (chpwz1(2:2).eq.' ') chpwz1(2:2)='0'
   if (ipwz1.ge.0) then
     chpwz1(1:1)='+'
   else
     chpwz1(1:1)='-'
   endif
   write(chpwx2,'(i3)') abs(ipwx2)
   if (chpwx2(2:2).eq.' ') chpwx2(2:2)='0'
   if (ipwx2.ge.0) then
     chpwx2(1:1)='+'
   else
     chpwx2(1:1)='-'
   endif
   write(chpwy2,'(i3)') abs(ipwy2)
   if (chpwy2(2:2).eq.' ') chpwy2(2:2)='0'
   if (ipwy2.ge.0) then
     chpwy2(1:1)='+'
   else
     chpwy2(1:1)='-'
   endif
   write(chpwz2,'(i3)') abs(ipwz2)
   if (chpwz2(2:2).eq.' ') chpwz2(2:2)='0'
   if (ipwz2.ge.0) then
     chpwz2(1:1)='+'
   else
     chpwz2(1:1)='-'
   endif
   strpw1=chpwx1//','//chpwy1//','//chpwz1
   strpw2=chpwx2//','//chpwy2//','//chpwz2
   fileoutr='pol'//strpw1//':'//strpw2//'.bin'
   fileout3=trim(bindir)//trim(fileoutr)
   open(unit=25,file=trim(fileout3),status='old',form='unformatted',iostat=ios)
   if (ios.ne.0.or.ipolcalc.ne.0) then
     igglf=ipwx(:,ipw2)-ipwx(:,ipw1)
!     write(6,'(a,3i3,4x,a,3i3)') "calculating polarization, G = ",ipwx(:,ipw1)," G' = ",ipwx(:,ipw2)
     write(6,'(a,2x,2i3,2x,i4,a,i4)') "calculating polarization "//trim(fileoutr),ipw1,ipw2,iipw," of ",ipolpw(2)
     !===============================
     !simplest possible OMP parallelization of this loop.
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iqpt,ii,jj,iks,iqq,iqp,qq,qs,qps,qq2,qp2,vq,iqs,iw,qp,iiq)
     !!$OMP DO SCHEDULE(GUIDED)
     !$OMP DO SCHEDULE(dynamic)
     do iqpt=ipolqpt(1),ipolqpt(2)
       !write(*,*) 'iqpt, ithread: ',iqpt,OMP_GET_THREAD_NUM()
       do ii=1,3
         iks(ii)=nint(qpt(ii,iqpt)*ngkpt(ii))
         iqq(ii)=iks(ii)+ipwx(ii,ipw1)*ngkpt(ii)
         iqp(ii)=iks(ii)+ipwx(ii,ipw2)*ngkpt(ii)
       enddo
       call fqq(shiftk,shiftk,iqq,ngkpt,qq,qs)
       call fqq(shiftk,shiftk,iqp,ngkpt,qp,qps)
       write(6,'(i3,a,i3,2(a,3(f6.3,a)))') iqpt," of ",ipolqpt(2), &
&        ",  q  = (",qq(1),",",qq(2),",",qq(3),"), ", &
&        "q' = (",qp(1),",",qp(2),",",qp(3),")"

       qq2=0.d0
       qp2=0.d0
       do ii=1,3
       do jj=1,3
         qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
         qp2=qp2+qp(ii)*bmet(ii,jj)*qp(jj)
       enddo
       enddo
       vq=4.d0*pi/sqrt(qq2*qp2)
       vqgrid(iqpt)=vq

       iqs=iks+ngkpt/2

!write(28,*) iipw
       if (iks(1).eq.0.and.iks(2).eq.0.and.iks(3).eq.0) then
!       if (.false.) then
         call mkpolj1(qq,ipwx(:,ipw1),ipwx(:,ipw2),iks,vq,pi, &
&         test_bands_pol, &
&         nsppol,shiftk,npw,igmx,igmn,ipaw,itetrahedron, &
&         nsym,nkpt, &
&         ncg,bantot,npwt,nlmnmax, &
&         nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
&         xred, &
&         projwf, &
&         kg, &
&         enrgy, eigen, &
&         cg, &
&         indxkpw, &
&         indxkbnd, &
&         indxkcg, &
&         npwarr, &
&         kpt, &
&         nsymk,symk,symrel,syminv, &
&         ihlf,lvtrans, &
&         ntypepaw, &
&         pwjmatel,tpwjmatel, &
&         pwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)), &
&         pwmatel(:,:,:,ipw2,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,ipw2,iqs(1),iqs(2),iqs(3)), &
&         igndx,ikndx,isymndx,isymg, &
&         nband,pola(:,nqpt+1:nqpt+9,iipw),polb(:,nqpt+1:nqpt+9,iipw))
         do iw=1,nwpt
           pola(iw,iqpt,iipw) = (0.d0,0.d0)
           polb(iw,iqpt,iipw) = (0.d0,0.d0)
           do ii=1,3
             iiq = ii+3*(ii-1)
             pola(iw,iqpt,iipw) = pola(iw,iqpt,iipw) + pola(iw,nqpt+iiq,iipw)
             polb(iw,iqpt,iipw) = polb(iw,iqpt,iipw) + polb(iw,nqpt+iiq,iipw)
           enddo
           pola(iw,iqpt,iipw) = pola(iw,iqpt,iipw)/3.d0
           polb(iw,iqpt,iipw) = polb(iw,iqpt,iipw)/3.d0
!write(28,'(i5,f9.2,9e12.4)') iw,(dimag(polb(iw,nqpt+iiq,iipw)),iiq=1,9)
         enddo
       else
         call mkpol1(qq,ipwx(:,ipw1),ipwx(:,ipw2),iks,vq,vol,pi, &
&         test_bands_pol, &
&         nsppol,shiftk,npw,igmx,igmn,ipaw,itetrahedron, &
&         nsym,nkpt, &
&         ncg,bantot,npwt,nlmnmax, &
&         nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom, &
&         xred, &
&         projwf, &
&         kg, &
&         enrgy, eigen, &
&         cg, &
&         indxkpw, &
&         indxkbnd, &
&         indxkcg, &
&         npwarr, &
&         kpt, &
&         nsymk,symk,symrel,syminv, &
&         ihlf,lvtrans, &
&         ntypepaw, &
&         pwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)), &
&         pwmatel(:,:,:,ipw2,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,ipw2,iqs(1),iqs(2),iqs(3)), &
&         igndx,ikndx,isymndx,isymg, &
&         nband,pola(:,iqpt,iipw),polb(:,iqpt,iipw))
       endif
!write(6,*) 1.d0-pola(1,iqpt,iipw)-polb(1,iqpt,iipw)
  
       write(chpw1,'(i2)') ipw1
       if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
       write(chpw2,'(i2)') ipw2
       if (chpw2(1:1).eq.' ') chpw2(1:1)='0'

!       open(unit=20,file="losfn.dat",status="unknown")
!       open(unit=21,file="dielf.dat",status="unknown")
!       open(unit=22,file="pola.dat",status="unknown")
!       open(unit=23,file="polb.dat",status="unknown")
!       do iw=1,nwpt
!         ww=dble(iw)*dw
!!         write(20,1000) ww*hart,1.d0/(1.d0-pola(iw,iqpt,iipw)-polb(iw,iqpt,iipw))
!         write(21,1000) ww*hart,1.d0-pola(iw,iqpt,iipw)-polb(iw,iqpt,iipw)
!         write(22,1000) ww*hart,pola(iw,iqpt,iipw)
!         write(23,1000) ww*hart,polb(iw,iqpt,iipw)
!!1000     format(f6.2,3x,2f20.6)
!1000     format(f6.2,3x,2es20.6)
!       enddo
!       close(20)
!       close(21)
!       close(22)
!       close(23)

     enddo
     !$OMP END DO 
     !$OMP END PARALLEL
     !===============================
     open(unit=25,file=trim(fileout3),form='unformatted',status='unknown',iostat=ios)
     if (ios.ne.0) then
       write(6,*) "Error - cannot open file ",trim(fileout3)
       stop
     endif
     write(25) pola(:,:,iipw),polb(:,:,iipw)
     close(25)
   endif
   close(25)
 enddo
 !mpi regropu here. we need full polarizablility matrix here.
 call par_barrier
 do iipw=ipolpw(1),ipolpw(2)
   ipw1=ipwndx(1,iipw)
   ipw2=ipwndx(2,iipw)
   write(chpw1,'(i2)') ipw1
   if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
   write(chpw2,'(i2)') ipw2
   if (chpw2(1:1).eq.' ') chpw2(1:1)='0'
   ipwx1=ipwx(1,ipwndx(1,iipw))
   ipwy1=ipwx(2,ipwndx(1,iipw))
   ipwz1=ipwx(3,ipwndx(1,iipw))
   ipwx2=ipwx(1,ipwndx(2,iipw))
   ipwy2=ipwx(2,ipwndx(2,iipw))
   ipwz2=ipwx(3,ipwndx(2,iipw))
   !write(*,*) 'Top of loop 1'
   !write(*,*) ipwx1,ipwy1,ipwz1,ipwx2,ipwy2,ipwz2
   !call system('date')
   !write(*,*) 'my_rank: ',my_rank
   write(chpwx1,'(i3)') abs(ipwx1)
   if (chpwx1(2:2).eq.' ') chpwx1(2:2)='0'
   if (ipwx1.ge.0) then
     chpwx1(1:1)='+'
   else
     chpwx1(1:1)='-'
   endif
   write(chpwy1,'(i3)') abs(ipwy1)
   if (chpwy1(2:2).eq.' ') chpwy1(2:2)='0'
   if (ipwy1.ge.0) then
     chpwy1(1:1)='+'
   else
     chpwy1(1:1)='-'
   endif
   write(chpwz1,'(i3)') abs(ipwz1)
   if (chpwz1(2:2).eq.' ') chpwz1(2:2)='0'
   if (ipwz1.ge.0) then
     chpwz1(1:1)='+'
   else
     chpwz1(1:1)='-'
   endif
   write(chpwx2,'(i3)') abs(ipwx2)
   if (chpwx2(2:2).eq.' ') chpwx2(2:2)='0'
   if (ipwx2.ge.0) then
     chpwx2(1:1)='+'
   else
     chpwx2(1:1)='-'
   endif
   write(chpwy2,'(i3)') abs(ipwy2)
   if (chpwy2(2:2).eq.' ') chpwy2(2:2)='0'
   if (ipwy2.ge.0) then
     chpwy2(1:1)='+'
   else
     chpwy2(1:1)='-'
   endif
   write(chpwz2,'(i3)') abs(ipwz2)
   if (chpwz2(2:2).eq.' ') chpwz2(2:2)='0'
   if (ipwz2.ge.0) then
     chpwz2(1:1)='+'
   else
     chpwz2(1:1)='-'
   endif
   strpw1=chpwx1//','//chpwy1//','//chpwz1
   strpw2=chpwx2//','//chpwy2//','//chpwz2
   fileoutr='pol'//strpw1//':'//strpw2//'.bin'
   fileout3=trim(bindir)//trim(fileoutr)
   open(unit=25,file=trim(fileout3),status='old',form='unformatted',iostat=ios)
   !write(6,'(a,2x,2i3,2x,i4)') "reading polarization "//trim(fileoutr),ipw1,ipw2,iipw
   read(25) pola(:,:,iipw),polb(:,:,iipw)
   close(25)

 enddo
 write(6,*) "Polarizability function found."
 if (iprtdielf(1).ne.0.and.master) then 
   write(6,*) "Printing dielectric function."
   call prtdielf(iprtdielf,npwc,npwlf,ipwndx,napwndx,nqpt,qpt,qtensor,ngkpt,nwpt,npwndx,pola,polb,dw,ipwx,bmet)
 end if
 if (ipolpw(1).ne.1.or.ipolpw(2).ne.npwndx) stop
  
   !if(master) write(*,*) 'Done interpolating polarizability to q=0'
 !==================================
 write(6,*) "Polarizability function found."
 if (iprtdielf(1).ne.0.and.master) write(6,*) "Printing dielectric function."
 if (iprtdielf(1).ne.0.and.master) &
 call prtdielf(iprtdielf,npwc,npwlf,ipwndx,napwndx,nqpt,qpt,qtensor,ngkpt,nwpt,npwndx,pola,polb,dw,ipwx,bmet)
 if (ipolpw(1).ne.1.or.ipolpw(2).ne.npwndx) stop

 if (master) then
   allocate(fsum(npwc,nqpt))
   open(unit=44,file='sumrule.dat',status='unknown')
   write(44,*) "Plane wave    q-point"
   do ipw=1,npwc
     iipw=invpw2ndx(ipw,ipw)
     do iqpt=1,nqpt
       fsum(ipw,iqpt)=0.d0
       gsum=0.d0
       ikk=nint((qpt(:,ikpt)+0.5d0)*dble(ngkpt)-shiftk)
       ikpt=ikndx(ikk(1),ikk(2),ikk(3))
       do iw=1,nwpt
         ww=dble(iw)*dw
         fsum(ipw,iqpt)=fsum(ipw,iqpt)+2.d0*ww*dw*dimag(1.d0/(1.d0-pola(iw,iqpt,iipw)-polb(iw,iqpt,iipw)))/(pi*omegap2)
         gsum=gsum-2.d0*ww*dw*dimag(1.d0-pola(iw,iqpt,iipw)-polb(iw,iqpt,iipw))/(pi*omegap2)
       enddo
       write(44,'(3i3,3x,3f6.2,3x,"fsum = ",f10.6,", gsum = ",f10.6)') ipwx(:,ipw),qpt(:,iqpt),fsum(ipw,iqpt),gsum
     enddo
   enddo
   close(44)
 end if

if (icore.ne.0) call CoreCorrectPol(wmax,nwpt,nqpt,npwndx,napwndx,ipwndx, &
&                    pola,polb,vol,pi,hart)

!!!!!! Calculate loss function
 !first compute which q-points belong to this mpi node; we divide the q-points
 !up between non-master nodes if any exist. 
 if(.not.master) then
   call MPE_DECOMP1D(nqpt,numprocs-1,my_rank-1,iqpt_start,iqpt_end)
   write(*,*) 'mpi proc ',my_rank,' will do q-points from ',iqpt_start&
   ,' to ',iqpt_end,' out of ',nqpt,'.'
 else
   iqpt_start=1
   iqpt_end=nqpt
 end if
 !if we need a local array to stor partial calculations of loss function,
 !allocate it
 if (.not.master) then
   allocate(lossfn_loc(nwpt,iqpt_start:iqpt_end,ntpwndx))
 else
   allocate(lossfn_loc(nwpt,nqpt+9,ntpwndx))
 end if
 !write(*,*) 'size of lossfn_loc: ', size(lossfn_loc)
 !!if the master is the only proc, his loss function array needs to be able to
 !!hold all the q-points and the tensor part
 !if (numprocs==1) then
 !  allocate(lossfn_loc(nwpt,nqpt+9,ntpwndx))
 !end if
 allocate(dielf(npwndx),dielf2(npwc,npwc),lossfn2(npwc,npwc))

 !allocate full loss function which will not be fully filled until the lossfn
 !file is read for worker procs
! q-points greater than nqpt describe the 
! q->(0,0,0) part of the loss function tensor
 allocate(lossfn(nwpt,nqpt+9,ntpwndx))
!write(*,*) 'allocating lossfn', nwpt,nqpt+9,ntpwndx,size(lossfn)

 !master checks to see if we want to reuse loss function that is already saved
 !in a file; if so master checks for existence of a suitable file. if all looks
 !good, we can use it; else recalculate. after figureing all this out, master
 !broadcasts her decision to everybody.
 if (master) then
   if (irdloss(1).eq.0) then
     iinv=1
   else
     iinv=0
     open(unit=26,file=trim(bindir)//'losfn.bin',status='old',form='unformatted',iostat=ios)
     if (ios.eq.0) then
       read(26) npwct,npwlft
       if (npwct.ne.npwc) iinv=1
       if (npwlft.ne.npwlf) iinv=1
       if (iinv.ne.0) close(26)
     else
       iinv=1
     endif
   endif
 end if
 call par_bcast_int(iinv,1,0)

! debug
!stop
! debug
 if (iinv.ne.0) then
   write(6,'(a)') 'Inverting dielectric function'
   !openmp loop over iw (frequency); mpi loop over iqpt (q-points)

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iw,iqpt,dielf,dielf2,lossfn2)
   !$OMP DO SCHEDULE(dynamic)
   do iw=1,nwpt
     if(mod(iw,100)==0) write(*,*) 'iw: ',iw
     !do iqpt=1,nqpt
     if ((.not.master).or.(numprocs==1)) then
       !non-master nodes do the q-points if there are any; else master does this
       !part
       do iqpt=iqpt_start,iqpt_end
         call invertdielf(pola(iw,iqpt,:),polb(iw,iqpt,:),ipwndx(:,1:ntpwndx), &
              npwndx,ntpwndx,ipwx(:,1:npwc),npwc,nwpt,nqpt,&
              lossfn_loc(iw,iqpt,:)) 
       enddo
     end if
     if (master) then
       !master node does tensor part in any event
       call invertdielfT(pola(iw,nqpt+1:nqpt+9,:),polb(iw,nqpt+1:nqpt+9,:), &
            ipwndx(:,1:ntpwndx), &
            npwndx,ntpwndx,ipwx(:,1:npwc),npwc,nwpt,nqpt,&
            lossfn_loc(iw,nqpt+1:nqpt+9,:))
            !lossfn(iw,nqpt+1:nqpt+9,:))
     end if
   enddo
   !$OMP END DO 
   !$OMP END PARALLEL
   !at this point each mpi node has computed the loss function for its
   !q-points; collect full loss function here
   if (.not.master) then
     !worker sends her work to master
     call par_send_dc(lossfn_loc,size(lossfn_loc),master,my_rank)
   else
     if (numprocs/=1) then
       !master collects contributions from all workers
       do iproc=1,numprocs-1
         call MPE_DECOMP1D(nqpt,numprocs-1,iproc-1,iqpt_start,iqpt_end)
         call par_recv_dc&
         (lossfn(:,iqpt_start:iqpt_end,:),size(lossfn(:,iqpt_start:iqpt_end,:)),iproc,iproc)
       end do
     else
       write(*,*) &
       'before loading global loss function, its size is ', size(lossfn)
       !there were no other workers, so master had to do everything herself
       lossfn=lossfn_loc
       write(*,*) &
       'global loss function loaded from local one. size: ', size(lossfn)
     end if
   end if
   if(allocated(lossfn_loc)) deallocate(lossfn_loc)
   !master now has full loss function; write it to disk
   if (master) then
     !write(*,*) 'Writing loss funciton', size(lossfn)
     open(unit=26,file=trim(bindir)//'losfn.bin',form='unformatted',status='unknown')
     write(26) npwc,npwlf
     write(26) lossfn
     close(26)
   endif
 endif
 !everybody wait until the loss function is written
 call par_barrier
 !now we can read the full loss function
 write(6,fmt='(a,i0,a)') 'MPI proc ',my_rank,' Reading loss function'
 inquire(file=trim(bindir)//'losfn.bin',number=unit_num)
 if (unit_num>=0) close(unit_num)
 open(unit=26,file=trim(bindir)//'losfn.bin',form='unformatted',status='unknown')
 read(26) npwc,npwlf
 read(26) lossfn
 close(26)
! iw=54
! iqpt=1
! call invertdielf(pola(iw,iqpt,:),polb(iw,iqpt,:),ipwndx(:,1:ntpwndx),npwndx,ntpwndx,ipwx(1:npwc),npwc,nwpt,nqpt,dielf,lossfn(iw,iqpt,:))
!stop
   write(6,fmt='(a,i0,a)') 'MPI proc ',my_rank,": loss function found."
   if (iprtlosfn(1).ne.0.and.master) then
     write(6,*) "Printing loss function."
     call prtlosfn(iprtlosfn,npwc,npwlf,ipwndx,napwndx,nqpt,qpt,qtensor,ngkpt,nwpt,ntpwndx,lossfn,dw)
   end if
   if(master) write(6,*) "dielectric properties printed."
!   if (iprtdielf(1).ne.0.or.iprtlosfn(1).ne.0) then
!     call par_end
!     stop
!   end if

   call par_barrier
   if (calc_xtdielf.gt.0.or.calc_pwploss.gt.0.or.calc_pwse.gt.0) then
     nxtqpt = 0
     do iqpt=1,nqpt
       do ii=1,3
         ikk(ii) = int(qpt(ii,iqpt)*2*ngkpt(ii)+0.5)
       enddo
       if (((ikk(1).eq.ngkpt(1)).or.(ikk(1).eq.0)).and. &
&          ((ikk(2).eq.ngkpt(2)).or.(ikk(2).eq.0)).and. &
&          ((ikk(3).eq.ngkpt(3)).or.(ikk(3).eq.0))) then
         nxtqpt = nxtqpt+1
       endif
     enddo
     if (iprtxtdielf(5).eq.0) iprtxtdielf(5)=nxtqpt
     if (iprtxtlosfn(5).eq.0) iprtxtlosfn(5)=nxtqpt
     allocate (xtqpt(3,nxtqpt),xtqpt2qpt(nxtqpt),qpt2xtqpt(nqpt),qpt2xtg(nqpt,3))
     jj = 0
     do iqpt=1,nqpt
       do ii=1,3
         ikk(ii) = int(qpt(ii,iqpt)*2*ngkpt(ii)+0.5)
       enddo
       if (((ikk(1).eq.ngkpt(1)).or.(ikk(1).eq.0)).and. &
&          ((ikk(2).eq.ngkpt(2)).or.(ikk(2).eq.0)).and. &
&          ((ikk(3).eq.ngkpt(3)).or.(ikk(3).eq.0))) then
         jj = jj + 1
         xtqpt(:,jj) = qpt(:,iqpt)
         xtqpt2qpt(jj) = iqpt
         qpt2xtqpt(iqpt) = jj
write(6,'(1x,2i4,3f10.6,3i3)') iqpt,nqpt,qpt(:,iqpt),ikk
       endif
     enddo
     do ix=1,nxtgkpt(1)
     do iy=1,nxtgkpt(2)
     do iz=1,nxtgkpt(3)
       ikx=(/ix,iy,iz/)
       ikk = ikx*ngkpt/nxtgkpt
       iqpt = iqndx(ikk(1),ikk(2),ikk(3))
       jj = qpt2xtqpt(iqpt)
       ixtqndx(ikx(1),ikx(2),ikx(3)) = jj
!write(6,'(1x,3i3,2x,i3,2x,i3,2x,i3)') ikx,iqpt,ixtqndx(ikx(1),ikx(2),ikx(3)),qpt2xtqpt(iqpt)
     enddo
     enddo
     enddo
     do iqpt=1,nqpt
       do ii=1,3
         ikk(ii) = int(qpt(ii,iqpt)*2*ngkpt(ii)+0.5)
       enddo
       ikx = (((ikk+ngkpt)/ngkpt)*ngkpt)/2
       do ii=1,3
         qpt2xtg(iqpt,ii) = 0
         if (ikx(ii).le.0) then
           ikx(ii) = ikx(ii) + ngkpt(ii)
           qpt2xtg(iqpt,ii) = -1
         endif
       enddo
       qpt2xtqpt(iqpt) = ixtqndx((2*ikx(1))/ngkpt(1),(2*ikx(2))/ngkpt(2),(2*ikx(3))/ngkpt(3))
!write(6,'(1x,2i4,3f10.6,3i3,2x,3i3,2x,i3)') iqpt,nqpt,qpt(:,iqpt),ikk,ikx,qpt2xtqpt(iqpt)
     enddo
     allocate(xtpola(nwpt,nxtqpt,npwx),xtpolb(nwpt,nxtqpt,npwx),xtlossfn(nwpt,nxtqpt,npwx))
     do ipw1=1,npwx
!do ipw1=1,3
       write(chpw1,'(i2)') ipw1
       if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
       ipwx1=ipwx(1,ipw1)
       ipwy1=ipwx(2,ipw1)
       ipwz1=ipwx(3,ipw1)
       write(chpwx1,'(i3)') abs(ipwx1)
       if (chpwx1(2:2).eq.' ') chpwx1(2:2)='0'
       if (ipwx1.ge.0) then
         chpwx1(1:1)='+'
       else
         chpwx1(1:1)='-'
       endif
       write(chpwy1,'(i3)') abs(ipwy1)
       if (chpwy1(2:2).eq.' ') chpwy1(2:2)='0'
       if (ipwy1.ge.0) then
         chpwy1(1:1)='+'
       else
         chpwy1(1:1)='-'
       endif
       write(chpwz1,'(i3)') abs(ipwz1)
       if (chpwz1(2:2).eq.' ') chpwz1(2:2)='0'
       if (ipwz1.ge.0) then
         chpwz1(1:1)='+'
       else
         chpwz1(1:1)='-'
       endif
       strpw1=chpwx1//','//chpwy1//','//chpwz1
       fileoutr='xtpol'//strpw1//'.bin'
       fileout3=trim(bindir)//trim(fileoutr)
       open(unit=25,file=trim(fileout3),status='old',form='unformatted',iostat=ios)
       if (ios.ne.0.or.ipolcalc.ne.0) then
         write(6,'(a,2x,i4,a,i4)') "calculating extended polarization "//trim(fileoutr),ipw1," of ",npwx
         do iqpt=1,nxtqpt
           do ii=1,3
             iks(ii)=nint(xtqpt(ii,iqpt)*ngkpt(ii))
             iqq(ii)=iks(ii)+ipwx(ii,ipw1)*ngkpt(ii)
             qq(ii)=dble(iqq(ii))/dble(ngkpt(ii))
           enddo
           call fqq(shiftk,shiftk,iqq,ngkpt,qq,qs)
           write(6,'(i3,a,i3,a,3(f6.3,a))') iqpt," of ",nxtqpt, &
&            ",  q  = (",qq(1),",",qq(2),",",qq(3),")"
  
           qq2=0.d0
           do ii=1,3
           do jj=1,3
             qq2=qq2+qq(ii)*bmet(ii,jj)*qq(jj)
           enddo
           enddo
           vq=4.d0*pi/qq2
  
           iqs=iks+ngkpt/2
  
!           call mkxtpol(iqpt,qq,ipwx(:,ipw1),iks,vq,vol,pi,bmet,blat, &
!&           pwse_mass, &
!&           test_bands_pol, &
!&           nsppol,shiftk,npw,igmx,igmn,ipaw,itetrahedron, &
!&           nsym,nkpt, &
!&           ncg,bantot,npwt,nlmnmax, &
!&           nxtwpt,xtwmax,nbcore,nbocc,ncband,ngkpt,natom, &
!&           xred, &
!&           projwf, &
!&           kg, &
!&           enrgy, eigen, &
!&           cg, &
!&           indxkpw, &
!&           indxkbnd, &
!&           indxkcg, &
!&           npwarr, &
!&           kpt, &
!&           nsymk,symk,symrel,syminv, &
!&           ihlf,lvtrans, &
!&           ntypepaw, &
!&           pwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)), &
!&           igndx,ikndx,isymndx,isymg, &
!&           nband,xtpola(:,iqpt,ipw1),xtpolb(:,iqpt,ipw1))
           call mkpol1(qq,ipwx(:,ipw1),ipwx(:,ipw1),iks,vq,vol,pi, &
&           test_bands_pol, &
&           nsppol,shiftk,npw,igmx,igmn,ipaw,itetrahedron, &
&           nsym,nkpt, &
&           ncg,bantot,npwt,nlmnmax, &
&           nxtwpt,xtwmax,nbcore,nbocc,ncband,ngkpt,natom, &
&           xred, &
&           projwf, &
&           kg, &
&           enrgy, eigen, &
&           cg, &
&           indxkpw, &
&           indxkbnd, &
&           indxkcg, &
&           npwarr, &
&           kpt, &
&           nsymk,symk,symrel,syminv, &
&           ihlf,lvtrans, &
&           ntypepaw, &
&           pwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,ipw1,iqs(1),iqs(2),iqs(3)), &
&           pwmatel(:,:,:,ipw2,iqs(1),iqs(2),iqs(3)),tpwmatel(:,:,:,ipw2,iqs(1),iqs(2),iqs(3)), &
&           igndx,ikndx,isymndx,isymg, &
&           nband,xtpola(:,iqpt,ipw1),xtpolb(:,iqpt,ipw1))
         enddo
         open(unit=25,file=trim(fileout3),form='unformatted',status='unknown',iostat=ios)
         if (ios.ne.0) then
           write(6,*) "Error - cannot open file ",trim(fileout3)," for writing"
           stop
         endif
         write(25) xtpola(:,:,ipw1),xtpolb(:,:,ipw1)
         close(25)
       endif
       close(25)
     enddo
     call par_barrier
     do ipw1=1,npwx
       write(chpw1,'(i2)') ipw1
       if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
       ipwx1=ipwx(1,ipw1)
       ipwy1=ipwx(2,ipw1)
       ipwz1=ipwx(3,ipw1)
       write(chpwx1,'(i3)') abs(ipwx1)
       if (chpwx1(2:2).eq.' ') chpwx1(2:2)='0'
       if (ipwx1.ge.0) then
         chpwx1(1:1)='+'
       else
         chpwx1(1:1)='-'
       endif
       write(chpwy1,'(i3)') abs(ipwy1)
       if (chpwy1(2:2).eq.' ') chpwy1(2:2)='0'
       if (ipwy1.ge.0) then
         chpwy1(1:1)='+'
       else
         chpwy1(1:1)='-'
       endif
       write(chpwz1,'(i3)') abs(ipwz1)
       if (chpwz1(2:2).eq.' ') chpwz1(2:2)='0'
       if (ipwz1.ge.0) then
         chpwz1(1:1)='+'
       else
         chpwz1(1:1)='-'
       endif
       strpw1=chpwx1//','//chpwy1//','//chpwz1
       fileoutr='xtpol'//strpw1//'.bin'
       fileout3=trim(bindir)//trim(fileoutr)
!write(6,*) fileout3
       open(unit=25,file=trim(fileout3),status='old',form='unformatted',iostat=ios)
       if (ios.ne.0) then
         write(6,*) "Error - cannot open file ",trim(fileout3)," for reading"
         stop
       endif
       read(25) xtpola(:,:,ipw1),xtpolb(:,:,ipw1)
       close(25)
       do iqpt=1,nxtqpt
         do iw=0,nwpt
           xtlossfn(iw,iqpt,ipw1) = 1.d0/(1.d0-dble(xtpola(iw,iqpt,ipw1))-dimag(xtpolb(iw,iqpt,ipw1))*(0.d0,1.d0))
         enddo
       enddo
     enddo
     write(6,*) "Extended polarizability function found."
     if (iprtxtdielf(1).ne.0.and.master) then 
       write(6,*) "Printing extended dielectric function."
       call prtxtdielf(iprtxtdielf,ipwndx,napwndx,nxtqpt,xtqpt,ngkpt,nxtwpt,npwndx,xtpola,xtpolb,dxtw)
     end if
     if (iprtxtlosfn(1).ne.0.and.master) then 
       write(6,*) "Printing extended loss function."
       call prtxtlosfn(iprtxtlosfn,ipwndx,napwndx,nxtqpt,xtqpt,ngkpt,nxtwpt,npwndx,xtlossfn,dxtw)
     end if
   endif

   allocate(cse(ncband,nqpt),xse(ncband,nqpt),se(ncband,nqpt),zz(ncband,nqpt),isecalculated(ncband,nqpt))
   if (calc_se.gt.0.or.calc_ploss.gt.0) then
!!!!! Calculate self energy and rate of scattering
     fileout3=trim(bindir)//"selfenergy.bin"
     open(unit=26,file=trim(fileout3),status='old',form='unformatted',iostat=ios)
     if (ios.ne.0) then
       tcalc=1
     else
       read(26) tncband,tnqpt
       read(26) t_test_bands_se(1),t_test_bands_se(2)
       if (tncband.ne.ncband.or.tnqpt.ne.nqpt.or. &
&          t_test_bands_se(1).ne.test_bands_se(1).or. &
&          t_test_bands_se(2).ne.test_bands_se(2)) then
         tcalc=1
       else
         tcalc=0
       endif
     endif
     if (tcalc.eq.0) then
       read(26) isecalculated
       read(26) cse,xse,zz
     else
       do ikpt=1,nqpt
         do iband=1,ncband
           isecalculated(iband,ikpt)=0
           cse(iband,ikpt)=(0.d0,0.d0)
           xse(iband,ikpt)=0.d0
           zz(iband,ikpt)=1.d0
         enddo
       enddo
     endif
     write(6,'(a)') 'calculating correlation self energy'
     open (unit=27,file='se.dat',status='unknown')
     write(27,'(2a)') " band  kpt        energy    eigenvalue          cse", &
& "                 xse         Z                E_xc        SE"
     open (unit=28,file='imcse.dat',status='unknown')
     write(28,'(a)') " band  kpt        energy                     self energy"
     do ikpt=sekptlo,sekpthi
!     do ikpt=1,nqpt
!     do ikpt=1,1
!     do ikpt=8,8
       do ii=1,3
         ikk(ii)=nint((qpt(ii,ikpt)+0.5d0)*dble(ngkpt(ii))-shiftk(ii))
       enddo
       ikpt2=ikndx(ikk(1),ikk(2),ikk(3))
       write(6,'(a,3f10.4)') '        k-point ',qpt(:,ikpt)
       !call MPE_DECOMP1D(sebandhi-sebandlo+1,numprocs,my_rank,iqpt_start,iqpt_end)
       !iqpt_start=iqpt_start+sebandlo-1
       !iqpt_end=iqpt_end+sebandlo-1
       !do iband=iqpt_start,iqpt_end
       do iband=sebandlo,sebandhi
!       do iband=1,20
!write(6,*) "twenak ",iband,ikpt,isecalculated(iband,ikpt)
         if (isecalculated(iband,ikpt).eq.0) then
           call mkse(iband,ikk,lossfn, &
&               vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmnmax, &
&               test_bands_se, &
&               ntypepaw, &
&               pwmatel,tpwmatel,pwjmatel,tpwjmatel, &
&               kg,eigen,cg,npwt,bantot,ncg, &
&               indxkpw,indxkbnd,indxkcg,npwarr, &
&               kpt,nkpt,nqpt,nsymk,symk,nsym,symrel,syminv, &
&               ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
&               ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
&               npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
&               igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
&               nband,nsppol,shiftk,zz(iband,ikpt),cse(iband,ikpt),xse(iband,ikpt))
           isecalculated(iband,ikpt)=1.d0
         endif
         !write(6,'(i3,4f16.8)') iband,xse(iband,ikpt)*hart,cse(iband,ikpt)*hart
!  write(6,*)
         if (ivxc.eq.1) then
!           call mkxclda(iband,ikpt,ngfft,vxc,vol,ihlf,nkpt,igndx,igmx,igmn,natom, &
!    &           exc(iband,ikpt))
!           write(6,'(i3,4f16.8)') iband,exc(iband,ikpt)*hart
           write(27,'(2i5,2f14.6,3x,7f10.5)') iband,ikpt, &
&     enrgy(indxkbnd(ikpt2)+iband)*hart, &
&     eigen(indxkbnd(ikpt2)+iband)*hart, &
&     cse(iband,ikpt)*hart,xse(iband,ikpt)*hart, &
&     zz(iband,ikpt), &
&     exc(iband,ikpt)*hart, &
&     (xse(iband,ikpt)+dble(cse(iband,ikpt))*dble(zz(iband,ikpt))-exc(iband,ikpt))*hart
         endif
         if (ivxc.eq.1) then
           se(iband,ikpt) = dble(cse(iband,ikpt))*dble(zz(iband,ikpt))+xse(iband,ikpt)-exc(iband,ikpt)
           write(28,'(2i5,f14.6,3x,2f20.12)') iband,ikpt, &
&     enrgy(indxkbnd(ikpt2)+iband)*hart, &
&     (dble(cse(iband,ikpt))*dble(zz(iband,ikpt))+xse(iband,ikpt)-exc(iband,ikpt)) &
&     *hart, &
&     dimag(cse(iband,ikpt))*dble(zz(iband,ikpt))*hart
         else
           se(iband,ikpt) = dble(cse(iband,ikpt))*dble(zz(iband,ikpt))+xse(iband,ikpt)
           write(28,'(2i5,f14.6,3x,2f20.12)') iband,ikpt, &
&     enrgy(indxkbnd(ikpt2)+iband)*hart, &
&     (dble(cse(iband,ikpt))*dble(zz(iband,ikpt))+xse(iband,ikpt))*hart, &
&     dimag(cse(iband,ikpt))*dble(zz(iband,ikpt))*hart
         endif
!         write(6,*)
!!         rate=2.d0*abs(dimag(cse(iband,ikpt)))
!!       write(6,'(a,"(",es12.4,",",es12.4,")",a)') '      correlation self energy ',cse(iband,ikpt)*hart,' eV'
!!       write(6,'(a,es12.4,a)') '      exchange self energy ',xse(iband,ikpt)*hart,' eV'
!!       write(6,'(a,"(",es12.4,",",es12.4,")",a)') '      total self energy ',(xse(iband,ikpt)+cse(iband,ikpt))*hart,' eV'
!!       write(6,'(a,es12.4,a)') '      DFT-LDA exchange-correlation energy ',exc(iband,ikpt)*hart,' eV'
!!       write(6,'(a,"(",es12.4,",",es12.4,")",a)') '      energy difference ',((xse(iband,ikpt)+cse(iband,ikpt))-exc(iband,ikpt))*hart,' eV'
!!       write(6,'(a,"(",es12.4,",",es12.4,")",a)') '      energy correction ',((xse(iband,ikpt)+cse(iband,ikpt))-exc(iband,ikpt))*zz(iband,ikpt)*hart,' eV'
!!         write(6,'(a,es12.4,a)') '      interaction rate ',rate/atomictime,' s^(-1)'
!!         ie=indxkbnd(ikpt)+iband
!!         speed=0.d0
!!         do ii=1,3
!!           speed=speed+egrad(ii,ie)**2
!!         enddo
!!         speed=sqrt(speed)
!!         mfp=speed/rate
!!         xsect=1.d0/(mfp*atomdens)
!!         write(6,'(a,es12.4,a)') '      mean free path ',mfp*bohr,' m'
!!         write(6,'(a,es12.4,a)') '      cross section ',xsect*bohr**2,' m^2'
       enddo
     enddo
     close(27)
     close(28)
     close(26)

     open(unit=26,file=trim(fileout3),status='unknown',form='unformatted',iostat=ios)
     write(26) ncband,nqpt
     write(26) test_bands_se(1),test_bands_se(2)
     write(26) isecalculated
     write(26) cse,xse,zz
     close(26)
   endif
   if (calc_ploss.gt.0) then
!!!!! Calculate rate of energy loss
     allocate(ploss(ncband,nqpt))
     fileout3=trim(bindir)//"ploss.bin"
     open(unit=26,file=trim(fileout3),status='old',form='unformatted',iostat=ios)
     if (ios.ne.0) then
       tcalc=1
     else
       read(26) tncband,tnqpt
       read(26) t_test_bands_se(1),t_test_bands_se(2)
       if (tncband.ne.ncband.or.tnqpt.ne.nqpt.or. &
&          t_test_bands_se(1).ne.test_bands_se(1).or. &
&          t_test_bands_se(2).ne.test_bands_se(2)) then
         tcalc=1
       else
         tcalc=0
       endif
     endif
     if (tcalc.eq.0) then
       read(26) isecalculated
       read(26) ploss
     else
       do ikpt=1,nqpt
         do iband=1,ncband
           isecalculated(iband,ikpt)=0
           ploss(iband,ikpt)=0.d0
         enddo
       enddo
     endif

     write(6,'(a)') 'calculating power loss'
     open (unit=27,file='ploss.dat',status='unknown')
     do ikpt=sekptlo,sekpthi
!     do ikpt=1,nqpt
!     do ikpt=1,1
!     do ikpt=8,8
       do ii=1,3
         ikk(ii)=nint((qpt(ii,ikpt)+0.5d0)*dble(ngkpt(ii))-shiftk(ii))
       enddo
       ikpt2=ikndx(ikk(1),ikk(2),ikk(3))
       write(6,'(a,3f10.4)') '        k-point ',qpt(:,ikpt)
       do iband=sebandlo,sebandhi
!       do iband=1,20
!  ! power loss calculated in mkploss.f90
!write(6,*) "qopf ",se(iband,ikpt),ploss(iband,ikpt)
         if (isecalculated(iband,ikpt).eq.0) then
           call mkploss(iband,ikk,lossfn, &
&             vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmnmax, &
&             test_bands_se, &
&             ntypepaw, &
&             pwmatel,tpwmatel,pwjmatel,tpwjmatel, &
&             kg,eigen,cg,npwt,bantot,ncg, &
&             indxkpw,indxkbnd,indxkcg,npwarr, &
&             kpt,nkpt,nqpt,nsymk,symk,nsym,symrel,syminv, &
&             ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
&             ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
&             npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
&             igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
&             nband,nsppol,shiftk,se(iband,ikpt),ploss(iband,ikpt))
           isecalculated(iband,ikpt)=1.d0
         endif
!write(6,*) "kwls ",se(iband,ikpt),ploss(iband,ikpt)
         write(6,*) "ploss: ",ploss(iband,ikpt)," au"
!         write(6,*) "ploss: ",ploss(iband,ikpt)*hart/atomictime," eV/second"
         write(27,'(2i5,f14.6,3x,e20.6)') iband,ikpt, &
&     enrgy(indxkbnd(ikpt2)+iband)*hart, &
&     ploss(iband,ikpt) !, ploss(iband,ikpt)*hart/atomictime
       enddo
     enddo
     close(27)
     close(26)

     open(unit=26,file=trim(fileout3),status='unknown',form='unformatted',iostat=ios)
     write(26) ncband,nqpt
     write(26) test_bands_se(1),test_bands_se(2)
     write(26) isecalculated
     write(26) ploss
     close(26)
     deallocate(ploss)
   endif

   if (calc_pwse.gt.0.or.calc_pwploss.gt.0) then
!!!!! Calculate self energy and rate of scattering in plane wave approximation
     mode = 1
     allocate(pwcse(npwx,nqpt),pwxse(npwx,nqpt),pwse(npwx,nqpt),pwzz(npwx,nqpt),ipwsecalculated(npwx,nqpt))
     fileout3=trim(bindir)//"pwselfenergy.bin"
     open(unit=26,file=trim(fileout3),status='old',form='unformatted',iostat=ios)
     if (ios.ne.0) then
       tcalc=1
     else
       read(26) tnpwx,tnqpt
       if (tnpwx.ne.npwx.or.tnqpt.ne.nqpt) then
         tcalc=1
       else
         tcalc=0
       endif
     endif
     if (tcalc.eq.0) then
       read(26) ipwsecalculated
       read(26) pwcse,pwxse,pwzz
     else
       do ikpt=1,nqpt
         do ipw=1,npwx
           ipwsecalculated(ipw,ikpt)=0
           pwcse(ipw,ikpt)=(0.d0,0.d0)
           pwxse(ipw,ikpt)=0.d0
           pwzz(ipw,ikpt)=1.d0
         enddo
       enddo
     endif
     write(6,'(a)') 'calculating correlation self energy in plane wave approximation'
     open (unit=27,file='pwse.dat',status='unknown')
     open (unit=28,file='imcpwse.dat',status='unknown')
     if (pwse_nmax.ge.pwse_nmin) then
       do lkpt=pwse_nmin,pwse_nmax
         kv = pwse_kstp*lkpt
         igg = nint(kv)
         kr = kv-igg
         ikk = nint((kr+0.5d0)*dble(ngkpt))
         ek=0.d0
         do ii=1,3
         do jj=1,3
           ek=ek+kv(ii)*bmet(ii,jj)*kv(jj)/(2.d0*pwse_mass)
         enddo
         enddo
         call mkpwse(igg,ikk,mode,lossfn, &
&                 pwse_mass, &
&                 vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmnmax, &
&                 test_bands_se, &
&                 ntypepaw, &
&                 pwmatel,tpwmatel,pwjmatel,tpwjmatel, &
&                 kg,eigen,cg,npwt,bantot,ncg, &
&                 indxkpw,indxkbnd,indxkcg,npwarr, &
&                 kpt,nkpt,nqpt,nsymk,symk,nsym,symrel,syminv, &
&                 ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
&                 ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
&                 npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
&                 igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
&                 nband,nsppol,shiftk,pwzz_x,pwdcse_x,pwcse_x,pwxse_x)
         pwlo = npwc+1
         pwhi = npwx
!xtpwxse_x = 1.2345e6
!write(6,*) pwcse_x,hart
         call mkxtpwse(igg,ikk,mode,xtlossfn, &
&                 pwse_mass, &
&                 vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmnmax, &
&                 test_bands_se, &
&                 ntypepaw, &
&                 pwmatel,tpwmatel,pwjmatel,tpwjmatel, &
&                 kg,eigen,cg,npwt,bantot,ncg, &
&                 indxkpw,indxkbnd,indxkcg,npwarr, &
&                 kpt,nkpt,nqpt,nxtqpt,nsymk,symk,nsym,symrel,syminv, &
&                 ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
&                 ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
&                 pwlo,pwhi,npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
&                 igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
&                 nxtgkpt,xtqpt,xtqpt2qpt,qpt2xtqpt,qpt2xtg,ixtqndx, &
&                 nband,nsppol,shiftk,xtpwzz_x,xtpwdcse_x,xtpwcse_x,xtpwxse_x)
!write(6,*) pwcse_x,hart
         pwzz_x = 1.d0/(1.d0 - pwdcse_x - xtpwdcse_x)
         pwse_x = dble(pwcse_x + xtpwcse_x)*dble(pwzz_x)
         write(29,'(i3,a,i3,5f16.8)') lkpt," of ",pwse_nmax,ek*hart,pwcse_x*hart,xtpwcse_x*hart ! for convergence tests
         write(6,'(i3,a,i3,4f16.8)') lkpt," of ",pwse_nmax,ek*hart,(pwcse_x + xtpwcse_x)*dble(pwzz_x)*hart
         write(28,'(2i5,f14.6,3x,3f20.12)') lkpt,ikpt, &
&           (ek+pwse_x)*hart, &
&           (pwcse_x + xtpwcse_x)*dble(pwzz_x)*hart, &
&           dble(pwzz_x)  ! the useful data
         write(27,'(2i5,2f14.6,3x,7f10.5)') lkpt,ikpt, &
&           (ek+pwse_x)*hart, &
&           ek*hart, &
&           (pwcse_x + xtpwcse_x)*hart, &
&           pwzz_x, &
&           (pwcse_x + xtpwcse_x)*dble(pwzz_x)*hart
!if (lkpt.ge.2) stop
       enddo
     else
       do ikpt=sekptlo,sekpthi
         do ii=1,3
           ikk(ii)=nint((qpt(ii,ikpt)+0.5d0)*dble(ngkpt(ii))-shiftk(ii))
         enddo
         ikpt2=ikndx(ikk(1),ikk(2),ikk(3))
         write(6,'(a,3f10.4)') '        k-point ',qpt(:,ikpt)
         do ipw=pwse_pwlo,pwse_pwhi
           igg = ipwx(:,ipw)
           ek=0.d0
           do ii=1,3
           do jj=1,3
             ek=ek+(qpt(ii,ikpt)+igg(ii))*bmet(ii,jj)*(qpt(jj,ikpt)+igg(jj))/(2.d0*pwse_mass)
           enddo
           enddo
           if (ipwsecalculated(ipw,ikpt).eq.0) then
             call mkpwse(igg,ikk,mode,lossfn, &
&                 pwse_mass, &
&                 vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmnmax, &
&                 test_bands_se, &
&                 ntypepaw, &
&                 pwmatel,tpwmatel,pwjmatel,tpwjmatel, &
&                 kg,eigen,cg,npwt,bantot,ncg, &
&                 indxkpw,indxkbnd,indxkcg,npwarr, &
&                 kpt,nkpt,nqpt,nsymk,symk,nsym,symrel,syminv, &
&                 ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
&                 ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
&                 npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
&                 igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
&                 nband,nsppol,shiftk,pwzz(ipw,ikpt),pwdcse,pwcse(ipw,ikpt),pwxse(ipw,ikpt))
             ipwsecalculated(ipw,ikpt)=1.d0
           endif
           write(6,'(i3,a,i3,4f16.8)') ipw," of ",pwse_pwhi,ek*hart,pwcse(ipw,ikpt)*hart
           pwse(ipw,ikpt) = dble(pwcse(ipw,ikpt))*dble(pwzz(ipw,ikpt))
           write(28,'(2i5,f14.6,3x,2f20.12)') ipw,ikpt, &
&       (ek+pwse(ipw,ikpt))*hart, &
&       pwse(ipw,ikpt)*hart, &
&       dimag(cse(iband,ikpt))*dble(zz(iband,ikpt))*hart
           write(27,'(2i5,2f14.6,3x,7f10.5)') ipw,ikpt, &
&       (ek+pwse(ipw,ikpt))*hart, &
&       ek*hart, &
&       pwcse(ipw,ikpt)*hart, &
&       pwzz(ipw,ikpt), &
&       pwse(ipw,ikpt)
!        enddo
         enddo
         if (ipwsecalculated(ipw,ikpt).eq.0) then
           call mkpwse(ipw,ikk,lossfn, &
&               vol,pi,nwpt,wmax,nbcore,nbocc,ncband,ngkpt,natom,xred,projwf,nlmnmax, &
&               test_bands_se, &
&               ntypepaw, &
&               pwmatel,tpwmatel,pwjmatel,tpwjmatel, &
&               kg,eigen,cg,npwt,bantot,ncg, &
&               indxkpw,indxkbnd,indxkcg,npwarr, &
&               kpt,nkpt,nqpt,nsymk,symk,nsym,symrel,syminv, &
&               ihlf,lvtrans,bmet,blat,ipaw,itetrahedron, &
&               ipwx,ipwndx,npwndx,ntpwndx,napwndx, &
&               npwc,npwx,invpw2ndx,pwsymndx,iqsymndx, &
&               igmx,igmn,igndx,ikndx,iqndx,isymndx,isymg,npw, &
&               nband,nsppol,shiftk,pwzz(ipw,ikpt),pwdcse,pwcse(ipw,ikpt),pwxse(ipw,ikpt))
           ipwsecalculated(ipw,ikpt)=1.d0
         endif
         write(6,'(i3,a,i3,4f16.8)') ipw," of ",npwx,pwxse(ipw,ikpt)*hart,pwcse(ipw,ikpt)*hart
         pwse(ipw,ikpt) = dble(pwcse(ipw,ikpt))*dble(pwzz(ipw,ikpt))+pwxse(ipw,ikpt)
         write(28,'(2i5,f14.6,3x,2f20.12)') ipw,ikpt, &
&     (ek+pwse(ipw,ikpt))*hart, &
&     pwse(ipw,ikpt)*hart, &
&     dimag(cse(iband,ikpt))*dble(zz(iband,ikpt))*hart
         write(27,'(2i5,2f14.6,3x,7f10.5)') ipw,ikpt, &
&     (ek+pwse(ipw,ikpt))*hart, &
&     ek*hart, &
&     pwcse(ipw,ikpt)*hart,pwxse(ipw,ikpt)*hart, &
&     pwzz(ipw,ikpt), &
&     pwse(ipw,ikpt)
       enddo
     endif
     close(27)
     close(28)
     close(26)

     open(unit=26,file=trim(fileout3),status='unknown',form='unformatted',iostat=ios)
     write(26) npwx,nqpt
     write(26) ipwsecalculated
     write(26) pwcse,pwxse,pwzz
     close(26)
   endif
 endif

call par_barrier
call par_end
stop
end program mkprm

!**************************************************************************

subroutine rdvxc(fname,ipaw,ngfft,npsp,nspden,vxc)
implicit none
character(80) :: fname
integer :: ipaw,ngfft(3),npsp,nspden
double precision :: vxc(ngfft(1)*ngfft(2)*ngfft(3))
integer :: ii,ispden,ir,ipsp
character*6 :: cdum
integer :: idum1,idum2,idum3,idum4,idum5,idum6,idum7,iadum(3)

 open(unit=10,file=fname,status='old',form='unformatted')
 read(10) ! codvsn,headform,fform

 read(10) ! bantot,date,intxc,ixc,natom,ngfft(1:3),&
!& nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
!& ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear

 read(10) !istwfk(1:nkpt),nband(1:nkpt*nsppol),&
!& npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym), &
!& symrel(1:3,1:3,1:nsym), &
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

 do ispden=1,nspden
   read(10) (vxc(ir),ir=1,ngfft(1)*ngfft(2)*ngfft(3))
 enddo

 close(10)

end subroutine rdvxc

!**************************************************************************

subroutine rdgwvxc(fname,ncband,nqpt,exc,gwbandmax)
implicit none
character(80) :: fname
character(80) :: line
integer :: ncband,nqpt
double precision :: exc(ncband,nqpt)
integer :: ii,ios,iend,iwrd,iband,iqpt,gwbandmax
double precision :: E0,VxcLDA,SigX,SigC,ZZ,dSigC,Sig,EmE0,EE

 iqpt=0
 open(unit=9,file=fname,status='old',iostat=ios)
 do
   read(9,'(a)',iostat=ios) line
   if (ios.ne.0) exit
   iend=len(line)
   iwrd=index(line(:iend),'<VxcLDA>')
   if (iwrd.ne.0) then
     iqpt=iqpt+1
     gwbandmax=0
     do
!       read(line(iwrd+6:iend),*) nbocc
       read(9,'(a)') line
       if (len(line).lt.70) exit
       if (gwbandmax.ge.ncband) exit
       read(line,*,iostat=ios) iband,E0,VxcLDA,SigX,SigC,ZZ,dSigC,Sig,EmE0,EE
       if (ios.ne.0) exit
       exc(iband,iqpt)=VxcLDA
       gwbandmax=iband
     enddo
     if (gwbandmax.lt.ncband) then
       do iband=gwbandmax+1,ncband
         exc(iband,iqpt)=exc(iband-1,iqpt)
       enddo
     endif
   endif
 enddo

 close(10)

end subroutine rdgwvxc

!**************************************************************************

subroutine rdden(fname,ipaw,ngfft,npsp,nspden,rhor)
implicit none
character(80) :: fname
integer :: ngfft(3),npsp,nspden,ipaw
double precision :: rhor(ngfft(1)*ngfft(2)*ngfft(3))
integer :: ii,ispden,ir,ipsp

 open(unit=10,file=fname,status='old',form='unformatted')
 read(10) ! codvsn,headform,fform

 read(10) ! bantot,date,intxc,ixc,natom,ngfft(1:3),&
!& nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
!& ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear

 read(10) !istwfk(1:nkpt),nband(1:nkpt*nsppol),&
!& npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym), &
!& symrel(1:3,1:3,1:nsym), &
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

 do ispden=1,nspden
   read(10) (rhor(ir),ir=1,ngfft(1)*ngfft(2)*ngfft(3))
 enddo

 close(10)

end subroutine rdden

!**************************************************************************

subroutine rdwfk(fnamewf1,fnameinp1,ipaw,ncband)
use wfkvars
use geometry
implicit none
character(80) :: fnamewf1,fnameinp1
integer :: ipaw,ncband
character(20) :: dir
character*80 :: title1,title2,line
double precision :: hart,ryd
parameter (ryd=13.6056981,hart=2.d0*ryd)
integer :: ik,ib,ips,il,ipwt,ish,is,jtypat,ieof,ikn,ikwfm,jps,kat,in,ioffmx, &
& isym,ioff,jloop,nln,iloop,iwrite,i,j,k,ig,iban,ikpt,iat,ipsp,ibwfm, &
& ict,ibndt,iatom,ispden,cplex,nspden2,selmax
double precision :: wfmax
integer :: ii,npw2,nspinor2,nbnd2,ibnd,isppol,ios,lv,lend
integer, allocatable :: occ2(:)
integer, allocatable :: nrhoijsel(:,:),rhoijselect(:,:,:)
double complex, allocatable :: rhoijp(:,:,:)
integer :: idummy

 open(unit=10,file=trim(fnamewf1),status='old',form='unformatted',iostat=ios)
 if (ios.ne.0) then
   write(6,*) "wave function file ",trim(fnamewf1)," cannot be found."
   stop
 endif
 read(10) codvsn,headform,fform

 read(10) bantot,date,intxc,ixc,natom,ngfft(1:3),&
& nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
& ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear

 allocate (istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),so_typat(ntypat),&
& symafm(nsym),symrel(3,3,nsym),typat(natom))
 allocate (kpt(3,nkpt),occ(bantot),tnons(3,nsym),znucltypat(ntypat), &
& xred(3,natom))
 read(10) istwfk(1:nkpt),nband(1:nkpt*nsppol),&
& npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym),symrel(1:3,1:3,1:nsym), &
& typat(1:natom),kpt(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym), &
& znucltypat(1:ntypat)

 nln=6
 do ipsp=1,npsp
! (npsp lines, 1 for each pseudopotential ; npsp=ntypat, except if alchemical pseudo-atoms)
  read(10) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc
 enddo

!(final record: residm, coordinates, total energy, Fermi energy)
 read(10) residm,xred(1:3,1:natom),etotal,fermie
!write(6,*) residm
!write(6,*) xred(1:3,1:natom)
!write(6,*) etotal
!write(6,*) fermie

 if (ipaw.ne.0) read(10) idummy
 if (ipaw.ne.0) read(10) idummy

 npwt=sum(npwarr)
 ncg=0
 do ikpt=1,nkpt
   if (ncband.gt.nband(ikpt)) then
     write(6,*) "Not enough bands in wave function input file"
     write(6,*) "Bands in file ",fnamewf1,": ",nband(ikpt)
     write(6,*) "Bands in computaion: ",ncband
     stop
   endif
   ncg=ncg+npwarr(ikpt)*nband(ikpt)
 enddo
 npw=maxval(npwarr)
 nbnd=maxval(nband)

 allocate (kg(3,npwt),eigen(bantot),cg(ncg),occ2(bantot))
 allocate (indxkbnd(nkpt),indxkpw(nkpt),indxkcg(nkpt))
 ibndt=0
 ipwt=0
 ict=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     indxkbnd(ikpt)=ibndt
     indxkpw(ikpt)=ipwt
     indxkcg(ikpt)=ict
     read(10) npw2,nspinor2,nbnd2
!write(6,*) npw2,nspinor2,nbnd2
     read(10) kg(1:3,1+ipwt:npw2+ipwt)
! kg(ikpt,ipw)=kg(indxkpw(ikpt)+ipw)
     read(10) eigen(1+ibndt:nbnd2+ibndt),occ2(1+ibndt:nbnd2+ibndt)
! en(ikpt,ibnd)=eigen(indxkbnd(ikpt)+ibnd)
     do ibnd=1,nbnd2
       read(10) (cg(ii+ict),ii=1,nspinor2*npw2)
       ict=ict+nspinor2*npw2
! cg(ikpt,ibnd,ipw)=cg(indxkcg(ikpt)+(ibnd-1)*npwarr(ikpt)+ipw)
     enddo
     ibndt=ibndt+nbnd2
     ipwt=ipwt+npw2
   enddo
 enddo
 deallocate(occ2)

 close(10)

 open(unit=10,file=fnameinp1,status='old',iostat=ios)
 if (ios.ne.0) then
   write(6,*) "input file ",fnameinp1," cannot be found."
   stop
 endif
 do
   read(10,'(a)',iostat=ios) line
   if (ios.ne.0) exit
   lend=index(line,'#')
   if (lend.eq.0) lend=len(line)
   lv=index(line(1:lend),'ngkpt')
   if (lv.ne.0) then
     read(line(lv+5:lend),*) ngkpt(1:3)
   endif
   lv=index(line(1:lend),'kptopt')
   if (lv.ne.0) then
     read(line(lv+6:lend),*) kptopt
   endif
   if (index(line(1:lend),'nshiftk').ne.0) cycle
   lv=index(line(1:lend),'shiftk')
   if (lv.ne.0) then
     read(line(lv+6:lend),*) shiftk(1:3)
   endif
 enddo
 close(10)

end subroutine rdwfk

!***************************************************************************
! q(ii)=k(1)*blat(1,ii)+k(2)*blat(2,ii)+k(3)*blat(3,ii)

subroutine mkgeom(pi)
use geometry
implicit none
double precision :: pi
integer :: ic,jc,iwrite,ik
double precision :: tempmat(3,3)

 vol=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(2,3)*rprimd(3,2)) &
& -rprimd(1,2)*(rprimd(2,1)*rprimd(3,3)-rprimd(2,3)*rprimd(3,1)) &
& +rprimd(1,3)*(rprimd(2,1)*rprimd(3,2)-rprimd(2,2)*rprimd(3,1))
 vbz=(2.d0*pi)**3/vol
 blat(1,1)= (rprimd(2,2)*rprimd(3,3)-rprimd(2,3)*rprimd(3,2))*(2.d0*pi)/vol
 blat(1,2)=-(rprimd(2,1)*rprimd(3,3)-rprimd(2,3)*rprimd(3,1))*(2.d0*pi)/vol
 blat(1,3)= (rprimd(2,1)*rprimd(3,2)-rprimd(2,2)*rprimd(3,1))*(2.d0*pi)/vol
 blat(2,1)= (rprimd(3,2)*rprimd(1,3)-rprimd(3,3)*rprimd(1,2))*(2.d0*pi)/vol
 blat(2,2)=-(rprimd(3,1)*rprimd(1,3)-rprimd(3,3)*rprimd(1,1))*(2.d0*pi)/vol
 blat(2,3)= (rprimd(3,1)*rprimd(1,2)-rprimd(3,2)*rprimd(1,1))*(2.d0*pi)/vol
 blat(3,1)= (rprimd(1,2)*rprimd(2,3)-rprimd(1,3)*rprimd(2,2))*(2.d0*pi)/vol
 blat(3,2)=-(rprimd(1,1)*rprimd(2,3)-rprimd(1,3)*rprimd(2,1))*(2.d0*pi)/vol
 blat(3,3)= (rprimd(1,1)*rprimd(2,2)-rprimd(1,2)*rprimd(2,1))*(2.d0*pi)/vol

 do ic=1,3
   do jc=1,3
     rmet(ic,jc)=rprimd(ic,1)*rprimd(jc,1)+rprimd(ic,2)*rprimd(jc,2)+rprimd(ic,3)*rprimd(jc,3)
     bmet(ic,jc)=blat(ic,1)*blat(jc,1)+blat(ic,2)*blat(jc,2)+blat(ic,3)*blat(jc,3)
   enddo
 enddo
 tempmat = bmet
 call inverse(tempmat,bmetinv,3,3)

return
end subroutine mkgeom

!***************************************************************************

subroutine mkiptind(shiftk,nsym,nkpt,kpt,nkc,ist,ind,symrel,iptsym,ipttr,idg,iptind)
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
! do ix=4,4
   xkc1(1)=(dble(ix)+shiftk(1))/dble(nkc)
   do iy=ist,ind
!   do iy=-4,-4
     xkc1(2)=(dble(iy)+shiftk(2))/dble(nkc)
     do iz=ist,ind
!     do iz=5,5
       xkc1(3)=(dble(iz)+shiftk(3))/dble(nkc)
       iptind(ix,iy,iz)=0
l1:    do ik=1,nkpt
         do isym1=1,nsym*2
           if (isym1.le.nsym) then
             isym=isym1
             do ic=1,3
               xkc(ic)=xkc1(ic)
!               if (xkc(ic).gt.(0.5-tol)) xkc(ic)=xkc(ic)-1.0
!               if (xkc(ic).le.-(0.5-tol)) xkc(ic)=xkc(ic)+1.0
               if (nint((xkc(ic)-shiftk(ic))*nkc).gt.nint(0.5*nkc)) xkc(ic)=xkc(ic)-1.0
               if (nint((xkc(ic)-shiftk(ic))*nkc).le.nint(-0.5*nkc)) xkc(ic)=xkc(ic)+1.0
             enddo
           else
             isym=isym1-nsym
             do ic=1,3
               xkc(ic)=-xkc1(ic)
!               if (xkc(ic).gt.(0.5-tol)) xkc(ic)=xkc(ic)-1.0
!               if (xkc(ic).le.-(0.5-tol)) xkc(ic)=xkc(ic)+1.0
               if (nint((xkc(ic)-shiftk(ic))*nkc).gt.nint(0.5*nkc)) xkc(ic)=xkc(ic)-1.0
               if (nint((xkc(ic)-shiftk(ic))*nkc).le.nint(-0.5*nkc)) xkc(ic)=xkc(ic)+1.0
             enddo
           endif
           do ic=1,3
             ykc(ic)=0.d0
             do jc=1,3
               ykc(ic)=ykc(ic)+symrel(jc,ic,isym)*kpt(jc,ik)
             enddo
!             if (ykc(ic).gt.(0.5-tol)) then
             if (nint((ykc(ic)-shiftk(ic))*nkc).gt.nint(0.5*nkc)) then
               ykc(ic)=ykc(ic)-1.0
               idg(ic,ix,iy,iz)=-1
!             elseif (ykc(ic).le.-(0.5-tol)) then
             elseif (nint((ykc(ic)-shiftk(ic))*nkc).le.nint(-0.5*nkc)) then
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
!           write(6,'(2i3,3f6.2,3x,3f6.2,3x,f10.6)') ik,isym,ykc,xkc,diff
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

subroutine rdinput(calc_dielf, calc_se, calc_ploss, calc_xtdielf, calc_pwse, calc_pwploss, calc_scspec, &
& nbcore,nbocc,nbsev,nbsec,nwpt,wmax,nxtwpt,xtwmax,ncband,epwx,epwc,epwlf,epwg, &
& test_bands_se,test_bands_pol, &
& ngrid,nvalel,nprefsym,iprefsym, &
& fnamewf1,fnameinp1,fnamekss1,fnamekss2,fnameq,fnamevxc,fnamegwvxc,fnamecore,bindir, &
& iscissors,scissors,igw,sebandlo,sebandhi,sekptlo,sekpthi,pwse_pwlo,pwse_pwhi, &
& pwse_mass, pwse_nmin, pwse_nmax, pwse_kstp, &
& iprtsym,iprtdielf,iprtlosfn,iprtxtdielf,iprtxtlosfn,ipaw,ipolcalc,ipolpw,ipolqpt, &
& iprtpolpw,irdloss,ivxc,igwvxc,icore,itetrahedron, &
& sskpt, sspw, ssnrg, ssdnrg)
use pawvars
implicit none
character(280) :: line,line2
character(80) :: fnamewf1,fnameinp1,fnameq,bindir
character(80) :: fnamekss1,fnamekss2,fnamegwvxc,fnamevxc,fnamecore
integer :: calc_dielf, calc_se, calc_ploss, calc_xtdielf, calc_pwse, calc_pwploss, calc_scspec
integer :: ncband,nbocc,nbcore,nbsev,nbsec,nwpt,ngcut,ngrid(3),nvalel,natom,nprefsym,iprefsym(20)
integer :: nxtwpt
integer :: iscissors,igw,sebandlo,sebandhi,sekptlo,sekpthi,iprtsym, &
& iprtdielf(5),iprtlosfn(5),iprtxtdielf(5),iprtxtlosfn(5),irdloss(2), &
& ivxc,igwvxc,icore,itetrahedron
integer :: pwse_pwlo, pwse_pwhi
double precision ::  pwse_mass
integer :: pwse_nmin, pwse_nmax
double precision ::  pwse_kstp(3)
integer :: ipaw,ipolcalc,ipolpw(2),ipolqpt(2),iprtpolpw
integer :: test_bands_se(2),test_bands_pol(4)
integer :: iwrd,iend,ii,ios
double precision :: wmax,xtwmax,epwx,epwlf,epwc,epwg,xtemp,scissors,hart
 double precision :: ssnrg, ssdnrg
 integer :: sskpt, sspw, start

 hart=27.2113845d0            ! atomic unit of energy = 27.2 eV
 nwpt=300
 wmax=30
 nxtwpt=300
 xtwmax=300.d0
 ncband=20
 nbcore=0
 nbsev=0
 nbsec=0
 test_bands_se=(/0,0/)
 test_bands_pol=(/0,0,0,0/)
 ngrid=(/10,10,10/)
 nvalel=1
 natom=1
 epwlf=0.d0
 epwc=0.d0
 epwx=0.d0
 epwg=0.d0
 nprefsym=0
 iprefsym=(/(0*ii,ii=1,20)/)
 bindir='./'
 igw=0
 iscissors=0
 scissors=0.d0
 sebandlo=0
 sebandhi=0
 sekptlo=0
 sekpthi=0
 pwse_pwlo=0
 pwse_pwhi=0
 pwse_mass=-1.0
 pwse_nmin=0
 pwse_nmax=-1
 pwse_kstp=(/1.0,0.0,0.0/)
 iprtsym=1
 iprtdielf=(/0,0,0,0,0/)
 iprtlosfn=(/0,0,0,0,0/)
 iprtxtdielf=(/0,0,0,0,0/)
 iprtxtlosfn=(/0,0,0,0,0/)
 irdloss=(/1,1/) ! read from lossfn.bin if it exists, write it if it does not.
 ipaw=0
 ipolcalc=0
 ipolpw=(/1,0/)
 ipolqpt=(/1,0/)
 iprtpolpw=0
 ivxc=0 ! find LDA exchange correlation potential?
 igwvxc=0 ! read LDA exchange correlation potential from ABINIT GW calcuation?
 itetrahedron=1 ! use tetrahedron integration
 icore=0 ! perform core polarization corrections?
 ssnrg=-1.d0
 ssdnrg=1.d0
 sskpt=-1
 sspw=-1
 calc_dielf = 0
 calc_se = 0
 calc_ploss = 0
 calc_pwse = 0
 calc_xtdielf = 0
 calc_pwploss = 0
 calc_scspec = 0
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
   iwrd=index(line(:iend),'calc_dielf')
   if (iwrd.ne.0) then
     read(line(iwrd+11:iend),*) calc_dielf
   endif
   iwrd=index(line(:iend),'calc_se')
   if (iwrd.ne.0) then
     read(line(iwrd+8:iend),*) calc_se
   endif
   iwrd=index(line(:iend),'calc_ploss')
   if (iwrd.ne.0) then
     read(line(iwrd+11:iend),*) calc_ploss
   endif
   iwrd=index(line(:iend),'calc_xtdielf')
   if (iwrd.ne.0) then
     read(line(iwrd+12:iend),*) calc_xtdielf
   endif
   iwrd=index(line(:iend),'calc_pwse')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) calc_pwse
   endif
   iwrd=index(line(:iend),'calc_pwploss')
   if (iwrd.ne.0) then
     read(line(iwrd+13:iend),*) calc_pwploss
   endif
   iwrd=index(line(:iend),'calc_scspec')
   if (iwrd.ne.0) then
     read(line(iwrd+12:iend),*) calc_scspec
   endif
   write(*,*) trim(line(1:iend))
   iwrd=index(line(:iend),'nbocc')
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) nbocc
   endif
   iwrd=index(line(:iend),'nbcore')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),*) nbcore
   endif
   iwrd=index(line(:iend),'nbsev')  ! number of valence bands to include in BSE
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) nbsev
   endif
   iwrd=index(line(:iend),'nbsec')  ! number of conduction bands in BSE
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) nbsec
   endif
   iwrd=index(line(:iend),'nwpt')
   if (iwrd.ne.0) then
     read(line(iwrd+5:iend),*) nwpt
   endif
   iwrd=index(line(:iend),'wmax')
   if (iwrd.ne.0) then
     read(line(iwrd+5:iend),*) wmax
   endif
   iwrd=index(line(:iend),'nxtwpt')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),*) nxtwpt
   endif
   iwrd=index(line(:iend),'xtwmax')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),*) xtwmax
   endif
   iwrd=index(line(:iend),'ncband')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),*) ncband
   endif
   iwrd=index(line(:iend),'epwlf')
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) epwlf
   endif
   iwrd=index(line(:iend),'epwc')
   if (iwrd.ne.0) then
     read(line(iwrd+5:iend),*) epwc
   endif
   iwrd=index(line(:iend),'epwg')
   if (iwrd.ne.0) then
     read(line(iwrd+5:iend),*) epwg
   endif
   iwrd=index(line(:iend),'epwx')
   if (iwrd.ne.0) then
     read(line(iwrd+5:iend),*) epwx
   endif
   iwrd=index(line(:iend),'test_bands_se')
   if (iwrd.ne.0) then
     read(line(iwrd+14:iend),*) test_bands_se(1:2)
   endif
   iwrd=index(line(:iend),'test_bands_pol')
   if (iwrd.ne.0) then
     read(line(iwrd+15:iend),*) test_bands_pol(1:4)
   endif
   iwrd=index(line(:iend),'ngrid')
   iwrd=index(line(:iend),'ngrid')
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) ngrid(1:3)
   endif
   iwrd=index(line(:iend),'nvalel')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),*) nvalel
   endif
   iwrd=index(line(:iend),'natom')
   if (iwrd.ne.0) then
     read(line(iwrd+6:iend),*) natom
   endif
   iwrd=index(line(:iend),'nprefsym')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),*) nprefsym
   endif
   iwrd=index(line(:iend),'iprefsym')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),*) iprefsym(1:nprefsym)
   endif
   iwrd=index(line(:iend),'fnamewf1')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),'(a)') fnamewf1
   endif
   iwrd=index(line(:iend),'fnameinp1')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),'(a)') fnameinp1
   endif
   iwrd=index(line(:iend),'fnamekss1')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),'(a)') fnamekss1
   endif
   iwrd=index(line(:iend),'fnamekss2')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),'(a)') fnamekss2
   endif
   iwrd=index(line(:iend),'fnameq')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),'(a)') fnameq
   endif
   iwrd=index(line(:iend),'fnamepaw')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),*) ntypepaw
     ipaw=1
     allocate(fnamepaw(ntypepaw))
     do ii=1,ntypepaw
       read(9,'(a)') line2
       start=1
       do 
         if (line2(start:start).eq.' ') then
           start=start+1
         else 
           exit
         endif
       enddo
       fnamepaw(ii) = line2(start:)
     enddo
   endif
   iwrd=index(line(:iend),'fnamevxc')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),'(a)') fnamevxc
     ivxc=1
   endif
   iwrd=index(line(:iend),'fnamegwvxc')
   if (iwrd.ne.0) then
     read(line(iwrd+11:iend),'(a)') fnamegwvxc
     igwvxc=1
   endif
   iwrd=index(line(:iend),'fnamecore')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),'(a)') fnamecore
     icore=1
   endif
   iwrd=index(line(:iend),'bindir')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),'(a)') bindir
     ii=len_trim(bindir)
     if (bindir(ii:ii).ne.'/') then
       bindir(ii+1:ii+1)='/'
     endif
   endif
   iwrd=index(line(:iend),'scissors')
   if (iwrd.ne.0) then
     iscissors=1
     read(line(iwrd+9:iend),*) scissors
     if (index(line,'eV').ne.0) scissors=scissors/hart
   endif
   iwrd=index(line(:iend),'igw')
   if (iwrd.ne.0) then
     read(line(iwrd+4:iend),*) igw
   endif
   iwrd=index(line(:iend),'itetrahedron')
   if (iwrd.ne.0) then
     read(line(iwrd+13:iend),*) itetrahedron
   endif
   iwrd=index(line(:iend),'sebandlo')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),*) sebandlo
   endif
   iwrd=index(line(:iend),'sebandhi')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),*) sebandhi
   endif
   iwrd=index(line(:iend),'sekptlo')
   if (iwrd.ne.0) then
     read(line(iwrd+8:iend),*) sekptlo
   endif
   iwrd=index(line(:iend),'sekpthi')
   if (iwrd.ne.0) then
     read(line(iwrd+8:iend),*) sekpthi
   endif
   iwrd=index(line(:iend),'pwse_pwlo')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) pwse_pwlo
   endif
   iwrd=index(line(:iend),'pwse_pwhi')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) pwse_pwhi
   endif
   iwrd=index(line(:iend),'pwse_mass')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) pwse_mass
   endif
   iwrd=index(line(:iend),'pwse_kstp')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) pwse_kstp(1:3)
   endif
   iwrd=index(line(:iend),'pwse_nmin')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) pwse_nmin
   endif
   iwrd=index(line(:iend),'pwse_nmax')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) pwse_nmax
   endif
   iwrd=index(line(:iend),'iprtsym')
   if (iwrd.ne.0) then
     read(line(iwrd+8:iend),*) iprtsym
   endif
   iwrd=index(line(:iend),'iprtsym')
   if (iwrd.ne.0) then
     read(line(iwrd+8:iend),*) iprtsym
   endif
   iwrd=index(line(:iend),'iprtdielf')
   if (iwrd.ne.0) then
     iprtdielf(1)=1
     read(line(iwrd+10:iend),*) iprtdielf(2:5)
   endif
   iwrd=index(line(:iend),'iprtlosfn')
   if (iwrd.ne.0) then
     iprtlosfn(1)=1
     read(line(iwrd+10:iend),*) iprtlosfn(2:5)
   endif
   iwrd=index(line(:iend),'iprtxtdielf')
   if (iwrd.ne.0) then
     iprtxtdielf(1)=1
     read(line(iwrd+12:iend),*) iprtxtdielf(2:5)
   endif
   iwrd=index(line(:iend),'iprtxtlosfn')
   if (iwrd.ne.0) then
     iprtxtlosfn(1)=1
     read(line(iwrd+12:iend),*) iprtxtlosfn(2:5)
   endif
   iwrd=index(line(:iend),'irdloss')
   if (iwrd.ne.0) then
     read(line(iwrd+8:iend),*) irdloss
   endif
   iwrd=index(line(:iend),'ipolcalc')
   if (iwrd.ne.0) then
     read(line(iwrd+9:iend),*) ipolcalc
   endif
   iwrd=index(line(:iend),'ipolpw')
   if (iwrd.ne.0) then
     read(line(iwrd+7:iend),*) ipolpw
   endif
   iwrd=index(line(:iend),'ipolqpt')
   if (iwrd.ne.0) then
     read(line(iwrd+8:iend),*) ipolqpt
   endif
   iwrd=index(line(:iend),'iprtpolpw')
   if (iwrd.ne.0) then
     read(line(iwrd+10:iend),*) iprtpolpw
   endif
   iwrd=index(line(:iend),'secondaryspectrum')
   if (iwrd.ne.0) then
     read(line(iwrd+18:iend),*) sskpt, sspw, ssnrg, ssdnrg
     ssnrg=ssnrg/hart
     ssdnrg=ssdnrg/hart
   endif
 enddo
 close(9)

 if (epwc.lt.epwlf) then
   epwc=epwlf
 endif
 if (epwx.lt.epwc) then
   epwx=epwc
 endif

 if (sebandlo.eq.0) sebandlo=nbcore+1
 if (sebandhi.eq.0) sebandhi=ncband
 if (pwse_pwlo.eq.0) pwse_pwlo=1
 if (test_bands_se(1).eq.0) test_bands_se(1)=nbcore+1
 if (test_bands_se(2).eq.0) test_bands_se(2)=ncband
 if (test_bands_pol(1).eq.0) test_bands_pol(1)=nbcore+1
 if (test_bands_pol(2).eq.0) test_bands_pol(2)=nbocc
 if (test_bands_pol(3).eq.0) test_bands_pol(3)=nbocc+1
 if (test_bands_pol(4).eq.0) test_bands_pol(4)=ncband
 if (pwse_mass.gt.0) then 
   pwse_mass = pwse_mass * 1822.8884845
 else 
   pwse_mass = 1.0
 endif

return
end subroutine rdinput

!**************************************************************************

subroutine fndsymk(symrel,ipsymndx,kpt,nkpt,nsym,ngkpt,shiftk,kptopt,nsymk,symk,ikndx,isymndx,isymg)
! for each k-point ik, nsymk(ik) is the number of distinguishable symmetry 
! operations of the lattice, and symk(ik,isym) gives the list of distinguishable
! symmetries such that symrel(:,:,symk(ik,isym)) enumerates all possible 
! k-points obtainable by symmetry operations exactly once for each 
! isym=1,nsymk(ik)
implicit none
integer :: nkpt,nsym,kptopt
integer :: symrel(3,3,nsym),ngkpt(3),ipsymndx(nsym)
double precision :: kpt(3,nkpt),shiftk(3)
integer :: nsymk(nkpt),symk(nkpt,nsym*2),ika(3),ikb(3)
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymg(3,nkpt,nsym)
integer :: ik,isym,jsym,iisym,jjsym,ii,jj,kk,iset,itrsign,iwrite,ix,iy,iz
double precision :: ksym(3,nsym*2),ksymsh(3),dk(3),dk2,tol,resk(3),ksym2(3,nsym*2),ksym0(3),ksymt(3)

iwrite=0

do ix=1,ngkpt(1)
do iy=1,ngkpt(2)
do iz=1,ngkpt(3)
  ikndx(ix,iy,iz)=0
  isymndx(ix,iy,iz)=0
enddo
enddo
enddo

tol=1.d-8
if (kptopt.eq.1) then
  do ik=1,nkpt
!  do ik=46,46
    nsymk(ik)=0
    do isym=1,nsym*2
      ksym(:,isym)=(/0.d0,0.d0,0.d0/)
      if (isym.gt.nsym) then
        itrsign=-1   ! time reversal symmetry
        iisym=ipsymndx(isym-nsym)
!        iisym=isym-nsym
      else
        itrsign=1
        iisym=ipsymndx(isym)
!        iisym=isym
      endif
      iset=1
      do ii=1,3
        do jj=1,3
          ksym(ii,isym)=ksym(ii,isym)+symrel(jj,ii,iisym)*kpt(jj,ik)*itrsign
!          ksym(ii,isym)=ksym(ii,isym)+symrel(ii,jj,iisym)*kpt(jj,ik)*itrsign
        enddo
        if (ksym(ii,isym).gt.0.5d0+1.d0/(2.d0*dble(ngkpt(ii)))) then
          ksym(ii,isym)=ksym(ii,isym)-int(ksym(ii,isym)+0.5d0)
        elseif (ksym(ii,isym).le.-0.5d0+1.d0/(2.d0*dble(ngkpt(ii)))) then
          ksym(ii,isym)=-abs(ksym(ii,isym))+int(abs(ksym(ii,isym))+0.5d0)
        endif
        ksymsh(ii)=ksym(ii,isym)-shiftk(ii)/ngkpt(ii)
        ika(ii)=nint((ksymsh(ii)+0.5d0)*ngkpt(ii))
        if (ika(ii).eq.0) then
          ksym(ii,isym)=ksym(ii,isym)+1.d0
          ksymsh(ii)=ksymsh(ii)+1.d0
          ika(ii)=ika(ii)+ngkpt(ii)
        endif
        resk(ii)=ksym(ii,isym)-((dble(ika(ii))+shiftk(ii))/dble(ngkpt(ii))-0.5d0)
        ksym2(ii,iisym+nsym*((1-itrsign)/2))=ksym(ii,isym)
      enddo
!      write(6,'(4i4,3f5.1,2x,3i3)') ik,isym,iisym,itrsign,ksym(1:3,isym),ika
!      write(6,'(8x,3i3)') ngkpt(1:3)
!      write(6,'(8x,3f10.6)') shiftk(1:3)
!      write(6,'(8x,3f10.6)') ksymsh(1:3)
!      write(6,'(8x,3i3)') ika(1:3)
!      write(6,'(8x,3f10.6)') resk(1:3)
      if ((resk(1)**2+resk(2)**2+resk(3)**2).gt.tol) then
        iset=0
        cycle
      endif
      do jsym=1,isym-1
        do ii=1,3
          dk(ii)=ksym(ii,isym)-ksym(ii,jsym)
        enddo
        dk2=dk(1)**2+dk(2)**2+dk(3)**2
!        write(6,'(3i3,3(3f6.2,2x),f10.6)') ik,isym,jsym,ksym(1:3,isym),ksym(1:3,jsym),dk(1:3),dk2
        if (dk2.le.tol) then
          iset=0
          exit
        endif
      enddo
      if (iset.eq.1) then
        nsymk(ik)=nsymk(ik)+1
!        write(6,'(3i8)') ik,isym,nsymk(ik)
        symk(ik,nsymk(ik))=iisym+nsym*((1-itrsign)/2)
        do ii=1,3
          ika(ii)=mod(ika(ii)+ngkpt(ii)-1,ngkpt(ii))+1
        enddo
!        write(6,'(4i8,5x,3i3)') ik,isym,iisym+nsym*((1-itrsign)/2),nsymk(ik),ika
!        write(6,*) ika
        ikndx(ika(1),ika(2),ika(3))=ik
        isymndx(ika(1),ika(2),ika(3))=iisym+nsym*((1-itrsign)/2)
        ksym0=(/0.d0,0.d0,0.d0/)
        do ii=1,3
          do jj=1,3
            ksym0(ii)=ksym0(ii)+symrel(jj,ii,isym)*kpt(jj,ik)
          enddo
          ksymt(ii)=dble(ika(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
          ksymt(ii)=mod(ksymt(ii),1.d0)
          if (ksymt(ii).lt.0.d0) ksymt(ii)=ksymt(ii)+1.d0
          ksymt(ii)=1.d0-mod(1.d0-ksymt(ii),1.d0)
          ksymt(ii)=ksymt(ii)-0.5d0
        enddo
        isymg(:,ik,isym)=nint(ksymt(:)-ksym0(:))
!        write(6,'(2i3,3x,3i3,3x,3f7.3,3x,3f7.3)') ik,isym,isymg(:,ik,isym),ksym0,ksymt
      endif
    enddo
    if (iwrite.eq.1) then
      write(47,*) ik
      do isym=1,nsymk(ik)
        write(47,'(3i3,2x,3f9.5,2x,3i3)') symk(ik,isym), &
&          mod(symk(ik,isym)-1,nsym)+1, &
&          -2*((symk(ik,isym)-1)/nsym)+1,ksym2(1:3,symk(ik,isym)), &
&          isymg(1:3,ik,symk(ik,isym))
      enddo
      write(47,*)
    endif
  enddo
elseif (kptopt.eq.3) then
  do ik=1,nkpt
    nsymk(ik)=1
    symk(ik,1)=1
    do ii=1,3
      ika(ii)=nint((kpt(ii,ik)+0.5d0)*ngkpt(ii)-shiftk(ii))
      ikb(ii)=mod(ika(ii)+ngkpt(ii)-1,ngkpt(ii))+1
    enddo
    ikndx(ikb(1),ikb(2),ikb(3))=ik
    isymndx(ikb(1),ikb(2),ikb(3))=1
    do isym=1,nsym
      isymg(:,ik,isym)=(/0,0,0/)
    enddo
    if (iwrite.eq.1) then
      write(47,'(i4,4x,3i4,3x,3i4,3f9.5)') ik,ika(1:3),ikb(1:3),kpt(1:3,ik)
    endif
  enddo
else
  write(6,'(a)') 'unsupported kptopt in abinit calculations'
  write(6,'(a)') 'change kptopt in abinit input file to 1 or 3 and try again'
  stop
endif
close(46)

return
end subroutine fndsymk

!**************************************************************************

subroutine figmxn(kg,npwt,igmx,igmn)
implicit none
integer :: npwt
integer :: kg(3,npwt)
integer :: igmx(3),igmn(3)
integer :: ii,jj

igmx=(/0,0,0/)
igmn=(/0,0,0/)
do ii=1,npwt
  do jj=1,3
    igmx(jj)=max(igmx(jj),kg(jj,ii))
    igmn(jj)=min(igmn(jj),kg(jj,ii))
  enddo
enddo

return
end subroutine figmxn

!**************************************************************************

subroutine mkigndx(kg,indxkpw,npwarr,igmx,igmn,nkpt,npwt,igndx)
implicit none
integer :: nkpt,npwt,igmx(3),igmn(3)
integer :: kg(3,npwt),indxkpw(nkpt),npwarr(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ix,iy,iz,ik,ipw,ipwt

do ix=igmn(1),igmx(1)
do iy=igmn(2),igmx(2)
do iz=igmn(3),igmx(3)
  do ik=1,nkpt
    igndx(ik,ix,iy,iz)=0
  enddo
enddo
enddo
enddo

do ik=1,nkpt
  do ipw=1,npwarr(ik)
    ipwt=indxkpw(ik)+ipw
    igndx(ik,kg(1,ipwt),kg(2,ipwt),kg(3,ipwt))=ipw
  enddo
enddo

return
end subroutine mkigndx

!**************************************************************************

subroutine mkigkndx(gbig,npwkss,igmx,igmn,igkndx)
implicit none
integer :: igmx(3),igmn(3),npwkss
integer :: gbig(3,npwkss)
integer :: igkndx(igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: ix,iy,iz,ik,ipw

do ix=igmn(1),igmx(1)
do iy=igmn(2),igmx(2)
do iz=igmn(3),igmx(3)
  igkndx(ix,iy,iz)=0
enddo
enddo
enddo

do ipw=1,npwkss
  igkndx(gbig(1,ipw),gbig(2,ipw),gbig(3,ipw))=ipw
enddo

return
end subroutine mkigkndx

!**************************************************************************

subroutine fndwf(ikk,rr,ib,ngkpt,wfb,iwfg,npwarr,nkpt,ikndx,npw,npwt,blat,wf)
implicit none
integer :: ikk(3),ib,nkpt,npw,npwt,ngkpt(3)
double precision :: rr(3),xk(3),xg(3),blat(3,3),xkgdr
integer :: npwarr(nkpt),iwfg(3,npw,ngkpt(1),ngkpt(2),ngkpt(3)),ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: wfb(npw,ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: wf
integer :: ii,jj,ig,ikpt

wf=(0.d0,0.d0)
ikpt=ikndx(ikk(1),ikk(2),ikk(3))
do ii=1,3
  xk(ii)=0.d0
  do jj=1,3
    xk(ii)=xk(ii)+(dble(ikk(jj))/dble(ngkpt(jj))-0.5d0)*blat(ii,jj)
  enddo
enddo
!write(6,*) ikpt,npwarr(ikpt)
do ig=1,npwarr(ikpt)
  do ii=1,3
    xg(ii)=0.d0
    do jj=1,3
      xg(ii)=xg(ii)+iwfg(jj,ig,ikk(1),ikk(2),ikk(3))*blat(ii,jj)
    enddo
  enddo
  do ii=1,3
    xkgdr=xkgdr+(xk(ii)+xg(ii))*rr(ii)
  enddo
  wf=wf+wfb(ig,ikk(1),ikk(2),ikk(3))*exp((0.d0,1.d0)*xkgdr)
!  write(6,*) wf
enddo

return
end subroutine fndwf

!**************************************************************************

subroutine fqq(shiftk,shiftkq,iks,ngkpt,qq,qs)
implicit none
integer :: ngkpt(3),iks(3)
double precision :: shiftk(3),shiftkq(3)
double precision :: qq(3),qs(3)
integer :: ii

do ii=1,3
  qs(ii)=(shiftk(ii)-shiftkq(ii))/ngkpt(ii)
  qq(ii)=(shiftk(ii)-shiftkq(ii)+iks(ii))/ngkpt(ii)
enddo

return
end

!**************************************************************************

subroutine mkpwx(epwx,bmet)
use pwx
implicit none
double precision :: epwx,bmet(3,3)
double precision :: eg
integer :: igg(3),npwxold
integer :: is,ix,iy,ic,id,ipw
integer, allocatable :: indx(:),itemp(:,:)
double precision, allocatable :: xtemp(:)

is=0
npwx=1
allocate (ipwx(3,npwx))
ipwx(1:3,1)=(/0,0,0/)

do
  is=is+1
  npwxold=npwx
  igg=(/is,is,is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  igg=(/-is,is,is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  igg=(/is,-is,is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  igg=(/-is,-is,is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  igg=(/is,is,-is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  igg=(/-is,is,-is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  igg=(/is,-is,-is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  igg=(/-is,-is,-is/)
  call  fpwe(igg,bmet,eg)
  if (eg.le.epwx) then
    call expandpwx(igg)
  endif
  do ix=-is+1,is-1
    igg=(/ix,is,is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/ix,-is,is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/ix,is,-is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/ix,-is,-is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/is,ix,is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/-is,ix,is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/is,ix,-is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/-is,ix,-is/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/is,is,ix/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/-is,is,ix/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/is,-is,ix/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    igg=(/-is,-is,ix/)
    call  fpwe(igg,bmet,eg)
    if (eg.le.epwx) then
      call expandpwx(igg)
    endif
    do iy=-is+1,is-1
      igg=(/ix,iy,is/)
      call  fpwe(igg,bmet,eg)
      if (eg.le.epwx) then
        call expandpwx(igg)
      endif
      igg=(/ix,iy,-is/)
      call  fpwe(igg,bmet,eg)
      if (eg.le.epwx) then
        call expandpwx(igg)
      endif
      igg=(/ix,is,iy/)
      call  fpwe(igg,bmet,eg)
      if (eg.le.epwx) then
        call expandpwx(igg)
      endif
      igg=(/ix,-is,iy/)
      call  fpwe(igg,bmet,eg)
      if (eg.le.epwx) then
        call expandpwx(igg)
      endif
      igg=(/is,ix,iy/)
      call  fpwe(igg,bmet,eg)
      if (eg.le.epwx) then
        call expandpwx(igg)
      endif
      igg=(/-is,ix,iy/)
      call  fpwe(igg,bmet,eg)
      if (eg.le.epwx) then
        call expandpwx(igg)
      endif
    enddo
  enddo
  if (npwx.eq.npwxold) exit
enddo

allocate(engpw(npwx),indx(npwx),itemp(3,npwx),xtemp(npwx))
do ipw=1,npwx
  call fpwe(ipwx(:,ipw),bmet,engpw(ipw))
enddo
call indxhpsort(npwx,npwx,engpw,indx)
do ipw=1,npwx
  itemp(:,ipw)=ipwx(:,indx(ipw))
  xtemp(ipw)=engpw(indx(ipw))
enddo
ipwx=itemp
engpw=xtemp

! write(6,*) npwx
! do ipw=1,npwx
!   write(6,'(i4,3x,3i3,3x,2f14.6,3x,3i3)') ipw,ipwx(:,ipw),engpw(ipw),xtemp(ipw),itemp(:,ipw)
! enddo

!deallocate(engpw,indx,itemp,xtemp)
deallocate(indx,itemp,xtemp)

return
end subroutine mkpwx

!**************************************************************************

subroutine fpwe(igg,bmet,eg)
implicit none
integer :: igg(3)
double precision :: bmet(3,3)
double precision :: eg
integer :: ic,id

eg=0.d0
do ic=1,3
do id=1,3
  eg=eg+igg(ic)*bmet(ic,id)*igg(id)
enddo
enddo
!write(6,'(3i3,f10.3)') igg,eg

return
end subroutine fpwe

!**************************************************************************

subroutine expandpwx(igg)
use pwx
implicit none
integer :: igg(3)
integer, allocatable :: itemp(:,:)

npwx=npwx+1
allocate (itemp(3,npwx))
itemp(1:3,1:npwx-1)=ipwx(1:3,1:npwx-1)
itemp(1:3,npwx)=igg(1:3)
deallocate (ipwx)
allocate (ipwx(3,npwx))
ipwx=itemp
deallocate (itemp)

return
end subroutine expandpwx

!**************************************************************************

subroutine fhlf(kpt,nkpt,ihlf)
implicit none
integer :: nkpt
double precision :: kpt(3,nkpt)
integer :: ihlf(nkpt)
integer :: ikpt,ii

do ikpt=1,nkpt
  ihlf(ikpt)=1
  do ii=1,3
    if (mod(2.d0*kpt(ii,ikpt),0.5d0).ne.0.d0) ihlf(ikpt)=0
  enddo
!  write(6,'(i5,3f10.6,5x,i3)') ikpt,kpt(1:3,ikpt),ihlf(ikpt)
enddo

return
end subroutine fhlf

!**************************************************************************

subroutine invsym(nsym,symrel,syminv)
implicit none
integer :: nsym
integer :: symrel(3,3,nsym)
integer :: syminv(3,3,nsym)
integer :: ii,jj,isym
double precision :: xmat(3,3),ymat(3,3)

do isym=1,nsym
  do ii=1,3
    do jj=1,3
      xmat(ii,jj)=dble(symrel(ii,jj,isym))
    enddo
  enddo
  call inverse(xmat,ymat,3,3)
  do ii=1,3
    do jj=1,3
      syminv(ii,jj,isym)=nint(ymat(ii,jj))
    enddo
  enddo
enddo

return
end subroutine invsym

!**************************************************************************

subroutine flvtrans(symrel,shiftk,kpt,ngkpt,nsym,nkpt,nsymk,symk,ikndx,isymndx,lvtrans)
implicit none
integer :: nsym,ngkpt(3),nkpt
integer :: symrel(3,3,nsym)
double precision :: kpt(3,nkpt),shiftk(3)
integer :: nsymk(nkpt),symk(nkpt,nsym*2)
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: lvtrans(3,ngkpt(1),ngkpt(2),ngkpt(3))
integer :: ii,jj,kk,ix,iy,iz,ikk(3),iisym,isign,ikpt
double precision :: xck(3),xck2(3)

do ix=1,ngkpt(1)
do iy=1,ngkpt(2)
do iz=1,ngkpt(3)
  ikk=(/ix,iy,iz/)
  do ii=1,3
    xck(ii)=dble(ikk(ii))/dble(ngkpt(ii))+shiftk(ii)/dble(ngkpt(ii))
    xck(ii)=mod(xck(ii),1.d0)
    if (xck(ii).lt.0.d0) xck(ii)=xck(ii)+1
    xck(ii)=1.d0-mod(1.d0-xck(ii),1.d0)
    xck(ii)=xck(ii)-0.5d0
  enddo
  xck2=(/0.d0,0.d0,0.d0/)
  if (isymndx(ikk(1),ikk(2),ikk(3)).gt.nsym) then
    iisym=isymndx(ikk(1),ikk(2),ikk(3))-nsym
    isign=-1
  else
    iisym=isymndx(ikk(1),ikk(2),ikk(3))
    isign=1
  endif
  ikpt=ikndx(ikk(1),ikk(2),ikk(3))
  do ii=1,3
  do jj=1,3
    xck2(ii)=xck2(ii)+symrel(jj,ii,iisym)*isign*kpt(jj,ikpt)
  enddo
  enddo
  do ii=1,3
    lvtrans(ii,ikk(1),ikk(2),ikk(3))=nint(xck(ii)-xck2(ii))
  enddo
!  write(6,'(3f10.6)') kpt(1:3,ikpt)
!  do ii=1,3
!    write(6,'(3i3)') symrel(ii,1:3,iisym)
!  enddo
!  write(6,'(3f10.6)') xck
!  write(6,'(3f10.6)') xck2
!  write(6,'(3i10)') lvtrans(1:3,ikk(1),ikk(2),ikk(3))
!  stop
enddo
enddo
enddo

return
end subroutine flvtrans

!**************************************************************************

subroutine mkipsymndx(nprefsym,iprefsym,nsym,ipsymndx)
implicit none
integer :: nprefsym,iprefsym(20),nsym
integer :: ipsymndx(nsym)
integer :: ii,jj,isym,jsym,iadd

do isym=1,nprefsym
  ipsymndx(isym)=iprefsym(isym)
enddo
jsym=nprefsym
do ii=1,nsym
  iadd=1
  do isym=1,nprefsym
    if (ii.eq.ipsymndx(isym)) iadd=0
  enddo
  if (iadd.eq.1) then
    jsym=jsym+1
    ipsymndx(jsym)=ii
  endif
enddo

return
end subroutine mkipsymndx

!**************************************************************************

subroutine rdqpts(fname)
use qpoints
implicit none
character(80) :: line
character(80) :: fname
integer :: ipos,iq,ios

open(unit=10,file=fname,status='old')
do
  read(10,'(a)',iostat=ios) line
  if (ios.ne.0) exit
  ipos=index(line,'nkpt')
  if (ipos.ne.0) then
    read(line(ipos+6:ipos+15),*) nqpt
    allocate(qpt(3,nqpt))
    exit
  endif
enddo
l1: do
  read(10,'(a)',iostat=ios) line
  if (ios.ne.0) exit
  if (index(line,'outvars:').ne.0) then
    do
      read(10,'(a)',iostat=ios) line
      if (ios.ne.0) exit
      ipos=index(line,' kpt ')
      if (ipos.ne.0) then
        read(line(ipos+4:),*) qpt(1:3,1)
        do iq=2,nqpt
          read(10,'(a)',iostat=ios) line
          if (index(line,'do not print more k-points')/=0) then
            write(*,*) 'k-point list truncated in abinit output.'//&
            ' Try adjusting prtvol upwards.'
            stop
          end if
          if (ios.ne.0) exit
          read(line,*) qpt(1:3,iq)
        enddo
        exit l1
      endif
    enddo
  endif
enddo l1
close(10)

return
end

!**************************************************************************

subroutine mkpwndx(npwx,npwc,npwlf,ipwndx,npwndx,ntpwndx,napwndx)
implicit none
integer :: npwx,npwc,npwlf
integer :: ipwndx(2,napwndx)
integer :: npwndx,ntpwndx,napwndx,ipw1,ipw2,iipw

  iipw=0
  do ipw1=1,npwlf
  do ipw2=1,ipw1
    iipw=iipw+1
    ipwndx(1,iipw)=ipw1
    ipwndx(2,iipw)=ipw2
!write(6,'(a,2x,2i4,2x,i3,2x,3i5)') 'a ',ipw1,ipw2,iipw,npwndx,ntpwndx,napwndx
  enddo
  enddo

  do ipw1=npwlf+1,npwc
    iipw=iipw+1
    ipwndx(1,iipw)=ipw1
    ipwndx(2,iipw)=ipw1
!write(6,'(a,2x,2i4,2x,i3,2x,3i5)') 'b ',ipw1,ipw1,iipw,npwndx,ntpwndx,napwndx
  enddo

  do ipw1=1,npwlf
  do ipw2=ipw1+1,npwlf
    iipw=iipw+1
    ipwndx(1,iipw)=ipw1
    ipwndx(2,iipw)=ipw2
!write(6,'(a,2x,2i4,2x,i3,2x,3i5)') 'c ',ipw1,ipw2,iipw,npwndx,ntpwndx,napwndx
  enddo
  enddo

  do ipw1=npwc+1,npwx
    iipw=iipw+1
    ipwndx(1,iipw)=ipw1
    ipwndx(2,iipw)=ipw1
!write(6,'(a,2x,2i4,2x,i3,2x,3i5)') 'd ',ipw1,ipw1,iipw,npwndx,ntpwndx,napwndx
  enddo

return
end

!**************************************************************************

subroutine readpaw()
use pawvars
implicit none
integer :: itpaw,ios
character*160 :: title,line,line2
character*4 :: pspfmt
double precision :: zatom,zion,r2well,rcut,rshape
integer :: pspdat,pspcod,pspxc,mmax,bshape
integer :: ii,ix,ibasis,imesh,iset,id,ndm,ndl,idsm,idsp,irad,mm,iorb
integer :: illmax,illoc,inbasis,inlmn,inmesh
integer :: tlorbital,tgridtype,tgridsize
integer :: linlog
double precision :: trstep,tlstep

 do itpaw=1,ntypepaw
   open (unit=10,file=trim(fnamepaw(itpaw)),status='old',iostat=ios)
   if (ios.ne.0) then
     write(6,*) "Error - cannot open file ",trim(fnamepaw(itpaw))
     stop
   endif
   read(10,'(a)') title
   read(10,*) zatom,zion,pspdat
   read(10,*) pspcod,pspxc,illmax,illoc,mmax,r2well
   llmaxmax=max(llmaxmax,illmax)
   llocmax=max(llocmax,illoc)
   read(10,*) pspfmt
   read(10,*) inbasis,inlmn
   nbasismax=max(nbasismax,inbasis)
   nlmnmax=max(nlmnmax,inlmn)
   read(10,*) tlorbital
   read(10,*) inmesh
   nmeshmax=max(nmeshmax,inmesh)
   do ii=1,inmesh
     read(10,*) ix,tgridtype,tgridsize
     gridsizemax=max(gridsizemax,tgridsize)
   enddo
   close(10)
 enddo

 ntlmn=(nlmnmax*(nlmnmax+1))/2
 allocate(nmesh(ntypepaw),nbasis(ntypepaw),nlmn(ntypepaw),llmax(ntypepaw),lloc(ntypepaw))
 allocate(iphigrid(ntypepaw),iprojgrid(ntypepaw),icoregrid(ntypepaw),ivgrid(ntypepaw))
 allocate(lorbital(ntypepaw,nbasismax))
 allocate(gridtype(ntypepaw,nmeshmax),gridsize(ntypepaw,nmeshmax), &
&         rstep(ntypepaw,nmeshmax),lstep(ntypepaw,nmeshmax))
 allocate(rphigrid(ntypepaw,gridsizemax))
 allocate(phi(ntypepaw,nbasismax,gridsizemax))
 allocate(tphi(ntypepaw,nbasismax,gridsizemax))
 allocate(rprojgrid(ntypepaw,gridsizemax))
 allocate(tproj(ntypepaw,nbasismax,gridsizemax))
 allocate(rcoregrid(ntypepaw,gridsizemax))
 allocate(coredens(ntypepaw,gridsizemax))
 allocate(tcoredens(ntypepaw,gridsizemax))
 allocate(dij0(ntypepaw,ntlmn),rhoij0(ntypepaw,ntlmn))
 allocate(rvgrid(ntypepaw,gridsizemax))
 allocate(vhntzc(ntypepaw,gridsizemax))
 allocate(iorbno(itpaw,nlmnmax),llorb(itpaw,nlmnmax),mmorb(itpaw,nlmnmax))

 do itpaw=1,ntypepaw
   open (unit=10,file=trim(fnamepaw(itpaw)),status='old')
   read(10,'(a)') title
   read(10,*) zatom,zion,pspdat
   read(10,*) pspcod,pspxc,llmax(itpaw),lloc(itpaw),mmax,r2well
   read(10,*) pspfmt
   read(10,*) nbasis(itpaw),nlmn(itpaw)
   read(10,*) lorbital(itpaw,1:nbasis(itpaw))
   read(10,*) nmesh(itpaw)
   do ii=1,nmesh(itpaw)
     read(10,'(a)') line2
     read(line2,*,iostat=ios) ix,gridtype(itpaw,ii),gridsize(itpaw,ii),rstep(itpaw,ii),lstep(itpaw,ii)
     if (ios.ne.0) then
       read(line2(1:),*,iostat=ios) ix,gridtype(itpaw,ii),gridsize(itpaw,ii),rstep(itpaw,ii)
       lstep(itpaw,ii) = -1.d0
       linlog=1
     else
       linlog=0
     endif
   enddo
   read(10,*) rcut
   read(10,*) bshape,rshape

   read(10,*)
   read(10,*) imesh
   iphigrid(itpaw)=imesh
   do irad=1,gridsize(itpaw,iphigrid(itpaw))
     if (linlog.eq.0) then
       rphigrid(itpaw,irad)=rstep(itpaw,iphigrid(itpaw)) &
&              *(exp(dble(irad-1)*lstep(itpaw,iphigrid(itpaw)))-1.d0)
     else
       rphigrid(itpaw,irad)=rstep(itpaw,iphigrid(itpaw))*dble(irad-1)
     endif
   enddo
   do ibasis=1,nbasis(itpaw)
     ndl=gridsize(itpaw,iphigrid(itpaw))/3
     ndm=mod(gridsize(itpaw,iphigrid(itpaw)),3)
     do id=1,ndl
       read(10,*) phi(itpaw,ibasis,3*id-2:3*id)
     enddo
     if (ndm.ne.0) read(10,*) phi(itpaw,ibasis,3*ndl+1:3*ndl+ndm)
     if (ibasis.lt.nbasis(itpaw)) read(10,*)
     if (ibasis.lt.nbasis(itpaw)) read(10,*) imesh
   enddo

   do ibasis=1,nbasis(itpaw)
     read(10,*)
     read(10,*) imesh
     ndl=gridsize(itpaw,iphigrid(itpaw))/3
     ndm=mod(gridsize(itpaw,iphigrid(itpaw)),3)
     do id=1,ndl
       read(10,*) tphi(itpaw,ibasis,3*id-2:3*id)
     enddo
     if (ndm.ne.0) read(10,*) tphi(itpaw,ibasis,3*ndl+1:3*ndl+ndm)
   enddo

   read(10,*)
   read(10,*) imesh
   iprojgrid(itpaw)=imesh
   do irad=1,gridsize(itpaw,iprojgrid(itpaw))
     rprojgrid(itpaw,irad)=rstep(itpaw,iprojgrid(itpaw)) &
&            *(exp(dble(irad-1)*lstep(itpaw,iprojgrid(itpaw)))-1.d0)
   enddo
   do ibasis=1,nbasis(itpaw)
     ndl=gridsize(itpaw,iprojgrid(itpaw))/3
     ndm=mod(gridsize(itpaw,iprojgrid(itpaw)),3)
     do id=1,ndl
       read(10,*) tproj(itpaw,ibasis,3*id-2:3*id)
     enddo
     if (ndm.ne.0) read(10,*) tproj(itpaw,ibasis,3*ndl+1:3*ndl+ndm)
     if (ibasis.lt.nbasis(itpaw)) read(10,*)
     if (ibasis.lt.nbasis(itpaw)) read(10,*) imesh
   enddo

   read(10,*)
   read(10,*) imesh
   icoregrid(itpaw)=imesh
   do irad=1,gridsize(itpaw,icoregrid(itpaw))
     rcoregrid(itpaw,irad)=rstep(itpaw,icoregrid(itpaw)) &
&            *(exp(dble(irad-1)*lstep(itpaw,icoregrid(itpaw)))-1.d0)
   enddo
   ndl=gridsize(itpaw,icoregrid(itpaw))/3
   ndm=mod(gridsize(itpaw,icoregrid(itpaw)),3)
   do id=1,ndl
     read(10,*) coredens(itpaw,3*id-2:3*id)
   enddo
   if (ndm.ne.0) read(10,*) coredens(itpaw,3*ndl+1:3*ndl+ndm)

   read(10,*)
   read(10,*) imesh
   ndl=gridsize(itpaw,icoregrid(itpaw))/3
   ndm=mod(gridsize(itpaw,icoregrid(itpaw)),3)
   do id=1,ndl
     read(10,*) tcoredens(itpaw,3*id-2:3*id)
   enddo
   if (ndm.ne.0) read(10,*) tcoredens(itpaw,3*ndl+1:3*ndl+ndm)

   read(10,*)
   do id=1,nlmn(itpaw)
     idsp=(id*(id+1))/2
     idsm=(id*(id-1))/2+1
     read(10,*) dij0(itpaw,idsm:idsp)
   enddo

   read(10,*)
   do id=1,nlmn(itpaw)
     idsp=(id*(id+1))/2
     idsm=(id*(id-1))/2+1
     read(10,*) rhoij0(itpaw,idsm:idsp)
   enddo

   read(10,*)
   read(10,*) imesh
   ivgrid(itpaw)=imesh
   do irad=1,gridsize(itpaw,ivgrid(itpaw))
     rvgrid(itpaw,irad)=rstep(itpaw,ivgrid(itpaw)) & 
&         *(exp(dble(irad-1)*lstep(itpaw,ivgrid(itpaw)))-1.d0)
   enddo
   ndl=gridsize(itpaw,ivgrid(itpaw))/3
   ndm=mod(gridsize(itpaw,ivgrid(itpaw)),3)
   do id=1,ndl
     read(10,*) vhntzc(itpaw,3*id-2:3*id)
   enddo
   if (ndm.ne.0) read(10,*) vhntzc(itpaw,3*ndl+1:3*ndl+ndm)

   close(10)

   iorb=1
   do ibasis=1,nbasis(itpaw)
     do mm=-lorbital(itpaw,ibasis),lorbital(itpaw,ibasis)
       iorbno(itpaw,iorb)=ibasis
       llorb(itpaw,iorb)=lorbital(itpaw,ibasis)
       mmorb(itpaw,iorb)=mm
!write(6,'(i3,2x,3i3)') iorb,iorbno(iorb),llorb(iorb),mmorb(iorb)
       iorb=iorb+1
     enddo
   enddo

 enddo

return
end subroutine readpaw

!**************************************************************************

subroutine prtsym(nsym,symrel,syminv,tnons,ngkpt,isymndx,ikndx)
implicit none
integer :: nsym,ngkpt(3)
integer :: symrel(3,3,nsym),syminv(3,3,nsym),tnons(3,nsym)
integer :: isymndx(ngkpt(1),ngkpt(2),ngkpt(3)),ikndx(ngkpt(1),ngkpt(2),ngkpt(3))
integer :: nln,ilsym,ioffmx,jloop,ii,jj,kk,ioff

 open(unit=10,file='symmetry',status='unknown')
 nln=6
 write(10,*) "Symmetry operations"
 do ilsym=1,(nsym+nln-1)/nln
   ioffmx=-max(0,nln*ilsym-nsym)
   write(10,'(7(i2,10x))') (nln*ilsym+ioff,ioff=-(nln-1),ioffmx)
   do jloop=1,3
      write(10,'(7(3i3,3x))') (symrel(1:3,jloop,nln*ilsym+ioff),ioff=-(nln-1),ioffmx)
   enddo
   write(10,*)
 enddo
 write(10,*) "Inverse symmetry operations"
 do ilsym=1,(nsym+nln-1)/nln
   ioffmx=-max(0,nln*ilsym-nsym)
   write(10,'(7(i2,10x))') (nln*ilsym+ioff,ioff=-(nln-1),ioffmx)
   do jloop=1,3
      write(10,'(7(3i3,3x))') (syminv(1:3,jloop,nln*ilsym+ioff),ioff=-(nln-1),ioffmx)
   enddo
   write(10,*)
 enddo
 write(10,*) "Translation non-symorphic vectors"
 do ilsym=1,(nsym+nln-1)/nln
   ioffmx=-max(0,nln*ilsym-nsym)
   write(10,'(7(i2,10x))') (nln*ilsym+ioff,ioff=-(nln-1),ioffmx)
   write(10,'(7(3i3,3x))') (tnons(1:3,nln*ilsym+ioff),ioff=-(nln-1),ioffmx)
   write(10,*)
 enddo
 write(10,*) "Symmetry indices on k-point grid"
 do ii=1,ngkpt(1)
 write(10,*) ii
 do jj=1,ngkpt(2)
   write(10,'(10i3,5x,10i3)') (mod(isymndx(ii,jj,kk)-1,nsym)+1,kk=1,ngkpt(3)), &
&   (-(2*((isymndx(ii,jj,kk)-1)/nsym)-1),kk=1,ngkpt(3))
!   write(10,'(16i3)') isymndx(ii,jj,1:ngkpt(3))
!   write(10,'(16i5)') ikndx(ii,jj,1:ngkpt(3))
 enddo
 write(10,*)
 enddo
 write(10,*) nsym
 close(10)

return
end subroutine prtsym

!**************************************************************************

subroutine chkwfsym(nsym,syminv,symrel,natom,ihlf,nkpt,igndx,igmx,igmn)
implicit none
integer :: nsym,natom,igmn(3),igmx(3),nkpt
integer :: syminv(3,3,nsym),symrel(3,3,nsym),ihlf(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: iband,ikpt,iatom,isym,ii,jj
double precision :: rr(3),rrs(3),orbden
double complex :: wfval,twfval,pwfval,pwfval1,ptwf,ppwf,ptwf1,ppwf1

 open(unit=10,file='chkwfsym',status='unknown')
 rr=(/0.12d0,0.12d0,0.12d0/)

 iband=16
 ikpt=5

 call pwrwf(rr,iband,ikpt,twfval,ihlf,igndx,igmx,igmn)
 pwfval=0.d0
 ptwf=0.d0
 ppwf=0.d0
 do iatom=1,natom
   call prjrwf(rr,iatom,iband,ikpt,pwfval1)
   call prjrtwf(rr,iatom,iband,ikpt,ptwf1)
   call prjrpwf(rr,iatom,iband,ikpt,ppwf1)
   pwfval=pwfval+pwfval1
   ptwf=ptwf+ptwf1
   ppwf=ppwf+ppwf1
 enddo
 wfval=twfval+pwfval
 write(10,'(a,2es15.6)') 'plane wave expansion          ',twfval
 write(10,'(a,2es15.6)') 'auxiliary orbital expansion   ',ptwf
 write(10,'(a,2es15.6)') 'atomic orbital expansion      ',ppwf
 write(10,'(a,2es15.6)') 'PAW correction                ',wfval
! wfval=twfval
 orbden=dble(wfval*conjg(wfval))
 write(10,'(i4,es15.6)') 0,orbden
 write(10,*)

return
end subroutine chkwfsym

!**************************************************************************

subroutine prtdielf(iprtdielf,npwc,npwlf,ipwndx,napwndx,nqpt,qpt,qtensor,ngkpt,nwpt,npwndx,pola,polb,dw)
use geometry
use pwx
implicit none
integer :: iprtdielf(5),npwc,npwlf,nqpt,ngkpt(3),nwpt,npwndx,iks(3)
integer :: iw,iqpt,qtensor,ii,jj,kk,ll,mm,nn,iiq,iipw,ipw,ipw1,ipw2,napwndx
integer :: ipwndx(2,napwndx)
double precision :: qpt(3,nqpt),dw,ww,hart
double complex :: pola(nwpt,nqpt+9,npwndx),polb(nwpt,nqpt+9,npwndx),dielf(npwndx)
character(80) :: fileout1,fileout2,dirout
character(2) :: chq1,chq2,chq3
character(3) :: chpw1,chpw2
double precision :: qh1(3),qh2(3),qh1m,qh2m
double complex :: t_comp,t_d(3,3)

 hart=27.2113845d0            ! atomic unit of energy = 27.2 eV
 !call MPE_DECOMP1D(iprtdielf(2)-iprtdielf(1)+1,numprocs,my_rank,iat_start,iat_end)
 !iat_start=iat_start+iprtdielf(1)-1
 !iat_end=iat_end+iprtdielf(1)-1
 !write(*,*) 'proc ',my_rank,' will print dielf for plane waves from ',iat_start&
 !,' to ',iat_end,'.'
 !write(*,*) 'Global plane wave bounds: ',ipolpw(1),ipolpw(2)
 !do iipw=iat_start,iat_end
 do iipw=iprtdielf(2),iprtdielf(3)
   ipw1=ipwndx(1,iipw)
   ipw2=ipwndx(2,iipw)
write(28,*) iipw
   write(chpw1,'(i3)') ipw1
   if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
   if (chpw1(2:2).eq.' ') chpw1(2:2)='0'
   write(chpw2,'(i3)') ipw2
   if (chpw2(1:1).eq.' ') chpw2(1:1)='0'
   if (chpw2(2:2).eq.' ') chpw2(2:2)='0'
   do iqpt=iprtdielf(4),iprtdielf(5)
     do ii=1,3
       iks(ii)=nint(qpt(ii,iqpt)*ngkpt(ii))
     enddo
     if (qtensor.eq.iqpt) then
do iw=1,nwpt
write(28,'(i5,f9.2,9e12.4)') iw,(dimag(polb(iw,nqpt+iiq,iipw)),iiq=1,9)
enddo
       do ii=1,3
       do jj=1,3
!        q in cartesian direction ii, q' in cartesian direction jj
         write(chq1,'(i2)') ii
         if (chq1(1:1).eq.' ') chq1(1:1)='0'
         write(chq2,'(i2)') jj
         if (chq2(1:1).eq.' ') chq2(1:1)='0'
         dirout='Data/'
         fileout1='dielf-t_'//chq1//','//chq2//'-'//chpw1//','//chpw2//'.dat'
         fileout2='alpha-r_'//chq1//','//chq2//'-'//chpw1//','//chpw2//'.dat'
         open(unit=20,file=trim(dirout)//trim(fileout1),status="unknown")
         open(unit=21,file=trim(dirout)//trim(fileout2),status="unknown")
         do iw=1,nwpt
           ww=dble(iw)*dw
!          find D_{k,l} for 1 + D = dielctric tensor in basis of reciprocal lattice vectors
           do kk=1,3
           do ll=1,3
             iiq = kk + 3*(ll-1)
             t_d(kk,ll)=-pola(iw,nqpt+iiq,iipw)-polb(iw,nqpt+iiq,iipw)
           enddo
           enddo
!          find T = q_k * D_{k,l} * q'_l for component of (dielectric tensor - 1) in specified cartesian coordinates
           t_comp = (0.d0,0.d0)
           do kk=1,3
           do ll=1,3
             t_comp = t_comp + blat(kk,ii)*t_d(kk,ll)*blat(ll,jj)
           enddo
           enddo
           if (ipw1.eq.ipw2.and.ii.eq.jj) then
             write(20,'(f6.2,3x,2f20.6,3x,2f20.6)') &
&                  ww*hart,1.d0+t_comp
           else
             write(20,'(f6.2,3x,2f20.6)') &
&                  ww*hart,t_comp
           endif
           write(21,'(f6.2,3x,2f20.6)') &
&                ww*hart,t_d(ii,jj)
         enddo
         close(21)
         close(20)
       enddo
       enddo
     else 
       write(chq1,'(i2)') iks(1)
       if (chq1(1:1).eq.' ') chq1(1:1)='0'
       write(chq2,'(i2)') iks(2)
       if (chq2(1:1).eq.' ') chq2(1:1)='0'
       write(chq3,'(i2)') iks(3)
       if (chq3(1:1).eq.' ') chq3(1:1)='0'
       dirout='Data/'
       fileout1='dielf'//chq1//','//chq2//','//chq3//'-'//chpw1//','//chpw2//'.dat'
       open(unit=20,file=trim(dirout)//trim(fileout1),status="unknown")
       do iw=1,nwpt
         ww=dble(iw)*dw
         if (ipw1.eq.ipw2) then
           write(20,'(f6.2,3x,2f20.6,3x,2f20.6)') &
&                  ww*hart,1.d0-pola(iw,iqpt,iipw)-polb(iw,iqpt,iipw)
         else
           write(20,'(f6.2,3x,2f20.6)') &
&                  ww*hart,-pola(iw,iqpt,iipw)-polb(iw,iqpt,iipw)
         endif
       enddo
       close(20)
     endif
   enddo
 enddo

return
end subroutine prtdielf

!**************************************************************************

subroutine prtlosfn(iprtlosfn,npwc,npwlf,ipwndx,napwndx,nqpt,qpt,qtensor,ngkpt,nwpt,ntpwndx,lossfn,dw)
use geometry
implicit none
integer :: iprtlosfn(5),npwc,npwlf,nqpt,ngkpt(3),nwpt,ntpwndx,iks(3)
integer :: iw,iqpt,qtensor,ii,jj,kk,ll,iiq,iipw,ipw1,ipw2,napwndx
integer :: ipwndx(2,napwndx)
double precision :: qpt(3,nqpt),dw,ww,hart
double complex :: lossfn(nwpt,nqpt+9,ntpwndx)
character(80) :: fileout1,fileout2,dirout
character(2) :: chq1,chq2,chq3
character(3) :: chpw1,chpw2
double complex :: t_comp,t_d(3,3)

 hart=27.2113845d0            ! atomic unit of energy = 27.2 eV
 do iipw=iprtlosfn(2),iprtlosfn(3)
   ipw1=ipwndx(1,iipw)
   ipw2=ipwndx(2,iipw)
   write(chpw1,'(i3)') ipw1
   if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
   if (chpw1(2:2).eq.' ') chpw1(2:2)='0'
   write(chpw2,'(i3)') ipw2
   if (chpw2(1:1).eq.' ') chpw2(1:1)='0'
   if (chpw2(2:2).eq.' ') chpw2(2:2)='0'
   do iqpt=iprtlosfn(4),iprtlosfn(5)
     do ii=1,3
       iks(ii)=nint(qpt(ii,iqpt)*ngkpt(ii))
     enddo
     if (qtensor.eq.iqpt) then
       do ii=1,3
       do jj=1,3
         iiq = ii + 3*(jj-1)
         write(chq1,'(i2)') ii
         if (chq1(1:1).eq.' ') chq1(1:1)='0'
         write(chq2,'(i2)') jj
         if (chq2(1:1).eq.' ') chq2(1:1)='0'
         dirout='Data/'
         fileout2='losfn-t_'//chq1//','//chq2//'-'//chpw1//','//chpw2//'.dat'
         open(unit=20,file=trim(dirout)//trim(fileout2),status="unknown")
         do iw=1,nwpt
           ww=dble(iw)*dw
!          for lossfn = loss function tensor in basis of reciprocal lattice vectors
           do kk=1,3
           do ll=1,3
             iiq = kk + 3*(ll-1)
             t_d(kk,ll)=lossfn(iw,nqpt+iiq,iipw)
           enddo
           enddo
!          find T = q_k * lossfn_{k,l} * q'_l for component of loss function tensor in specified cartesian coordinates
           t_comp = (0.d0,0.d0)
           do kk=1,3
           do ll=1,3
             t_comp = t_comp + blat(kk,ii)*t_d(kk,ll)*blat(ll,jj)
           enddo
           enddo
           write(20,'(f6.2,3x,2f20.6)') ww*hart,t_comp
         enddo
         close(20)
       enddo
       enddo
     else 
       write(chq1,'(i2)') iks(1)
       if (chq1(1:1).eq.' ') chq1(1:1)='0'
       write(chq2,'(i2)') iks(2)
       if (chq2(1:1).eq.' ') chq2(1:1)='0'
       write(chq3,'(i2)') iks(3)
       if (chq3(1:1).eq.' ') chq3(1:1)='0'
       dirout='Data/'
       fileout2='losfn'//chq1//','//chq2//','//chq3//'-'//chpw1//','//chpw2//'.dat'
       open(unit=20,file=trim(dirout)//trim(fileout2),status="unknown")
       do iw=1,nwpt
         ww=dble(iw)*dw
         write(20,'(f6.2,3x,2f20.6)') ww*hart,lossfn(iw,iqpt,iipw)
       enddo
       close(20)
     endif
   enddo
 enddo

return
end subroutine prtlosfn

!**************************************************************************

subroutine prtxtdielf(iprtxtdielf,ipwndx,napwndx,nxtqpt,xtqpt,ngkpt,nwpt,npwndx,xtpola,xtpolb,dw)
use geometry
use pwx
implicit none
integer :: iprtxtdielf(5),nxtqpt,ngkpt(3),nwpt,npwndx,iks(3)
integer :: iw,iqpt,ii,jj,kk,ll,mm,nn,iiq,iipw,ipw,ipw1,ipw2,napwndx
integer :: ipwndx(2,napwndx)
double precision :: xtqpt(3,nxtqpt),dw,ww,hart
double complex :: xtpola(nwpt,nxtqpt,npwndx),xtpolb(nwpt,nxtqpt,npwndx),dielf(npwndx)
character(80) :: fileout1,fileout2,dirout
character(2) :: chq1,chq2,chq3
character(3) :: chpw1,chpw2
double precision :: qh1(3),qh2(3),qh1m,qh2m
double complex :: t_comp,t_d(3,3)

 hart=27.2113845d0            ! atomic unit of energy = 27.2 eV
 do ipw1=iprtxtdielf(2),iprtxtdielf(3)
   write(chpw1,'(i3)') ipw1
   if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
   if (chpw1(2:2).eq.' ') chpw1(2:2)='0'
   do iqpt=iprtxtdielf(4),iprtxtdielf(5)
     do ii=1,3
       iks(ii)=nint(xtqpt(ii,iqpt)*ngkpt(ii))
     enddo
     write(chq1,'(i2)') iks(1)
     if (chq1(1:1).eq.' ') chq1(1:1)='0'
     write(chq2,'(i2)') iks(2)
     if (chq2(1:1).eq.' ') chq2(1:1)='0'
     write(chq3,'(i2)') iks(3)
     if (chq3(1:1).eq.' ') chq3(1:1)='0'
     dirout='Data/'
     fileout1='xt_dielf'//chq1//','//chq2//','//chq3//'-'//chpw1//'.dat'
     open(unit=20,file=trim(dirout)//trim(fileout1),status="unknown")
     do iw=1,nwpt
       ww=dble(iw)*dw
       write(20,'(f6.2,3x,2f20.6,3x,2f20.6)') &
&                ww*hart,1.d0-xtpola(iw,iqpt,ipw1)-xtpolb(iw,iqpt,ipw1)
     enddo
     close(20)
   enddo
 enddo

return
end subroutine prtxtdielf

!**************************************************************************

subroutine prtxtlosfn(iprtxtlosfn,ipwndx,napwndx,nxtqpt,xtqpt,ngkpt,nwpt,npwndx,xtlossfn,dw)
use geometry
use pwx
implicit none
integer :: iprtxtlosfn(5),nxtqpt,ngkpt(3),nwpt,npwndx,iks(3)
integer :: iw,iqpt,ii,jj,kk,ll,mm,nn,iiq,iipw,ipw,ipw1,ipw2,napwndx
integer :: ipwndx(2,napwndx)
double precision :: xtqpt(3,nxtqpt),dw,ww,hart
double complex :: xtlossfn(nwpt,nxtqpt,npwndx)
character(80) :: fileout1,fileout2,dirout
character(2) :: chq1,chq2,chq3
character(3) :: chpw1,chpw2
double precision :: qh1(3),qh2(3),qh1m,qh2m
double complex :: t_comp,t_d(3,3)

 hart=27.2113845d0            ! atomic unit of energy = 27.2 eV
 do ipw1=iprtxtlosfn(2),iprtxtlosfn(3)
   write(chpw1,'(i3)') ipw1
   if (chpw1(1:1).eq.' ') chpw1(1:1)='0'
   if (chpw1(2:2).eq.' ') chpw1(2:2)='0'
   do iqpt=iprtxtlosfn(4),iprtxtlosfn(5)
     do ii=1,3
       iks(ii)=nint(xtqpt(ii,iqpt)*ngkpt(ii))
     enddo
     write(chq1,'(i2)') iks(1)
     if (chq1(1:1).eq.' ') chq1(1:1)='0'
     write(chq2,'(i2)') iks(2)
     if (chq2(1:1).eq.' ') chq2(1:1)='0'
     write(chq3,'(i2)') iks(3)
     if (chq3(1:1).eq.' ') chq3(1:1)='0'
     dirout='Data/'
     fileout1='xt_losfn'//chq1//','//chq2//','//chq3//'-'//chpw1//'.dat'
     open(unit=20,file=trim(dirout)//trim(fileout1),status="unknown")
     do iw=1,nwpt
       ww=dble(iw)*dw
       write(20,'(f6.2,3x,2f20.6)') ww*hart,xtlossfn(iw,iqpt,ipw1)
     enddo
     close(20)
   enddo
 enddo

return
end subroutine prtxtlosfn

!**************************************************************************

subroutine fovlp(ikpt,iband,jkpt,jband,ipaw,nkpt,ncg,ntypepaw,typat,nlmn, &
& natom,ncband,cg,indxkcg,npwarr,projwf,aoovlp,taoovlp,pwovlp,pawovlp,ovlp)
implicit none
integer :: ikpt,jkpt,iband,jband,ipaw,nkpt,ncg,nlmn,natom,ncband,itpaw,ntypepaw
integer :: icg,jcg,ipw,iatom,ilmn,jlmn
integer :: indxkcg(nkpt),npwarr(nkpt)
integer :: typat(natom)
double complex :: cg(ncg),projwf(natom,nlmn,nkpt,ncband)
double complex :: aoovlp(ntypepaw,nlmn,nlmn),taoovlp(ntypepaw,nlmn,nlmn)
double complex :: pwovlp,pawovlp,ovlp

  icg=indxkcg(ikpt)+(iband-1)*npwarr(ikpt)
  jcg=indxkcg(jkpt)+(jband-1)*npwarr(jkpt)
  pwovlp=(0.d0,0.d0)
  do ipw=1,npwarr(ikpt)
    pwovlp=pwovlp+cg(icg+ipw)*conjg(cg(jcg+ipw))
  enddo
  pawovlp=(0.d0,0.d0)
  if (ipaw.ne.0) then
    do iatom=1,natom
      itpaw=typat(iatom)
      do ilmn=1,nlmn
      do jlmn=1,nlmn
        pawovlp=pawovlp+projwf(iatom,ilmn,ikpt,iband) &
&               *conjg(projwf(iatom,jlmn,jkpt,jband)) &
&               *(aoovlp(itpaw,ilmn,jlmn)-taoovlp(itpaw,ilmn,jlmn))
      enddo
      enddo
    enddo
  endif
  ovlp=pwovlp+pawovlp

return
end subroutine fovlp

!**************************************************************************

subroutine fsymovlp(ikpt,iband,jband,isym,nkpt,ncg,npwt,ncband,ngkpt, &
& nsym,igmn,igmx,igndx,ikndx,symrel,cg,kg,indxkcg,indxkpw,npwarr,kpt,ovlp)
implicit none
integer :: ikpt,ikpts,iband,jband,isym,nkpt,ncg,ncband,npwt,nsym
integer :: igmn(3),igmx(3),ngkpt(3),symrel(3,3,nsym)
integer :: icg,jcg,ipw,jpw,igg(3),iggs(3),ikk(3),ikks(3),ii,jj
integer :: ikndx(ngkpt(1),ngkpt(2),ngkpt(3)),indxkpw(nkpt)
integer :: igndx(nkpt,igmn(1):igmx(1),igmn(2):igmx(2),igmn(3):igmx(3))
integer :: indxkcg(nkpt),npwarr(nkpt),kg(3,npwt)
double complex :: cg(ncg),ovlp
double precision :: kpt(3,nkpt)

  ikk=nint(kpt(:,ikpt)*ngkpt)
  ikks=(/0,0,0/)
  do ii=1,3
  do jj=1,3
    ikks(jj)=ikks(jj)+ikk(ii)*symrel(ii,jj,isym)
  enddo
  enddo
  ikks=ikks+ngkpt/2
  ikpts=ikndx(ikks(1),ikks(2),ikks(3))
  icg=indxkcg(ikpts)+(iband-1)*npwarr(ikpts)
  jcg=indxkcg(ikpt)+(jband-1)*npwarr(ikpt)
  ovlp=(0.d0,0.d0)
  do ipw=1,npwarr(ikpt)
    igg=kg(:,indxkpw(ikpt)+ipw)
    iggs=(/0,0,0/)
    do ii=1,3
    do jj=1,3
      iggs(jj)=iggs(jj)+igg(ii)*symrel(ii,jj,isym)
    enddo
    enddo
    jpw=igndx(ikpt,iggs(1),iggs(2),iggs(3))
    ovlp=ovlp+cg(icg+jpw)*conjg(cg(jcg+ipw))
  enddo

return
end subroutine fsymovlp






























