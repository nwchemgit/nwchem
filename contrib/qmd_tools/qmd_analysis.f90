!
! qmd_analysis
!
! Written by Sean A. Fischer
!
! Analysis code written specifically for the
!   qmd module of NWChem.
!
! The vibrational density of states and the IR 
!   spectrum are computed.
!
! The vibrational density of states is broken down
!   in contributions from each element
!
! The IR spectrum is output as-is and with the
!   harmonic approximation quantum correction
!   see J. Chem. Phys. 121, 3973 (2004)
!
! Multiple trajectories can be concatenated together
!   in order to compute the average spectrum.
!
! User has the option of aligning all the molecules
!   to the first frame of the trajectory, thereby
!   eliminating contributions from rotations.
!
! Linking to LAPACK REQUIRED.
!
program qmd_analysis
implicit none

integer :: i, j, is, ig
integer :: pn
!Number of atoms
integer :: nat
!Number of data points
integer :: nsteps
!Offset: Data points in input file >= nsteps+skip
integer :: skip
integer :: cstep, cstep_prev, cla, ioerror
integer,dimension(:,:),allocatable :: atgrp
logical :: step_flag, ts_flag, xyz_flag
logical :: ignoreH, align_flag, peratom
logical,dimension(:),allocatable :: mask
character(2) :: elem
character(2),dimension(:),allocatable :: grptags,altags
character(255) :: arg
character*(255) :: xyzname
character*(255) :: outnameIR, outnameVDOS, outnameAL
double precision,parameter :: c=2.99792458d-5 !cm/fs
double precision,parameter :: au2fs=2.418884326505d-2
double precision,parameter :: Pi=4.d0*atan(1.d0)
double precision,parameter :: kb=3.1668114d-6
double precision,parameter :: hbar=1.d0
double precision :: ts, specmax, mean
double precision :: width, damp, temper, qc
double precision :: freq_au, freq_cm, energy
double precision,dimension(3) :: scrdip, com
double precision,dimension(:,:),allocatable :: coord, ref
double precision,dimension(:),allocatable :: atmass
double complex,dimension(:,:),allocatable :: vdos_fft
double complex,dimension(:),allocatable :: ir_fft
double complex,dimension(:),allocatable :: v
double complex,dimension(:,:),allocatable :: vdos_acf
double complex,dimension(:),allocatable :: ir_acf
double precision,dimension(:,:),allocatable :: dip
double precision,dimension(:,:,:),allocatable :: vel
double precision,dimension(:),allocatable :: wsave
character*6 vdos_tag
integer jfirst,jlast

width=10.d0
skip=0
specmax=5.d3
temper=298.15d0
ts_flag=.false.
xyz_flag=.false.
step_flag=.false.
align_flag=.false.
ignoreH=.false.
peratom=.false.
jfirst=1
jlast=3
i=1
cla=command_argument_count()
if (cla==0) then
 call usage_mess
 stop
end if
do
 if (i>cla) exit
 call get_command_argument(i,arg)
 if (trim(arg)=='-help') then
  call usage_mess
  stop
 else if (trim(arg)=='-xyz') then
  call get_command_argument(i+1,arg,j,ioerror)
  if (ioerror.ne.0) then
   write(*,*) 'IO error: Problem with input file name'
   stop
  end if
  xyzname=trim(arg)
  if (xyzname(j-3:j).ne.'.xyz') xyzname=trim(arg)//'.xyz'
  outnameIR=xyzname(1:j-4)//'_IR.dat'
  outnameVDOS=xyzname(1:j-4)//'_VDOS.dat'
  xyz_flag=.true.
  i=i+2
 else if (trim(arg)=='-width') then
  call get_command_argument(i+1,arg)
  read(arg,*) width
  i=i+2
 else if (trim(arg)=='-steps') then
  call get_command_argument(i+1,arg)
  read(arg,*) nsteps
  step_flag=.true.
  i=i+2
 else if (trim(arg)=='-skip') then
  call get_command_argument(i+1,arg)
  read(arg,*) skip
  i=i+2
 else if (trim(arg)=='-ts') then
  call get_command_argument(i+1,arg)
  read(arg,*) ts
  ts_flag=.true.
  i=i+2
 else if (trim(arg)=='-temp') then
  call get_command_argument(i+1,arg)
  read(arg,*) temper
  i=i+2
 else if (trim(arg)=='-smax') then
  call get_command_argument(i+1,arg)
  read(arg,*) specmax
  i=i+2
 else if (trim(arg)=='-align') then
  align_flag=.true.
  outnameAL=xyzname(1:j-4)//'_aligned.xyz'
  i=i+1
 else if (trim(arg)=='-ignH') then
  ignoreH=.true.
  i=i+1
 else if (trim(arg)=='-peratom') then
  peratom=.true.
  i=i+1
 else if (trim(arg)=='-peratomx') then
    peratom=.true.
    jfirst=1
    jlast=1
  i=i+1
 else if (trim(arg)=='-peratomy') then
    peratom=.true.
    jfirst=2
    jlast=2
  i=i+1
 else if (trim(arg)=='-peratomz') then
    peratom=.true.
    jfirst=3
    jlast=3
  i=i+1
 else
  write(*,*)
  write(*,*) 'Unrecognized option:', TRIM(arg)
  call usage_mess
  stop
 end if
end do
if (.not.xyz_flag) then
 write(*,*)
 write(*,*) 'No input specified'
 call usage_mess
 stop
end if
if (.not.ts_flag) then
 write(*,*)
 write(*,*) 'No time step specified'
 call usage_mess
 stop
end if
if (.not.step_flag) then
 write(*,*)
 write(*,*) 'Number of steps not specified'
 call usage_mess
 stop
end if
open(unit=9,file=xyzname,position='append')
open(unit=10,file=outnameIR)
open(unit=20,file=outnameVDOS)
if (align_flag) open(unit=30,file=outnameAL)
write(*,*)
write(*,'(A,A)') ' Input file:                    ', trim(xyzname)
write(*,'(A,A)') ' Ouput files:                   ', trim(outnameVDOS)
write(*,'(A,A)') '                                ', trim(outnameIR)
if (align_flag) then
 write(*,'(A,A)') '                                ', trim(outnameAL)
end if
write(*,'(A,I8)') ' Number of steps to analyze:    ', nsteps
write(*,'(A,I8)') ' Number of steps to be skipped: ', skip
write(*,'(A,F8.2,A)') ' Temperature:                   ', temper, ' K'
write(*,'(A,F8.2,A)') ' Time step of simulation:       ', ts, ' a.u.'
write(*,'(A,F8.2,A)') ' FWHM of resulting spectrum:    ', width, ' cm^{-1}'
write(*,*)
damp=2.d0/(width*c*au2fs)
!Set number of points for FFT
!Need at least 2*nsteps so ACF will be done correctly
!Routine is most efficient with a power of 2
pn=2**ceiling(log(dble(2*nsteps))/log(2.d0))
allocate(v(pn))
allocate(ir_acf(pn))
allocate(wsave(4*pn+15))
call zffti(pn,wsave)
allocate(dip(3,nsteps))
com=0.d0
ig=1
!Make sure we are at beginning of xyz file
rewind(9)
!Get number of atoms
read(9,*,iostat=ioerror) nat
if (ioerror.ne.0) then
 write(*,*) 'Error reading number of atoms'
 stop
end if
allocate(mask(nat))
allocate(coord(3,nat))
allocate(ref(3,nat))
allocate(grptags(nat))
allocate(altags(nat))
allocate(atgrp(nat,2))
allocate(vel(3,nat,nsteps))
allocate(atmass(nat))
mask=.true.
atgrp=0
!Find out if comment line is as expected
read(9,*,iostat=ioerror) cstep, energy, scrdip
if (ioerror.ne.0) then
 write(*,*) 'Error reading comment line'
 write(*,*) 'Expecting step number, energy, dipole'
 stop
end if
do i=1,nat
!Only concerned with the element tags right now
 read(9,*,iostat=ioerror) elem, ref(:,i), vel(:,i,1)
 if (ioerror.ne.0) then
  write(*,*) 'Error reading molecular data'
  write(*,*) 'Expecting atomic symbol, coordinates, velocities'
  stop
 end if
 call get_mass(elem,atmass(i),atgrp(i,1))
 com(:)=com(:)+atmass(i)*ref(:,i)
 altags(i)=elem
 if (ignoreH .and. elem.eq.'H ') then
  mask(i)=.false.
 end if
!doesn't enter loop on i==1
 if(peratom) then
   atgrp(ig,2)=ig
   grptags(ig)=elem
   ig=ig+1
 else
 do j=1,i-1
  if (atgrp(i,1).eq.atgrp(j,1)) then
   atgrp(i,2)=atgrp(j,2)
   exit
  end if
 end do
 if (atgrp(i,2).eq.0) then
  atgrp(i,2)=ig
  grptags(ig)=elem
  ig=ig+1
 end if
endif
end do
write(*,*) 'Found ',ig-1,' unique elements in trajectory file'
!Rewind xyz file
rewind(9)
com=com/sum(atmass)
do i=1,nat
 ref(:,i)=ref(:,i)-com(:)
end do
allocate(vdos_acf(ig,pn))
vdos_acf=(0.d0,0.d0)
ir_acf=(0.d0,0.d0)
cstep_prev=0
is=0
do
 cstep_prev=0
!If end of file, break out of loop.
!If some other error, abort.
!If no error, back up one line.
!Above we checked the first frame
! for whether the data was as expected.
!If there is a read error here we are
! not specific about the problem.
 read(9,*,iostat=ioerror)
 if (ioerror.lt.0) then
  exit
 else if (ioerror.gt.0) then
  write(*,*) ioerror,'READ ERROR'
  stop
 else
  backspace(9)
 end if
 do
  read(9,*,iostat=ioerror)
  if (ioerror.lt.0) then
     backspace(9)
   exit
  else if (ioerror.gt.0) then
   write(*,*) ioerror,'READ ERROR'
   stop
  end if
  read(9,*) cstep, energy, scrdip
  if (cstep.lt.cstep_prev) then
   cstep=cstep_prev
   backspace(9);backspace(9)
   exit
  else if (cstep.gt.skip .and. cstep.le.nsteps+skip) then
   dip(:,cstep-skip)=scrdip
   do i=1,nat
    read(9,*) elem, coord(:,i), vel(:,i,cstep-skip)
    vel(:,i,cstep-skip)=sqrt(atmass(i))*vel(:,i,cstep-skip)
   end do
   if (align_flag) then
    call align(ref,nat,coord,vel(:,:,cstep-skip),dip(:,cstep-skip),atmass,mask)
    write(30,*) nat
    write(30,'(I10,F20.10,3ES15.6)') cstep, energy, dip(:,cstep-skip)
    do i=1,nat
     write(30,'(1X,A2,10X,3F15.8,2X,3F15.8)') altags(i),coord(:,i),vel(:,i,cstep-skip)
    end do
   end if
  else
   do i=1,nat
    read(9,*)
   end do
  end if
  cstep_prev=cstep
 end do
 if (cstep.lt.nsteps+skip) then
  write(*,*) 'End of file before reading requested number of steps'
  write(*,*) 'Make sure steps+skip is less than or equal to number of steps in trajectory'
  write(*,*) 'Aborting calculation'
  stop
 end if
!IR
 do j=jfirst,jlast
  mean=sum(dip(j,:))/dble(nsteps)
  dip(j,:)=dip(j,:)-mean
  v=(0.d0,0.d0)
  v(1:nsteps)=dip(j,:)
  call zfftf(pn,v,wsave)
  v(:)=v(:)/dble(pn)
  v(:)=v(:)*conjg(v(:))
  call zfftb(pn,v,wsave)
  v=v*(dble(pn)/dble(nsteps))
  ir_acf(:)=ir_acf(:)+v(:)
 end do
!VDOS
 do i=1,nat
!  do j=1,3
  do j=jfirst,jlast
   v=(0.d0,0.d0)
   v(1:nsteps)=vel(j,i,:)
   call zfftf(pn,v,wsave)
   v(:)=v(:)/dble(pn)
   v(:)=v(:)*conjg(v(:))
   call zfftb(pn,v,wsave)
   v=v*(dble(pn)/dble(nsteps))
   vdos_acf(atgrp(i,2),:)=vdos_acf(atgrp(i,2),:)+v(:)
   vdos_acf(ig,:)=vdos_acf(ig,:)+v(:)
  end do
 end do
 is=is+1
end do
write(*,*) 'Averaging results over ',is,' trajectories'
write(*,*)
ir_acf=ir_acf/dble(is)
vdos_acf=vdos_acf/dble(is)
deallocate(vel,dip,v)
allocate(ir_fft(pn))
allocate(vdos_fft(ig,pn))
!Rearrange ACF so negative lags come before positive lags
ir_fft=(0.d0,0.d0)
ir_fft(:)=ir_acf(:)
ir_acf(1:pn/2)=ir_fft(pn/2+1:pn)
ir_acf(pn/2+1:pn)=ir_fft(1:pn/2)
ir_fft(:)=ir_acf(:)
vdos_fft=(0.d0,0.d0)
vdos_fft(:,:)=vdos_acf(:,:)
vdos_acf(:,1:pn/2)=vdos_fft(:,pn/2+1:pn)
vdos_acf(:,pn/2+1:pn)=vdos_fft(:,1:pn/2)
vdos_fft(:,:)=vdos_acf(:,:)
!Apply exponential window to ACF before FFT
do j=1,pn
 ir_fft(j)=ir_fft(j)*exp(-abs((dble(j-1)-dble(pn/2))*ts)*Pi/damp)
 vdos_fft(:,j)=vdos_fft(:,j)*exp(-abs((dble(j-1)-dble(pn/2))*ts)*Pi/damp)
end do
call zfftf(pn,ir_fft,wsave)
do i=1,ig
 call zfftf(pn,vdos_fft(i,:),wsave)
end do
write(10,'(A17,2A20)') 'Freq. (cm^{-1})', 'IR Inten.', 'QC-IR Inten.'
write(20,'((A17),$)') 'Freq. (cm^{-1})'
vdos_tag='VDOS-'
if(jlast.eq.1) then
   vdos_tag='VDOSx-'
elseif (jlast.eq.2) then
   vdos_tag='VDOSy-'
elseif (jfirst.eq.3) then
   vdos_tag='VDOSz-'
endif
do i=1,ig-1
 write(20,'((A20),$)') vdos_tag//grptags(i)
end do
write(20,'(A20)') vdos_tag//'all'
do i=2,pn
 freq_au=dble(i-1)/(dble(pn)*ts)
 freq_cm=(dble(i-1)/(dble(pn)*ts*au2fs))/c
 if (freq_cm.gt.specmax) exit
 qc=(hbar*freq_au)/(kb*temper*(1.d0-exp(-hbar*freq_au/(kb*temper))))
 write(10,'(f17.4,2es20.6)') freq_cm, abs(ir_fft(i)), abs(ir_fft(i))*qc
 write(20,'((f17.4),$)') freq_cm
 do j=1,ig
  write(20,'((es20.6),$)') abs(vdos_fft(j,i))
 end do
 write(20,*)
end do 
close(10)
close(20)
if (align_flag) close(30)
deallocate(ir_acf,ir_fft)
deallocate(vdos_acf,vdos_fft,wsave)

end program qmd_analysis

subroutine usage_mess
implicit none

 write(*,*)
 write(*,*) 'USAGE: qmd_analysis -xyz <xyz file> -steps <integer> -ts <real>'
 write(*,*) ' -xyz <file>:   xyz trajectory file, can concatenate multiple trajectories to average them'
 write(*,*) ' -steps <int>:  number of time steps to use in analysis'
 write(*,*) ' -ts <real>:    time step of the simulation in atomic units'
 write(*,*) ' Optional arguments'
 write(*,*) ' -align:        align molecules in trajectory to orientation of first frame (default off)'
 write(*,*) ' -ignH:         do NOT use hydrogen atoms for alignment (default use hydrogens)'
 write(*,*) ' -width <real>: broadening for the spectrum in wavenumbers, optional (default 10 cm^{-1})'
 write(*,*) ' -skip <int>:   number of time steps to skip in the beginning, optional (default 0)'
 write(*,*) ' -smax <real>:  maximum frequency for output in wavenumbers, optional (default 5000 cm^{-1})'
 write(*,*) ' -temp <real>:  temperature of the simulation in Kelvin (default 298.15 K)'
 write(*,*) '                 used for the quantum correction to the correlation function'
 write(*,*) ' -peratom:      write VDOS for each atom (instead of atom type)'
 write(*,*) ' -peratomx:     write x-component of VDOS for each atom '
 write(*,*) ' -peratomy:     write y-component of VDOS for each atom'
 write(*,*) ' -peratomz:     write z-component of VDOS for each atom'
 write(*,*) ' -help:         print this message and quit'
 write(*,*)

end subroutine usage_mess

subroutine get_mass(elem,mass,atgrp)
implicit none
character(2), intent(in) :: elem
character(2),dimension(109) :: pertab
integer :: i
integer,intent(out) :: atgrp
double precision,dimension(109) :: atommass
double precision,intent(out) :: mass

pertab=(/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',&
'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',&
'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',&
'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',&
'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',&
'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',&
'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',&
'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',&
'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',&
'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',&
'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt'/)

atommass=(/1.007825d0, 4.0026d0,    7.016d0,    9.01218d0, 11.00931d0,&
 12.0d0,     14.00307d0,  15.99491d0, 18.9984d0,  19.99244d0,&
 22.9898d0,  23.98504d0,  26.98154d0, 27.97693d0, 30.97376d0,&
 31.97207d0, 34.96885d0,  39.9624d0,  38.96371d0, 39.96259d0,&
 44.95592d0, 45.948d0,    50.9440d0,  51.9405d0,  54.9381d0,&
 55.9349d0,  58.9332d0,   57.9353d0,  62.9298d0,  63.9291d0,&
 68.9257d0,  73.9219d0,   74.9216d0,  78.9183d0,  79.9165d0,&
 83.912d0,   84.9117d0,   87.9056d0,  88.9054d0,  89.9043d0,&
 92.9060d0,  97.9055d0,   97.9072d0, 101.9037d0, 102.9048d0,&
105.9032d0, 106.90509d0, 113.9036d0, 114.9041d0, 117.9018d0,&
 120.9038d0, 129.9067d0, 126.9004d0, 131.9042d0, 132.9051d0,&
 137.9050d0, 138.9061d0, 139.9053d0, 140.9074d0, 143.9099d0,&
 144.9128d0, 151.9195d0, 152.9209d0, 157.9241d0, 159.9250d0,&
 163.9288d0, 164.9303d0, 165.9304d0, 168.9344d0, 173.9390d0,&
 174.9409d0, 179.9468d0, 180.948d0,  183.9510d0, 186.9560d0,&
 189.9586d0, 192.9633d0, 194.9648d0, 196.9666d0, 201.9706d0,&
 204.9745d0, 207.9766d0, 208.9804d0, 209.9829d0, 210.9875d0,&
 222.0175d0, 223.0198d0, 226.0254d0, 227.0278d0, 232.0382d0,&
 231.0359d0, 238.0508d0, 237.0482d0, 244.0642d0, 243.0614d0,&
 247.0704d0, 247.0703d0, 251.0796d0, 252.0829d0, 257.0950d0,&
 258.0986d0, 259.1009d0, 262.1100d0, 261.1087d0, 262.1138d0,&
 266.1219d0, 262.1229d0, 267.1318d0, 268.1388d0 /)

mass=0.d0
do i=1,109
 if (elem==pertab(i)) then
  mass=atommass(i)
  atgrp=i
  exit
 end if
end do
if (mass==0.d0) then
 write(*,*) 'Element not found'
 stop
end if
end subroutine get_mass

subroutine align(ref,nat,coord,vel,dip,mass,mask)
implicit none

!Input/output variables
integer,intent(in) :: nat
logical,dimension(nat),intent(in) :: mask
double precision,dimension(nat),intent(in) :: mass
double precision,dimension(3,nat),intent(in) :: ref
double precision,dimension(3),intent(inout) :: dip
double precision,dimension(3,nat),intent(inout) :: coord, vel
!Local variables
integer :: i, j
double precision,dimension(3) :: rdip
double precision,dimension(3) :: ro_coord, com, ro_vel
double precision,dimension(3,3) :: R, U 
double precision,dimension(4) :: val, q
double precision,dimension(4,4) :: F
!LAPACK
INTEGER :: LDA=4
INTEGER :: LWORK=3*4-1
INTEGER :: INFO
CHARACTER, PARAMETER :: JOBZ = 'V'
CHARACTER, PARAMETER :: UPLO = 'U'
DOUBLE PRECISION,DIMENSION(3*4-1) :: WORK

com=0.d0
do i=1,nat
 com(:)=com(:)+mass(i)*coord(:,i)
end do
com=com/sum(mass)
do i=1,nat
 coord(:,i)=coord(:,i)-com(:)
end do
do i=1,3
 do j=1,3
  R(i,j)=sum(coord(i,:)*ref(j,:),mask)
 end do
end do
F=0.d0
F(1,1)=R(1,1)+R(2,2)+R(3,3)
F(2,1)=R(2,3)-R(3,2);F(1,2)=F(2,1)
F(3,1)=R(3,1)-R(1,3);F(1,3)=F(3,1)
F(4,1)=R(1,2)-R(2,1);F(1,4)=F(4,1)
F(2,2)=R(1,1)-R(2,2)-R(3,3)
F(3,2)=R(1,2)+R(2,1);F(2,3)=F(3,2)
F(4,2)=R(1,3)+R(3,1);F(2,4)=F(4,2)
F(3,3)=R(2,2)-R(1,1)-R(3,3)
F(4,3)=R(2,3)+R(3,2);F(3,4)=F(4,3)
F(4,4)=R(3,3)-R(1,1)-R(2,2)
call DSYEV(JOBZ,UPLO,4,F,LDA,val,WORK,LWORK,INFO)
q(:)=F(:,4)
U(1,1)=q(1)**2+q(2)**2-q(3)**2-q(4)**2
U(2,1)=2.d0*(q(2)*q(3)+q(1)*q(4))
U(3,1)=2.d0*(q(2)*q(4)-q(1)*q(3))
U(1,2)=2.d0*(q(2)*q(3)-q(1)*q(4))
U(2,2)=q(1)**2-q(2)**2+q(3)**2-q(4)**2
U(3,2)=2.d0*(q(3)*q(4)+q(1)*q(2))
U(1,3)=2.d0*(q(2)*q(4)+q(1)*q(3))
U(2,3)=2.d0*(q(3)*q(4)-q(1)*q(2))
U(3,3)=q(1)**2-q(2)**2-q(3)**2+q(4)**2
rdip(:)=matmul(U(:,:),dip(:))
dip(:)=rdip(:)
do i=1,nat
 ro_coord(:)=matmul(U(:,:),coord(:,i))
 coord(:,i)=ro_coord(:)
 ro_vel(:)=matmul(U(:,:),vel(:,i))
 vel(:,i)=ro_vel(:)
end do

end subroutine align


!The following were taken from
!                      DFFTPACK V1.0
!*****************************************************************
!
!        A Double precision clone by Hugh C. Pumphrey  of:
!
!                      FFTPACK
!               version 4  april 1985
!
!     a package of fortran subprograms for the fast fourier
!      transform of periodic and other symmetric sequences
!
!                         by
!
!                  paul n swarztrauber
!
!  national center for atmospheric research  boulder,colorado 80307
!
!   which is sponsored by the national science foundation
!
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ZFFTI (N,WSAVE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

      SUBROUTINE ZFFTF (N,C,WSAVE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       C(1)       ,WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

      SUBROUTINE ZFFTB (N,C,WSAVE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       C(1)       ,WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

      SUBROUTINE CFFTI1 (N,WA,IFAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI =  6.28318530717958647692D0
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.0D0
            WA(I) = 0.0D0
            LD = LD+L1
            FI = 0.0D0
            ARGLD = FLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END

      SUBROUTINE CFFTF1 (N,C,CH,WA,IFAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

      SUBROUTINE CFFTB1 (N,C,CH,WA,IFAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

      SUBROUTINE PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,&
                     C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP),&
                     CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END

      SUBROUTINE PASSF2 (IDO,L1,CC,CH,WA1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,&
                     WA1(1)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSF3 (IDO,L1,CC,CH,WA1,WA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,&
                     WA1(1)     ,WA2(1)
!     *** TAUI IS -SQRT(3)/2 ***
      DATA TAUR,TAUI /-0.5D0,-0.86602540378443864676D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,&
                     WA1(1)     ,WA2(1)     ,WA3(1)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,&
                     WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
!     *** TR11=COS(2*PI/5), TI11=-SIN(2*PI/5)
!     *** TR12=-COS(4*PI/5), TI12=-SIN(4*PI/5)  
      DATA TR11,TI11,TR12,TI12 /0.3090169943749474241D0,&
          -0.95105651629515357212D0,&
          -0.8090169943749474241D0, -0.58778525229247312917D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,&
                     C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP),&
                     CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END

      SUBROUTINE PASSB2 (IDO,L1,CC,CH,WA1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,&
                     WA1(1)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,&
                     WA1(1)     ,WA2(1)
!     *** TAUI IS SQRT(3)/2 *** 
      DATA TAUR,TAUI /-0.5D0,0.86602540378443864676D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,&
                     WA1(1)     ,WA2(1)     ,WA3(1)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,&
                     WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
!     *** TR11=COS(2*PI/5), TI11=SIN(2*PI/5)
!     *** TR12=COS(4*PI/5), TI12=SIN(4*PI/5)      
      DATA TR11,TI11,TR12,TI12 /0.3090169943749474241D0,&
          0.95105651629515357212D0,&
          -0.8090169943749474241D0,0.58778525229247312917D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
