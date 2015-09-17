      subroutine vib_vib(rtdb,hess_file,in_file,hess_ma,in_ma,
     &    hess_ga,in_ga,natomin)
* $Id$
*:: utah messkit code adapted to nwchem by Ricky A. Kendall July 1995
*utah:C
*utah:C   This program was written to mass weight and diagonalize the force
*utah:C        constant matrix to get the harmonic frequencies.  The program
*utah:C        will read the Updated hessian from tape 11 or the analytic
*utah:C        hessian (force constant matrix) from tape10.  You may use
*utah:C        default masses or enter in a mass for the atom.
*utah:C   You need tape10 for the "equilibrium" geometry and tape11 if the
*utah:C        Updated Hessian matrix is to be used.
*utah:C
*utah:C
*utah:C   Written by R.A. Kendall
*utah:C   start date: 18-JAN-88
*utah:C   mutilated by H. Taylor (TO MAKE IT WORK!!!!)
*utah:C   This version has the option for the zero point energy calculation. xiao
*utah:C   Internal coordinate representation of Hessian and Normal modes added
*utah:C   by J. Anchell (12-DEC-1990).  Subroutines which generate the Cartesian
*utah:C   to internal coordinate transformation matrix was lifted from GAMESS.
*utah:C
*
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
#include "errquit.fh"
      LOGICAL PROJEC,ZEROPE,HESOUT,INTERN
*      CHARACTER*7 INPFIL
      INTEGER NATOM, NAT3, NHESS, NHESST
      COMMON /cvib_HESS/ NATOM,NAT3,NHESS,NHESST    ! hessian information
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c:: passed
      integer rtdb             ! [input] rtdb handle
      character*(*) hess_file  ! [input] name of file storing lower triangular packed hessian
      integer hess_ma          ! [input] MA handle to square hessian
      integer hess_ga          ! [input] GA handle to square hessian
      logical in_file          ! [input] hessian is in file get it there
      logical in_ma            ! [input] hessian is in MA array
      logical in_ga            ! [input] hessian is in GA array
      integer natomin          ! [input] number of atoms
c
      logical status
c
      integer i_core, h_core, iii, ioldlabs, ivc, itot
      integer nels, npri, ihess, icoord, ihesst, ihesstcp,
     &    ihessp, iegval, iegvec, iddpol, iddpolq, intense
      integer imass, iscr, i10, i20, i30, i40,i_w1,l_w1,i_w2,l_w2
      integer isystype
      double precision dbl_tmp
      logical first_pass
      character*255 dipole_file
      logical dipole_file_exists
      logical animation_on
      logical iscrystal

      integer  vib_isystype
      external vib_isystype
c
      first_pass = .true.
*... check input logic 
      status =           (in_file.and.(in_ma.or.in_ga))
      status = status.or.(in_ma.and.(in_file.or.in_ga))
      status = status.or.(in_ga.and.(in_file.or.in_ma))
      if (status) then
        write(luout,*)' ERROR: more than one source for hessian '
        write(luout,*)' in_file :',in_file
        write(luout,*)' in_ma   :',in_ma
        write(luout,*)' in_ga   :',in_ga
        call errquit(' vib_vib: error ',911, UNKNOWN_ERR)
      endif
      if (in_ga)
     &    call errquit
     &    ('vib_vib: ga access to hessian not implemented yet',911,
     &       CAPMIS_ERR)
C
C Zero core
C
      isystype = vib_isystype(rtdb)
      iscrystal = (isystype.ne.0)

      call vib_setup ! subroutine to set up some constants
      NATOM  =  natomin ! number of atoms in species.
      IF (NATOM.LT.1) THEN      ! check for incorrect number of atoms
          WRITE(6,*)' You want to calculate the vibrational ',
     +              'frequencies for ',NATOM,' atoms?'
          WRITE(6,*)' Unfortunately this is not possible '
          CALL errquit('vib_vib: bomb',911, INPUT_ERR)
      ENDIF
      NAT3   =  NATOM*3         ! 3-N (as in degrees of freedom)
      NHESS  =  NAT3*NAT3       ! dimension of hessian
      NHESST =  NAT3*(NAT3+1)/2 ! dimension of lower triangular hessian
      NELS   =  7*MAX(3*NATOM-6,1)
      NPRI = 0
C
C Calculate pointers
C
      IHESS    =  1                ! square hessian
      IHESST   =  IHESS  + NHESS   ! lower-tri Hessian
      ihesstcp =  IHESST + NHESST  ! copy of lower-tri hessian
      ICOORD   =  IHESSTcp + NHESST  ! geometrical coordinates
      IMASS    =  ICOORD + NAT3    ! mass of each atom
      IEGVAL   =  IMASS  + NATOM   ! eigenvalues from Hessian matrix
      IEGVEC   =  IEGVAL + NAT3    ! eigenvectors from Hessian matrix
      ISCR     =  IEGVEC + NHESS   ! dynamic bottom of core array
      IHESSP   =  ISCR   + 8*NAT3  ! addition of scratch space needed
C--------The following are pointers for GAMESS internal coordinate subroutines.
C
      I10      =  IHESSP + NAT3*NAT3 ! space for zmat
      I20      =  I10    + NELS
      I30      =  I20    + NAT3*NAT3 ! Space to represent internal coord. Hessian
      I40      =  I30    + NAT3*NAT3 !
      Iddpol   =  I40    + 8*NAT3    ! derivative dipole in cartesians
      Iddpolq  =  Iddpol + 3*NAT3    ! derivative dipole in normal modes
      Intense  =  Iddpolq+ 3*NAT3    ! intensities
      ITOT     =  Intense + 3*natom*4
      itot = itot + 2*natom+1 + 6*nat3 ! extra for call to rdinp
c
      if (.not.ma_push_get
     &    (MT_DBL,itot,' core for vib ',h_core, i_core))
     &    call errquit('vib_vib: ma_push_get failed ',911, MA_ERR)
C
C Reset pointers for MA array
C
      IHESS    =  i_core           ! square hessian
      IHESST   =  IHESS  + NHESS   ! lower-tri Hessian
      ihesstcp =  IHESST + NHESST  ! copy lower-tri Hessian
      ICOORD   =  IHESSTcp + NHESST  ! geometrical coordinates
      IMASS    =  ICOORD + NAT3    ! mass of each atom
      IEGVAL   =  IMASS  + NATOM   ! eigenvalues from Hessian matrix
      IEGVEC   =  IEGVAL + NAT3    ! eigenvectors from Hessian matrix
      ISCR     =  IEGVEC + NHESS   ! dynamic bottom of core array
      IHESSP   =  ISCR   + 8*NAT3  ! addition of scratch space needed
C--------The following are pointers for GAMESS internal coordinate subroutines.
C
      I10      =  IHESSP + NAT3*NAT3 ! space for zmat
      I20      =  I10    + NELS
      I30      =  I20    + NAT3*NAT3 ! Space to represent internal coord. Hessian
      I40      =  I30    + NAT3*NAT3 !
      ioldlabs =  I40    + 8*NAT3    !
      ivc      =  ioldlabs + 2*natom + 1
      iddpol   =  ivc + 6*nat3
      iddpolq  =  iddpol + 3*nat3
      intense  =  iddpolq + 3*nat3
      Itot     =  intense + 3*natom*4
c
c read/load hessian and form triangle/square as needed
c
      if (in_ma) then
        ihess = hess_ma   ! simply reset ptr to dbl_mb
*        form triangle
        call vib_DTRNGL(dbl_mb(ihess),dbl_mb(ihesst),nat3,nat3)
      endif
      if (in_file) then
        open(unit=69,file=hess_file,form='formatted',status='old',
     &      err=99900,access='sequential')
        do iii = 0,(nhesst-1)
          read(69,*,err=99901,end=99902)dbl_tmp
          dbl_mb(ihesst+iii) = dbl_tmp
        enddo
        close(unit=69,status='keep')
        call vib_DSQUAR(dbl_mb(ihesst),dbl_mb(ihess),nat3,nat3)
      endif
      call util_file_name('fd_ddipole',.false., .false.,dipole_file)
      dipole_file_exists = .false.
      inquire(file=dipole_file,exist=dipole_file_exists)
      if (dipole_file_exists) then
        open(unit=70,file=dipole_file,form='formatted',status='old',
     &      err=89900,access='sequential')
        do iii = 0,((3*nat3)-1)
          read(70,*,err=89901,end=89902) dbl_tmp
          dbl_mb(iddpol+iii) = dbl_tmp
        enddo
        close(unit=70,status='keep')
      endif
00001 continue
      if (.not.iscrystal) then
      write(luout,*)
      write(luout,*)
      if (first_pass) then
        WRITE(luout,*)' Vibrational analysis via the FX method '
        WRITE(luout,*)
        WRITE(luout,*)' See chapter 2 in "Molecular Vibrations" by ',
     +      'Wilson, Decius and Cross'
        WRITE(luout,*)
      else
        WRITE(luout,*)
     &      '       Vibrational analysis via the FX method '
        write(luout,*)
     &      ' --- with translations and rotations projected out ---'
        write(luout,*)
     &      ' --- via the Eckart algorithm                      ---'
      endif
      end if

      if (first_pass) then
c
c save a copy of hesst
c
        call dcopy(nhesst,dbl_mb(ihesst),1,dbl_mb(ihesstcp),1)
      else
c
c restore copy of hesst and hess
c
        call dcopy(nhesst,dbl_mb(ihesstcp),1,dbl_mb(ihesst),1)
        call vib_DSQUAR(dbl_mb(ihesst),dbl_mb(ihess),nat3,nat3)
      endif

c
c set the hessian to zero if atom and not crystal
c

C
C Read in user input and tape10 arrays.  NO INPUT REQUIRED NOW
C      Note: ! scratch pointer for atom charges (real NAT words)
C              and atom lables (real 2*nat words)
      call vib_RDINP(
     &    dbl_mb(ihess),dbl_mb(ihesst),dbl_mb(icoord),
     &    dbl_mb(imass),dbl_mb(iscr),  dbl_mb(ioldlabs),
     &    dbl_mb(i10),nels,projec,zerope,hesout,intern,
     &    rtdb,first_pass)

*     **** call vib_phonon spectra if crystal and then exit ****
      if (iscrystal) then
         call vib_phonon(rtdb,NATOM,dbl_mb(icoord),
     >                   NAT3,dbl_mb(ihess),dbl_mb(imass))
         goto 9999
      end if
*     *** otherwise continuue on with implied loop insantity ***

      if (projec) then
        if (.not.ma_push_get
     &    (MT_DBL,nat3*nat3,' w1 ',l_w1, i_w1))
     &    call errquit('vib_vib: ma_push_get failed ',911, MA_ERR)
        if (.not.ma_push_get
     &    (MT_DBL,nat3*nat3,' w2 ',l_w2, i_w2))
     &    call errquit('vib_vib: ma_push_get failed ',911, MA_ERR)

        call vib_eckart(dbl_mb(ihess),dbl_mb(ihessp),dbl_mb(ihesst),
     &         dbl_mb(icoord),dbl_mb(ivc),dbl_mb(i_w1),dbl_mb(i_w2))
           
        if(.not.ma_chop_stack(l_w1)) call errquit(
     '       ' vib_vib: machopstack failed',1, MA_ERR)
      end if
* rak dfill
      CALL Dfill(NAT3,0.0d00,DBL_MB(ISCR),1) ! zero scratch used
C
      CALL vib_HMASS(DBL_MB(IHESST),DBL_MB(IMASS),first_pass) ! mass weight and scale hessian
C
C Diagonalize mass-weighted, scaled hessian matrix
C     Note: ! scratch pointer for givens (real 5*NAT3 words)
c use hessp as scratch now calling rsg
C

      if (NATOM.gt.1) then
         CALL vib_CALLG(DBL_MB(IHESSt),nhesst,DBL_MB(IHESSP),
     &               dbl_mb(iscr),dbl_mb(iscr+nat3),DBL_MB(IEGVAL),
     &               DBL_MB(IEGVEC), NAT3,NAT3)
         call vib_vecphase(dbl_mb(icoord),dbl_mb(iegvec),nat3)
      else
         call dcopy(3,0.0d0,0,dbl_mb(IEGVAL),1)
         call dcopy(3*3,0.0d0,0,dbl_mb(IEGVEC),1)
      end if
      call dfill(5*nat3,0.0d00,dbl_mb(iscr),1)    ! zero scratch used
C
         CALL vib_WRTFREQ(rtdb,DBL_MB(IEGVAL),NAT3,ZEROPE,NPRI) ! Write out the zero-point energy

C
C  ---Added by Bojana September 11,2008 
        if(.not. rtdb_put(rtdb,'mc_data:eigenvalues',mt_dbl,nat3,
     &      dbl_mb(iegval)))
     &      call errquit('vib_vib:rtdb_put of eigenvalues failed',
     &      555, RTDB_ERR)

C      call printeig(nat3,nhess,dbl_mb(iegval),dbl_mb(iegvec))
C  ---End of addition by Bojana, September 11, 2008 
C
      CALL vib_NMASS(DBL_MB(IEGVEC),DBL_MB(IMASS)) ! "unmass" weight the normal modes.
C  ---Added by Bojana September 30,2008
        if(.not. rtdb_put(rtdb,'mc_data:eigenvectors',mt_dbl,nhess,
     &      dbl_mb(iegvec)))
     &      call errquit('vib_vib:rtdb_put of eigenvectors failed',
     &      555, RTDB_ERR)
c      call printeig(nat3,nhess,dbl_mb(iegval),dbl_mb(iegvec))
C  ---End of addition by Bojana, September 30, 2008

*      call vib_vecnormal(dbl_mb(iegvec),nat3)
C
C *** Note: DBL_MB(IHESST) now destroyed if needed reinitialize from DBL_MB(IHESS)
C
* rak dfill
c      call dfill(5*nat3,0.0d00,dbl_mb(iscr),1)    ! zero scratch used
C
c      CALL vib_WRTFREQ(rtdb,DBL_MB(IEGVAL),NAT3,ZEROPE,NPRI) ! Write out the zero-point energy
c      call printeig(nat3,nhess,dbl_mb(iegval),dbl_mb(iegvec))
C
      CALL vib_CLEAN(DBL_MB(IEGVEC),NAT3*NAT3,1.0D-27) ! CLEAN eigenvectors
      if (first_pass) then
        WRITE(6,9000)
      else
        WRITE(6,9001)
      endif
      CALL vib_PRDATA(rtdb,DBL_MB(IEGVEC),DBL_MB(IEGVAL),
     &    NAT3,NAT3,NAT3,NAT3,.TRUE.,NPRI,
     &    first_pass)     ! Write eigenvalues and eigenvectors
      if (dipole_file_exists) then
        call vib_intense(rtdb,dbl_mb(iegvec),dbl_mb(iegval),natom,
     &      dbl_mb(iddpol),dbl_mb(iddpolq),dbl_mb(intense),
     &      first_pass)
      endif
C
C----------Call GAMESS subroutines which express the HESSIAN and the normal 
C----------modes in terms of internal coordinates.
C
      IF(INTERN)CALL vib_BCALC(DBL_MB(I10),DBL_MB(I20),DBL_MB(I30),
     &    DBL_MB(I40),DBL_MB(ITOT),
     &    DBL_MB(IHESS),DBL_MB(IEGVEC),DBL_MB(IEGVAL))
c
      if (.not.first_pass) then
* if any negative eigenvalues print out steps in their direction
        call vib_istep(
     &          rtdb,nat3,natom,
     &          dbl_mb(iegvec),dbl_mb(iegval),
     &          dbl_mb(icoord),dbl_mb(iscr),
     &          dbl_mb(iscr+nat3),dbl_mb(iscr+(2*nat3)))
        if(.not.rtdb_get(rtdb,'vib:animate',mt_log,1,animation_on)) then
          animation_on = .false. ! no animation
        endif
        write(luout,*)'vib:animation ', animation_on
          
        if (animation_on)
     &        call vib_modestep(
     &        rtdb,nat3,natom,
     &        dbl_mb(iegvec),dbl_mb(iegval),
     &        dbl_mb(icoord),dbl_mb(iscr),
     &        dbl_mb(iscr+nat3),dbl_mb(iscr+(2*nat3)))
* compute thermodynamic properties from frequencies
        call vib_thermo(rtdb,natom,nat3,
     &        dbl_mb(iegval),
     &        dbl_mb(imass),
     &        dbl_mb(icoord),
     &        dbl_mb(iscr))            
      endif
c
      if (first_pass) then
        first_pass = .false.
        goto 00001
      endif
C
 9000 FORMAT(///,10X,'-------------------------------------------------'
     &        ,/,10X,'NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES'
     &        ,/,10X,'-------------------------------------------------'
     &        ,/,10X,'       (Frequencies expressed in cm-1)')
 9001 FORMAT(///,10X,'-------------------------------------------------'
     &        ,/,10X,'NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES'
     &        ,/,10X,'-------------------------------------------------'
     &        ,/,10X,'   (Projected Frequencies expressed in cm-1)')
 9999 if (.not.ma_pop_stack(h_core)) call errquit
     &    ('vib_rdinp ma_pop failed',911, MA_ERR)
      return
89900 continue
      write(luout,*)'dipole_file => ',dipole_file
      call errquit('vib_vib: error opening file: "dipole_file"',811,
     &       DISK_ERR)
89901 continue
      write(luout,*)'dipole_file => ',dipole_file
      call errquit('vib_vib: error reading file: "dipole_file"',811,
     &       DISK_ERR)
89902 continue
      write(luout,*)'dipole_file => ',dipole_file
      call errquit
     & ('vib_vib: unexpected EOF when reading file: "dipole_file"',811,
     &       DISK_ERR)
99900 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit('vib_vib: error opening file: "hess_file"',911,
     &       DISK_ERR)
99901 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit('vib_vib: error reading file: "hess_file"',911,
     &       DISK_ERR)
99902 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit
     & ('vib_vib: unexpected EOF when reading file: "hess_file"',911,
     &       DISK_ERR)
      END


C--- Added by Bojana
      Subroutine printeig(nat3,nhess,eigval,eigvec)

       integer nat3, nhess, i 
       double precision eigval(nat3),eigvec(nhess)

       write(*,*) 'Printing Eigenvalues in vib_vib:',nat3 
       do i=1, nat3 
        write(*,*) eigval(i)
       enddo  
       write(*,*) 'Printing Eigenvvectors in vib_vib:',nhess
       do i=1, nat3
         write(*,*) 'Vector: ',i
         do j=1,nat3
          write(*,*) eigvec((i-1)*nat3+j)
         enddo
       enddo


      END 

C--- Added by Bojana  
