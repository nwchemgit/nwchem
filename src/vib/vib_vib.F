      subroutine vib_vib(rtdb,hess_file,in_file,hess_ma,in_ma,
     &    hess_ga,in_ga,natomin)
* $Id: vib_vib.F,v 1.2 1995-07-31 06:09:04 d3e129 Exp $
*:: utah messkit code adapted to nwchem by Ricky A. Kendall July 1995n      
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
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL PROJEC,ZEROPE,HESOUT,INTERN
*      CHARACTER*7 INPFIL
      COMMON /HESS/ NATOM,NAT3,NHESS,NHESST    ! hessian information
#include "stdio.fh"
#include "mafdecls.fh"
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
      WRITE(luout,*)' Vibrational analysis via the FX method '
      WRITE(luout,*)
      WRITE(luout,*)' See chapter 2 in "Molecular Vibrations" by ',
     +          'Wilson, Decius and Cross'
      WRITE(luout,*)
*... check input logic 
      status =           (in_file.and.(in_ma.or.in_ga))
      status = status.or.(in_ma.and.(in_file.or.in_ga))
      status = status.or.(in_ga.and.(in_file.or.in_ma))
      if (status) then
        write(luout,*)' ERROR: more than one source for hessian '
        write(luout,*)' in_file :',in_file
        write(luout,*)' in_ma   :',in_ma
        write(luout,*)' in_ga   :',in_ga
        call errquit(' vib_vib: error ',911)
      endif
      if (in_ga)
     &    call errquit
     &    ('vib_vib: ga access to hessian not implemented yet',911)
C
C Zero core
C
      call vib_setup ! subroutine to set up some constants
      NATOM  =  natomin ! number of atoms in species.
      IF (NATOM.LE.1) THEN      ! check for incorrect number of atoms
          WRITE(6,*)' You want to calculate the vibrational ',
     +              'frequencies for ',NATOM,' atoms?'
          WRITE(6,*)' Unfortunately this is not possible '
          CALL errquit('vib_vib: bomb',911)
      ENDIF
      NAT3   =  NATOM*3         ! 3-N (as in degrees of freedom)
      NHESS  =  NAT3*NAT3       ! dimension of hessian
      NHESST =  NAT3*(NAT3+1)/2 ! dimension of lower triangular hessian
      NELS   =  7*MAX(3*NAT-6,1)
      NPRI = 0
C
C Calculate pointers
C
      IHESS    =  1                ! square hessian
      IHESST   =  IHESS  + NHESS   ! lower-tri Hessian
      ICOORD   =  IHESST + NHESST  ! geometrical coordinates
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
      ITOT     =  I40    + 8*NAT3    !
      itot = itot + 2*natom+1 + 6*nat3 ! extra for call to rdinp
c
      if (.not.ma_push_get
     &    (MT_DBL,itot,' core for vib ',h_core, i_core))
     &    call errquit('vib_vib: ma_push_get failed ',911)
C
C Reset pointers for MA array
C
      IHESS    =  i_core           ! square hessian
      IHESST   =  IHESS  + NHESS   ! lower-tri Hessian
      ICOORD   =  IHESST + NHESST  ! geometrical coordinates
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
      ITOT     =  I40    + 8*NAT3    !
c
c read/load hessian and form triangle/square as needed
c
      if (in_ma) then
        ihess = hess_ma   ! simply reset ptr to dbl_mb
*        form triangle
        call vib_dtrngl(dbl_mb(ihess),dbl_mb(ihesst),nat3,nat3)
      endif
      if (in_file) then
        open(unit=69,file=hess_file,form='formatted',status='old',
     &      err=99900,access='sequential')
        do iii = 0,(nhesst-1)
          read(69,*,err=99901,end=99902)dbl_tmp
          dbl_mb(ihesst+iii) = dbl_tmp
        enddo
        close(unit=69,status='keep')
        call vib_dsquar(dbl_mb(ihesst),dbl_mb(ihess),nat3,nat3)
      endif
C
C Read in user input and tape10 arrays.
C      Note: ! scratch pointer for atom charges (real NAT words)
C              and atom lables (real 2*nat words)
      IDUM=ITOT+1
      IVC=ITOT+2*NATOM+1 ! actually an integer, but real length used
      ITOT=IVC+6*NAT3
      call vib_rdinp(
     &    dbl_mb(ihess),dbl_mb(ihesst),dbl_mb(icoord),
     &    dbl_mb(imass),dbl_mb(iscr),  dbl_mb(idum),
     &    dbl_mb(i10),nels,projec,zerope,hesout,intern,
     &    rtdb)
      if (projec) then
        call vib_eckart( dbl_mb(ihess), dbl_mb(ihessp), dbl_mb(ihesst),
     &      dbl_mb(icoord),  dbl_mb(ivc) )
      end if
*rak:c****
*rak:c**** write projected Hessian to tape 10
*rak:C****
*rak:      IF (HESOUT) THEN
*rak:        IPH = IP10(57) ! pointer to hessian
*rak:C
*rak:        CALL WWRITW(ITAP10,DBL_MB(IHESS),INTOWP(NAT3*NAT3),
*rak:     &    IPH,IJUNK) ! Write it
*rak:        CALL MATOUT(DBL_MB(IHESS),NAT3,NAT3,NAT3)
*rak:      END IF
*rak:C
* rak dfill
      CALL Dfill(NAT3,0.0d00,DBL_MB(ISCR),1) ! zero scratch used
C
      CALL vib_hmass(DBL_MB(IHESST),DBL_MB(IMASS))         ! mass weight and scale hessian
C
C Diagonalize mass-weighted, scaled hessian matrix
C     Note: ! scratch pointer for givens (real 5*NAT3 words)
c use hessp as scratch now calling rsg
C
      CALL vib_CALLG(DBL_MB(IHESSt),nhesst,DBL_MB(IHESSP),
     &    dbl_mb(iscr),dbl_mb(iscr+nat3),DBL_MB(IEGVAL),
     &    DBL_MB(IEGVEC), NAT3,NAT3)
C
C
C
      CALL vib_NMASS(DBL_MB(IEGVEC),DBL_MB(IMASS)) ! "unmass" weight the normal modes.
C
C *** Note: DBL_MB(IHESST) now destroyed if needed reinitialize from DBL_MB(IHESS)
C
* rak dfill
      call dfill(5*nat3,0.0d00,dbl_mb(iscr),1)    ! zero scratch used
C
      CALL vib_WRTFREQ(DBL_MB(IEGVAL),NAT3,ZEROPE,NPRI) ! Write out the zero-point energy
C
      CALL vib_CLEAN(DBL_MB(IEGVEC),NAT3*NAT3,1.0D-07) ! CLEAN eigenvectors
      WRITE(6,9000)
      CALL vib_PRDATA(DBL_MB(IEGVEC),DBL_MB(IEGVAL),
     &    NAT3,NAT3,NAT3,NAT3,.TRUE.,NPRI) ! Write eigenvalues and eigenvectors
C
C----------Call GAMESS subroutines which express the HESSIAN and the normal 
C----------modes in terms of internal coordinates.
C
      IF(INTERN)CALL vib_BCALC(DBL_MB(I10),DBL_MB(I20),DBL_MB(I30),
     &    DBL_MB(I40),DBL_MB(ITOT),
     &    DBL_MB(IHESS),DBL_MB(IEGVEC),DBL_MB(IEGVAL))
C
 9000 FORMAT(///,10X,'-------------------------------------------------'
     &        ,/,10X,'NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES'
     &        ,/,10X,'-------------------------------------------------'
     &        ,/,10X,'       (Freqencies expressed in cm-1)')
      if (.not.ma_pop_stack(h_core)) call errquit
     &    ('vib_rdinp ma_pop failed',911)
      return
99900 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit('vib_vib: error opening file: "hess_file"',911)
99901 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit('vib_vib: error reading file: "hess_file"',911)
99902 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit
     & ('vib_vib: unexpected EOF when reading file: "hess_file"',911)
      END
