      SUBROUTINE vib_RDINP(
     &    HESS,HESST,COORD,ATMASS,ZAN,OLDLBLAT,IZMAT,NELS,
     &    PROJEC,ZEROPE,HESOUT,INTERN,
     &    rtdb,first_pass)
* $Id: vib_rdinp.F,v 1.7 1999-07-08 19:11:02 d3e129 Exp $
C
C This routine reads the user input and the hessian matrix from
C      the appropriate place.  The default is to use the "default"
C      masses and to read the analytic hessian from tape 10.
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"
*      CHARACTER*7 UNITNM,NAME
*      CHARACTER*6 USEROP         ! user input option
      CHARACTER*10 ULABEL        ! label to locate program input
      INTEGER NELS
      LOGICAL UPDATE,SETMAS,PROJEC,ZEROPE,HESOUT,INTERN
      logical dum_log
      integer natom, nat3, nhess, nhesst
      logical first_pass
      DOUBLE PRECISION HESS(NAT3,NAT3)  ! Hessian matrix
      DOUBLE PRECISION HESST(NHESST)    ! Lower triangle of hessian matrix
      DOUBLE PRECISION COORD(3,NATOM)   ! Atom x,y,z coordinates
      DOUBLE PRECISION ATMASS(NATOM)    ! Atoms mass
      DOUBLE PRECISION ZAN(NATOM)       ! Nuclear charge (i.e. atomic number)
      INTEGER OLDLBLAT(2,NATOM)     ! Atom character labels
      INTEGER IZMAT(NELS)
      character*16 tag
      character*22 element
      character*2 symbol
      double precision q
      double precision xyz(3)
      COMMON /cvib_HESS/ NATOM,NAT3,NHESS,NHESST   ! Hessian information
      integer numans
      double precision ams, wave
      COMMON /cvib_SETCON/ AMS(36),WAVE,NUMANS     ! setup parameters
c
      integer h_lblb, i_lblb, geom
      integer il_cnt, iii, ijunk
c
      integer rtdb
*rak:      DATA ITAP11 /11/
      DATA ULABEL/'U:VIB  '/   ! program input locator label
      DATA SETMAS/.TRUE./         ! logical USE default masses
      UPDATE = .FALSE.            ! logical READ update Hessian
      HESOUT = .FALSE.            ! do not write Hessian to tape 10
*rak:C
*rak:      CALL UTOPEN (5,'OLD')
*rak:C
*rak:      IERROR = 0                   ! set no errors
*rak:      CALL LOCATE(5,ULABEL,IERROR) ! set pointer to input for this module
*rak:      WRITE(6,*)'IERROR:',IERROR
*rak:*rak input options
*rak:      IF( IERROR .NE. 0) THEN      ! no input found for this module
*rak:          WRITE(6,*)' Default input used '
*rak:          UPDATE = .FALSE.          ! analytic Hessian used
*rak:          SETMAS = .TRUE.          ! default masses used
*rak:          PROJEC = .FALSE.
*rak:          ZEROPE = .FALSE.
*rak:          HESOUT = .FALSE.
*rak:          INTERN = .FALSE.
*rak:      ELSE ! input found for this module
*rak:* rak change to rtdb reads        
*rak:00001     CONTINUE
*rak:          READ(5,'(A6)',ERR=00003,END=00002)USEROP  ! read option
*rak:          WRITE(6,'('' option label '',A6,'' found'')')USEROP  !
*rak:          IF(USEROP.EQ.'UPDATE')THEN       ! read updated Hessian
*rak:             UPDATE = .TRUE.
*rak:             GOTO 00001
*rak:          ELSEIF(USEROP.EQ.'ANALYT')THEN   ! read analytic Hessian
*rak:             UPDATE = .FALSE.
*rak:             GOTO 00001
*rak:          ELSEIF(USEROP(1:4).EQ.'MASS')THEN ! use input masses
*rak:             WRITE(6,*)' atomic number and masses read from unit 5'
*rak:             DO 00100 III = 1,NATOM     ! loop over atoms
*rak:                READ(5,*,ERR=99900,END=99900)IREAD,ATMASS(III)
*rak:                IF(ATMASS(III).EQ.0.0D00) THEN  ! zero mass => use default
*rak:                   IF(IREAD.GT.NUMANS)THEN ! check if default mass available
*rak:                      WRITE(6,*)' ATOMIC NUMBER OUT OF RANGE'
*rak:                      WRITE(6,*)' MASSES FOR ATOMIC NUMBERS .LE. ',
*rak:     +                            NUMANS,' ONLY'
*rak:                      WRITE(6,*)' YOU MUST INPUT THE MASS '
*rak:*
*rak:                      CALL errquit('error in rdinp',911) ! if not then kill program execution
*rak:                   ENDIF
*rak:                   ATMASS(III) = AMS(IREAD) ! use default mass
*rak:                ENDIF
*rak:00100        CONTINUE
*rak:             SETMAS = .FALSE.  ! do not reset masses
*rak:             GOTO 00001
*rak:          ELSEIF(USEROP.EQ.'PROJEC')THEN  ! project out translations,rotations
*rak:             PROJEC = .TRUE.
*rak:             GOTO 00001
*rak:          ELSEIF(USEROP.EQ.'ZEROPE')THEN  ! Calculate Zero Point Energy
*rak:             PROJEC = .TRUE.
*rak:             ZEROPE = .TRUE.
*rak:             GOTO 00001
*rak:          ELSEIF(USEROP.EQ.'HESOUT')THEN  ! write out Hessian to tape 10
*rak:             HESOUT = .TRUE.
*rak:             GOTO 00001
*rak:          ELSEIF(USEROP.EQ.'INTERN')THEN  ! Calculate Hessian and normal modes
*rak:             INTERN = .TRUE.              ! in internal coordinates.
*rak:             CALL ZMATIN(IZMAT,NELS)
*rak:          ELSEIF(USEROP(1:2).EQ.'U:')THEN  ! new module input so leave
*rak:             GOTO 00002
*rak:          ELSE
*rak:             WRITE(6,'('' Unrecognizable user input ==>  '',A6)')
*rak:     &             USEROP
*rak:          ENDIF
*rak:      ENDIF
*rak:      GOTO 00002
*rak:00003 CONTINUE   ! error reading input file use default settings
*rak:      UPDATE = .FALSE.
*rak:      SETMAS = .TRUE.
*rak:00002 CONTINUE
*rak:C
*rak:      IPCHG =  IP10( 1)         ! tape10 pointer to atom charges
*rak:      IPLBL =  IP10(19)         ! tape10 pointer to atom labels
*rak:      IPC   =  IP10(42)         ! tape10 pointer to atom coordinates
*rak:C
*rak:*rak: zan not really needed map integer to atomic number instead of charges?
*rak:      CALL WREADW(ITAP10,ZAN,INTOWP(NATOM),IPCHG,IJUNK)  ! read charges
*rak:      CALL WREADW(ITAP10,OLDLBLAT,2*NATOM,IPLBL,IJUNK)      ! read labels
*rak:      CALL WREADW(ITAP10,COORD,INTOWP(NAT3),IPC,IJUNK)   ! read coordinates
      if (first_pass) WRITE(luout,*)' Vib: Default input used '
      UPDATE = .FALSE. ! analytic Hessian used
      SETMAS = .TRUE. ! default masses used
      PROJEC = .FALSE.
      ZEROPE = .FALSE.
      HESOUT = .FALSE.
      INTERN = .FALSE.
      if (rtdb_get(rtdb,'vib:project',MT_LOG,1,dum_log)) then
        projec = dum_log
        write(luout,*)' vib:project option set to',projec
      endif
      if (rtdb_get(rtdb,'vib:zero point energy',MT_LOG,1,dum_log)) then
        zerope = dum_log
        write(luout,*)' vib:zero point energy option set to',zerope
      endif
      if (zerope) projec = .true.
*
* over ride all input to do normal and projected analysis
*
      if (first_pass) then
        projec = .false.
        zerope = .false.
      else
        projec = .true.
        zerope = .true.
      endif
*.. hessian read in vib_vib hess and hesst set there
*rak:      IF (UPDATE) THEN  ! read updated Hessian
*rak:        WRITE(luout,*)' Updated Hessian read from unit 11 '
*rak:        call errquit('vib_rdinp should never be here',911)
*rak:          NAME=UNITNM(11)
*rak:          OPEN (UNIT=ITAP11, FORM='FORMATTED', ACCESS='SEQUENTIAL',
*rak:     +          STATUS='OLD',FILE=NAME)
*rak:          REWIND ITAP11
*rak:          READ (ITAP11,*,END=99901,ERR=99901) ! Skip record #1
*rak:          READ (ITAP11,*,END=99901,ERR=99901) ! Skip record #2
*rak:          READ (ITAP11,*,END=99901,ERR=99901) ! Skip record #3
*rak:          READ (ITAP11,*,END=99901,ERR=99901) ! Skip record #4
*rak:          READ (ITAP11,*,END=99901,ERR=99901) (HESST(I),I=1,NHESST)
*rak:          CLOSE(UNIT=ITAP11,STATUS='KEEP')
*rak:          CALL DSQUAR(HESST,HESS,NAT3,NAT3)   ! make the square Hessian
*rak:C          WRITE(6,'(//A/)')
*rak:C     +    '  ---   Nuclear Hessian (Hartree/Bohr/Bohr)  ---'
*rak:          WRITE(6,9000)
*rak:          CALL PRNTPD(HESST,NHESST,NAT3,6)
*rak:      ELSE ! read analytic Hessian
*rak:          WRITE(6,*)' Analytic Hessian read from tape 10'
*rak:          IPH = IP10(57) ! pointer to Hessian
*rak:          CALL WREADW(ITAP10,HESS,INTOWP(NAT3*NAT3),IPH,IJUNK) ! read it
*rak:          CALL CHKHESS(HESS,NAT3) ! check to make sure it is symmetric
*rak:          CALL DTRNGL(HESS,HESST,NAT3,NAT3) ! for[Am lower-triangle
*rak:          WRITE(6,9000)
*rak:c          CALL PRDATA(HESS,DUM,NAT3,NAT3,NAT3,NAT3,.FALSE.,NPRI)
*rak:          CALL PRNTPD(HESST,NHESST,NAT3,6)
*rak: ENDIF
      call vib_chkhess(hess,nat3,first_pass)
* set labels, zan and coords information
c
      if (.not.ma_push_get(MT_BYTE,(2*natom),' labels for vib',
     &    h_lblb,i_lblb))
     &    call errquit('vib_rdinp ma_get for labels failed',911)
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('vib_rdinp: error creating geometry',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('vib_rdinp: error loading geometry',911)
      il_cnt = 0
      do iii = 1,natom
        if(.not.geom_cent_get(geom,iii,tag,xyz,q))
     &      call errquit('vib_rdinp: geom_cent_get failed',911)
        zan(iii) = q
        coord(1,iii) = xyz(1)
        coord(2,iii) = xyz(2)
        coord(3,iii) = xyz(3)
        if(.not.geom_tag_to_element(tag,symbol,element,ijunk)) then
          if (symbol.ne.'bq')
     &        call errquit('vib_rdinp: tag2elem fail',911)
        endif
        byte_mb(il_cnt+i_lblb) = symbol(1:1)
        il_cnt = il_cnt+1
        byte_mb(il_cnt+i_lblb) = symbol(2:2)
        il_cnt = il_cnt+1
*        oldlblat(1,iii) = symbol(1:1)
*        oldlblat(2,iii) = symbol(2:2)
      enddo
*
*rak:      IF (SETMAS) THEN ! if default masses for all atoms set up here
*rak:        WRITE(6,'(/A/)')
*rak:     +      ' Default masses from atomic number (charge) used'
*rak:        DO 00200 III = 1,NATOM ! loop over atoms
*rak:             IDUM        = INT(ZAN(III))  ! atomic charge =?= atomic number
*rak:             if (idum.gt.36.or.idum.lt.1)
*rak:     &           call errquit
*rak:     &           ('vib_rdinp: error setting masses',911)
*rak:             ATMASS(III) = AMS(IDUM)      ! get mass of atom
*rak:00200     CONTINUE
*rak:      ENDIF
      SETMAS = .true.   ! only place to set mass is in geometry object
      setmas = geom_masses_get(geom,natom,atmass)
      if (.not.setmas)
     &    call errquit
     &    ('vib_rdinp: geom_masses_get failed',911)
      if (.not.geom_destroy(geom)) call errquit
     &    ('vib_rdinp: geom_destroy failed',911)
C
C Write out atom information read here
C
      if (first_pass) then
        WRITE(6,10000)
C
        il_cnt = 0
        DO 00300 III = 1,NATOM
          WRITE(6,10001)
     &        byte_mb(il_cnt+i_lblb),
     &        byte_mb(il_cnt+i_lblb+1),
     &        III,
     &        COORD(1,III),COORD(2,III),COORD(3,III),
     &        ATMASS(III)
          il_cnt = il_cnt + 2
00300   CONTINUE
        WRITE(6,10002)
      endif
      if (.not.ma_pop_stack(h_lblb)) call errquit
     &    ('vib_rdinp ma_pop failed',911)
      RETURN ! leave routine
10001 FORMAT(4x,a1,a1,3x,I5,4(1PD15.7))
10002 FORMAT(1X,74('-')///)
C
 9000 FORMAT(/,10X,'-----------------------------------',/,10X,
     &             'NUCLEAR HESSIAN (Hartree/Bohr/Bohr)',/,10X,
     &             '-----------------------------------',//)
10000 FORMAT(//,1X,28('-'),' Atom information ',28('-'),/,
     +5X,'atom',4X,'#',8X,'X',14X,'Y',14X,'Z',12X,'mass',/,
     +1X,74('-'))
99900 WRITE(6,*)' ERROR READING MASSES FROM UNIT 5'
      CALL errquit('rdinp error',911)
99901 WRITE(6,*)' ERROR READING TAPE 11 '
      CALL errquit('rdinp error',911)
      END
