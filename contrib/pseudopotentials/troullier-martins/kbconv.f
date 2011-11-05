      PROGRAM KBCONV
c
c ********************************************************************
c *                                                                  *
c *  Program converts a nonlocal or spin-orbit pseudopotential       *
c *  to the Kleinman and Bylander form and does the `Fourier'        *
c *  like Bessel transforms.  See Phys. Rev. Lett. v.48 no.20        *
c *  page 1425 as ref. The local potential indicator                 *
c *  are stored at the end of the                                    *
c *  pseudokb.dat data file.  Note that the projector is stored      *
c *  as Phi(r)deltaV(r), not as rPhi(r)deltaV(r).                    *
c *                                                                  *
c *  Written by Norm Troullier while at the U. of MN                 *
c *  Copyright Norm Troullier and Jose Luis Martins                  *
c *  Version 1.25 Dated Nov. 12, 1990                                *
c *                                                                  *
c ********************************************************************
c
c
      implicit double precision (a-h,o-z)
      PARAMETER (NRMAX=2000,LMAX=5)
c
      PARAMETER (ZERO=0.D0,SMALL=0.5D-5,ONE=1.D0,TWO=2.D0)
c
      DIMENSION R(NRMAX),RAB(NRMAX),NO(2*LMAX),LO(2*LMAX),SO(LMAX),
     1 VIOD(LMAX,NRMAX),VIOU(LMAX,NRMAX),VID(NRMAX),CDD(NRMAX),
     2 CDC(NRMAX),EV(LMAX),Y(NRMAX),YP(NRMAX),YPP(NRMAX),W(NRMAX*3),
     3 CDU(NRMAX),VIU(NRMAX),S1(NRMAX),S2(NRMAX),ETOT(10),EVL(2),
     4 AR(NRMAX,LMAX),ANORM(LMAX),INORM(LMAX),VQL(NRMAX),ARD(NRMAX),
     5 EVI(5),EVD(5)
c
      CHARACTER*1 ISPP
      CHARACTER*2 ICORR,NAMEAT
      CHARACTER*3 IREL
      CHARACTER*4 NICORE
      CHARACTER*10 IRAY(6),ITITLE(7),DATED
c
      DATA SO/5*ZERO/
      DATA EVI/5*ZERO/
      DATA EVD/5*ZERO/
c      CALL DROPFILE(0)
      DO 1 I=1,NRMAX
        CDU(I)=ZERO
 1    CONTINUE
c
c  Open the `kb.dat' input file - unit=2, the `kb.out'
C  output file - unit=3, and the `kbplot.dat' plotting
c  file - unit=4.
c  Read in the angular momentum quantum number
c  for the local potential part, the number of q-space
c  points and the seperation of the q-space points(a.u.).
c
      OPEN (UNIT=2,FILE='kb.dat',FORM='FORMATTED',
     1 STATUS='OLD')
      OPEN (UNIT=6,FILE='kb.out',FORM='FORMATTED',
     1 STATUS='NEW')
      OPEN (UNIT=4,FILE='kbplot.dat',FORM='FORMATTED',
     1 STATUS='NEW')
      WRITE(6,9301)
 9301 FORMAT(5X,'UNIT2 INPUT DATA',/)
      READ(2,2000)LLOCAL,IGHOST,COEF,ALP
      WRITE(6,9302)LLOCAL,IGHOST
 9302 FORMAT(5X,'LLOCAL=',I10,2X,'IGHOST=',I5,/)
      IF(IGHOST.GT.0)WRITE(6,9937)COEF,ALP
 9937 FORMAT(2X,'COEF=',F10.5,2X,'ALP=',F10.5)
      READ(2,2001)NQL,DELQL
      WRITE(6,9303)NQL,DELQL
 9303 FORMAT(2X,'NQL=',I10,2X,'DELQL=',F10.6)
 2000 FORMAT(5X,I2,I5,2F10.5)
 2001 FORMAT(10X,I4,10X,F10.8)
      WRITE(6,9305)
 9305 FORMAT(/,4X,'J',5X,'EVI(J)')
      DO 2 J=1,5
        READ(2,4001)EVI(J)
      WRITE(6,9304)J,EVI(J)
 9304 FORMAT(2X,I5,2X,E12.5)
 2    CONTINUE
 4001 FORMAT(5X,F10.8)
      CLOSE(UNIT=2)
c
c  Print out heading to `kb.out' file.
c
      CALL ZEDATE(DATED)
      WRITE(6,2002)DATED
 2002 FORMAT(2X,'Kleinman and Bylander pseudopotential conversion',
     1 /,'  and Bessel/Fourier transform program. Version 1.25',//,
     2 '  Ran on ',A10)
      WRITE(6,2003)LLOCAL
 2003 FORMAT(//,'  The l=',I2,' pseudopotentail is treated as the',
     1 /,'  local potential')
      WRITE(6,2004)NQL,DELQL
 2004 FORMAT(/,'  NQL=',I4,' and DELQL=',F10.8)
c
c  Open and read in data from file `pseudo.dat'.
c  Open and output some of the same data to file `pseudokb.dat'
c  for use in the atomkb program.  Open and output some of
c  the same data to `potfourkb.dat' for use in the planewave
c  program.
c
c ** On the Cray the RECL must be removed
c
      OPEN (UNIT=7,FILE='pseudo.dat',FORM='UNFORMATTED')
      READ(7) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1 (ITITLE(I),I=1,7),NORB,NUMU,NR,A,B,ZNUC
      NR=NR+1
      IZNUC=INT(ZNUC+0.1)
      READ(7) (R(I),I=2,NR)
      DO 5 I=1,NORB
        READ(7) LO(I),(VIOD(LO(I)+1,J),J=2,NR)
 5    CONTINUE
      DO 10 I=1,NUMU
        READ(7) LO(I+NORB),(VIOU(LO(I+NORB)+1,J),J=2,NR)
 10   CONTINUE
      READ(7) (CDC(I),I=2,NR)
      READ(7) (CDD(I),I=2,NR)
      CLOSE (UNIT=7)
c
c  Set up RAB integration grid
c
      R(1)=ZERO
      DO 15 I=1,NR
        RAB(I) = (R(I)+A)*B
 15   CONTINUE
c
c   Find the total amount of charge in the valence charge.
c   Fit cdd to splines and integrate.
c
      IF (NICORE .EQ. 'nc  ') THEN
        ITYPE=0
      ELSEIF(NICORE .EQ. 'pcec' .OR. NICORE .EQ. 'fcec') THEN
        ITYPE=1
      ELSE
        ITYPE=2
      ENDIF
      Y(1) = ZERO
      DO 20 I=2,NR
        Y(I) = CDD(I)
 20   CONTINUE
      ISX = 0
      A1 = ZERO
      AN = ZERO
      B1 = ZERO
      BN = ZERO
      NRM = NR
      CALL SPLIFT(R,Y,YP,YPP,NRM,W,IERR,ISX,A1,B1,AN,BN)
      IF ( IERR .NE. 1) THEN
        WRITE(6,2005)IERR
        STOP 'SPLIFT'
      ENDIF
 2005 FORMAT(1X,'****** Error in splift ierr =',I2)
      XLO = ZERO
      NUP = 1
      CALL SPLIQ(R,Y,YP,YPP,NRM,XLO,R(NR),NUP,TOTVEL,IERR)
      IF ( IERR .NE. 1 ) THEN
        WRITE(6,2006)IERR
        STOP 'SPLIQ'
      ENDIF
 2006 FORMAT(1X,'****** Error in spliq ierr =',I2)
c
c  Write out `pseudo.dat' info to `kb.out'.
c
      WRITE(6,2007)NAMEAT,ICORR,IREL,NICORE,TOTVEL,
     1 (IRAY(J),J=1,6),
     1 (ITITLE(J),J=1,7),A,B
 2007 FORMAT(/,A2,2X,A2,2X,A3,2X,A4,' pseudopotential with ',
     1 F10.6,' valence electrons',//,1X,6A10,4X,//,1X,7A10,//,
     2 ' a = ',F9.7,' b = ',F9.7,//)
c
c  Find el-el potential
c
      CALL VELECT(0,0,ICORR,' ',ITYPE,NR,R,RAB,TOTVEL,CDD,
     1 CDU,CDC,VID,VIU,ETOT,Y,YP,YPP,S1,S2,W)
c
c  Set up potentials
c
      C2 = -ONE/B**2
      C1 = -2*C2 + ONE/4
      DO 25 I=1,NORB
        LP=LO(I)+1
        NO(I)=LP
        LLP=LP*(LP-1)
c
c  Set up hamiltonian matrix for kinetic energy,
c  only the diagonal depends on the potential.
c
        Y(1)  = C1 / (R(2)+A)**2
        YP(1)  = ZERO
        YPP(1) = ZERO
        DO 30 K=3,NR
          Y(K-1)  = C1 / (R(K)+A)**2
          YP(K-1)  = C2 / ((R(K)+A)*(R(K-1)+A))
          YPP(K-1) = YP(K-1)**2
 30     CONTINUE
        DO 35 K=2,NR
          Y(K-1)=Y(K-1)+(VIOD(LP,K)+LLP/R(K))/R(K)+VID(K)
 35     CONTINUE
c
c  Diagonalize and find wave function and store in AR().
c
        EPS = -ONE
        CALL TRIDIB(NR-1,EPS,Y,YP,YPP,BL,BU,1,1,E,IND,IERR,
     1   S1,S2)
        EV(I)=E
        DO 40 J=2,NR
          Y(J)=(VIOD(LP,J)+LLP/R(J))/R(J)+VID(J)
 40     CONTINUE
        IFLAG=0
        CALL DIFNRL(1,I,Y,AR(1,LP),YP,LMAX,NR,A,B,R,RAB,NORB,NO,
     1   LO,SO,ZNUC,VIOD,VIOU,VID,VIU,EV,IFLAG,YPP,S1,S2,EVI)
 25   CONTINUE
c
c  Create deltaV(r) for nonlocal parts, and find local
c  part of potential.
c
      IF (LLOCAL .GE. 0 ) THEN
        DO 48 J=2,NR
          VQL(J) = VIOD(LLOCAL+1,J)
      IF(IGHOST.EQ.0)GO TO 48
       FAC=1.0-COEF*EXP(-ALP*R(J)*R(J))
      VQL(J)=VQL(J)*FAC
 48     CONTINUE
      ELSEIF(LLOCAL .EQ. -1) THEN
        DO 49 J=2,NR
          VQL(J) = (VIOD(1,J)/R(J)+VIOD(2,J)/R(J))/TWO*R(J)
 49     CONTINUE
      ENDIF
      DO 50 I=1,NORB
        LP = LO(I) + 1
        ANORM(LP) = ZERO
        A2NORM = ZERO
        DO 54 J=2,30
          VIOD(LP,J)=VIOD(LP,J)/R(J)-VQL(J)/R(J)
 54     CONTINUE
        DO 55 J=31,NR
          VIOD(LP,J)=(VIOD(LP,J)-VQL(J))/R(J)
 55     CONTINUE
c
c  Calculate the normalizing coeffcient for nonlocal parts,
c  uses Bode's rule for integration.
c
        AR(1,LP) = ZERO
        AR(NR+1,LP) = ZERO
        AR(NR+2,LP) = ZERO
        AR(NR+3,LP) = ZERO
        AR(NR+4,LP) = ZERO
        W(J+1) = ZERO
        W(J+2) = ZERO
        W(J+3) = ZERO
        W(J+4) = ZERO
        DO 56 J=1,NR
          W(J) = AR(J,LP)*AR(J,LP)*RAB(J)*VIOD(LP,J)
 56     CONTINUE
        DO 60 J=1,NR,4
          ANORM(LP)=ANORM(LP)+7*(W(J)+W(J+4))+
     1     32*(W(J+1)+W(J+3))+12*W(J+2)
 60     CONTINUE
        ANORM(LP)=2*ANORM(LP)/45
        IF ( ANORM(LP) .LT. ZERO) THEN
          INORM(LP) = -1
          ANORM(LP) = SQRT(-ANORM(LP))
        ELSEIF(ANORM(LP) .GT. ZERO) THEN
          INORM(LP) = 1
          ANORM(LP) = SQRT(ANORM(LP))
        ENDIF
        IF (LLOCAL .EQ. LO(I)) THEN
          INORM(LP) = 0
          ANORM(LP) = ONE
        ENDIF
        IF (INORM(LP) .NE. 0) THEN
          WMAX = ZERO
          DO 59 J=1,NR
            WMAX=MAX(WMAX,ABS(W(J)))
 59       CONTINUE
          DO 61 J=1,NR
            Y(J) = W(J)/WMAX
 61       CONTINUE
        ENDIF
        CALL PLOTKB(NR,Y,R,INORM(I),LO(I),'t')
        DO 57 J=1,NR
          W(J) = AR(J,LP)*AR(J,LP)*RAB(J)*ABS(VIOD(LP,J))
 57     CONTINUE
        DO 58 J=1,NR,4
          A2NORM=A2NORM+7*(W(J)+W(J+4))+32*(W(J+1)+W(J+3))+12*W(J+2)
 58     CONTINUE
        A2NORM=2*A2NORM/45
        IF ( A2NORM .NE. ZERO) THEN
          RATIO = ANORM(LP)*ANORM(LP)/A2NORM
          WRITE(6,2021)LP-1,INORM(LP)*RATIO
        ENDIF
 2021 FORMAT(1X,'Ratio for l=',I1,' is ',F10.7)
c
c  Calculate projector, NOTE stored as Phi(r)deltaV(r).
c
        IF (ANORM(LP) .EQ. ZERO) THEN
          DO 64 J=2,NR
            VIOD(LP,J)=ZERO
 64       CONTINUE
        ELSE
          DO 65 J=2,NR
            VIOD(LP,J)=(AR(J,LP)/R(J))*VIOD(LP,J)/ANORM(LP)
 65       CONTINUE
        ENDIF
 50   CONTINUE
c
      DO 67 J=2,NR
        VQL(J) = VQL(J)/R(J)
 67   CONTINUE
c
c   Find the 1st and 2nd eigenvalues for all angular
c   momentum using just the local potential.  Print
c   out results, this is used to find if any ghost
c   states exist for the potential.  See (preprint)paper
c   of Gonze, Kackell, and Scheffler, Phy. Rev. B.
c
      WRITE(6,301)
 301  FORMAT(//,'  Ghost State Existence Data',/
     1 ' l     0,node-eigen   1,node-eigen    True-eigen   inorm',/)
      DO 300 I=1,NORB
c
c  Set up the potential.
c
        IF (LLOCAL .EQ. LO(I)) GOTO 300
        EVT=EV(I)
        IF (EV(I) .GE. 0.0) EV(I)=-SMALL
        LP=LO(I)+1
        LLP = LP * (LP-1)
        STORE = VIOD(I,2)
        VIOD(I,2) = VQL(2)*R(2)
        DO 318 J=2,NR
          Y(J) = VQL(J) + LLP/R(J)**2 + VID(J)
 318    CONTINUE
c
c  Call the integration routine.
c
        NOT=NO(I)
        DO 345 JJ=0,1
          DO 317 J=1,NR
            ARD(J)=ZERO
 317      CONTINUE
          NO(I)=NOT+JJ
          IFLAG=0
          CALL DIFNRL(1,I,Y,ARD,YP,LMAX,NR,A,B,R,RAB,NORB,NO,
     1     LO,SO,ZNUC,VIOD,VIOU,VID,VIU,EV,IFLAG,YPP,S1,S2,EVD)
           EVL(JJ+1)=EV(I)
 345    CONTINUE
c
        VIOD(I,2) = STORE
        EV(I) = EVT
        NO(I) = NOT
c
        WRITE(6,355)LO(I),EVL(1),EVL(2),EV(I),INORM(I)
 355  FORMAT(1X,I1,5X,3(F12.9,4X),I2)
        IF (INORM(I) .LT. 0 ) THEN
          IF (EV(I) .GT. EVL(1)) THEN
            WRITE(6,302)
            IF (EVI(I) .NE. ZERO) WRITE(6,903)
          ENDIF
        ELSEIF (INORM(I) .GT. 0 ) THEN
          IF (EV(I) .GT. EVL(2)) THEN
            WRITE(6,302)
            IF (EVI(I) .NE. ZERO) WRITE(6,903)
          ENDIF
        ENDIF
 302  FORMAT(1X,' WARNING GHOST STATE WILL BE PRESENT!!!!!!')
 903  FORMAT(1X,' WARNING: GHOST STATE MAY BE DUE TO NON-EIGENVALUE')
 300  CONTINUE
c
c  Open and write heading to `pseudokb.dat' file.
c  Store data back into file `pseudokb.dat'.
c
c     call lin yang interface subroutine
c
      OPEN (UNIT=8,FILE='pseudokb.dat',FORM='UNFORMATTED')
      WRITE(8) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1 (ITITLE(I),I=1,7),NORB,NUMU,NR-1,A,B,ZNUC
      WRITE(8) (R(I),I=2,NR)
      DO 70 I=1,NORB
        WRITE(8) LO(I),(VIOD(I,J),J=2,NR)
 70   CONTINUE
      DO 75 I=1,NUMU
        WRITE(8) LO(I+NORB),(VIOU(LO(I+NORB)+1,J),J=2,NR)
 75   CONTINUE
      WRITE(8) (CDC(I),I=2,NR)
      WRITE(8) (CDD(I),I=2,NR)
      WRITE(8) (VQL(J),J=2,NR)
      WRITE(8) NORB
      DO 80 I=1,NORB
        WRITE(8)INORM(LO(I)+1),ANORM(LO(I)+1)
 80   CONTINUE
      CLOSE(UNIT=8)
c
c   Send operator data to plotkb for writing of data to
c   the `kbplot.dat' file.
c
      DO 85 I=1,NORB
        DO 90 J=2,NR
          Y(J) = VIOD(LO(I)+1,J)
 90     CONTINUE
       CALL PLOTKB(NR,Y,R,INORM(LO(I)+1),LO(I),'r')
 85   CONTINUE
c
c  Print out the local and K & B potentials
c
      WRITE(6,2008)NORB,NR
 2008 FORMAT(//,' Pseudo potentials in real space',//,
     1 ' number of potentials =',I2,/,
     2 ' number of radial grid points =',I4,//)
      IF (NORB .EQ. 2) THEN
        WRITE(6,2010)
      ELSEIF (NORB .EQ. 3) THEN
        WRITE(6,2011)
      ELSEIF (NORB .EQ. 4) THEN
        WRITE(6,2012)
      ELSE
        WRITE(6,3013)
      ENDIF
 2010 FORMAT('Index   r(I)    Local    S pot.   P pot.   Core c',
     1 '  Val c',/)
 2011 FORMAT('Index   r(I)    Local    S pot.   P pot.   D pot.',
     1 '  Core c  Val c',/)
 2012 FORMAT('Index   r(I)    Local    S pot.   P pot.   D pot.',
     1 '  F pot.  Core c  Val c',/)
 3013 FORMAT('Index   r(I)    Local    S pot.   P pot.   D pot.',
     1 '   F pot.  G pot.  Core c  Val c',/)
      DO 95 J=2,NR
        WRITE(6,2013) J,R(J),VQL(J),(VIOD(I,J),I=1,NORB),CDC(J),CDD(J)
 95   CONTINUE
      CALL KBLY(VQL,VIOD,CDC,CDD,R,INORM,NR, LLOCAL, NAMEAT)
 2013 FORMAT(1X,I3,1X,F8.5,6F9.4,2F7.4)
c
c  Fourier transform local potential, first for q=0 and
c  then for rest q>0.  Use Bode's rule to integrate.
c
      PI4 = 16*ATAN(ONE)
      VQL0 = ZERO
      W(1) =ZERO
      VT = 2*ZNUC
      LOC = LLOCAL + 1
      DO 100 K=2,NR
        W(K) = RAB(K)*R(K)*(R(K)*VQL(K)+VT)
 100  CONTINUE
      DO 105 K=NR+1,NR+4
        W(K) = ZERO
 105  CONTINUE
      DO 110 K=1,NR,4
        VQL0 = VQL0 + 7*W(K) + 32*W(K+1) + 12*W(K+2) +
     1   32*W(K+3) + 7*W(K+4)
 110  CONTINUE
      VQL0 = 2 * VQL0 * PI4 / 45
c
      Y(1) = ZERO
      Y(NR+1) = ZERO
      Y(NR+2) = ZERO
      Y(NR+3) = ZERO
      Y(NR+4) = ZERO
      DO 111 K=2,NR
        W(K) = W(K)/R(K)
 111  CONTINUE
      YP(1) = ZERO
      DO 112 J=2,NQL+1
        YP(J) = DELQL * (J-1)
 112  CONTINUE
      DO 115 J=1,NQL
        VQL(J) = ZERO
        DO 120 K=2,NR
          Y(K) = SIN(YP(J+1)*R(K))*W(K)
 120    CONTINUE
        DO 125 K=1,NR,4
          VQL(J) = VQL(J) + 7*Y(K) + 32*Y(K+1) + 12*Y(K+2) +
     1     32*Y(K+3) + 7*Y(K+4)
 125    CONTINUE
        VQL(J) = PI4 * (2 * VQL(J) / 45 - VT/YP(J+1))/YP(J+1)
 115  CONTINUE
c
c  Do loop over K&B projector potentials
c
      DO 130 I=1,NORB
        DO 131 K=NR,2,-1
          IF (VIOD(I,K) .EQ. ZERO) THEN
            NRM=K
          ELSE
            GOTO 134
          ENDIF
 131    CONTINUE
 134    DO 135 K=2,NRM
          W(K) = RAB(K)*R(K)*R(K)*VIOD(I,K)
 135    CONTINUE
        Y(1) = ZERO
        DO 136 K=NRM+1,NRM+7
          Y(K) = ZERO
 136    CONTINUE
        CRNORM = 7*SQRT((2*LO(I)+1)*PI4)/17280
        DO 145 J=1,NQL+1
          DO 150 K=2,NRM
            Y(K) = SBESSJ(I-1,YP(J)*R(K))*W(K)
 150      CONTINUE
c
c  Due to the high number of occilations in the intagrand,
c  an eight point Newton-Cotes intagration method is used.
c  See  Abramowitz and Stegun Eq. 25.4.17
c
          AR(J,I) = ZERO
          DO 155 K=1,NRM,7
            AR(J,I) = AR(J,I)+751*(Y(K)+Y(K+7))+3577*(Y(K+1)+Y(K+6))+
     1       1323*(Y(K+2)+Y(K+5))+2989*(Y(K+3)+Y(K+4))
 155      CONTINUE
          AR(J,I) = CRNORM * AR(J,I)
 145    CONTINUE
 130  CONTINUE
c
c  Fourtier transform core charge density, ignore if
c  itype is equal to 0.
c
      IF (ITYPE .EQ. 0) THEN
        DO 160 J=2,NQL+1
          CDC(J) = ZERO
 160    CONTINUE
      ELSE
        DO 165 J=2,NR
          W(J) = RAB(J)*CDC(J)/R(J)
 165    CONTINUE
        DO 170 K=2,NQL+1
          CDC(K) = ZERO
          DO 175 J=1,NR
            Y(J) = W(J)*SIN(YP(K)*R(J))
 175      CONTINUE
          DO 180 J=1,NR,4
            CDC(K) = CDC(K) + 7*Y(J) + 32*Y(J+1) + 12*Y(J+2) +
     1       32*Y(J+3) + 7*Y(J+4)
 180      CONTINUE
          CDC(K) = 2*CDC(K)/45/YP(K)
 170    CONTINUE
      ENDIF
c
c  Fourier transform the valence charge density.
c
      DO 185 J=2,NR
        W(J) = RAB(J)*CDD(J)/R(J)
 185  CONTINUE
      DO 190 K=2,NQL+1
        CDD(K) = ZERO
        DO 195 J=1,NR
          Y(J) = W(J)*SIN(YP(K)*R(J))
 195    CONTINUE
        DO 200 J=1,NR,4
          CDD(K) = CDD(K) + 7*Y(J) + 32*Y(J+1) + 12*Y(J+2) +
     1     32*Y(J+3) + 7*Y(J+4)
 200    CONTINUE
        CDD(K) = 2*CDD(K)/45/YP(K)
 190  CONTINUE
c
c  Send K&B potentials to plotkb
c
      DO 205 J=1,NORB
        IF (J .NE. LOC) THEN
          CALL PLOTKB(NQL,AR(1,J),YP,INORM(J),LO(J),'q')
        ENDIF
 205  CONTINUE
      CLOSE(UNIT=4)
c
c  Open and write heading to `potfourkb.dat' file.
c  Write the transforms to file `potfourkb.dat'
c
      OPEN (UNIT=9,FILE='potfourkb.dat',FORM='UNFORMATTED')
      WRITE(9) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1 (ITITLE(I),I=1,7)
      WRITE(9) IZNUC,NQL,DELQL,VQL0
      WRITE(9) NORB,(INORM(I),I=1,NORB)
      WRITE(9) (VQL(J),J=1,NQL)
      DO 210 I=1,NORB
        WRITE(9) (AR(J,I),J=1,NQL+1)
 210  CONTINUE
      WRITE(9) (CDC(J),J=2,NQL+1)
      WRITE(9) (CDD(J),J=2,NQL+1)
      CLOSE(UNIT=9)
c
c  Print out info to `kb.out' file.
c
      WRITE(6,2014)NQL,VQL0
 2014 FORMAT(//,' Pseudo potentials in fourier space',//,
     1 ' Number of q-space points = ',I4,/,
     2 ' V(q=0) = ',f10.5,/)
      WRITE(6,2020)NORB,(INORM(I),I=1,NORB)
 2020 FORMAT(1X,'Number of potentials = ',I1,/
     1 ' Normalization indexs ',5(I2,5X))
      IF (NORB .EQ. 2) THEN
        WRITE(6,2016)
      ELSEIF (NORB .EQ. 3) THEN
        WRITE(6,2017)
      ELSEIF (NORB .EQ. 4) THEN
        WRITE(6,2018)
      ELSE
        WRITE(6,3018)
      ENDIF
 2016 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   Core c ',
     1 '   Val c',/)
 2017 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   D pot. ',
     1 '   Core c   Val c',/)
 2018 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   D pot. ',
     1 '   F pot.   Core c   Val c',/)
 3018 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   D pot. ',
     1 '    F pot.   G pot.   Core c   Val c',/)
      WRITE(6,2030)(AR(1,I),I=1,NORB)
 2030 FORMAT(2X,'0',3X,'0.000',14X,5F9.5)
      DO 215 J=1,NQL
        WRITE(6,2019)J,YP(J+1),VQL(J),(AR(J+1,I),I=1,NORB),
     1   CDC(J+1),CDD(J+1)
 215  CONTINUE
 2019 FORMAT(I4,F7.3,F14.5,7F9.5)
      CLOSE(UNIT=6)
c
      CALL EXIT
      END

      subroutine plotkb(nr,vd,r,inorm,lo,rq)
       implicit double precision (a-h, o-z)
c
c ***********************************************************
c *                                                         *
c *  This is a plotting routine for kbconv.  Prints         *
c *  out the K & B operator or transform.                   *
c *                                                         *
c ***********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter statement
c  ###    to single precision statement.
c  ###  Cray conversions
c  njtj
c
      parameter (zero=0.D0,pzf=0.05D0,small=0.0001D0)
Cray      parameter (zero=0.0,pzf=0.05)
c
c
      character*1 rq
      character*3 marker
c
      dimension vd(nr),r(nr)
c
c  Step size of 0.05 is adjustable as seen fit to give
c  a reasonalble plot.
c
      if (inorm .eq. 0) then
        return
      elseif (rq .eq. 'r' .or. rq .eq. 't') then
        nrm = nr - 120
        nst = 2
        step=small
      elseif (rq .eq. 'q') then
        nrm = nr
        nst = 1
        step=zero
      endif
      do 150,j=nst,nrm
        if (r(j) .ge. step) then
          write(4,6000)r(j),vd(j)
          step=step+pzf
        endif
 150  continue
      if (rq .eq. 'r') then
        if (lo .eq. 0) then
          marker='psr'
        elseif (lo .eq. 1) then
          marker='ppr'
        elseif (lo .eq. 2) then
          marker='pdr'
        elseif (lo .eq. 3) then
          marker='pfr'
        elseif (lo .eq. 4) then
          marker='pgr'
        endif
      elseif (rq .eq. 'q') then
        if (lo .eq. 0) then
          marker='psq'
        elseif (lo .eq. 1) then
          marker='ppq'
        elseif (lo .eq. 2) then
          marker='pdq'
        elseif (lo .eq. 3) then
          marker='pfq'
        elseif (lo .eq. 4) then
          marker='pgq'
        endif
      else
        if (lo .eq. 0) then
          marker='pst'
        elseif (lo .eq. 1) then
          marker='ppt'
        elseif (lo .eq. 2) then
          marker='pdt'
        elseif (lo .eq. 3) then
          marker='pft'
        elseif (lo .eq. 4) then
          marker='pgt'
        endif
      endif
      write(4,6001)marker
      return
c
c  Format statements
c
 6000 format(1x,f7.4,3x,f10.6)
 6001 format(1x,'marker ',a3)
      end
C       SUBROUTINE ZEDATE(BDATE)
C       implicit double precision (a-h, o-z)
CC
CC    GETS THE DATE (DAY-MONTH-YEAR)
CC    CRAY-2 VERSION
CCCC
C       CHARACTER*10 BDATE
C       CHARACTER*8 ADATE
C       CHARACTER*3 MONTH(12)
C       CHARACTER*1 DASH,DUM1,DUM2
C       DATA DASH/'-'/
C       DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
C     2  'OCT','NOV','DEC'/
CC
C       WRITE(ADATE,100) DATE()
C       READ(ADATE,101) LMONTH,DUM1,LDAY,DUM2,LYEAR
C       WRITE(BDATE,102) LDAY,DASH,MONTH(LMONTH),DASH,LYEAR
C  100  FORMAT(A8)
C  101  FORMAT(I2,A1,I2,A1,I2)
C  102  FORMAT(I2,A1,A3,A1,I2,' ')
C       RETURN
C       END
CC
CC  *****************Cray end***********************
       SUBROUTINE ZEDATE(BDATE)
C
C   GETS THE DATE (DAY-MONTH-YEAR)
C   Sun version
C
       CHARACTER*1 BDATE(10)
       CHARACTER*1 LOCTIM(24)
       CALL FDATE(LOCTIM)
       DO 101 I = 11, 20
       II = I - 10
       BDATE(II) = LOCTIM(I)
101    CONTINUE
       RETURN
       END
      subroutine difnrl(iter,iorb,v,ar,br,lmax,
     1 nr,a,b,r,rab,norb,no,lo,so,znuc,viod,viou,
     2 vid,viu,ev,iflag,rab2,fa,fb,evi)
      implicit double precision (a-h,o-z)
c
c    difnrl integrates the Schroedinger equation
c    if finds the eigenvalue ev, the wavefunction ar
c    and the derivative br = d(ar)/dr
c
c  njtj  ***  modifications  ***
c    This routine has had major modifications.  Some
c    of the data used inside the main loop has been
c    calculated outside the main loop to reduce the number
c    of operations(uses extra array space to gain speed)
c    and are passed as work arrays form the main.
c    The predictor-corrector functions have been put
c    into a array.
c    The iflag variable was added to indicate nonconvergence
c    for other programs.  It has no use in the atom program
c    and can be removed by the user.
c    All output from the routine is compatible to
c    the Berkeley/Sverre Froyen version.
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch the double precision parameter statements
c  ###      to single precision parameter statements.
c  ###  Cray conversions
c  njtj
c
c  njtj
c  &&&  Machine dependent Parameter
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for
c  &&&    the machine in use from the list, or find
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c
       parameter(zero=0.0,pnine=0.9,two=2.0,etol=-1.E-7)
c
c  Tolerence
c
       parameter(tol=1.E-10,five=5.0)
c
c  Integration coefficients
c
       parameter(abc1=190.1/72,abc2=-138.7/36,abc3=10.9/3,
     1 abc4=-63.7/36,abc5=25.1/72,amc0=25.1/72,amc1=32.3/36,
     2 amc2=-1.1/3,amc3=5.3/36,amc4=-1.9/72)
c
c
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr),no(norb),
     1 lo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),evi(norb)
c
c  njtj  *** start modification  ***
c    Arrays added to gain speed.
c
      dimension rabrlo(5),rlp(5),rab2(nr),fa(nr),fb(nr)
c
c  njtj  ***  end modification  ***
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
cApollo      expzer = 3.7D2
C SUN:
      expzer = 3.7D2
cVax      expzer = 44.D0
C       expzer =  2.8E3
c
c  njtj  *** major modification start  ***
c    Loop data calculated outside loop to gain speed.
c
      itmax = 100
      iflag = 0
      lp = lo(iorb)+1
      ar(1) = zero
      if (lo(iorb) .eq. 0) then
        br(1) = b*a
      else
        br(1) = zero
      endif
      do 1 j=2,nr
        ar(j) = zero
 1    continue
      do 2 j=2,nr
        br(j) =zero
 2    continue
      do 4 j=2,5
        rlp(j)=r(j)**lp
 4    continue
      do 5 j=2,5
        rabrlo(j)=rab(j)*r(j)**lo(iorb)
 5    continue
      do 6 j=1,nr
        rab2(j)=rab(j)*rab(j)
 6    continue
c
c   set underflow trap, error from Berkeley version,
c   fixed by Troy Barbee sqrt(expzer) should be expzer/2
c   4/17/90
c
      juflow=1
      do 42 j=2,nr
        if (lp*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue
c
c  njtj  *** end major modification  ***
c
c   determine effective charge and vzero for startup of
c   outward integration
c   ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
c   aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
c
      zeff = zero
      if (so(iorb) .lt. 0.1 .and. viod(lp,2) .lt. -0.1) zeff=znuc
      if (so(iorb) .gt. 0.1 .and. viou(lp,2) .lt. -0.1) zeff=znuc
      aa = -zeff/lp
      vzero = -2*zeff*aa
      if (zeff .eq. zero) then
        if (so(iorb) .lt. 0.1 ) then
          vzero=vzero+viod(lp,2)/r(2)
        else
          vzero=vzero+viou(lp,2)/r(2)
        endif
      endif
      if (so(iorb) .lt. 0.1) then
        vzero=vzero+vid(2)
      else
        vzero=vzero+viu(2)
      endif
      var0 = zero
      if (lo(iorb) .eq. 0) var0=-2*zeff
      if (lo(iorb) .eq. 1) var0=two
c
      emax = zero
      emin = -two*100000
      if (ev(iorb) .gt. emax) ev(iorb) = emax
 10   if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15   format(' iorb =',i3,' iter =',i3,' ev =',e18.10,' nodes =',i2)
      if (itmax .eq. 0) then
        iflag =1
        return
      endif
      if (ev(iorb) .gt. zero) then
        write(6,1000)iorb
        call ext(620+iorb)
      endif
 1000 format(//,' error in difnrl - ev(',i2,
     1 ') greater then v(infinty)')
c
c   find practical infinity ninf and classical turning
c   point nctp for orbital
c
      icount=0
 20   icount=icount+1
      do 22 j=nr,2,-1
        temp = v(j) -ev(iorb)
        if (temp .lt. zero) temp = zero
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
 22   continue
 23   ninf=j
      nctp = ninf - 5
      do 25 j=2,ninf-5
        if (v(j) .lt. ev(iorb)) nctp = j
 25   continue
      if (ev(iorb) .ge. etol*10) nctp=ninf-5
      if (ev(iorb) .ge. etol) ev(iorb)=zero
      if (evi(iorb) .ne. zero) then
        ev(iorb) = evi(iorb)
        do 26 j=1,nr
          if (r(j) .lt. five) nctp=j
 26     continue
      endif
      if (nctp .le. 6) then
        ev(iorb) = pnine*ev(iorb)
        if (icount .gt. 200) then
          write(6,1010)iorb
          call ext(650+iorb)
        endif
        goto 20
      endif
 1010 format(//,'error in difnrl - cannot find the classical '
     1 ,/' turning point for orbital ',i2)
c
c   outward integration from 1 to nctp
c   startup
c
      bb = (vzero-ev(iorb))/(4*lp+2)
      do 35 j=2,5
        ar(j) = rlp(j) * (1+(aa+bb*r(j))*r(j))
        br(j) = rabrlo(j) * (lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
 35   continue
c
c  njtj  ***  start major modification  ***
c    Predictor-corrector array added.
c
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-ev(iorb))*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-ev(iorb))*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-ev(iorb))*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-ev(iorb))*ar(5)
c
c   intergration loop
c
      nodes = 0
      do 40 j=6,nctp
c
c   predictor (Adams-Bashforth)
c
        j1=j-1
        j2=j-2
        j3=j-3
        j4=j-4
        j5=j-5
        vev=v(j)-ev(iorb)
        arp = ar(j1) + abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5)
        brp = br(j1) + abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5)
        fb1 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) + amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4)
        brc = br(j1) + amc0*fb1+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4)
        fb0 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc + amc0*(brc-brp)
        br(j) = brc + amc0*(fb0-fb1)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
c
c   count nodes - if no underflow
c
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.zero)nodes=nodes+1
 40   continue
c
c  njtj  ***  end major modification  ***
c
      arctp = ar(nctp)
      brctp = br(nctp)
c
c   end outward integration
c
c   if number of nodes correct, start inward integration
c   else modify energy stepwise and try again
c
      if (evi(iorb) .ne. zero) goto 111
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
        if (nodes .lt. no(iorb)-lo(iorb)-1) then
c
c  too few nodes; increase ev
c
          if (ev(iorb) .gt. emin) emin = ev(iorb)
          ev(iorb) = ev(iorb) - ev(iorb)/10
        else
c
c  too many nodes; decrease ev
c
          if (ev(iorb) .lt. emax) emax = ev(iorb)
          ev(iorb) = ev(iorb) + ev(iorb)/10
        endif
        itmax = itmax-1
        goto 10
      endif
c
c   inward integration from ninf to nctp
c   startup
c
      do 71 j=ninf,ninf-4,-1
        alf = v(j) - ev(iorb)
        if (alf .lt. zero) alf = zero
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = -rab(j)*alf*ar(j)
 71   continue
c
c  njtj  ***  start major modification  ***
c    Array for predictor-corrector added.
c
      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*
     1 (v(ninf)-ev(iorb))*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*
     1       (v(ninf1)-ev(iorb))*ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*
     1       (v(ninf2)-ev(iorb))*ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*
     1       (v(ninf3)-ev(iorb))*ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*
     1       (v(ninf4)-ev(iorb))*ar(ninf4)
c
c   integration loop
c
      istop = ninf - nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
c
c   predictor (Adams-Bashforth)
c
        j1 = j + 1
        j2 = j + 2
        j3 = j + 3
        j4 = j + 4
        j5 = j + 5
        vev = v(j)-ev(iorb)
        arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5))
        brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5))
        fb0 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4))
        brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4))
c
        fb1 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc - amc0*(brc-brp)
        br(j) = brc - amc0*(fb1-fb0)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
 80   continue
c
c   end inward integration
c
c  njtj  *** end major modification  ***
c
c   rescale ar and br outside nctp to match ar(nctp) from
c   outward integration
c
  222 factor = arctp/ar(nctp)
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
c
c   find normalizing factor
c
      factor = zero
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*ar(j)*ar(j)*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
c
c   modify eigenvalue ev
c
      dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev(iorb)) dev=dsign(ev(iorb),dev)/5
      itmax = itmax-1
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) / 2
      if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) / 2
      if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c   normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
c
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j) / rab(j)
 110  continue
 111  continue
      if (evi(iorb) .ne. zero) then
        factor = zero
        ll = 4
        do 112 j=2,nctp
          factor = factor + ll*ar(j)*ar(j)*rab(j)
          ll = 6 - ll
 112    continue
        factor = factor / 3
        factor = 1 / sqrt(factor)
        do 113 j=1,nctp
          ar(j) = factor*ar(j)
          br(j) = factor*br(j) / rab(j)
 113    continue
      endif
      return
      end
C
C
      subroutine velect(iter,iconv,icorr,ispp,ifcore,
     1 nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,y,yp,
     2 ypp,s1,s2,w)
c
c    velect generates the electronic output potential from
c    the electron charge density.  The ionic part is
c    added in dsolv1/dsolv2.
c
c  njtj  ***  modifications  ***
c    The only major modiication is that the constants for the
c    ceperly-alder 'ca' method are placed in parameter
c    statements, this was done so non-opt compiliers
c    would minimize the number of calculations.
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter statements
c  ###      to single precision parameter statements.
c  ###  Cray conversions
c  njtj
c
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 icorr
c
c  njtj  *** modification start  ***
c
       parameter (zero=0.D0,one=1.D0,pfive=.5D0,opf=1.5D0,pnn=.99D0)
       parameter (pthree=0.3D0,psevf=0.75D0,c0504=0.0504D0)
       parameter (c0254=0.0254D0,c014=0.014D0,c0406=0.0406D0)
       parameter (c15p9=15.9D0,c0666=0.0666D0,c11p4=11.4D0)
       parameter (c045=0.045D0,c7p8=7.8D0,c88=0.88D0,c20p592=20.592D0)
       parameter (c3p52=3.52D0,c0311=0.0311D0,c0014=0.0014D0)
       parameter (c0538=0.0538D0,c0096=0.0096D0,c096=0.096D0)
       parameter (c0622=0.0622D0,c004=0.004D0,c0232=0.0232D0)
       parameter (c1686=0.1686D0,c1p3981=1.3981D0,c2611=0.2611D0)
       parameter (c2846=0.2846D0,c1p0529=1.0529D0,c3334=0.3334D0)
Cray       parameter (zero=0.0,one=1.0,pfive=0.5,opf=1.5,pnn=0.99)
Cray       parameter (pthree=0.3,psevf=0.75,c0504=0.0504)
Cray       parameter (c0254=0.0254,c014=0.014,c0406=0.0406)
Cray       parameter (c15p9=15.9,c0666=0.0666,c11p4=11.4)
Cray       parameter (c045=0.045,c7p8=7.8,c88=0.88,c20p592=20.592)
Cray       parameter (c3p52=3.52,c0311=0.0311,c0014=0.0014)
Cray       parameter (c0538=0.0538,c0096=0.0096,c096=0.096)
Cray       parameter (c0622=0.0622,c004=0.004,c0232=0.0232)
Cray       parameter (c1686=0.1686,c1p3981=1.3981,c2611=0.2611)
Cray       parameter (c2846=0.2846,c1p0529=1.0529,c3334=0.3334)
c
c    Ceperly-Alder 'ca' constants
c
       parameter (con1=1.D0/6, con2=0.008D0/3, con3=0.3502D0/3)
       parameter (con4=0.0504D0/3, con5=0.0028D0/3, con6=0.1925D0/3)
       parameter (con7=0.0206D0/3, con8=9.7867D0/6, con9=1.0444D0/3)
       parameter (con10=7.3703D0/6, con11=1.3336D0/3)
Cray       parameter (con1=1.0/6, con2=0.008/3, con3=0.3502/3)
Cray       parameter (con4=0.0504/3, con5=0.0028/3, con6=0.1925/3)
Cray       parameter (con7=0.0206/3, con8=9.7867/6, con9=1.0444/3)
Cray       parameter (con10=7.3703/6, con11=1.3336/3)
c
c  njtj  ***  modification end  ***
c
      dimension r(nr),rab(nr),cdd(nr),cdu(nr),cdc(nr),
     1 vod(nr),vou(nr),etot(10),y(nr),yp(nr),ypp(nr),
     2 s1(nr),s2(nr),w(3*nr)
c
       pi=4*atan(one)
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
c      fit cd/r by splines
c
       y(1) = zero
       do 10 i=2,nr
         y(i) = (cdd(i)+cdu(i))/r(i)
 10    continue
       if (ifcore .eq. 2) then
         do 11 i=2,nr
           y(i) = y(i) + cdc(i)/r(i)
 11      continue
       endif
       isx = 0
       a1 = zero
       an = zero
       b1 = zero
       bn = zero
       nrm=nr
       call splift(r,y,yp,ypp,nrm,w,ierr,isx,a1,b1,an,bn)
       if(ierr.ne.1) then
         write(6,20000)ierr
         call ext(420+ierr)
       endif
20000  format(1x,'****** Error in splift ierr =',i2)
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
       xlo = zero
       call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s2,ierr)
       if(ierr.ne.1) then
         write(6,20001)ierr
         call ext(440+ierr)
       endif
20001  format(1x,'****** Error in spliq ierr =',i2)
       do 20 i=1,nr
         ypp(i) = r(i)*ypp(i) + 2*yp(i)
         yp(i)  = r(i)*yp(i)  + y(i)
         y(i)   = r(i)*y(i)
 20    continue
       call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s1,ierr)
       if(ierr.ne.1) then
         write(6,20002)ierr
         call ext(460+ierr)
       endif
20002  format(1x,'****** Error in spliq ierr =',i2)
c
c      check normalization
c
       xnorm = zero
       if (ifcore .eq. 2 .and. iter .eq. 0 ) zel=s1(nr)
       if (zel .ne. zero) xnorm = zel/s1(nr)
       if (iter .gt. 3 .and. abs(zel-s1(nr)) .gt. 0.01) then
         if (zel .lt. s1(nr)+1.0 ) then
           write(6,24) iter,xnorm
 24    format(/,' warning *** charge density rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         else
           xnorm=pnn*xnorm
           write(6,25) iter,xnorm
 25    format(/,' warning *** charge density partially rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         endif
       endif
c
c      compute new hartree potential
c      renormalize the charge density
c
       do 30 i=2,nr
         vod(i) = 2 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
         vou(i) = vod(i)
         cdd(i) = xnorm*cdd(i)
         cdu(i) = xnorm*cdu(i)
 30    continue
c
c      compute hartree contribution to total energy
c
       if (iconv .eq. 1) then
         ehart = zero
         ll = 4
         do 40 i=2,nr
           ehart = ehart+ll*(cdd(i)+cdu(i))*vod(i)*rab(i)
           ll = 6 - ll
 40      continue
         ehart = ehart / 6
       endif
c
c      add exchange and correlation
c
       trd = one/3
       ftrd = 4*trd
       tftm = 2**ftrd-2
       a0 = (4/(9*pi))**trd
c
c      set x-alpha
c
       alp = one
       if (icorr .ne. 'xa') alp = 2 * trd
       vxc = zero
       vc  = zero
       exc = zero
       ec  = zero
c
c      start loop
c
       ll = 4
       do 210 i=2,nr
         cdsum = cdd(i) + cdu(i)
         if (ifcore .ge. 1) cdsum=cdsum+cdc(i)
         if (cdsum .le. zero) goto 210
c
c  Vax bug fix.  Troy Barbee - 4/17/90
c
         if (log(3*r(i)**2/cdsum) .gt. 2*expzer) goto 210
         rs = (3*r(i)**2/cdsum)**trd
         z = zero
         fz = zero
         fzp = zero
         if (ispp .eq. 's') then
           z = (cdd(i)-cdu(i)) / cdsum
           fz = ((1+z)**ftrd+(1-z)**ftrd-2)/tftm
           fzp = ftrd*((1+z)**trd-(1-z)**trd)/tftm
         endif
c
c      exchange (only use (xa))
c
         vxp = -3*alp/(pi*a0*rs)
         exp = 3*vxp/4
         if (ispp .eq. 'r') then
           beta = c014/rs
           sb = sqrt(1+beta*beta)
           alb = log(beta+sb)
           vxp = vxp * (-pfive + opf * alb / (beta*sb))
           exp = exp *(one-opf*((beta*sb-alb)/beta**2)**2)
         endif
 65      vxf = 2**trd*vxp
         exf = 2**trd*exp
         vcp = zero
         ecp = zero
         vcf = zero
         ecf = zero
         if (icorr .eq. 'ca') then
c          ceperly-alder (ca)
c          The Perdew-Zunger parameterization is used.
c          See Phys. Rev. B 23 5075 (1981).
           if (rs .gt. one) then
             sqrs=sqrt(rs)
             te = one+con10*sqrs+con11*rs
             be = one+c1p0529*sqrs+c3334*rs
             ecp = -c2846/be
             vcp = ecp*te/be
             te = one+con8*sqrs+con9*rs
             be = one+c1p3981*sqrs+c2611*rs
             ecf = -c1686/be
             vcf = ecf*te/be
           else
             rslog=log(rs)
             ecp=(c0622+c004*rs)*rslog-c096-c0232*rs
             vcp=(c0622+con2*rs)*rslog-con3-con4*rs
             ecf=(c0311+c0014*rs)*rslog-c0538-c0096*rs
             vcf=(c0311+con5*rs)*rslog-con6-con7*rs
           endif
         elseif (icorr .eq. 'xa') then
c          correlation
         elseif (icorr .eq. 'wi') then
c          wigner (wi)
           vcp = -(c3p52*rs+c20p592)/(3*(rs+c7p8)**2)
           ecp = -c88/(rs+c7p8)
         elseif (icorr .eq. 'hl') then
c          hedin-lundqvist (hl)
           x = rs/21
           aln = log(1+1/x)
           vcp = -c045*aln
           ecp = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecp=((con1/x-pthree)/x+psevf)/x
           ecp = -c045*ecp
         elseif (icorr .eq. 'gl') then
c          gunnarson-lundqvist-wilkins (gl)
           x = rs/c11p4
           aln = log(1+1/x)
           vcp = -c0666*aln
           ecp = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecp=((con1/x-pthree)/x+psevf)/x
           ecp = -c0666*ecp
           x = rs/c15p9
           aln = log(1+1/x)
           vcf = -c0406*aln
           ecf = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecf=((con1/x-pthree)/x+psevf)/x
           ecf = -c0406*ecf
         elseif (icorr .eq. 'bh') then
c          von barth - hedin (bh)
           x = rs/30
           aln = log(1+1/x)
           vcp = -c0504*aln
           ecp = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecp=((con1/x-pthree)/x+psevf)/x
           ecp = -c0504*ecp
           x = rs/75
           aln = log(1+1/x)
           vcf = -c0254*aln
           ecf = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecf=((con1/x-pthree)/x+psevf)/x
           ecf = -c0254*ecf
         else
           write(6,70) icorr
           call ext(400)
         endif
 70   format('error in velect - icorr =',a2,' not implemented')
         vxcp = vxp + vcp
         vxcf = vxf + vcf
         vxcd = vxcp
         vxcu = vxcp
         excp = exp + ecp
         excf = exf + ecf
         vcd = vcp
         vcu = vcp
         exct = excp
         ect = ecp
         if (z .ne. zero) then
           vxcd = vxcd + fz*(vxcf-vxcp) + (1-z)*fzp*(excf-excp)
           vxcu = vxcu + fz*(vxcf-vxcp) - (1+z)*fzp*(excf-excp)
           vcd = vcd + fz*(vcf-vcp) + (1-z)*fzp*(ecf-ecp)
           vcu = vcu + fz*(vcf-vcp) - (1+z)*fzp*(ecf-ecp)
           exct = exct + fz*(excf-excp)
           ect = ect + fz*(ecf-ecp)
         endif
         vod(i) = vod(i) + vxcd
         vou(i) = vou(i) + vxcu
         vxc = vxc + ll * (cdd(i)*vxcd + cdu(i)*vxcu) * rab(i)
         vc  = vc  + ll * (cdd(i)*vcd  + cdu(i)*vcu ) * rab(i)
         exc = exc + ll * cdsum * exct * rab(i)
         ec  = ec  + ll * cdsum * ect  * rab(i)
         ll = 6 - ll
 210   continue
       etot(4) = ehart
       etot(5) = vxc / 3
       etot(6) = (3*vc - 4*ec) / 3
       etot(7) = exc / 3
       vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
       vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))
       return
       end
C
C
       SUBROUTINE SPLIQ(X,Y,YP,YPP,N,XLO,XUP,NUP,ANS,IERR)
       implicit double precision (a-h, o-z)
C
C
C  NJTJ
C  ###  CRAY CONVERSIONS
C  ###    1)Comment out implicit double precision.
C  ###  CRAY CONVERSIONS
C  NJTJ
C
       DIMENSION X(N),Y(N),YP(N),YPP(N),XUP(NUP),ANS(NUP)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     THIS ROUTINE WAS WRITTEN BY M. K. GORDON
C
C     ABSTRACT
C
C     SUBROUTINE SPLIQ INTEGRATES A CUBIC SPLINE (GENERATED BY
C     SPLIFT, SMOO, ETC.) ON THE INTERVALS (XLO,XUP(I)), WHERE XUP
C     IS A SEQUENCE OF UPPER LIMITS ON THE INTERVALS OF INTEGRATION.
C     THE ONLY RESTRICTIONS ON XLO AND XUP(*) ARE
C                XLO .LT. XUP(1),
C                XUP(I) .LE. XUP(I+1)   FOR EACH I .
C     ENDPOINTS BEYOND THE SPAN OF ABSCISSAS ARE ALLOWED.
C     THE SPLINE OVER THE INTERVAL (X(I),X(I+1)) IS REGARDED
C     AS A CUBIC POLYNOMIAL EXPANDED ABOUT X(I) AND IS INTEGRATED
C     ANALYTICALLY.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YP(N), YPP(N), XUP(NUP), ANS(NUP)
C
C      --INPUT--
C
C        X    - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C        Y    - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C        YP   - ARRAY OF FIRST DERIVATIVES OF THE SPLINE AT ABSCISSAS.
C               USUALLY YP IS THE SAME AS YP IN SPLIFT OR R1 IN SMOO.
C        YPP  - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C        N    - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C        XLO  - LEFT ENDPOINT OF INTEGRATION INTERVALS.
C        XUP  - RIGHT ENDPOINT OR ARRAY OF RIGHT ENDPOINTS OF
C               INTEGRATION INTERVALS IN ASCENDING ORDER.
C        NUP  - THE NUMBER OF RIGHT ENDPOINTS.  IF NUP IS GREATER THAN
C               1, THEN XUP AND ANS MUST BE DIMENSIONED AT LEAST NUP.
C
C      --OUTPUT--
C
C        ANS -- ARRAY OF INTEGRAL VALUES, THAT IS,
C               ANS(I) = INTEGRAL FROM XLO TO XUP(I)
C        IERR -- ERROR STATUS
C                = 1 INTEGRATION SUCCESSFUL
C                = 2 IMPROPER INPUT - N.LT.4 OR NUP.LT.1
C                = 3 IMPROPER INPUT - ABSCISSAS NOT IN
C                        STRICTLY ASCENDING ORDER
C                = 4 IMPROPER INPUT - RIGHT ENDPOINTS XUP NOT
C                        IN ASCENDING ORDER
C                = 5 IMPROPER INPUT - XLO.GT.XUP(1)
C                = 6 INTEGRATION SUCCESSFUL BUT AT LEAST ONE ENDPOINT
C                        NOT WITHIN SPAN OF ABSCISSAS
C              ** NOTE.  ERRCHK PROCESSES DIAGNOSTICS FOR CODES 2,3,4,5
C
C   CHECK FOR IMPROPER INPUT
C
       IERR = 2
       IF(N .LT. 4  .OR.  NUP .LT. 1) THEN
         RETURN
       ENDIF
       NM1 = N-1
       NM2 = N-2
       IERR = 3
       DO 2 I = 1,NM1
         IF(X(I) .GE. X(I+1)) THEN
           RETURN
         ENDIF
 2     CONTINUE
       IF(NUP .NE. 1) THEN
         IERR = 4
         DO 3 I = 2,NUP
           IF(XUP(I-1) .GT. XUP(I)) THEN
             RETURN
           ENDIF
 3       CONTINUE
       ENDIF
       IERR = 5
       IF(XLO .GT. XUP(1)) THEN
         RETURN
       ENDIF
       IERR = 1
       IF(XLO .LT. X(1)  .OR.  XUP(NUP) .GT. X(N)) IERR = 6
C
C   LOCATE XLO IN INTERVAL (X(I),X(I+1))
C
       DO 10 I = 1,NM2
         IF(XLO .LT. X(I+1)) GO TO 20
 10      CONTINUE
       I = NM1
 20    HLO = XLO-X(I)
       HLO2 = HLO*HLO
       HI = X(I+1)-X(I)
       HI2 = HI*HI
       DO 30 J = 1,NUP
         IF(XUP(J) .GT. X(I+1)  .AND.  XLO .LT. X(NM1)) GO TO 40
C
C   COMPUTE SPECIAL CASES OF XUP IN INTERVAL WITH XLO
C
         HUP = XUP(J)-X(I)
         HSUM = HUP+HLO
         HDIFF = HUP-HLO
         HUP2 = HUP*HUP
         SUM = (YPP(I+1)-YPP(I))*HSUM*HDIFF*(HUP2+HLO2)/(24*HI)
         SUM = SUM + YPP(I)*HDIFF*(HUP2+HLO*HUP+HLO2)/6
         SUM = SUM + YP(I)*HDIFF*HSUM/2
         SUM = SUM + Y(I)*HDIFF
 30    ANS(J) = SUM
       RETURN
C
C   COMPUTE INTEGRAL BETWEEN XLO AND X(I+1) AS FOUR TERMS IN TAYLOR
C   POLYNOMIAL AND ADVANCE I TO I+1
C
 40    HDIFF = HI-HLO
       HSUM = HI+HLO
       SUM0 = Y(I)*HDIFF
       SUM1 = YP(I)*HDIFF*HSUM
       SUM2 = YPP(I)*HDIFF*(HI2+HI*HLO+HLO2)
       SUM3 = (YPP(I+1)-YPP(I))*HDIFF*HSUM*(HI2+HLO2)/HI
       I = I+1
C
C   LOCATE EACH XUP(M) IN INTERVAL (X(I),X(I+1))
C
       DO 80 M = J,NUP
 50      IF(XUP(M) .LT. X(I+1)  .OR.  I .EQ. NM1) GO TO 60
C
C   AUGMENT INTEGRAL BETWEEN ABSCISSAS TO INCLUDE INTERVAL
C   (X(I),X(I+1)) AND ADVANCE I TO I+1
C
         HI = X(I+1)-X(I)
         HI2 = HI*HI
         HI3 = HI2*HI
         SUM0 = SUM0 + Y(I)*HI
         SUM1 = SUM1 + YP(I)*HI2
         SUM2 = SUM2 + YPP(I)*HI3
         SUM3 = SUM3 + (YPP(I+1)-YPP(I))*HI3
         I = I+1
         GO TO 50
C
C   INTEGRAL BETWEEN X(I) AND XUP(M) IS ZERO
C
 60      IF(XUP(M) .NE. X(I)) THEN
C
C   COMPUTE INTEGRAL BETWEEN X(I) AND XUP(M) AND EVALUATE
C   TAYLOR POLYNOMIAL IN REVERSE ORDER
C
           HUP = XUP(M)-X(I)
           HUP2 = HUP*HUP
           HUP3 = HUP2*HUP
           HUP4 = HUP3*HUP
           HI = X(I+1)-X(I)
           PSUM0 = Y(I)*HUP
           PSUM1 = YP(I)*HUP2
           PSUM2 = YPP(I)*HUP3
           PSUM3 = (YPP(I+1)-YPP(I))*HUP4/HI
           SUM = (SUM3+PSUM3)/24 + (SUM2+PSUM2)/6
           SUM = SUM + (SUM1+PSUM1)/2
           SUM = SUM + (SUM0+PSUM0)
         ELSE
           SUM = ((SUM3/24 + SUM2/6) + SUM1/2) + SUM0
         ENDIF
 80    ANS(M) = SUM
       RETURN
       END
C
C
C
c
c  ********************************************************
c  *                                                      *
c  *   njtj                                               *
c  *     These are machine dependent routines.            *
c  *   Included are routine for Apollo, Sun,               *
c  *   Vax, and Cray systems.  The user must              *
c  *   1)compile with their systems lines uncommented     *
c  *   or 2)supply their own                              *
c  *   or 3)remove-comment out all references to          *
c  *   these calls in the program.                        *
c  *                                                      *
c  ********************************************************
c
c  ****************Apollo start***********************
c
C
C **************Cray start***********************
C
       DOUBLE PRECISION FUNCTION SBESSJ(N,X)

       implicit double precision (a-h, o-z)

       PARAMETER(ONE=1.D0,TWO=2.D0,THREE=3.D0,ZERO=0.D0)
       PARAMETER( FIVE = 5.0D0 , TEN = 10.0D0 , FOURTN = 14.0D0 )
C      SPHERICAL BESSEL FUNCTION OF THE FIRST KIND
C
       IF(ABS(X) .GT. 0.001) THEN
         SB0 = SIN(X)/X
       ELSE
         X2 = X*X/TWO
         SB0 = ONE - (X2/THREE)*(ONE - X2/TEN)
       ENDIF
       IF(N .EQ. 0) THEN
         SBESSJ = SB0
       ELSE
         IF(ABS(X) .GT. 0.001) THEN
           SB1 = (SIN(X)/X - COS(X)) / X
         ELSE
           X2 = X*X/TWO
           SB1 = (X/THREE)*(ONE - (X2/FIVE)*(1.0 - X2/FOURTN))
         ENDIF
         IF(N .EQ. 1) THEN
           SBESSJ = SB1
         ELSEIF(X .EQ. ZERO) THEN
           SBESSJ = ZERO
         ELSE
           BY = SB1
           BYM = SB0
           UX = ONE / X
           DO 10 J=1,N-1
             BYP = REAL(2*J+1)*UX*BY - BYM
             BYM = BY
             BY = BYP
 10        CONTINUE
           SBESSJ = BY
         ENDIF
       ENDIF
       RETURN
       END

       SUBROUTINE SPLIFT (X,Y,YP,YPP,N,W,IERR,ISX,A1,B1,AN,BN)
C
       implicit double precision (a-h, o-z)

       PARAMETER (FOUR=4.D0)
CRAY      PARAMETER (FOUR=4.0)
C
C  NJTJ
C  ###  CRAY CONVERSIONS
C  ###    1)Comment out the implicit double precision.
C  ###    2)Switch double precision parameter
C  ###      to single precision parameter
C  ###  CRAY CONVERSIONS
C  NJTJ
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
  237 FORMAT(F5.1,F5.1,2I5)
C
C     ABSTRACT
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
C
C       --INPUT--
C
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
C         Y    - ARRAY OF ORDINATES OF DATA
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
C                SET TO ONE FOR FASTER EXECUTION.
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
C                OF THE SPLINE AT THE END POINTS (SEE YPP).
C                THE END CONDITION CONSTRAINTS ARE
C                        YPP(1) = A1*YPP(2) + B1
C                AND
C                        YPP(N) = AN*YPP(N-1) + BN
C                WHERE
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
C
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
C
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
C                FOLLOWING CONDITIONS.
C
C                A1=-0.5
C
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
C
C                AN=-0.5
C
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
C
C       --OUTPUT--
C
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
C         IERR - A STATUS CODE
C              --NORMAL CODE
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
C              --ABNORMAL CODES
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.
C
C       --WORK--
C
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
       DIMENSION X(N),Y(N),YP(N),YPP(N),W(N,3)
C
       IF (N.LT.4) THEN
         IERR = 2
         RETURN
       ENDIF
       NM1  = N-1
       NM2  = N-2
       IF (ISX.GT.0) GO TO 40
       DO 5 I=2,N
         IF (X(I)-X(I-1) .LE. 0) THEN
           IERR = 3
           RETURN
         ENDIF
 5     CONTINUE
C
C     DEFINE THE TRIDIAGONAL MATRIX
C
       W(1,3) = X(2)-X(1)
       DO 10 I=2,NM1
         W(I,2) = W(I-1,3)
         W(I,3) = X(I+1)-X(I)
 10      W(I,1) = 2*(W(I,2)+W(I,3))
       W(1,1) = FOUR
       W(1,3) =-4*A1
       W(N,1) = FOUR
       W(N,2) =-4*AN
C
C     L U DECOMPOSITION
C
       DO 30 I=2,N
         W(I-1,3) = W(I-1,3)/W(I-1,1)
 30    W(I,1) = W(I,1) - W(I,2)*W(I-1,3)
C
C     DEFINE *CONSTANT* VECTOR
C
 40   YPP(1) = 4*B1
      DOLD = (Y(2)-Y(1))/W(2,2)
      DO 50 I=2,NM2
        DNEW   = (Y(I+1) - Y(I))/W(I+1,2)
        YPP(I) = 6*(DNEW - DOLD)
        YP(I)  = DOLD
 50   DOLD = DNEW
      DNEW = (Y(N)-Y(N-1))/(X(N)-X(N-1))
      YPP(NM1) = 6*(DNEW - DOLD)
      YPP(N) = 4*BN
      YP(NM1)= DOLD
      YP(N) = DNEW
C
C     FORWARD SUBSTITUTION
C
      YPP(1) = YPP(1)/W(1,1)
      DO 60 I=2,N
 60   YPP(I) = (YPP(I) - W(I,2)*YPP(I-1))/W(I,1)
C
C     BACKWARD SUBSTITUTION
C
       DO 70 J=1,NM1
         I = N-J
   70 YPP(I) = YPP(I) - W(I,3)*YPP(I+1)
C
C     COMPUTE FIRST DERIVATIVES
C
      YP(1) = (Y(2)-Y(1))/(X(2)-X(1)) - (X(2)-X(1))*(2*YPP(1)
     1  + YPP(2))/6
      DO 80 I=2,NM1
 80   YP(I) = YP(I) + W(I,2)*(YPP(I-1) + 2*YPP(I))/6
      YP(N) = YP(N) + (X(N)-X(NM1))*(YPP(NM1) + 2*YPP(N))/6
C
      IERR = 1
      RETURN
      END
C
C
C
      subroutine ext(i)
       implicit double precision (a-h, o-z)
c
c  Stops program in case of errors or completion.
c
c  i is a stop parameter
c   000-099 main (0 is normal exit)
c   100-199 input
c   200-299 charge
c   300-399 vionic
c   400-499 velect
c   500-599 dsolv1
c   600-699 dsolv2 (including difnrl and difrel)
c   700-799 etotal
c   800-899 pseudo, pseudk, pseudt and pseudv
c
      if (i .ne. 0) write(6,10) i
 10   format('stop parameter =',i3)
      close (unit=1)
      close (unit=3)
      close (unit=5)
      close (unit=6)
      call exit
      end
C
C
       SUBROUTINE SPLINT (X,Y,YPP,N,XI,YI,YPI,YPPI,NI,KERR)
       implicit double precision (a-h,o-z)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C
C         SPLINT EVALUATES A CUBIC SPLINE AND ITS FIRST AND SECOND
C         DERIVATIVES AT THE ABSCISSAS IN XI.  THE SPLINE (WHICH
C         IS DEFINED BY X, Y, AND YPP) MAY HAVE BEEN DETERMINED BY
C         SPLIFT OR SMOO OR ANY OTHER SPLINE FITTING ROUTINE THAT
C         PROVIDES SECOND DERIVATIVES.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YPP(N), XI(NI), YI(NI), YPI(NI), YPPI(NI)
C
C       --INPUT--
C
C         X   - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C         Y   - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C         YPP - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C         N   - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C               THE ARRAYS X, Y, AND YPP MUST BE DIMENSIONED AT LEAST N
C               N MUST BE GREATER THAN OR EQUAL TO 2.
C         XI  - THE ABSCISSA OR ARRAY OF ABSCISSAS (IN ARBITRARY ORDER)
C               AT WHICH THE SPLINE IS TO BE EVALUATED.
C               EACH XI(K) THAT LIES BETWEEN X(1) AND X(N) IS A CASE OF
C               INTERPOLATION.  EACH XI(K) THAT DOES NOT LIE BETWEEN
C               X(1) AND X(N) IS A CASE OF EXTRAPOLATION.  BOTH CASES
C               ARE ALLOWED.  SEE DESCRIPTION OF KERR.
C         NI  - THE NUMBER OF ABSCISSAS AT WHICH THE SPLINE IS TO BE
C               EVALUATED.  IF NI IS GREATER THAN 1, THEN XI, YI, YPI,
C               AND YPPI MUST BE ARRAYS DIMENSIONED AT LEAST NI.
C               NI MUST BE GREATER THAN OR EQUAL TO 1.
C
C       --OUTPUT--
C
C         YI  - ARRAY OF VALUES OF THE SPLINE (ORDINATES) AT XI.
C         YPI - ARRAY OF VALUES OF THE FIRST DERIVATIVE OF SPLINE AT XI
C         YPPI- ARRAY OF VALUES OF SECOND DERIVATIVES OF SPLINE AT XI.
C         KERR- A STATUS CODE
C             --NORMAL CODES
C                1 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI USING ONLY INTERPOLATION.
C                2 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI, BUT AT LEAST ONE EXTRAPOLATION WAS PERFORMED.
C             -- ABNORMAL CODE
C                3 MEANS THAT THE REQUESTED NUMBER OF EVALUATIONS, NI,
C                  WAS NOT POSITIVE.
C
       DIMENSION X(N),Y(N),YPP(N),XI(NI),YI(NI),YPI(NI),YPPI(NI)
C
C     CHECK INPUT
C
      IF (NI) 1,1,2
 1    CONTINUE
C    1 CALL ERRCHK(67,67HIN SPLINT,  THE REQUESTED NUMBER OF INTERPOLATI
C     1NS WAS NOT POSITIVE)
      KERR = 3
      RETURN
    2 KERR = 1
      NM1= N-1
C
C     K IS INDEX ON VALUE OF XI BEING WORKED ON.  XX IS THAT VALUE.
C     I IS CURRENT INDEX INTO X ARRAY.
C
       K  = 1
       XX = XI(1)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IL = 1
       IR = N
C
C     BISECTION SEARCH
C
   10 I  = (IL+IR)/2
       IF (I.EQ.IL) GO TO 100
       IF (XX-X(I)) 20,100,30
   20 IR = I
       GO TO 10
   30 IL = I
       GO TO 10
C
C     LINEAR FORWARD SEARCH
C
   50 IF (XX-X(I+1)) 100,100,60
   60 IF (I.GE.NM1) GO TO 80
       I  = I+1
       GO TO 50
C
C     EXTRAPOLATION
C
   80 KERR = 2
      I  = NM1
      GO TO 100
   90 KERR = 2
      I  = 1
C
C     INTERPOLATION
C
  100 H  = X(I+1) - X(I)
       H2 = H*H
       XR = (X(I+1)-XX)/H
       XR2= XR*XR
       XR3= XR*XR2
       XL = (XX-X(I))/H
       XL2= XL*XL
       XL3= XL*XL2
       YI(K) = Y(I)*XR + Y(I+1)*XL
     1       -H2*(YPP(I)*(XR-XR3) + YPP(I+1)*(XL-XL3))/6.0D0
       YPI(K) = (Y(I+1)-Y(I))/H
     1 +H*(YPP(I)*(1.0D0-3.0D0*XR2)-YPP(I+1)*(1.0D0-3.0D0*XL2))/6.0D0
       YPPI(K) = YPP(I)*XR + YPP(I+1)*XL
C
C     NEXT POINT
C
       IF (K.GE.NI) RETURN
       K = K+1
       XX = XI(K)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IF (XX-XI(K-1)) 110,100,50
  110 IL = 1
       IR = I+1
       GO TO 10
C
       END
C Added LLOCAL argument 10/29/91...
       SUBROUTINE KBLY(VQL,VIOD,CDC,CDD,R,INORM,NR,LLOCAL,NAMEAT)
       implicit double precision (a-h, o-z)
c
      PARAMETER (NRMAX=2000,LMAX=5)
c
      PARAMETER (ZERO=0.D0,SMALL=0.5D-5,ONE=1.D0,TWO=2.D0)
c
      DIMENSION R(NRMAX),RAB(NRMAX),NO(2*LMAX),LO(2*LMAX),SO(LMAX),
     1 VIOD(LMAX,NRMAX),VIOU(LMAX,NRMAX),VID(NRMAX),CDD(NRMAX),
     2 CDC(NRMAX),EV(LMAX),Y(NRMAX),YP(NRMAX),YPP(NRMAX),W(NRMAX*3),
     3 CDU(NRMAX),VIU(NRMAX),S1(NRMAX),S2(NRMAX),ETOT(10),EVL(2),
     4 AR(NRMAX,LMAX),ANORM(LMAX),INORM(LMAX),VQL(NRMAX),ARD(NRMAX),
     5 EVI(5),EVD(5)
       DIMENSION R_H(NRMAX),VNL_H(NRMAX,4),SNORM(10),VC_RHO(NRMAX)
     C,IND1(NRMAX),XI(NRMAX),YI(NRMAX),YPI(NRMAX),YPPI(NRMAX)
     C,RCR2(NRMAX),RVR2(NRMAX)
c
      CHARACTER*1 ISPP
      CHARACTER*2 ICORR,NAMEAT
      CHARACTER*3 IREL
      CHARACTER*4 NICORE
      CHARACTER*10 IRAY(6),ITITLE(7),DATED
c
      DATA SO/5*ZERO/
      DATA EVI/5*ZERO/
      DATA EVD/5*ZERO/
C
C     READ PARAMETERS FOR LIN YANG FORMAT OF PPOT
      OPEN(UNIT=11,FILE='kbly.d',FORM='FORMATTED')
C     CALL LINK("UNIT11=(kbly.d,open,text)//")
      OPEN(UNIT=13,FILE='chplot.d',FORM='FORMATTED')
C     READ(11,*)DRH,ZVH,LMAXLY,IHMAX
       READ(11,237)DRH,ZVH,LMAXLY,IHMAX,I_PCC
  237 FORMAT(F5.1,F5.1,3I5)
      WRITE(6,311)DRH,ZVH,LMAXLY,IHMAX, I_PCC
  311 FORMAT(/,2X,'KBLY,DRH,ZVH,LMAXLY,IHMAX,I_PCC=',2F6.3,3I5)
C     ONE PROJECTION
      NLMAXLY=1
      SQ2I=1./SQRT(2.)
C
      R_H(1)=0.0
      IHM2=2*IHMAX
      XI(1)=0.0
      DO 348 L=1,LMAXLY
      VNL_H(1,L)=0.0
  348 CONTINUE
      VNL_H(1,LMAXLY+1)=0.5*VQL(2)
      DO 100 I=2,IHM2
      XI(I)=(I-1)*DRH
  100 R_H(I)=(I-1)*DRH

C Changes made to the following loop, 10/29/91 to deal with LLOCAL.NE.LMAXLY
C Inserted:
      lvalue = 0
C Loop to LMAXLY+1 since one of the iterations will be skipped:
      DO 109 L=1,LMAXLY+1
C Following lines inserted 10/29/91
	if (l.eq.llocal + 1) goto 109
      lvalue = lvalue + 1
      ISX = 0
      A1 = ZERO
      AN = ZERO
      B1 = ZERO
      BN = ZERO
      NRM = NR
      DO 108 N=1,NR
      Y(N)=VIOD(L,N)
  108 CONTINUE
      CALL SPLIFT(R,Y,YP,YPP,NRM,W,IERR,ISX,A1,B1,AN,BN)
      CALL SPLINT(R,Y,YPP,NRM,XI,YI,YPI,YPPI,IHMAX,KERR)
      DO 110 I=2,IHMAX
C     CHANGE FROM RY TO HARTREE UNITS
      RLY=R_H(I)
      VNL_H(I,lvalue)=RLY*SQ2I*YI(I)
  110 CONTINUE
  109 CONTINUE
      ISX = 0
      A1 = ZERO
      AN = ZERO
      B1 = ZERO
      BN = ZERO
      NRM = NR
      DO 308 N=1,NR
      Y(N)=VQL(N)
  308 CONTINUE
      CALL SPLIFT(R,Y,YP,YPP,NRM,W,IERR,ISX,A1,B1,AN,BN)
      CALL SPLINT(R,Y,YPP,NRM,XI,YI,YPI,YPPI,IHMAX,KERR)
      DO 501 I=2,IHMAX
      VNL_H(I,LMAXLY+1)=0.5*YI(I)
  501 CONTINUE
      ISX = 0
      A1 = ZERO
      AN = ZERO
      B1 = ZERO
      BN = ZERO
      NRM = NR
      DO 608 N=1,NR
      Y(N)=CDC(N)
  608 CONTINUE
      CALL SPLIFT(R,Y,YP,YPP,NRM,W,IERR,ISX,A1,B1,AN,BN)
      CALL SPLINT(R,Y,YPP,NRM,XI,YI,YPI,YPPI,IHM2,KERR)
      PI=3.14159265358E0
      FAC=4.0*PI
      VC_RHO(1)=CDC(2)/(FAC*R(2)**2)
      DO 601 I=2,IHM2
      VC_RHO(I)=YI(I)/(FAC*R_H(I)**2)
      RCR2(I)=YI(I)
  601 CONTINUE
      DO 611 N=1,NR
  611 Y(N)=CDD(N)
      CALL SPLIFT(R,Y,YP,YPP,NRM,W,IERR,ISX,A1,B1,AN,BN)
      CALL SPLINT(R,Y,YPP,NRM,XI,YI,YPI,YPPI,IHM2,KERR)
      DO 613 I=2,IHM2
      RVR2(I)=YI(I)
  613 CONTINUE
      SUMC=0.0
      SUMV=0.0
      DO 615 I=2,IHM2
      SUMC=SUMC+RCR2(I)
      SUMV=SUMV+RVR2(I)
      WRITE(13,638)R_H(I),RCR2(I),RVR2(I)
  638 FORMAT(2X,F10.5,3X,E11.4,3X,E11.4)
  615 CONTINUE

      OPEN(UNIT=12,FILE='vchden',FORM='FORMATTED')
      WRITE(12, 9250) 1
      WRITE(12, 9250) IHM2-1
 9250 FORMAT(2X, I5)
      DO 9615 I=2,IHM2
      WRITE(12,638) R_H(I), RVR2(I)/(FAC*R_H(I)**2)
 9615 CONTINUE
      CLOSE(UNIT=12)

      PI=3.14159265358E0
      FAC= DRH
      SUMC=SUMC*FAC
      SUMV=SUMV*FAC
      WRITE(6,621)SUMC,SUMV
  621 FORMAT(/,2X,'INTEGRATED CORE AND VALENCE CH. KBLY',E12.5,2X,
     CE12.5,/)
C     FIND MAXIMUM ELEMENT OF KB ARRAYS
      lvalue = 0
      DO 421 L=1,LMAXLY+1
	if (l .eq. llocal + 1) goto 421
	lvalue = lvalue + 1
      VARMAX=0.0
      IMAX=1
      DO 422 I=2,IHMAX
      VAR=ABS(VNL_H(I,lvalue))
      IF(VAR.LT.VARMAX)GO TO 422
      VARMAX=VAR
      IMAX=I
  422 CONTINUE
      SNORM(lvalue)=INORM(L)*VARMAX*VARMAX
      DO 425 I=1,IHMAX
      VNL_H(I,lvalue)=VNL_H(I,lvalue)/VARMAX
  425 CONTINUE
  421 CONTINUE
      OPEN(UNIT=12,FILE='wavefunctions',FORM='FORMATTED')
      WRITE(12,9239) IHMAX, (SNORM(L),L=1,LMAXLY)
 9239 FORMAT(2X,I5, E14.7,2X,E14.7,2X,E14.7,2X,E14.7)
      WRITE(12,9239) IHMAX
      DO 9145 I=1,IHMAX
      WRITE(12,9241)R_H(I),VNL_H(I,LMAXLY+1)
 9145 CONTINUE
      DO 9150 I=1,IHMAX
      WRITE(12,9241)R_H(I),(VNL_H(I,L),L=1,LMAXLY)
 9241 FORMAT(2X,F6.3,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7,
     C2X,E14.7)
 9150 CONTINUE
      CLOSE(unit=12)

C     WRITE DATASET IN LY FORMAT
      OPEN(UNIT=12,FILE='kbly.o',FORM='FORMATTED')
C     CALL LINK("UNIT12=(kbly.o,create,text)//")
      IONE = 1
      WRITE(12,279) IONE,I_PCC,NAMEAT
  279 FORMAT(2I5, 3X, A2, A10)
      NCRHO = IHMAX
      WRITE(12,243) ZVH,LMAXLY,NLMAXLY,IHMAX, LLOCAL, NCRHO
  243 FORMAT(2X,F5.1,5I5)
      WRITE(6,339)(INORM(L),L=1,LMAXLY + 1)
  339 FORMAT(2X,'KBLY,INORM=',5I4)
      WRITE(12,239)(SNORM(L),L=1,LMAXLY)
  239 FORMAT(2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7)
      DO 150 I=1,IHMAX
      WRITE(12,241)R_H(I),(VNL_H(I,L),L=1,LMAXLY+1),VC_RHO(I)
      WRITE(6,241)R_H(I),(VNL_H(I,L),L=1,LMAXLY+1)
  241 FORMAT(2X,F6.3,2X,E14.7,2X,E14.7,2X,E14.7,2X,E14.7,
     C2X,E14.7)
  150 CONTINUE
      WRITE(12,493)
  493 FORMAT(2X,'TROUL.46,6/11/91')
      CLOSE(unit=12)
      RETURN
      END
c $Id$
