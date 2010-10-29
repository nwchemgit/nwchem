C$Id$
      Program averagehist_5 

      Implicit none

      double precision au_to_kcal_mol, au_to_angs
      double precision rmin, rmax, rcut
      integer natom, nmcchain

      integer nob
      parameter(nob=200) ! bins in grid

      integer mxatom
      parameter(mxatom=100)

      double precision Number(NOB) ! grid
      double precision Number2(NOB),Number3(NOB,(2*mxatom)+2)
      double precision Number4(NOB)
      double precision sum, sum2, sum3((2*mxatom)+2)
      double precision sum4
      double precision RB(NOB), NDA(NOB) ! labels and average grid
      double precision UU, RCONFIG ! Store information
      double precision UUinit
      character*256 text ! text holder
      character*256 filename ! filename for read and write files
      character*4 prcfil ! stores mc number 0000
      character*256 suffix ! cocat holder
      double precision Coord(mxatom,3) ! Store information
      integer ndistatom
      integer atomlabel(mxatom)
      double precision distance(mxatom*(mxatom-1)/2)
      double precision averagedistance
      double precision averageenergy((2*mxatom)+2)
      double precision RB2(NOB), NDA2(NOB)
      double precision emin, emax
      double precision ohmin, ohmax
      double precision RB3(NOB), NDA3(NOB,(2*mxatom)+2)
      double precision averageohdistance
      double precision ohdistance(mxatom*4),hodistance(mxatom*4)
      double precision ohangle(mxatom*4),PI
      double precision ohanglemax
      double precision RB4(NOB),NDA4(NOB)
      double precision majordistance(2)
      double precision averageohangle
      double precision averagehbonds
      integer hold
      integer numberhbonds
      double precision weight, kb, temp, tempinit
      double precision DIPOLE(3)
      double precision mass(mxatom),dipmax
      double precision RB5(NOB), number5(NOB,3), NDA5(NOB,3), sum5(3)
      double precision averageMOI(3), averageDipole(3),averageDipole2(3)
      double precision AtomDist(mxatom), MOI(3,3), DMOI(3), PMOI(3,3)
      double precision PDIPOLE(3),wk(40)
      integer ik(20), stat

      external dsyevd

      double precision dipvec(2,3),diplen(2),rvec(3),rlen,theta(3)
      double precision  averagedipint, dipint, dipinthold
      double precision number6(NOB),RB6(NOB),sum6,NDA6(NOB)

      integer i,j,l,n,o,IR
      integer R, M, P

      character*16 eletags(mxatom)
      integer idum
      double precision dntmc_RAN3    

      double precision mhodistance(mxatom*4), hoangle(mxatom*4)
      double precision averagehoangle, hoanglemax
      double precision RB7(NOB), number7(NOB), NDA7(NOB), sum7

      integer hbcounter(mxatom*2)
      
C     Input Parameters
      PI = 3.141592653589793D+00 ! Constant PI
      au_to_kcal_mol = 627.5095D+00 ! Conversion from au to kcal/mole
      au_to_angs = 0.52917720832D+00 ! Conversion from au to Angstroms
      kb = 1.987905D-03 ! Constant Boltzmann kcal/mole/K
C      emin = -20.00D+00
C      emax = 5.00D+00
C      emin = -95430.00D+00 ! min for energy bins
C      emax = -95405.00D+00 ! max for energy bins
      ohmin = 1.50D+00 ! min for hbond length (selection) and bins
      ohmax = 2.50D+00 ! max for hbond length (selection) and bins
      ohanglemax = 60.00D+00 ! max for hbond angle (selection) and bins, min is 0
      hoanglemax = 120.00D+00 ! max for hbond angle (2nd) and bins, min is 0
      idum = 7

      Write(*,*)'Please input file name prefix:'
      Write(*,*)'Tempinit and Temp' !inital and final temperature (K) for scaled distributions
      Write(*,*)'rmin and rmax' !grid min and max for O-O distance (Ang.)
      Write(*,*)'natoms, nmcchains' !number of atoms per configuration and markov chains (files)
      Write(*,*)'rcut' !cutoff of RConfig (All values .le.)
      Write(*,*)'emin and emax' !grid min and max for Energies (kcal/mole)
      Write(*,*)'maximum dipole moment' !grid max for Dipole Moments(au)
      Write(*,*)'number of monomers' !number of distinct water monomers
      Write(*,*)'Array of oxygen atom numbers' !atom numbers of oxygen atoms of each monomer
      Write(*,*)'Array of atomic masses' !Atomic masses (amu or g/mol)
      Read(*,*)filename
      Read(*,*)tempinit, temp
      Read(*,*)rmin,rmax
      Read(*,*)natom,nmcchain
      Read(*,*)rcut
      Read(*,*)emin, emax
      Read(*,*)dipmax
      Read(*,*)ndistatom
      Read(*,*)(atomlabel(IR),IR=1,ndistatom)
      Read(*,*)(mass(IR),IR=1,natom)

C     Determine number of possible h-bonds + 2
       if (ndistatom .gt. 5) then
         hold = (2*ndistatom) + 2
       else
         hold = ((ndistatom*(ndistatom-1)/2) + 2)
       endif

C     Open Input/Output Files
      prcfil(1:4)='0000'
      Do i = 20, 19+nmcchain
        write(prcfil(1:4), '(i4)') i-19
        if((i-19).lt. 10) then
          prcfil(1:3)='000'
        else if((i-19).lt.100) then
          prcfil(1:2)='00'
        else if((i-19).lt.1000) then
          prcfil(1:1)='0'
        endif

        write(suffix, '(a,a)')'.MCdata.', prcfil(1:4)
        Open(Unit = i,FILE = filename(1:len_trim(filename))
     &    //suffix(1:len_trim(suffix))
     &    ,STATUS = 'OLD')
       enddo

        write(suffix, '(a,i3,a)')'.fdist.average.',int(temp),'K'
        Open(Unit = 10,FILE = filename(1:len_trim(filename))
     &     //suffix(1:len_trim(suffix))
     &     ,STATUS = 'REPLACE')

      prcfil(1:4)='0000'
      Do i = 100, 99+hold
        write(prcfil(1:4), '(i4)') i-100
        if((i-100).lt. 10) then
          prcfil(1:3)='000'
        else if((i-100).lt.100) then
          prcfil(1:2)='00'
        else if((i-100).lt.1000) then
          prcfil(1:1)='0'
        endif

        write(suffix, '(a,i3,a,a)')'.fdist.average.',int(temp),'K.hb.',
     &                             prcfil(1:4)
        Open(Unit = i,FILE = filename(1:len_trim(filename))
     &     //suffix(1:len_trim(suffix))
     &     ,STATUS = 'REPLACE')
      enddo

C     Setup Grid  General
      RB(1) = rmin + (rmax - rmin)/DBLE(NOB)
      RB2(1) = emin + (emax -emin)/DBLE(NOB)
      RB3(1) = ohmin + (ohmax-ohmin)/DBLE(NOB)
      RB4(1) = (ohanglemax)/DBLE(NOB)
      RB6(1) = -0.100D+00 + (0.200D+00)/DBLE(NOB)
      RB5(1) = (dipmax)/DBLE(NOB)
      RB7(1) = (hoanglemax)/DBLE(NOB)
      DO i = 2, nob
        RB(i) = RB(i-1) + (rmax - rmin)/DBLE(NOB)
        RB2(i) = RB2(i-1) + (emax -emin)/DBLE(NOB)
        RB3(i) = RB3(i-1) + (ohmax-ohmin)/DBLE(NOB)
        RB4(i) = RB4(i-1) + (ohanglemax)/DBLE(NOB)
        RB6(i) = RB6(i-1) + (0.200D+00)/DBLE(NOB)
        RB5(i) = RB5(i-1) + (dipmax)/DBLE(NOB)
        RB7(i) = RB7(i-1) + (hoanglemax/DBLE(NOB))
      ENDDO
      Do i = 1, nob
        RB(i) = RB(i) - ((rmax - rmin)/DBLE(NOB)/2.0D+00)
        RB2(i) = RB2(i) - ((emax -emin)/DBLE(NOB)/2.0D+00)
        RB3(i) = RB3(i) - ((ohmax-ohmin)/DBLE(NOB)/2.0D+00)
        RB4(i) = RB4(i) - ((ohanglemax)/DBLE(NOB)/2.0D+00) 
        RB6(i) = RB6(i) - ((0.200D+00)/DBLE(NOB)/2.00D+00)
        RB5(i) = RB5(i) - ((dipmax)/DBLE(NOB)/2.0D+00)
        RB7(i) = RB7(i) - ((hoanglemax)/DBLE(NOB)/2.0D+00)
      EndDo

C     Initialize Bins General
      Do i = 1, NOB 
       number(i) = 0.00D+00 
       number2(i) = 0.00D+00
       number4(i) = 0.00D+00
       number6(i) = 0.00D+00
       number7(i) = 0.00D+00
      Do j = 1, (2*mxatom)+2
       number3(i,j) = 0.00D+00
      EndDO
      Do j = 1, 3
        number5(i,j) = 0.00D+00
      ENDDO
      EndDo
      averagedistance = 0.00D+00
      Do i = 1, (2*mxatom)+2
      averageenergy(i) = 0.00D+00
      EndDo
      averageohdistance = 0.00D+00
      averageohangle = 0.00D+00
      averagehoangle = 0.00D+00
      averagehbonds = 0.00D+00
      averagedipint = 0.00D+00
      UUinit = 0.00D+00
      Do i = 1, 3
      averageMOI(i) = 0.00D+00
      averageDipole(i) = 0.00D+00
      averageDipole2(i) = 0.00D+00
      EndDo

C     Begin a Do loop over records General
      R = 0
      Do i = 20, 19+nmcchain

C     Loop over items in record with while loop
5678    Continue
C     Data read
        Do j = 1, natom
C          Read(i,*,END=1234)(Coord(j,IR),IR=1,3)
          Read(i,*,END=1234)eletags(j),(Coord(j,IR),IR=1,3)
        EndDo
        Read(i,*,END=1234)text,text,UU
        Read(i,*,END=1234)text,text,(DIPOLE(IR),IR=1,3)
        Read(i,*,END=1234)text,text,RCONFIG
C     Check if within scope of average
        if(RCONFIG .gt. rcut) goto 5678 ! not in scope, new data
        R = R + 1 ! have usable data

C     Set initial Energy only once
        if (UUinit .eq. 0.00D+00) then
          UUinit = UU
        endif
C     Analze data General O-O distance
        M = 0
        Do j = 1, ndistatom - 1
          Do l = j+1, ndistatom
            M = M + 1
            distance(M) = 0.00D+00
            Do n = 1, 3
              distance(M) = distance(M) +
     &        (Coord(atomlabel(j),n) - Coord(atomlabel(l),n))**2.00D+00
            EndDo
            distance(M) = DSQRT(distance(M))
          EndDo
        EndDo

C     Analyze some more data Dipole Interaction Energy
        dipint = 0.00D+00
        Do j = 1, ndistatom - 1
        Do l = j+1, ndistatom
          Do n = 1, 3
            rvec(n) = Coord(atomlabel(l),n) - Coord(atomlabel(j),n)
            dipvec(1,n) =Coord(atomlabel(j)+1,n)+Coord(atomlabel(j)+2,n)
            dipvec(1,n) = dipvec(1,n)/2.00D+00
            dipvec(1,n) = dipvec(1,n) - Coord(atomlabel(j),n)
            dipvec(2,n) =Coord(atomlabel(l)+1,n)+Coord(atomlabel(l)+2,n)
            dipvec(2,n) = dipvec(2,n)/2.00D+00
            dipvec(2,n) = dipvec(2,n) - Coord(atomlabel(l),n)
          ENDDO
          rlen = 0.00D+00
          diplen(1) = 0.00D+00
          diplen(2) = 0.00D+00
          Do n = 1, 3
            rlen = rlen + (rvec(n)**2.00D+00)
            diplen(1) = diplen(1) + (dipvec(1,n)**2.00D+00)
            diplen(2) = diplen(2) + (dipvec(2,n)**2.00D+00)
          ENDDO
          rlen = DSQRT(rlen)
          diplen(1) = DSQRT(diplen(1))
          diplen(2) = DSQRT(diplen(2))
          Do n = 1, 3
            rvec(n) = rvec(n)/rlen
            dipvec(1,n) = dipvec(1,n)/diplen(1)
            dipvec(2,n) = dipvec(2,n)/diplen(2)
          ENDDO
          Do n = 1, 3
            theta(n) = 0.00D+00
          ENDDO
          Do n = 1, 3
            theta(1) = theta(1) + (rvec(n)*dipvec(1,n))
            theta(2) = theta(2) + (rvec(n)*dipvec(2,n))
            theta(3) = theta(3) + (dipvec(1,n)*dipvec(2,n))
          ENDDO
          Do n = 1, 3
            theta(n) = DACOS(theta(n))
          ENDDO
          dipinthold = -2.00D+00*DCOS(theta(1))*DCOS(theta(2))
          dipinthold = dipinthold + DCOS(theta(3))
          dipinthold = dipinthold * (rmin**3.00D+00) / (rlen**3.00D+00)
          dipint = dipint + dipinthold
        ENDDO
        ENDDO
             
C     Analyze some more data O--H distance and O--H-O angle
C     Added Bisect - O--H  angle and Bisect - O distance
        P = 0
        Do j = 1, ndistatom
          Do l = 1, ndistatom
            if (.not.(l .eq. j)) then
              Do n = 1, 2
                P = P + 1
                ohdistance(P) = 0.00D+00
                hodistance(P) = 0.00D+00
                mhodistance(P) = 0.00D+00
                ohangle(P) = 0.00D+00
                hoangle(P) = 0.00D+00
                Do o = 1, 3
                  ohdistance(P) = ohdistance(P) +
     &            (Coord(atomlabel(j),o) - Coord(atomlabel(l)+n,o))
     &            **2.00D+00
                  hodistance(P) = hodistance(P) +
     &            (Coord(atomlabel(l)+n,o) - Coord(atomlabel(l),o))
     &            **2.00D+00
                  mhodistance(P) = mhodistance(P) +
     &            (Coord(atomlabel(j),o) - ((Coord(atomlabel(j)+1,o) 
     &            + Coord(atomlabel(j)+2,o))/2.00D+00))**2.00D+00
                  ohangle(P) = ohangle(P) +
     &            (Coord(atomlabel(l)+n,o) - Coord(atomlabel(j),o))*
     &            (Coord(atomlabel(l),o) - Coord(atomlabel(l)+n,o))
                  hoangle(P) = hoangle(P) +
     &            (Coord(atomlabel(l)+n,o) - Coord(atomlabel(j),o))*
     &            (Coord(atomlabel(j),o) - ((Coord(atomlabel(j)+1,o)
     &            + Coord(atomlabel(j)+2,o))/2.00D+00))
                enddo
                ohdistance(P) = DSQRT(ohdistance(P))
                hodistance(P) = DSQRT(hodistance(P))
                mhodistance(P) = DSQRT(mhodistance(P))
                ohangle(P) = ohangle(P)/(ohdistance(P)*hodistance(P))
                ohangle(P) = DACOS(ohangle(P))
                hoangle(P) = hoangle(P)/(mhodistance(P)*ohdistance(P))
                hoangle(P) = DACOS(hoangle(P))
              enddo
            endif
          enddo
        enddo

C     Init counter
        Do j = 1, (ndistatom*2)
          hbcounter(j) = 0
        EndDo
C     Analyze number of hbonds per hydrogen
        P = 0
        Do j = 1, ndistatom
          Do l = 1, ndistatom
            if (.not. l .eq. j) then
              Do n = 1, 2
              P = P + 1
       if ((ohdistance(P) .ge. ohmin).and.(ohdistance(P).le.ohmax)) then
            if (.not. (ohangle(P) .gt. (ohanglemax*PI/180.00D+00))) then
            if (.not. (hoangle(P) .gt. (hoanglemax*PI/180.00D+00))) then
              hbcounter((2*l)-2+n) = hbcounter((2*l)-2+n) + 1
            endif
            endif
       endif
              EndDo
            endif
          EndDo
        EndDo
C     Count only one hydrogen bond per hydrogen
       numberhbonds = 0
       Do j = 1, (ndistatom*2)
         if (hbcounter(j) .ge. 1) then
           numberhbonds = numberhbonds + 1
         endif
       EndDo

CC     Analyze some more data number of hydrogen bonds
C       numberhbonds = 0
C        Do j = 1, (ndistatom*(ndistatom -1)*2)
C       if ((ohdistance(j) .ge. ohmin).and.(ohdistance(j).le.ohmax)) then
C            if (.not. (ohangle(j) .gt. (ohanglemax*PI/180.00D+00))) then
C            if (.not. (hoangle(j) .gt. (hoanglemax*PI/180.00D+00))) then
C               numberhbonds = numberhbonds + 1
C             endif
C             endif
C           endif
C        EndDo

C     Print coordinates (random) of some configurations of each # of hbonds
        if ((dntmc_RAN3(idum)*1.00D+05) .lt. 1.00D+00) then
        Do j = 1, hold
          if (j-1 .eq. numberhbonds) then
            write(j+99,*)natom
            write(j+99,*)numberhbonds, 'hydrogen bonds'
            Do l = 1, natom
              write(j+99,*)eletags(l),(Coord(l,IR),IR=1,3)
            EndDo
          endif
        ENDDO
        endif
          if (numberhbonds .gt. hold-2) then
            write(hold+99,*)natom
            write(hold+99,*)numberhbonds, 'hydrogen bonds'
            Do l = 1, natom
              write(hold+99,*)eletags(l),(Coord(l,IR),IR=1,3)
            EndDo
          endif

C     Analyze some more data, determine weight
        if (temp .eq. tempinit) then
          weight = 1.00D+00
        else
        weight = (UU - UUinit)*au_to_kcal_mol/kb*
     &           ((1.D+00/temp)-(1.D+00/tempinit))
        weight = DEXP(-weight)
        endif

C     Analyze some more data, MOI and Projected Dipole Moment
C     Distance to all atoms
      Do j = 1, natom
        AtomDist(j) = 0.00D+00
        Do l = 1, 3
          AtomDist(j) = AtomDist(j) +(Coord(j,l)**2.00D+00)
        ENDDO
      ENDDO
C     MOI Calculation
      Do j = 1, 3
        Do l = 1, 3
          MOI(j,l) = 0.00D+00
          PMOI(j,l) = 0.00D+00
          DMOI(j) = 0.00D+00
          Do n = 1, natom
            if(j .eq. l) then
              MOI(j,l) = MOI(j,l) +
     & (mass(n)*(AtomDist(n)-(Coord(n,j)*Coord(n,l))))
            else
              MOI(j,l) = MOI(j,l) +
     & (-mass(n)*Coord(n,j)*Coord(n,l))
            endif
          ENDDO
        ENDDO
      ENDDO
C     Copy MOI Matrix
      Do j = 1,3
        Do l = 1,3
         PMOI(j,l) = MOI(j,l)
        ENDDO
      ENDDO
C     Find PMOI and DMOI
      call dsyevd('V','U',3,PMOI(1,1),3,DMOI(1),WK(1),40,IK(1),20,stat)
      if (.not. (stat .eq. 0)) then
      WRITE(10,*)'ERROR with MOI Diagonalization'
      goto 1234
      endif
C     Reorient Dipole
      Do j = 1, 3
        Pdipole(j) = 0.00D+00
        Do l = 1, 3
          Pdipole(j) = Pdipole(j)+(PMOI(l,j)*DIPOLE(l))
        ENDDO
      ENDDO
C     Convert DMOI to au
      Do j = 1, 3
        DMOI(j) = DMOI(j)/(au_to_angs**2.00D+00)
      ENDDO
C     End of MOI, Projected Dipole Moment Calculation

C    Bin Data  O-O distance
        Do j = 1, (ndistatom*(ndistatom - 1)/2)
          l = INT((((distance(j)-rmin)/(rmax-rmin))*DBLE(NOB))) + 1
          if (.not.((l .lt. 1).or.(l .gt. NOB))) then
          number(l) = number(l) + weight
          averagedistance = averagedistance + (weight*distance(j))
          endif
        ENDDO

C    Bin Data O--H distance
       Do j = 1, (ndistatom*(ndistatom - 1)*2)
          l = INT((((ohdistance(j)-ohmin)/(ohmax-ohmin))*DBLE(NOB))) + 1
          if (.not.((l .lt. 1).or.(l .gt. NOB))) then
          number2(l) = number2(l) + weight
          averageohdistance = averageohdistance + (weight*ohdistance(j))
          endif
       EndDo

C    Bin Data Dipole interactions
       l = INT((((dipint+0.100D+00)/(0.200D+00))*DBLE(NOB))) + 1
       if (.not.((l .lt. 1).or.(l .gt. NOB))) then
         if (numberhbonds .eq. 0) then
           number6(l) = number6(l) + weight
           averagedipint = averagedipint + (dipint*weight)
         endif
       endif

C    Bin Data 0--H-O angle
       Do j = 1, (ndistatom*(ndistatom - 1)*2)
          l = INT((((ohangle(j)/PI*180.00D+00)/(ohanglemax))
     &        *DBLE(NOB))) + 1
          if (.not.((l .lt. 1).or.(l .gt. NOB))) then
          number4(l) = number4(l) + weight
      averageohangle=averageohangle+(weight*ohangle(j)/PI*180.00D+00)
          endif
       EndDo

C    Bin Data Bisect - O--H angle
       Do j = 1, (ndistatom*(ndistatom - 1)*2)
          l = INT((((hoangle(j)/PI*180.00D+00)/(hoanglemax))
     &        *DBLE(NOB))) + 1
          if (.not.((l .lt. 1).or.(l .gt. NOB))) then
          number7(l) = number7(l) + weight
      averagehoangle=averagehoangle+(weight*hoangle(j)/PI*180.00D+00)
          endif
       EndDo

C    Bin Data Energy Partition
       if (ndistatom .gt. 5) then
         hold = (2*ndistatom) + 2
       else
         hold = ((ndistatom*(ndistatom-1)/2) + 2)
       endif
       l = INT(((((UU*au_to_kcal_mol)-emin)/(emax-emin))*DBLE(NOB))) + 1
       if (.not.((l .lt. 1).or.(l .gt. NOB))) then
       Do j = 1, hold
       if ((j .gt. 1) .and. (numberhbonds+2 .eq. j)) then
       number3(l,j) = number3(l,j) + weight
       averageenergy(j) = averageenergy(j) + (UU*au_to_kcal_mol*weight)
       if (j .eq. hold) then
       averagehbonds = averagehbonds + (weight*DBLE(numberhbonds))
       endif
       else if (j .eq. 1) then
       number3(l,j) = number3(l,j) + weight
       averageenergy(j) = averageenergy(j) + (UU*au_to_kcal_mol*weight)
       else if ((j .eq. hold) .and. (numberhbonds+2 .gt. j)) then
       number3(l,j) = number3(l,j) + weight
       averageenergy(j) = averageenergy(j) + (UU*au_to_kcal_mol*weight)
       averagehbonds = averagehbonds + (weight*DBLE(numberhbonds))
       endif
       EndDo
       endif

C    Bin Data Projected Dipole Moments
      Do j = 1, 3
        l = INT(ABS(PDIPOLE(j))/dipmax*DBLE(NOB)) + 1
        if (.not.((l .lt. 1).or.(l .gt. NOB))) then
        number5(l,j) = number5(l,j) + weight
        averageMOI(j) = averageMOI(j) + (DMOI(j)*weight)
        averageDIPOLE(j) = averageDIPOLE(j) + ABS(PDIPOLE(j)*weight)
        averageDIPOLE2(j) = averageDIPOLE2(j) +
     &   (DSQRT((PDIPOLE(1)**2.0D+00)+(PDIPOLE(2)**2.0D+00)+
     &         (PDIPOLE(3)**2.0D+00)-(PDIPOLE(j)**2.0D+00))*weight)
        endif
      ENDDO

C     Continue with new data GENERAL
      CONTINUE
      goto 5678
C     EOF Condition GENERAL
1234  CONTINUE
      EndDO   

C     Close Read Files
      Do i = 20, 19+nmcchain
        Close(i)
      EndDo    
      Do i = 100, 99+hold
        Close(i)
      EndDo     

C     Calculate Final Distributions
C     Determine number of possible h-bonds + 2
       if (ndistatom .gt. 5) then
         hold = (2*ndistatom) + 2
       else
         hold = ((ndistatom*(ndistatom-1)/2) + 2)
       endif

C     Sum Distributions
        sum  = 0.00D+00
        sum2 = 0.00D+00
      Do i = 1, hold
        sum3(i) = 0.00D+00
      EndDo
        sum4 = 0.00D+00
        sum7 = 0.00D+00
        sum6 = 0.00D+00
      Do i = 1, 3
        sum5(i) = 0.00D+00
      ENDDO
      Do i = 1, nob
          sum  = sum + number(i)
          sum2 = sum2 + number2(i)
        Do j = 1, hold 
          sum3(j) = sum3(j) + number3(i,j)
        ENDDO
          sum4 = sum4 + number4(i)
          sum7 = sum7 + number7(i)
          sum6 = sum6 + number6(i)
        Do j = 1, 3
          sum5(j) = sum5(j) + number5(i,j)
        ENDDO
      EndDo 

C     Normalize Distributions
      Do i = 1, nob
         NDA(i) = number(i)/sum
         NDA2(i) = number2(i)/sum2
         NDA4(i) = number4(i)/sum4
         NDA7(i) = number7(i)/sum7
         NDA6(i) = number6(i)/sum6
         Do j = 1, hold
           NDA3(i,j) = number3(i,j)/sum3(1)
         EndDo
         Do j = 1, 3
           NDA5(i,j) = number5(i,j)/sum5(j)
         ENDDO
      EndDo

C     Find Averages
      averagedistance = averagedistance/sum
      averageohdistance = averageohdistance/sum2
      averageohangle = averageohangle/sum4
      averagehoangle = averagehoangle/sum7
      averagedipint = averagedipint/sum6
      Do i = 1, hold
        averageenergy(i) = averageenergy(i)/sum3(i)
      ENDDO
      averagehbonds = averagehbonds/sum3(hold)
      Do i = 2, hold
        sum3(i) = sum3(i)/sum3(1)
      ENDDO
      Do i = 1, 3
        averageMOI(i) = averageMOI(i)/sum5(i)
        averageDIPOLE(i) = averageDIPOLE(i)/sum5(i)
        averageDIPOLE2(i) = averageDIPOLE2(i)/sum5(i)
      ENDDO

C     Write Distributions
      WRITE(10,*)'***************************************'
      WRITE(10,*)'O-O Distances in Ang.'
      WRITE(10,*)'Thermal average is ',averagedistance
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB(i),NDA(i)
      ENDDO
      WRITE(10,*)'***************************************'
      WRITE(10,*)'***************************************'
      WRITE(10,*)'O--H Distances in Ang.'
      WRITE(10,*)'Thermal average is ',averageohdistance
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB3(i),NDA2(i)
      ENDDO
      WRITE(10,*)'***************************************'
      WRITE(10,*)'***************************************'
      WRITE(10,*)'O--H-O Angle in Degrees'
      WRITE(10,*)'Thermal average is ',averageohangle
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB4(i),NDA4(i)
      ENDDO
      WRITE(10,*)'***************************************'
      WRITE(10,*)'***************************************'
      WRITE(10,*)'Bisect - O--H Angle in Degrees'
      WRITE(10,*)'Thermal average is ',averagehoangle
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB7(i),NDA7(i)
      ENDDO
      WRITE(10,*)'***************************************'
      WRITE(10,*)'***************************************'
      WRITE(10,*)'Energy in kcal/mol'
      WRITE(10,*)'Thermal average is ',averageenergy(1)
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB2(i),NDA3(i,1)
      ENDDO
      WRITE(10,*)'***************************************'
      Do j = 2, hold
      WRITE(10,*)'***************************************'
      WRITE(10,*)'Energy in kcal/mol of configurations'
      WRITE(10,"(A,I3,A)")'with ',(j-2),' hydrogen bonds'
      WRITE(10,*)'Thermal average is',averageenergy(j)
      WRITE(10,*)'Fraction of Total is ',sum3(j)
      if (j .eq. hold) then
      WRITE(10,*)'Average hydrogen bonds is ',averagehbonds
      endif
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB2(i),NDA3(i,j)
      ENDDO
      WRITE(10,*)'***************************************'
      EndDo
      WRITE(10,*)'***************************************'
      WRITE(10,*)'Dipole Int. Energy in Arb. Units'
      WRITE(10,*)'For non bonded monomers'
      WRITE(10,*)'Thermal average is ',averagedipint
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB6(i),NDA6(i)
      ENDDO
      WRITE(10,*)'***************************************'

      Do j = 1, 3
      WRITE(10,*)'***************************************'
      WRITE(10,*)'Principal Dipole Moments in au'
      WRITE(10,*)'related to principal axis ',(j)
      WRITE(10,*)'with an average MOI (au) of ',averageMOI(j)
      WRITE(10,*)'Thermal average Dipole Moment is ',averageDIPOLE(j)
      WRITE(10,*)'Average Perp. Dipole Moment is ',averageDIPOLE2(j)
      WRITE(10,*)'***************************************'
      Do i = 1, nob
      WRITE(10,*)RB5(i),NDA5(i,j)
      ENDDO
      WRITE(10,*)'***************************************'
      EndDo


C     Close Distribution File
      Close(10)

      end

C***********************************************************************
      double precision FUNCTION dntmc_RAN3(IDUM)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
c     N.R. FUNCTION ran3(idum)
c
      implicit none

      integer idum  ! input

      integer mbig, mseed, ma
      integer iff, mj, mk, mz, inext, inextp
      integer i, ii, k
      double precision fac
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      dimension ma(55)
      data iff /0/
      save iff, mj, ma, mk, inext, inextp
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      dntmc_ran3=mj*fac
      return
      end
C***********************************************************************
C End of double precision FUNCTION dntmc_RAN3(IDUM)
C***********************************************************************
