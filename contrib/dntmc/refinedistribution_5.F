C$Id$
      Program refinedistribution_5

      Implicit none

      double precision kb, temp, tempinit
      double precision au_to_kcal_mol
      double precision rmin, rmax
      integer natom, nmcchain

      integer nob
      parameter(nob=100)

      integer mxatom
      parameter(mxatom=100)

      integer mxmcchain
      parameter(mxmcchain = 100)

      double precision Number(100,mxmcchain)
      double precision CVnumber(100,3)
      double precision RB(NOB), NDA(NOB), sum(mxmcchain)
      double precision UU(mxmcchain,2), RCONFIG(mxmcchain,2)
      double precision DIPOLE(mxmcchain,2,3)
      character*256 text
      character*256 filename
      character*4 prcfil
      character*256 suffix
      double precision Coord(mxmcchain,2,mxatom,3)

      integer i,j,IR
      integer R

      double precision TMP

C     Input Parameters
      kb = 1.987905D-03 ! Constant Boltzmann kcal/mole/K
      au_to_kcal_mol = 627.5095D+00 ! Conversion from au to kcal/mole

      Write(*,*)'Please input file name prefix:'
      Write(*,*)'Init. Temp  and Final Temp' !inital and final temperature (K) for scaled distributions
      Write(*,*)'rmin and rmax' !grid min and max for distribution (Ang.)
      Write(*,*)'natoms and nmcchains' !number of atoms per configuration and markov chains (files)
      Read(*,*)filename
      Read(*,*)tempinit, temp
      Read(*,*)rmin, rmax
      Read(*,*)natom, nmcchain

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

        write(suffix, '(a,i3,a)')'.fdist.output.',int(temp),'K'
        Open(Unit = 10,FILE = filename(1:len_trim(filename))
     &     //suffix(1:len_trim(suffix))
     &     ,STATUS = 'REPLACE')

C     Setup Grid
      RB(1) = rmin + (rmax - rmin)/DBLE(NOB)
      DO i = 2, nob
        RB(i) = RB(i-1) + (rmax - rmin)/DBLE(NOB)
      ENDDO        

C     Initialize Bins
      Do i = 1, NOB
        Do j = 1, nmcchain
       number(i,j) = 0.00D+00
        EndDo
      EndDo
     
C     Initialize MarkovChains
      Do i = 20, 19+nmcchain
        Do j = 1, natom 
          Read(i,*,END=1234)text,(Coord(i-19,1,j,IR),IR=1,3)
C          Read(i,*,END=1234)(Coord(i-19,1,j,IR),IR=1,3)
        EndDo
        Read(i,*,END=1234)text,text,UU(i-19,1)
        Read(i,*,END=1234)text,text,(DIPOLE(i-19,1,IR),IR=1,3)
        Read(i,*,END=1234)text,text,RCONFIG(i-19,1)
      EndDo

C     Bin Initial Values
      Do i = 1, nmcchain
        j = INT((((RCONFIG(i,1)-rmin)/(rmax-rmin))*DBLE(NOB))) + 1
        if (j .le. 0) j = 1
        if (j .gt. nob) goto 5678
        TMP = (UU(i,1)-UU(i,1))*au_to_kcal_mol/kb*
     &        ((1.D+00/temp)-(1.D+00/tempinit))
        TMP = DEXP(-TMP)
        number(j,i) = number(j,i) + TMP
      End Do

C     Begin a While Loop
      R = 1
9012  Continue

C     Read In Data New
      Do i = 20, 19+nmcchain
          Do j = 1, natom
C            Read(i,*,END=1234)(Coord(i-19,2,j,IR),IR=1,3)
            Read(i,*,END=1234)text,(Coord(i-19,2,j,IR),IR=1,3)
          EndDo
          Read(i,*,END=1234)text,text,UU(i-19,2)
          Read(i,*,END=1234)text,text,(DIPOLE(i-19,2,IR),IR=1,3)
          Read(i,*,END=1234)text,text,RCONFIG(i-19,2)
      EndDo

C     Bin Values
      Do i = 1, nmcchain
        j = INT((((RCONFIG(i,2)-rmin)/(rmax-rmin))*DBLE(NOB))) + 1
        if (j .le. 0) j = 1
        if (j .gt. nob) goto 5678
        TMP = (UU(i,2)-UU(i,1))*au_to_kcal_mol/kb*
     &        ((1.D+00/temp)-(1.D+00/tempinit))
        TMP = DEXP(-TMP)
        number(j,i) = number(j,i)+ TMP
      EndDo

C     Increment Counter
      R = R + 1

C     Continue Through While Loop
      goto 9012

C     ERROR termination
5678  CONTINUE
      WRITE(*,*)'Rconfig out of range. Check rmin and rmax'
C     EOF termination
1234  CONTINUE

C     Write final Distribution
C     Normalize Individual Markov Chains     
      Do j = 1, nmcchain
        sum(j) = 0.00D+00
      EndDo

      Do j = 1, nmcchain
        Do i = 1, nob
          sum(j) = sum(j) + number(i,j)
        EndDo
      EndDo
      
      Do j = 1, nmcchain
        Do i =1, nob
          number(i,j) = number(i,j)/sum(j)
        EndDo
      EndDo

C     Average Markov Chains into Distribution
      Do i = 1, nob
        NDA(i) = 0.00D+00
      EndDo

      Do i = 1, nob
        Do j = 1, nmcchain
          NDA(i) = NDA(i) + number(i,j)
        EndDo
      EndDo

C     Normalize Average
      sum(1) = 0.00D+00
      Do i = 1, nob
        sum(1) = sum(1) + NDA(i)
      EndDo

      Do i = 1, nob
        NDA(i) = NDA(i)/sum(1)
      EndDo

      WRITE(10,*)(R*nmcchain), 'Total Configurations at',temp
      Do i = 1, nob
        WRITE(10,*)RB(I),NDA(I),0.00D+00
      EndDo
      WRITE(10,*)'***************************************'

C     Close Files      
      Do i = 20, 19+nmcchain
        close(i)
      End do

C     Close Major Output
      close(10)

C     End Program
      end

