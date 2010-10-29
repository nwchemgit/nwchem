C$Id$
      Program Data_Analysis5

      implicit none

      integer nob
      parameter (nob = 100)
      double precision PI, kb, temp  
      integer i, j, z
      double precision sumx
      double precision mtot
      double precision flux(nob), sumerr(nob)
      double precision flux2(nob), sumerr2(nob) 
      double precision NDA(nob), ERR(nob)
      double precision RB(nob), ERRX
      double precision rflux(nob)
      double precision ratec(nob),errrate
      integer nc
      character*256 text
      logical EOF
      double precision intpress, intdens, press
      logical rconfig

C Input Parameters
      PI = 3.141592653589793D+00 ! Constant PI
      kb = 1.987905D-03 ! Constant Boltzmann kcal/mole/K 
      temp = 243.D+00 ! Temperature used in original data 
      mtot = 18.02D+00 ! Mass of monomers
      press = 1.7D-05 
      rconfig = .true. ! Method selection

C Get some input
      WRITE(*,*)'Input parameters'
      WRITE(*,*)'rconfig and temp'
      READ(*,*)rconfig, temp
      WRITE(*,"(A10, A10, A10, A20, A20,
     &          A20, A20)")
     &'N. Conf.','Rcut','Rcut Err.','Ev.Rate','Ev.Rate Err.','Int.Press'
     &,'Density'
      WRITE(*,"(A10, A10, A10, A20, A20,
     &          A20, A20)")
     &' ','(Ang)','(Ang)','(s^-1)','(s^-1)','(atm)'
     &,'(mol/ml)'



C Open Input File
      OPEN(UNIT = 20,FILE = 'data_analysis.inp' 
     +    ,STATUS = 'OLD')

      EOF = .false.
 1234 CONTINUE ! beginning of while loop

C Read Data
      READ(20,*,END=5678)nc
      Do i = 1, nob
        READ(20,*,END=5678)RB(i), NDA(i), ERR(i)
      END DO
      READ(20,*,END=5678)text

      if (rconfig) then
C Calculate Flux
        flux(1) = NDA(1)
        Do i=2, NOB
          flux(i) = NDA(i) + flux(i-1)
        End Do ! CPF in flux
        Do i =1, NOB
          flux2(i) = NDA(i)
        ENDDO ! Derivative CPF in flux2
        sumerr(1) = ERR(1) * ERR(1)
        Do i = 2, NOB
          sumerr(i) = (ERR(i)*ERR(i)) + sumerr(i-1)
        END DO
        Do i = 1, NOB
          sumerr(i) = DSQRT(sumerr(i))
        End DO !Error in CPF in sumerr
        DO i = 1, NOB
          sumerr2(i) = ERR(i)
        ENDDO !Error in Derivative CPF in sumerr2
      else
        Do i = 1, NOB
          flux(i) = NDA(i)/(4.D+00*PI*(RB(i)**2.D+00))
          flux(i) = flux(i)*DEXP((4.D+00/3.D+00*PI*(RB(i)**3.D+00)*press
     &              )/kb/temp)
        ENDDO ! CPF in flux
        Do i = 1, NOB
          if ((i .eq. 1).or.(i .eq. NOB)) then
            flux2(i) = 0.00D+00
          else
            flux2(i) = (flux(i+1) - flux(i-1))/2.00D+00
          endif
        ENDDO ! Derivative CPF in flux2
        Do i = 1, NOB
          sumerr(i) = ERR(i)/(4.D+00*PI*(RB(i)**2.D+00))
          sumerr(i) = sumerr(i)*DEXP((4.D+00/3.D+00*PI*(RB(i)**3.D+00)
     &                *press)/kb/temp)
        ENDDO !Error in CPF in sumerr
        DO i = 1, NOB
          if ((i .eq. 1).or.(i .eq. NOB)) then
            sumerr2(i) = 0.00D+00
          else
            sumerr2(i) = (sumerr(i+1)*sumerr(i+1))
     &                  +(sumerr(i-1)*sumerr(i-1))
            sumerr2(i) = DSQRT(sumerr2(i))/2.00D+00
          endif
        ENDDO !Error in Derivative CPF in sumerr2
      endif

C Calculate Error
      ERRX = 0.00D+00
      Do i=1, NOB - 1
        ERRX = ERRX + RB(i+1) - RB(i)
      End Do
      ERRX = ERRX/DBLE(NOB - 1) ! Error related to bin interval
C Reactive Flux
      Do i=1, NOB
        rflux(i) = flux2(i)*kb*temp/ERRX
      enddo ! Reactive flux in rflux
C Rate Constants
      Do i=1, NOB
        if(flux(i) .eq. 0.D+00) then
        ratec(i) = 1.D+99
        else
        ratec(i) = 1.D+00*rflux(i)/flux(i)*1.D+13/
     &  DSQRT(2.D+00*PI*kb*temp*mtot/4.184D+00)
        endif
      Enddo ! Rate Constants in ratec

C Find Maximum Flux
      sumx = rflux(1)
      j = 1
      Do i = 2, NOB
        IF (rflux(i) .gt. sumx) then
          sumx = rflux(i)
          j = i
        EndIf
      EndDo
      if (j .eq. NOB) j = NOB -1
      if (j .eq. 1) j = 2
      z = j

C Find Minimum Flux
      j = z
      if (j .gt. NOB/2) then

      DO i = z , 2, -1 
        IF ( (rflux(i) .lt. rflux(i-1)) .and.
     &       (rflux(i) .le. rflux(i+1)) .and.
     &       (rflux(i) .lt. sumx)       .and.
     &       (rflux(i+1) .gt. 0.00D+00) ) then
          sumx = rflux(i)
          j = i
        endif
      End Do

      else 

      DO i = z , NOB -1
        IF ( (rflux(i) .lt. rflux(i-1)) .and.
     &       (rflux(i) .le. rflux(i+1)) .and.
     &       (rflux(i) .lt. sumx)       .and.
     &       (rflux(i+1) .gt. 0.00D+00) ) then
          sumx = rflux(i)
          j = i
        endif
      End Do

      endif

C Minimum flux index in j

      errrate = ratec(j)*DSQRT((sumerr(j)/flux(j))**2.D+00 
     &          + (flux2(j)/sumerr2(j))**2.D+00)
      intpress = rflux(j)/flux(j)/4/PI/(RB(j)**2.00D+00)
      intpress = intpress*kb*temp
      intpress = intpress*4.184*1.00D+30/1.01325D+02/6.022142D+23
      intdens = 4.00D+00/3.00D+00*PI*(RB(j)**3.00D+00)
      intdens = 1.00D+00/intdens*1.00D+24/6.022142D+23

      WRITE(*,"(I10, F10.5, F10.5, ES20.10E3, ES20.10E3, 
     &          ES20.10E3, ES20.10E3)")
     & nc, RB(j), ERRX, ratec(j), errrate, intpress, intdens
    
      goto 9098 ! skip EOF
 5678 EOF = .true.
 9098 CONTINUE ! skipped EOF
      if (.not. EOF) goto 1234 ! begin while loop again
C End of while Loop
C close file
      close(20)

      end
