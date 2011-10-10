C***********************************************************************
      LOGICAL FUNCTION  DC_ENERGY(rtdb)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

C This code applies the dangchang potential to WATER ONLY.

      IMPLICIT NONE

#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "dc_common.fh"
#include "util.fh"

      DOUBLE PRECISION cau2kcal

      parameter(cau2kcal = 627.509552D00) !value from /src/drdy/drdyP.fh

      INTEGER rtdb

C      DOUBLE PRECISION R

      DOUBLE PRECISION R3,EFLDMAG

      DOUBLE PRECISION V,EF,EFWX,EFWY,EFWZ,EFIX,EFIY,EFIZ
      DOUBLE PRECISION EFUX,EFUY,EFUZ,EFX,EFY,EFZ,POL,UU,Upol1
      DOUBLE PRECISION RXIIJJ,RYIIJJ,RZIIJJ,RIIJJSQ,RIIJJ
      DOUBLE PRECISION DOTX,DOTY,DOTZ,DOT,UDOT1
      DOUBLE PRECISION Vh2o,DC_VLJONES,DC_VCOUL,DC_EFMAG ! FxNs
      DOUBLE PRECISION DCDIP(3) ! Dipole Moment holder
      DOUBLE PRECISION  NUX,NUY,NUZ,MUX,MUY,MUZ,MX,MY,MZ,MUTOL,MU,DMU
      
C      INTEGER MXMOLS,NMOLS,MXATOM
      INTEGER NSCF,SCF,I,II,J,JJ,COUNT,SUM

      PARAMETER (NSCF = 200)         

      DIMENSION COUNT(MXMOLS)
      DIMENSION V(MXMOLS,MXMOLS,MXATOM,MXATOM)
      DIMENSION EF(MXMOLS,MXMOLS,MXATOM,MXATOM,3)
      DIMENSION EFWX(MXMOLS),EFWY(MXMOLS),EFWZ(MXMOLS)
      DIMENSION EFIX(MXMOLS),EFIY(MXMOLS),EFIZ(MXMOLS)
      DIMENSION EFX(MXMOLS),EFY(MXMOLS),EFZ(MXMOLS)
      DIMENSION NUX(MXMOLS,MXMOLS),NUY(MXMOLS,MXMOLS),
     & NUZ(MXMOLS,MXMOLS)
      DIMENSION MUX(NSCF,MXMOLS),MUY(NSCF,MXMOLS),
     & MUZ(NSCF,MXMOLS),MU(NSCF,MXMOLS)
      DIMENSION MX(MXMOLS),MY(MXMOLS),MZ(MXMOLS)       
      DIMENSION Vh2o(MXATOM,MXATOM)
      DIMENSION R3(MXMOLS,MXMOLS)


      EXTERNAL DC_Initialize

      EXTERNAL dc_setgeom


C     INITIALIZE LOGICAL FxN
      dc_energy = .false.

      CALL DC_Initialize

      CALL dc_setgeom(rtdb)      

      UU = 0.D+00

      Upol1 = 0.D+00

      DO I = 1, NMOLS

         EFWX(I) = 0.D+00
         EFWY(I) = 0.D+00
         EFWZ(I) = 0.D+00

      ENDDO 


     
C THIS SECTION WILL BE ACCESSED FOR EVERY MONTE CARLO STEP
C DUE TO THE FACT THAT ALL MOLECULES ARE MOVED AT ONCE

C     The first part of 1st DO Loop calculates 
C     the initial electric field
C     for every atom


C     O atom = 1st point charge
C     H atom = 2nd-3rd point charges
C     M atom = 4th point charge added
C              to correct the dipole
C     I and J represent molecule of interest
C     II and JJ represent atom of interest

C     4's below on II and JJ loops
C     indicate 4 point charges on each H2O            

      DO I = 1, NMOLS - 1

        DO J = I + 1, NMOLS

            DO II = 1, 4

              DO JJ = 1, 4
           
              RXIIJJ = R(I,II,1) - R(J,JJ,1)
              RYIIJJ = R(I,II,2) - R(J,JJ,2)
              RZIIJJ = R(I,II,3) - R(J,JJ,3)
              RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2
              RIIJJ = DSQRT(RIIJJSQ)


C     if then applies to II being M atom
       
               IF(II.EQ.4)THEN
  
C               WRITE(6,*)'RIIJJ =',RIIJJ,' ANGSTROMS'

                EFUX = RXIIJJ/RIIJJ
               EFUY = RYIIJJ/RIIJJ
               EFUZ = RZIIJJ/RIIJJ
               
               EFLDMAG = DC_EFMAG(J,JJ,RIIJJSQ)
             
               EF(I,J,II,JJ,1) = EFLDMAG*EFUX
               EF(I,J,II,JJ,2) = EFLDMAG*EFUY
               EF(I,J,II,JJ,3) = EFLDMAG*EFUZ
               
               EFWX(I) = EFWX(I) + EF(I,J,II,JJ,1)
               EFWY(I) = EFWY(I) + EF(I,J,II,JJ,2)
               EFWZ(I) = EFWZ(I) + EF(I,J,II,JJ,3)
               
        
c     BELOW if then applies to II and JJ being the M atom                                                                   
                IF((II.EQ.4).AND.(JJ.EQ.4))THEN
                
                  EF(J,I,JJ,II,1) = -EF(I,J,II,JJ,1)
                  EF(J,I,JJ,II,2) = -EF(I,J,II,JJ,2)
                  EF(J,I,JJ,II,3) = -EF(I,J,II,JJ,3)
                  
                  EFWX(J) = EFWX(J) + EF(J,I,JJ,II,1)
                  EFWY(J) = EFWY(J) + EF(J,I,JJ,II,2)
                  EFWZ(J) = EFWZ(J) + EF(J,I,JJ,II,3)

                   NUX(I,J) = RXIIJJ/RIIJJ
                   NUY(I,J) = RYIIJJ/RIIJJ
                   NUZ(I,J) = RZIIJJ/RIIJJ
c234567                   
                   NUX(J,I) = -NUX(I,J)
                   NUY(J,I) = -NUY(I,J)
                   NUZ(J,I) = -NUZ(I,J)

                   R3(I,J) = RIIJJ*RIIJJ*RIIJJ*6.748332695D+00
                   R3(J,I) = R3(I,J)
                   
         ENDIF ! II being M atom

C     BELOW  Else if then applies to II NOT being M atom
C            and JJ being the M atom
        
                ELSEIF(JJ.EQ.4)THEN
                
                RXIIJJ = R(J,JJ,1) - R(I,II,1)
                RYIIJJ = R(J,JJ,2) - R(I,II,2)
                RZIIJJ = R(J,JJ,3) - R(I,II,3)
                RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2
                RIIJJ = DSQRT(RIIJJSQ)

                
C               WRITE(6,*)'RIIJJ =',RIIJJ,' ANGSTROMS'

                EFUX = RXIIJJ/RIIJJ
               EFUY = RYIIJJ/RIIJJ
               EFUZ = RZIIJJ/RIIJJ
               
               EFLDMAG = DC_EFMAG(I,II,RIIJJSQ)
             
               EF(J,I,JJ,II,1) = EFLDMAG*EFUX
               EF(J,I,JJ,II,2) = EFLDMAG*EFUY
               EF(J,I,JJ,II,3) = EFLDMAG*EFUZ
               
               EFWX(J) = EFWX(J) + EF(J,I,JJ,II,1)
               EFWY(J) = EFWY(J) + EF(J,I,JJ,II,2)
               EFWZ(J) = EFWZ(J) + EF(J,I,JJ,II,3)
               
        
               ENDIF ! II NOT being M atom 
                      ! JJ being M atom

C     The second part of 1st DO Loop calculates
C     the pairwise additive potential energy
                                                         
               IF(II.EQ.1.AND.JJ.EQ.1)THEN

                 Vh2o(II,JJ) = DC_Vljones(I,J,II,JJ,RIIJJSQ)                 
c                WRITE(12,*)'Vh2oOO =',Vh2o(II,JJ),' kcal/mole'              
                                               
              ELSEIF(II.EQ.4.AND.JJ.EQ.4)THEN
              
                 Vh2o(II,JJ) = DC_Vcoul(I,J,II,JJ,RIIJJ)
c                WRITE(12,*)'Vh2o =',Vh2o(II,JJ)
                 
               ELSEIF(II.EQ.1.OR.JJ.EQ.1)THEN                               

                 Vh2o(II,JJ) = 0.D+00               
                 
               ELSE
               
                 Vh2o(II,JJ) = DC_Vcoul(I,J,II,JJ,RIIJJ)
c                WRITE(12,*)'Vh2o =',DC_Vcoul(I,J,II,JJ,RIIJJ)              
                                                                            
              ENDIF


               V(I,J,II,JJ) = Vh2o(II,JJ)      

c234567                              
               UU = UU + V(I,J,II,JJ)
    
               END DO

            END DO       
 
         END DO

      END DO ! Ends all 4 loops



c      WRITE(12,*)'UU and DU before induced dipoles are computed'
c      WRITE(12,*)'UU =',UU

C      WRITE(12,*)'Done with pair potentials and point charge'
c      WRITE(12,*)'electric field calculation...'
C      WRITE(12,*)'...'

     
C This section calculates the self-consistent induced dipole moments
C for each molecule and once converged evaluates Upol


      POL = 9.74459258D+00
      MUTOL = 1.D-6

      DO I = 1, NMOLS

         COUNT(I) = 0
         MU(1,I) = 0.D+00
         MUX(1,I) = 0.D+00
         MUY(1,I) = 0.D+00
         MUZ(1,I) = 0.D+00
         EFX(I) = 0.D+00
         EFY(I) = 0.D+00
         EFZ(I) = 0.D+00

      ENDDO


      SUM = 0

C     First part of SCF DO Loop
C     calculates change in electric field
  
      DO SCF = 2, NSCF

        DO I = 1, NMOLS

          DO J = 1, I - 1

c         First dot product of induced dipole and Tensor
c         (MU) dot (unit vector between 2 atoms)

          DOTX = NUX(I,J)*MUX(SCF-1,J)
          DOTY = NUY(I,J)*MUY(SCF-1,J)
          DOTZ = NUZ(I,J)*MUZ(SCF-1,J)

          DOT = DOTX + DOTY + DOTZ

c         Second dot product of induced dipole and Tensor
c         (MU dot NU) dot (NU)

          EFIX(I)=(3.D+00*DOT*NUX(I,J)-MUX(SCF-1,J))/R3(I,J)
          EFIY(I)=(3.D+00*DOT*NUY(I,J)-MUY(SCF-1,J))/R3(I,J)
          EFIZ(I)=(3.D+00*DOT*NUZ(I,J)-MUZ(SCF-1,J))/R3(I,J)

C         Sum all Tensor-induced dipole dot products

          EFX(I) = EFX(I) + EFIX(I)
          EFY(I) = EFY(I) + EFIY(I)
          EFZ(I) = EFZ(I) + EFIZ(I)

          ENDDO ! J Loop

C         The next J DO Loop avoids self interactions

          DO J = I + 1, NMOLS

          DOTX = NUX(I,J)*MUX(SCF-1,J)
          DOTY = NUY(I,J)*MUY(SCF-1,J)
          DOTZ = NUZ(I,J)*MUZ(SCF-1,J)

          DOT = DOTX + DOTY + DOTZ

          EFIX(I)=(3.D+00*DOT*NUX(I,J)-MUX(SCF-1,J))/R3(I,J)
          EFIY(I)=(3.D+00*DOT*NUY(I,J)-MUY(SCF-1,J))/R3(I,J)
          EFIZ(I)=(3.D+00*DOT*NUZ(I,J)-MUZ(SCF-1,J))/R3(I,J)
          
          EFX(I) = EFX(I) + EFIX(I)
          EFY(I) = EFY(I) + EFIY(I)
          EFZ(I) = EFZ(I) + EFIZ(I)
          
          ENDDO ! 2nd J Loop

        ENDDO ! I Loop


C         Second part of SCF DO Loop 
C         Calculates and Saves
C         New Total Electric Field

        DO I = 1, NMOLS

          MUX(SCF,I)=POL*(EFWX(I) + EFX(I))
          MUY(SCF,I)=POL*(EFWY(I) + EFY(I))
          MUZ(SCF,I)=POL*(EFWZ(I) + EFZ(I))
          
          MU(SCF,I)=DSQRT(MUX(SCF,I)**2+MUY(SCF,I)**2+MUZ(SCF,I)**2)

          DMU=DABS(MU(SCF,I)-MU(SCF-1,I))

          IF(DMU.LT.MUTOL)THEN

C      MX(I) is the final self consistent induced dipole moment
C      from the final electric field times the polarizability

            MX(I)= MUX(SCF,I)
            MY(I)= MUY(SCF,I)
            MZ(I)= MUZ(SCF,I)

c            WRITE(12,*)'hey dipole',I,' converged'
c            WRITE(12,*)'after ',SCF,' steps'
c            WRITE(12,*)'************************'
c            WRITE(12,*)'MX(',I,') =',MX(I)
c            WRITE(12,*)'MY(',I,') =',MY(I)
c            WRITE(12,*)'MZ(',I,') =',MZ(I)
c            WRITE(12,*)'|MU| =',X
c            WRITE(12,*)'************************'

            COUNT(I) = 1 ! WHY WOULD HE LET IT EQUAL ONE AGAIN???

C            WRITE(12,*)'COUNT(',I,') =',COUNT(I)            

          ENDIF


        END DO ! 2nd I Loop
 
        DO J = 1, NMOLS

           SUM = SUM + COUNT(J)

        END DO ! J Loop

          
        IF(SUM.EQ.DBLE(NMOLS))THEN

C          WRITE(12,*)'all induced dipoles have converged!'
C          WRITE(12,*)'SUM =',SUM
          
          GOTO 817

        ENDIF 
           
        SUM = 0 ! Not Converged

        DO J = 1, NMOLS

           EFX(J) = 0.D+00 ! Reinitializing
           EFY(J) = 0.D+00
           EFZ(J) = 0.D+00
           COUNT(J) = 0.D+00

        ENDDO ! J Loop


C        WRITE(12,*)'finished with lower ith loop'
C        WRITE(12,*)'going back for another SCF iteration'


        IF(SCF.GE.NSCF)THEN
    
C        WRITE(12,*)'...there is a problem...'
C        WRITE(12,*)'the induced dipoles did not converge'
C        WRITE(12,*)'...stopping...sorry'


        ENDIF

      ENDDO ! SCF Loop

 817  CONTINUE

C******************************************************************
C This loop sums the dot products of the self consistent induced
C dipole moments with the total INITIAL electric field of the charges
C to give the polarization interaction


      DO I = 1, NMOLS
    
         DOTX=MX(I)*EFWX(I)
         DOTY=MY(I)*EFWY(I)
         DOTZ=MZ(I)*EFWZ(I)

         DOT = DOTX + DOTY + DOTZ

         UDOT1 = - 0.5D+00*DOT*627.5095D+00

         Upol1 = Upol1 + UDOT1


      ENDDO ! I Loop
      

      UU = UU + Upol1

      UU = UU / cau2kcal

C     Calculate Total Dipole Moment Components
      CALL DC_DIPOLEMOM(MX,MY,MZ,DCDIP)
      if (util_print('Information',print_low)) then
      WRITE(6,*)'Energy=', UU, 'Hartrees'
      WRITE(6,*)'Dipole Moments in (au)'
      WRITE(6,*)(DCDIP(I),I=1,3)
      endif

C     INITIALIZE LOGICAL FxN
      dc_energy = .true.

      if (.not. rtdb_put(rtdb,'dangchang:energy',mt_dbl,1,UU))
     & call errquit('dc_energy:  no energy in dangchang?', 0,
     & RTDB_ERR)

      if (.not. rtdb_put(rtdb,'dangchang:dipole',mt_dbl,3,DCDIP))
     & call errquit ('dc_energy:  no dipole in dangchang?', 0,
     & RTDB_ERR)

      RETURN
      END

C***********************************************************************
      FUNCTION DC_Vcoul(I,J,II,JJ,RIIJJ)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT NONE

      DOUBLE PRECISION kb,k,e,esqrd,PI,QSQRD,SGMSQ,EPSILON,EPS,Q
      DOUBLE PRECISION RIIJJ,DC_Vcoul
      INTEGER I,J,II,JJ

      COMMON/DCI/kb,k,e,esqrd,PI,Q(4),
     +QSQRD(4,4),SGMSQ(4,4),EPS(4,4)
                      
      DC_Vcoul = (k*QSQRD(II,JJ)*esqrd) / RIIJJ
      
      RETURN
      END 

C***********************************************************************
      FUNCTION DC_Vljones(I,J,II,JJ,RIIJJSQ)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT NONE

      DOUBLE PRECISION SR2,SR6,SR12,RIIJJSQ
      DOUBLE PRECISION DC_Vljones
      DOUBLE PRECISION kb,k,e,esqrd,PI,QSQRD,SGMSQ,EPS
      DOUBLE PRECISION EPSILON,Q
      INTEGER I,J,II,JJ

      COMMON/DCI/kb,k,e,esqrd,PI,Q(4),
     +QSQRD(4,4),SGMSQ(4,4),EPS(4,4)

      SR2 = SGMSQ(II,JJ)/RIIJJSQ
              
      SR6 = SR2*SR2*SR2
              
      SR12 = SR6*SR6
               
      DC_Vljones = 4.0D+00*EPS(II,JJ)*(SR12 - SR6)
     
      RETURN
      END
     
C***********************************************************************
      FUNCTION DC_EFMAG(I,II,RIIJJSQ)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT NONE
      DOUBLE PRECISION DC_EFMAG,RIIJJSQ,Q
      DOUBLE PRECISION  kb,k,e,esqrd,PI,QSQRD,SGMSQ,EPS
      DOUBLE PRECISION  EPSILON
      INTEGER I,J,II,JJ

      COMMON/DCI/kb,k,e,esqrd,PI,Q(4),
     +QSQRD(4,4),SGMSQ(4,4),EPS(4,4)

      DC_EFMAG= k*e*Q(II)/RIIJJSQ
      DC_EFMAG= DC_EFMAG/7.401250726D+21
      
      RETURN
      END         

C***********************************************************************
      SUBROUTINE DC_DIPOLEMOM(MX,MY,MZ,DCDIP)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      IMPLICIT NONE

#include "dc_common.fh"

      DOUBLE PRECISION kb,k,e,esqrd,PI,Q,QSQRD,SGMSQ,EPS
      INTEGER I, J
      DOUBLE PRECISION DCDIP(3)
      DOUBLE PRECISION MX(MXMOLS),MY(MXMOLS),MZ(MXMOLS)

      COMMON/DCI/kb,k,e,esqrd,PI,Q(4),
     +QSQRD(4,4),SGMSQ(4,4),EPS(4,4)

      DO I = 1, 3
          DCDIP(I) = 0.00D+00
      ENDDO
      DO I = 1, NMOLS
          DCDIP(1) = DCDIP(1) + MX(I)
          DCDIP(2) = DCDIP(2) + MY(I)
          DCDIP(3) = DCDIP(3) + MZ(I)
      DO J = 1, 4
          DCDIP(1) = DCDIP(1) + (Q(J)*R(I,J,1)/0.5291772083D+00)
          DCDIP(2) = DCDIP(2) + (Q(J)*R(I,J,2)/0.5291772083D+00)
          DCDIP(3) = DCDIP(3) + (Q(J)*R(I,J,3)/0.5291772083D+00)
      ENDDO
      ENDDO

      RETURN
      END
C***********************************************************************
C  End of SUBROUTINE DC_DIPOLEMOM(MX,MY,MZ,DCDIP)
C***********************************************************************
c $Id$
