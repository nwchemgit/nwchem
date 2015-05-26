      SUBROUTINE vib_setup
* $Id$
      IMPLICIT NONE! REAL*8 (A-H,O-Z)
c
#include "util_params.fh"
      double precision ams, wave
      integer numans
      COMMON /cvib_SETCON/ AMS(36),WAVE,NUMANS      ! setup parameters
      double precision d1, d2
      PARAMETER(D1=1.0D00,D2=2.0D00)
c
      double precision e, a0, c, e0, xl, pi, fact1
C
C This subroutine sets up constants and masses for use by the other
C      routines
C
C
C Set up constants
C
      E   =   1.6021892D00         ! unit charge
      A0  =   cau2ang              ! bohr radius
      C   =   2.9979245D00         ! speed of light
      E0  =   8.85418782D00        ! permittivity of vacuum
      XL  =   6.022045D00          ! Avagadro's number
      PI  =   D2*ACOS(0.0D00)      ! PI
C
      FACT1 = E*E*1.0D02/(D2*D2*PI*E0*A0)  ! HARTREE to ERG
      FACT1 = FACT1*XL/(A0*A0)             ! hartree/bohr**2/amu to sec**-2
      FACT1 = SQRT(FACT1)                  ! to sec**-1
      FACT1 = FACT1/(D2*PI*C)              ! to wavenumbers
      WAVE  = FACT1*1.0D04                 ! include all powers of ten

C
C Set up Masses
C
C From CRC Handbook of Chemistry and Physics 65th edition 1984-1985
C      Pages B-236 to B-252
C
C      The isotope with the most natural abundance is used.
C
      NUMANS    =  36                ! number of atoms with default masses
      AMS( 1)   =   1.007825D00      ! H
      AMS( 2)   =   4.00260D00       ! He
      AMS( 3)   =   7.01600D00       ! Li
      AMS( 4)   =   9.01218D00       ! Be
      AMS( 5)   =  11.00931D00       ! B
      AMS( 6)   =  12.00000D00       ! C
      AMS( 7)   =  14.00307D00       ! N
      AMS( 8)   =  15.99491D00       ! O
      AMS( 9)   =  18.99840D00       ! F
      AMS(10)   =  19.99244D00       ! Ne
      AMS(11)   =  22.9898D00        ! Na
      AMS(12)   =  23.98504D00       ! Mg
      AMS(13)   =  26.98153D00       ! Al
      AMS(14)   =  27.97693D00       ! Si
      AMS(15)   =  30.97376D00       ! P
      AMS(16)   =  31.97207D00       ! S
      AMS(17)   =  34.96885D00       ! Cl
      AMS(18)   =  39.948D00         ! Ar
      AMS(19)   =  38.96371D00       ! K
      AMS(20)   =  39.96259D00       ! Ca
      AMS(21)   =  44.95592D00       ! Sc
      AMS(22)   =  47.9D00           ! Ti
      AMS(23)   =  50.9440D00        ! V
      AMS(24)   =  51.9405D00        ! Cr
      AMS(25)   =  54.9380D00        ! Mn
      AMS(26)   =  55.9349D00        ! Fe
      AMS(27)   =  58.9332D00        ! Co
      AMS(28)   =  57.9353D00        ! Ni *
      AMS(29)   =  62.9296D00        ! Cu *
      AMS(30)   =  63.9291D00        ! Zn *
      AMS(31)   =  68.9257D00        ! Ga *
      AMS(32)   =  73.9219D00        ! Ge *
      AMS(33)   =  74.9216D00        ! As
      AMS(34)   =  79.9165D00        ! Se *
      AMS(35)   =  78.9183D00        ! Br *
      AMS(36)   =  83.80D00          ! Kr *
C
      END
