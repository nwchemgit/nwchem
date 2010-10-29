      SUBROUTINE stpr_setmass(AMS)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / MASS / TOTM, NUMAS
      DIMENSION AMS(NUMAS)
      AMUAU = (1.6605656D0/9.109534D0)*1.0D+04
C
C     Set up Masses.
C
C     From CRC Handbook of Chemistry and Physics 65th edition
C     1984-1985 Pages B-236 to B-252.
C
C     The isotope with the most natural abundance is used.
C
      AMS( 1)   =   1.007825D00
      AMS( 2)   =   4.00260D00
      AMS( 3)   =   7.01600D00
      AMS( 4)   =   9.01218D00
      AMS( 5)   =  11.00931D00
      AMS( 6)   =  12.00000D00
      AMS( 7)   =  14.00307D00
      AMS( 8)   =  15.99491D00
      AMS( 9)   =  18.99840D00
      AMS(10)   =  19.99244D00
      AMS(11)   =  22.9898D00
      AMS(12)   =  23.98504D00
      AMS(13)   =  26.98153D00
      AMS(14)   =  27.97693D00
      AMS(15)   =  30.97376D00
      AMS(16)   =  31.97207D00
      AMS(17)   =  34.96885D00
      AMS(18)   =  39.948D00
      AMS(19)   =  38.96371D00
      AMS(20)   =  39.96259D00
      AMS(21)   =  44.95592D00
      AMS(22)   =  47.9D00
      AMS(23)   =  50.9440D00
      AMS(24)   =  51.9405D00
      AMS(25)   =  54.9380D00
      AMS(26)   =  55.9349D00
      AMS(27)   =  58.9332D00
      AMS(28)   =  57.9353D00
      AMS(29)   =  62.9296D00
      AMS(30)   =  63.9291D00
      AMS(31)   =  68.9257D00
      AMS(32)   =  73.9219D00
      AMS(33)   =  74.9216D00
      AMS(34)   =  79.9165D00
      AMS(35)   =  78.9183D00
      AMS(36)   =  83.80D00
C
      DO 1 I=1,NUMAS
        AMS(I)= AMS(I)*AMUAU
    1 CONTINUE
      RETURN
      END
