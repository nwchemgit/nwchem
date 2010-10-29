      SUBROUTINE VSINT(M,N,X,XT,MDIMX,WSAVE)
*
* $Id$
*
      implicit double precision (a-h, o-z)
C***BEGIN PROLOGUE  VSINT
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
c***KEYWORDS  FAST FOURIER TRANSFORM, SINE TRANSFORM, MULTIPLE
C             SEQUENCES
C***AUTHOR  BOISVERT, R. F., (NIST)
C***PURPOSE  Sine transform of one or more real, odd sequences.
C***DESCRIPTION
C
C  Subroutine VSINT computes the discrete Fourier sine transform
C  of M odd sequences X(J,I), J=1,...,M.  The transform is defined
C  below at output parameter X.
C
C  The array WSAVE which is used by subroutine VSINT must be
C  initialized by calling subroutine VSINTI(N,WSAVE).
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is the product of small primes.
C
C  X       an array of size at least X(MDIMX,N+1) which contains the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of X.  Thus, the Jth sequence is stored in
C          X(J,I), I=1,..,N.  The extra column of X is used as work
C          storage.
C
C  XT      a work array of size at least XT(MDIMX,N+1).
C
C  MDIMX   the first dimension of the array X exactly as it appears in
C          the calling program.
C
C  WSAVE   a work array with dimension at least INT(2.5*N+15)
C          in the program that calls VSINT.  The WSAVE array must be
C          initialized by calling subroutine VSINTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.
C
C  Output Parameters
C
C  X       for I=1,...,N and J=1,...,M
C
C               X(J,I)= the sum from k=1 to k=N
C
C                    2*X(J,K)*SIN(K*I*PI/(N+1))/SQRT(2*(N+1))
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of VSINT.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VSINT followed immediately by another call
C           of VSINT will return the original sequences X.  Thus,
C           VSINT is the correctly normalized inverse of itself.
C
C  -----------------------------------------------------------------
C
C  VSINT is a straightforward extension of the subprogram SINT to
C  handle M simultaneous sequences.  The scaling of the sequences
C  computed by VSINT is different than that of SINT.  SINT was
C  originally developed by P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFFTF
C***END PROLOGUE  VSINT
      DIMENSION       X(MDIMX,*), XT(MDIMX,*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINT
      IF (M .LE. 0)  GO TO 900
      IF (N .LE. 1)  GO TO 900
      IF (N .GT. 2)  GO TO 300
C
C  CASE   N = 2
C
      SQRTH = SQRT(0.50d0)
      DO 201 J=1,M
         XH = SQRTH*(X(J,1)+X(J,2))
         X(J,2) = SQRTH*(X(J,1)-X(J,2))
         X(J,1) = XH
  201 CONTINUE
      GO TO 900
C
C  CASE   N .GT. 2
C
C     ... PREPROCESSING
C
  300 CONTINUE
      NP1 = N+1
      NS2 = N/2
      DO 301 J=1,M
         XT(J,1) = 0.0d0
  301 CONTINUE
      DO 310 K=1,NS2
         KC = NP1-K
         DO 310 J=1,M
            T1 = X(J,K)-X(J,KC)
            T2 = WSAVE(K)*(X(J,K)+X(J,KC))
            XT(J,K+1) = T1+T2
            XT(J,KC+1) = T2-T1
  310 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) THEN
         DO 320 J=1,M
            XT(J,NS2+2) = 4.0d0*X(J,NS2+1)
  320    CONTINUE
      ENDIF
C
C     ... REAL PERIODIC TRANSFORM
C
      NF = NS2+1
      CALL VRFFTF(M,NP1,XT,X,MDIMX,WSAVE(NF))
C
C     ... POSTPROCESSING
C
      DO 330 J=1,M
         X(J,1) = 0.5d0*XT(J,1)
  330 CONTINUE
      DO 350 I=3,N,2
         DO 340 J=1,M
            X(J,I-1) = -XT(J,I)
  340    CONTINUE
         DO 345 J=1,M
            X(J,I) = X(J,I-2)+XT(J,I-1)
  345    CONTINUE
  350 CONTINUE
      IF (MODN .EQ. 0) THEN
         DO 360 J=1,M
            X(J,N) = -XT(J,N+1)
  360    CONTINUE
      ENDIF
C
C     ... NORMALIZATION
C
      SCALE = SQRT(0.5d0)
      DO 370 I=1,N
         DO 370 J=1,M
            X(J,I) = SCALE*X(J,I)
  370 CONTINUE
C
C  EXIT
C
  900 CONTINUE
      RETURN
      END
      FUNCTION PIMACH(DUM)
      implicit double precision (a-h, o-z)
C***BEGIN PROLOGUE  PIMACH
C
C     This subprogram supplies the value of the constant PI correct to
C     machine precision where
C
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  PIMACH
C
C***FIRST EXECUTABLE STATEMENT  PIMACH
      PIMACH=3.1415926535897932d0
      RETURN
      END
      SUBROUTINE VRADF2 (MP,IDO,L1,CC,CH,MDIMC,WA1)
      implicit double precision (a-h, o-z)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION   CH(MDIMC,IDO,2,L1)  ,CC(MDIMC,IDO,L1,2)     ,
     1                WA1(IDO)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+CC(M,1,K,2)
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I,1,K) = CC(M,I,K,1)+(WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))
            CH(M,IC,2,K) = (WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-CC(M,I,K,1)
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)-(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         DO 1006 M=1,MP
         CH(M,1,2,K) = -CC(M,IDO,K,2)
         CH(M,IDO,1,K) = CC(M,IDO,K,1)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADF3 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2)
      implicit double precision (a-h, o-z)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION   CH(MDIMC,IDO,3,L1)  ,CC(MDIMC,IDO,L1,3)     ,
     1                WA1(IDO)     ,WA2(IDO)
      ARG=2.d0*PIMACH(1.0d0)/3.d0
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,2)+CC(M,1,K,3))
         CH(M,1,3,K) = TAUI*(CC(M,1,K,3)-CC(M,1,K,2))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TAUR*
     1      (CC(M,1,K,2)+CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))+(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))-(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))
            CH(M,IC,2,K) = (TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))-(CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADF4 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3)
      implicit double precision (a-h, o-z)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION    CC(MDIMC,IDO,L1,4)   ,CH(MDIMC,IDO,4,L1)     ,
     1                WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)
      HSQT2=SQRT(2.d0)/2.d0
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = (CC(M,1,K,2)+CC(M,1,K,4))
     1      +(CC(M,1,K,1)+CC(M,1,K,3))
         CH(M,IDO,4,K) = (CC(M,1,K,1)+CC(M,1,K,3))
     1      -(CC(M,1,K,2)+CC(M,1,K,4))
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,3)
         CH(M,1,3,K) = CC(M,1,K,4)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I-1,1,K) = ((WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4)))+(CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+
     1       WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,4,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))-(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I,3,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,2,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         DO 1006 M=1,MP
            CH(M,IDO,1,K) = (HSQT2*(CC(M,IDO,K,2)-CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,1)
            CH(M,IDO,3,K) = CC(M,IDO,K,1)-(HSQT2*(CC(M,IDO,K,2)-
     1       CC(M,IDO,K,4)))
            CH(M,1,2,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))-
     1       CC(M,IDO,K,3)
            CH(M,1,4,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,3)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADF5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
      implicit double precision (a-h, o-z)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,L1,5)    ,CH(MDIMC,IDO,5,L1)     ,
     1           WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.d0*PIMACH(1.0d0)/5.d0
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.d0*ARG)
      TI12=SIN(2.d0*ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,5)+CC(M,1,K,2))+
     1    (CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TR11*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR12*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,3,K) = TI11*(CC(M,1,K,5)-CC(M,1,K,2))+TI12*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
         CH(M,IDO,4,K) = CC(M,1,K,1)+TR12*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR11*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,5,K) = TI12*(CC(M,1,K,5)-CC(M,1,K,2))-TI11*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5)))+((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I-1,3,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1      ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1       +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))+TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4)))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1     ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1      +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))-(TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,2,K) = (TI11*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I-1,5,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))+(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))-(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I,5,K) = (CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,4,K) = (TI12*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADFG (MP,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,MDIMC,WA)
      implicit double precision (a-h, o-z)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION     CH(MDIMC,IDO,L1,IP)   ,CC(MDIMC,IDO,IP,L1)  ,
     1            C1(MDIMC,IDO,L1,IP)    ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMC,IDL1,IP)           ,WA(IDO)
      TPI=2.d0*PIMACH(1.0d0)
      ARG = TPI/DBLE(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         DO 1001 M=1,MP
         CH2(M,IK,1) = C2(M,IK,1)
 1001    CONTINUE
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            DO 1002 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               DO 1004 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1004          CONTINUE
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1008 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1008          CONTINUE
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               DO 1012 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1012          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               DO 1016 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1016          CONTINUE
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         DO 1020 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1020    CONTINUE
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            DO 1022 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)+CH(M,1,K,JC)
            C1(M,1,K,JC) = CH(M,1,K,JC)-CH(M,1,K,J)
 1022       CONTINUE
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.d0
      AI1 = 0.d0
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            DO 1024 M=1,MP
            CH2(M,IK,L) = C2(M,IK,1)+AR1*C2(M,IK,2)
            CH2(M,IK,LC) = AI1*C2(M,IK,IP)
 1024       CONTINUE
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               DO 1025 M=1,MP
               CH2(M,IK,L) = CH2(M,IK,L)+AR2*C2(M,IK,J)
               CH2(M,IK,LC) = CH2(M,IK,LC)+AI2*C2(M,IK,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            DO 1028 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+C2(M,IK,J)
 1028       CONTINUE
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            DO 1030 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1030       CONTINUE
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            DO 1033 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1033       CONTINUE
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            DO 1036 M=1,MP
            CC(M,IDO,J2-2,K) = CH(M,1,K,J)
            CC(M,1,J2-1,K) = CH(M,1,K,JC)
 1036       CONTINUE
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               DO 1038 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1038          CONTINUE
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               DO 1042 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1042          CONTINUE
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
      SUBROUTINE VRFFTF (M,N,R,RT,MDIMR,WSAVE)
      implicit double precision (a-h, o-z)
C***BEGIN PROLOGUE  VRFFTF
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, 
C             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Forward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTF computes the Fourier coefficients (forward 
C  transform) of a number of real periodic sequences.  Specifically,
C  for each sequence the subroutine claculates the independent
C  Fourier coefficients described below at output parameter R.
C
C  The array WSAVE which is used by subroutine VRFFTF must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes,
C          however n may be any positive integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of R.  Thus, the I-th sequence to be transformed,
C          X(I,J), J=0,1,...,N-1, is stored as
C
C               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly 
C          as they appear in the calling program.  This parameter is 
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by 
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTF and VRFFTB.
C
C  Output Parameters
C
C  R       contains the Fourier coefficients F(K) for each of the M 
C          input sequences.  Specifically, row I of R, R(I,J), 
C          J=1,2,..,N, contains the independent Fourier coefficients
C          F(I,K), for the I-th input sequence stored as
C
C             R(I,1) = DBLE( F(I,0) ),
C                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
C
C             R(I,2*K) = DBLE( F(I,K) )
C                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
C
C             R(I,2*K+1) = IMAG( F(I,K) )
C                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
C
C                   for K = 1, 2, . . . , M-1,
C
C              and, when N is even,
C
C              R(I,N) = DBLE( F(I,N/2) ).
C                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTF is a straightforward extension of the subprogram RFFTF to
C  handle M simultaneous sequences.  RFFTF was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTF1
C***END PROLOGUE  VRFFTF
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       R(MDIMR,N)  ,RT(MDIMR,N)    ,WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTF
      IF (N .EQ. 1) RETURN
      CALL VRFTF1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFTF1 (M,N,C,CH,MDIMC,WA,FAC)
      implicit double precision (a-h, o-z)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       CH(MDIMC,N) ,C(MDIMC,N)  ,WA(N)   ,FAC(15)
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADF4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL VRADF4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL VRADF2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 110
  103    CALL VRADF2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL VRADF3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  105    CALL VRADF3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
      CALL VRADF5(M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107 CALL VRADF5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL VRADFG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         NA = 1
         GO TO 110
  109    CALL VRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SCALE=SQRT(1.0d0/N)
      IF (NA .EQ. 1) GO TO 113
      DO 112 J=1,N
      DO 112 I=1,M
         C(I,J) = SCALE*CH(I,J)
  112 CONTINUE
      RETURN
  113 DO 114 J=1,N
      DO 114 I=1,M
         C(I,J)=SCALE*C(I,J)
  114 CONTINUE
      RETURN
      END


      SUBROUTINE VSINTI(N,WSAVE)
      implicit double precision (a-h, o-z)
C***BEGIN PROLOGUE  VSINTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A3
c***KEYWORDS  FAST FOURIER TRANSFORM, SINE TRANSFORM, MULTIPLE
C             SEQUENCES
C***AUTHOR  BOISVERT, R. F. (NIST)
C***PURPOSE  Initialize for VSINT.
C***DESCRIPTION
C
C  Subroutine VSINTI initializes the array WSAVE which is used in
C  subroutine SINT.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array with at least INT(2.5*N+15) locations.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of VSINT.
C
C  -----------------------------------------------------------------
C
C  VSINTI is a straightforward extension of the subprogram SINTI to
C  handle M simultaneous sequences.  SINTI was originally developed
C  P. N. Swarztrauber of NCAR.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFFTI
C***END PROLOGUE  VSINTI
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINTI
      PI = PIMACH(1.0d0)
      IF (N .LE. 1) RETURN
      NP1 = N+1
      NS2 = N/2
      DT = PI/DBLE(NP1)
      KS = 1
      KF = KS+NS2-1
      FK = 0.0d0
      DO 101 K=KS,KF
         FK = FK+1.d0
         WSAVE(K) = 2.d0*SIN(FK*DT)
  101 CONTINUE
      CALL VRFFTI (NP1,WSAVE(KF+1))
      RETURN
      END
      SUBROUTINE VRFFTI (N,WSAVE)
      implicit double precision (a-h, o-z)
C***BEGIN PROLOGUE  VRFFTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Initialization for VRFFTF and VRFFTB.
C***DESCRIPTION
C
C  Subroutine VRFFTI initializes the array WSAVE which is used in
C  both VRFFTF and VRFFTB.  The prime factorization of N together with
C  a tabulation of certain trigonometric functions are computed and
C  stored in the array WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  There is no
C          restriction on N.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least N+15.
C          The same work array can be used for both VRFFTF and VRFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of VRFFTF or VRFFTB.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTI1
C***END PROLOGUE  VRFFTI
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTI
      IF (N .EQ. 1) RETURN
      CALL VRFTI1 (N,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFTI1 (N,WA,FAC)
      implicit double precision (a-h, o-z)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       WA(N)      ,FAC(15)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      TPI = 2.d0*PIMACH(1.0d0)
      ARGH = TPI/DBLE(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = DBLE(LD)*ARGH
            FI = 0.d0
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.d0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END



