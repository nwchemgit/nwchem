      FUNCTION FAC_IN_BCH(NT,IT1,IT2,IT3,IT4)
*
* Determine factor for 1-4 operators in 
* BCH expansion of commutating operators. Permutation factors 
* are included, so only one term for combination is assumed 
* calculated  
*
* Jeppe Olsen, April 2003
*              Modified (bug for NT = 4) found Nov. 2004
*
      INCLUDE 'implicit.inc'
*
      IF(NT.EQ.0) THEN
        FAC = 1.0D0
      ELSE IF (NT.EQ.1) THEN
        FAC = 1.0D0
      ELSE IF (NT. EQ. 2) THEN
        IF(IT1.NE.IT2) THEN
          FAC = 1.0D0
        ELSE 
          FAC = 0.5D0
        END IF
      ELSE IF (NT. EQ. 3) THEN
        IF(IT1.NE.IT2.AND.IT2.NE.IT3.AND.IT1.NE.IT3) THEN
          FAC = 1.0D0
        ELSE IF (IT1.EQ.IT2.AND.IT2.EQ.IT3) THEN 
          FAC = 1.0D0/6.0D0
        ELSE 
          FAC = 1.0D0/2.0D0
        END IF
      ELSE IF (NT.EQ.4) THEN
        IF(IT1.NE.IT2.AND.IT1.NE.IT3.AND.IT1.NE.IT4.AND.
     &     IT2.NE.IT3.AND.IT2.NE.IT4.AND.IT3.NE.IT4     )THEN
*. Four different
          FAC = 1.0D0
        ELSE IF(IT1.EQ.IT2.AND.IT2.EQ.IT3.AND.IT3.EQ.IT4) THEN
*. Four identical 
          FAC = 1.0D0/24.0D0
        ELSE IF((IT1.EQ.IT2.AND.IT2.EQ.IT3) .OR. 
     &          (IT1.EQ.IT2.AND.IT2.EQ.IT4) .OR.
     &          (IT1.EQ.IT3.AND.IT3.EQ.IT4) .OR.
     &          (IT2.EQ.IT3.AND.IT3.EQ.IT4)     ) THEN
*. Three identical 
         FAC =  1.0D0/6.0D0
        ELSE IF ((IT1.EQ.IT2.AND.IT3.EQ.IT4) .OR. 
     &           (IT1.EQ.IT3.AND.IT2.EQ.IT4) .OR.
     &           (IT1.EQ.IT4.AND.IT2.EQ.IT3)    )  THEN
*. Two pairs of identical numbers pairs 
         FAC = 1.0D0/4.0D0
        ELSE 
*. Two identical numbers, the last two differs 
C        FAC = 1.0D0/6.0D0
         FAC = 1.0D0/2.0D0
        END IF
      END IF
*     ^ End of switch between different number of operators
      FAC_IN_BCH = FAC
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        IF(NT.EQ.0) THEN 
          WRITE(6,*) ' No T-ops, FAC == ', FAC
        ELSE IF (NT.EQ.1) THEN
          WRITE(6,*) ' One T-op, FAC = ', FAC
        ELSE IF (NT.EQ.2) THEN
          WRITE(6,*) ' T1, T2 , FAC = ', IT1,IT2, FAC
        ELSE IF (NT.EQ.3) THEN
          WRITE(6,*) 'T1, T2, T3, FAC = ', IT1,IT2,IT3, FAC
        ELSE IF (NT.EQ.4) THEN
          WRITE(6,*) ' T1, T2, T3, T4, FAC = ',  
     &                IT1,IT2,IT3,IT4,FAC
        END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE REO_CAABS(IN1,IN2,IN3,IN4,NIN,IOUT1,IOUT2,IOUT3,IOUT4,
     &                     INSM1,INSM2,INSM3,INSM4,
     &                     IOUTSM1,IOUTSM2,IOUTSM3,IOUTSM4,
     &                     IREO, NGAS)
*
* Reorder 0-4 CAAB operators
*
* Jeppe Olsen, April 2003
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER IN1(*),IN2(*),IN3(*),IN4(*)
*. New => old ( reordered => original)
      INTEGER IREO(*)
*. Output 
      INTEGER IOUT1(*),IOUT2(*),IOUT3(*),IOUT4(*)
*. Local scratch
      INTEGER ICAAB(4*MXPNGAS)
*
      ISMX = -1
      DO IOP = 1, NIN
        IF(IOP.EQ.1) THEN 
           CALL ICOPVE(IN1,ICAAB,4*NGAS)
           ISMX = INSM1
        END IF
        IF(IOP.EQ.2) THEN 
           CALL ICOPVE(IN2,ICAAB,4*NGAS)
           ISMX = INSM2
        END IF
        IF(IOP.EQ.3) THEN 
           CALL ICOPVE(IN3,ICAAB,4*NGAS)
           ISMX = INSM3
        END IF
        IF(IOP.EQ.4) THEN 
           CALL ICOPVE(IN4,ICAAB,4*NGAS)
           ISMX = INSM4
        END IF
*. What operator is IOP going into ?
*. Compiler warnings ...
        IOP_REO = -2206
        DO JOP = 1, NIN
          IF(IREO(JOP).EQ.IOP) IOP_REO = JOP
        END DO
        IF(IOP_REO.EQ.1) CALL ICOPVE(ICAAB,IOUT1,4*NGAS)
        IF(IOP_REO.EQ.2) CALL ICOPVE(ICAAB,IOUT2,4*NGAS)
        IF(IOP_REO.EQ.3) CALL ICOPVE(ICAAB,IOUT3,4*NGAS)
        IF(IOP_REO.EQ.4) CALL ICOPVE(ICAAB,IOUT4,4*NGAS)
*
        IF(IOP_REO.EQ.1) IOUTSM1 = ISMX
        IF(IOP_REO.EQ.2) IOUTSM2 = ISMX
        IF(IOP_REO.EQ.3) IOUTSM3 = ISMX
        IF(IOP_REO.EQ.4) IOUTSM4 = ISMX
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reordering of T-operators '
        WRITE(6,*) ' Reorder array, new => old '
        CALL IWRTMA(IREO,1,NIN,1,NIN)
        WRITE(6,*)
        WRITE(6,*) ' Input T-operators : '
        WRITE(6,*)
        IF(NIN .GE.1) THEN
          WRITE(6,*) ' T1_IN : '
          CALL WRT_SPOX_TP(IN1,1)
        END IF
        IF(NIN .GE.2) THEN
          WRITE(6,*) ' T2_IN : '
          CALL WRT_SPOX_TP(IN2,1)
        END IF
        IF(NIN .GE.3) THEN
          WRITE(6,*) ' T3_IN : '
          CALL WRT_SPOX_TP(IN3,1)
        END IF
        IF(NIN .GE.4) THEN
          WRITE(6,*) ' T4_IN : '
          CALL WRT_SPOX_TP(IN4,1)
        END IF
        WRITE(6,*)
        WRITE(6,*) ' Output T-operators : '
        WRITE(6,*)
        IF(NIN .GE.1) THEN
          WRITE(6,*) ' T1_OUT : '
          CALL WRT_SPOX_TP(IOUT1,1)
        END IF
        IF(NIN .GE.2) THEN
          WRITE(6,*) ' T2_OUT : '
          CALL WRT_SPOX_TP(IOUT2,1)
        END IF
        IF(NIN .GE.3) THEN
          WRITE(6,*) ' T3_OUT : '
          CALL WRT_SPOX_TP(IOUT3,1)
        END IF
        IF(NIN .GE.4) THEN
          WRITE(6,*) ' T4_OUT : '
          CALL WRT_SPOX_TP(IOUT4,1)
        END IF
        WRITE(6,*) ' Output symmetries : '
        IF(NIN.EQ.1) WRITE(6,*) IOUTSM1
        IF(NIN.EQ.2) WRITE(6,*) IOUTSM1, IOUTSM2
        IF(NIN.EQ.3) WRITE(6,*) IOUTSM1, IOUTSM2, IOUTSM3
        IF(NIN.EQ.4) WRITE(6,*) IOUTSM1, IOUTSM2, IOUTSM3, IOUTSM4
      END IF
*
      RETURN
      END
      SUBROUTINE HCT1234(IHCAAB,IHINDEX,NTOP,IT1,IT2,IT3,IT4,
     &                   IT1SM,IT2SM,IT3SM,IT4SM,IHSM,T1,T2,T3,T4,
     &                   HT1234,IONLY_CONN,FACX,LHT1234)
* 
* All possible contractions of a given type of the Hamiltonian 
* with upto four T-operators of given type. 
*
* IONLY_CONN .ne.0 => Only fully connected  terms are included  
*
* Jeppe Olsen, April 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
*
*=========
*. Input 
*=========
*
*.. H-operator in CAAB form
      INTEGER IHCAAB(4*NGAS),IHINDEX(4)
*.. T-operators : type and coefficients
      INTEGER IT1(4*NGAS),IT2(4*NGAS),IT3(4*NGAS),IT4(4*NGAS)
      DIMENSION T1(*), T2(*), T3(*), T4(*)
*
*. Local scratch
**. Assuming that H contains atmost 2 particle operator 
      PARAMETER(MXPCONN = 4*4*4*4)
      INTEGER ICONN(4,MXPCONN)
      INTEGER IHEXP(3,4)
      INTEGER ID1(4*MXPNGAS),ID2(4*MXPNGAS),ID3(4*MXPNGAS)
      INTEGER ID4(4*MXPNGAS), IEX(4*MXPNGAS)
      INTEGER IORD(4),IEXD1234(4)
      INTEGER IT1_REO(4*MXPNGAS),IT2_REO(4*MXPNGAS),IT3_REO(4*MXPNGAS)
      INTEGER IT4_REO(4*MXPNGAS)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN 
         WRITE(6,*) ' ==================='
         WRITE(6,*) ' Welcome to HCT1234 '
         WRITE(6,*) ' ==================='
         WRITE(6,*)
         WRITE(6,*) ' H-operator in action (CAAB)'
         CALL WRT_SPOX_TP(IHCAAB,1)
         WRITE(6,*) ' IT1SM, IT2SM, IT3SM, IT4SM ',
     &                IT1SM, IT2SM, IT3SM, IT4SM
         IF(NTOP.GE.1) THEN
           WRITE(6,*) ' T1 operator '
           CALL WRT_SPOX_TP(IT1,1)
         END IF
         IF(NTOP.GE.2) THEN
           WRITE(6,*) ' T2 operator '
           CALL WRT_SPOX_TP(IT2,1)
         END IF
         IF(NTOP.GE.3) THEN
           WRITE(6,*) ' T3 operator '
           CALL WRT_SPOX_TP(IT3,1)
         END IF
         IF(NTOP.GE.4) THEN
           WRITE(6,*) ' T4 operator '
           CALL WRT_SPOX_TP(IT4,1)
         END IF
         IF(IONLY_CONN.NE.0) THEN
           WRITE(6,*) ' Only connected terms are included '
         END IF
         IF (NTEST.GE.1000) THEN
           WRITE(6,*) ' HCT1234 : Input block of CC vector function '
           CALL WRTMAT(HT1234,1,LHT1234,1,LHT1234)
         END IF
      END IF
*. Reform H operator to expanded form - assumed in CONTR_POS
C     REFORM_HTYP(IHCAAB,IHOP,NOP,IWAY,NGAS,IDOREO,IREO)
C     CALL REFORM_HTYP(IHCAAB,IHEXP,NHOP,1,NGAS,1,IHINDEX)
*. Obtain in CAAB form 
      CALL REFORM_HTYP(IHCAAB,IHEXP,NHOP,1,NGAS,0,IHINDEX)
*. Find all possible contractions of H operator with the NTOP T operators
      CALL CONTR_POS(IHEXP,NHOP,NTOP,IT1,IT2,IT3,IT4,
     &               IONLY_CONN,NCONN,ICONN)
C     CONTR_POS(IHOP,NHOP,NTOP,IT1,IT2,IT3,IT4,
C    &                     IONLY_CONN,NCONN,ICONN)
      DO JCONN = 1, NCONN
*. Obtain contraction order, and the individual deexcitation operators
        CALL HTYPE_TO_ED(IHCAAB,IHINDEX,NHOP,ICONN(1,JCONN),NTOP,
     &                   ID1,ID2,ID3,ID4,IEX,IEXD1234,IORD,ISIGN_DE)
*. Reorder t-operators so they correspond to new order
C     REO_CAABS(IN1,IN2,IN3,IN4,NIN,IOUT1,IOUT2,IOUT3,IOUT4,
C    &                     INSM1,INSM2,INSM3,INSM4,
C    &                     IOUTSM1,IOUTSM2,IOUTSM3,IOUTSM4,
C    &                     IREO, NGAS)
*. I will use 4 operators, as the last (inactive have just been 
*. defined to be unit operator. This is done as OPCT1234 
*. currently requires four ops and four symmetries
        CALL REO_CAABS(IT1,IT2,IT3,IT4,4,
     &                 IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                 IT1SM,IT2SM,IT3SM,IT4SM,
     &                 IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                 IORD, NGAS)
*.   And then do the contraction 
*
* There are 24 different ways of making the contractions so, I must
* have 24 cases, and insert the appropriate T-operators, 
* If I do this without problems, I will give myself a beer ..
*
        IF(IORD(1).EQ.1.AND.IORD(2).EQ.2.AND.IORD(3).EQ.3
     &  .AND.IORD(4).EQ.4) THEN
*.1234
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T1,T2,T3,T4,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.1.AND.IORD(2).EQ.2.AND.IORD(3).EQ.4
     &  .AND.IORD(4).EQ.3)THEN
*.1243
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T1,T2,T4,T3,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.1.AND.IORD(2).EQ.3.AND.IORD(3).EQ.2
     &  .AND.IORD(4).EQ.4)THEN
*. 1324
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T1,T3,T2,T4,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.1.AND.IORD(2).EQ.3.AND.IORD(3).EQ.4.
     &  AND.IORD(4).EQ.2)THEN
*. 1342 
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T1,T3,T4,T2,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.1.AND.IORD(2).EQ.4.AND.IORD(3).EQ.2
     &  .AND.IORD(4).EQ.3)THEN
*. 1423
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T1,T4,T2,T3,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.1.AND.IORD(2).EQ.4.AND.IORD(3).EQ.3
     &  .AND.IORD(4).EQ.2)THEN
*. 1432
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T1,T4,T3,T2,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.2.AND.IORD(2).EQ.1.AND.IORD(3).EQ.3.AND
     &  .IORD(4).EQ.4)THEN
*.2134
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T2,T1,T3,T4,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.2.AND.IORD(2).EQ.1.AND.IORD(3).EQ.4.AND
     &  .IORD(4).EQ.3)THEN
*.2143
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T2,T1,T4,T3,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.2.AND.IORD(2).EQ.3.AND.IORD(3).EQ.1.AND
     &  .IORD(4).EQ.4)THEN
*.2314
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T2,T3,T1,T4,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.2.AND.IORD(2).EQ.3.AND.IORD(3).EQ.4.AND
     &  .IORD(4).EQ.1)THEN
*.2341
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T2,T3,T4,T1,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.2.AND.IORD(2).EQ.4.AND.IORD(3).EQ.1.AND
     &  .IORD(4).EQ.3)THEN
*.2413
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T2,T4,T1,T3,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.2.AND.IORD(2).EQ.4.AND.IORD(3).EQ.3.AND
     &  .IORD(4).EQ.1)THEN
*.2431
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T2,T4,T3,T1,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.3.AND.IORD(2).EQ.1.AND.IORD(3).EQ.2.AND.
     &  IORD(4).EQ.4)THEN
*.3124
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T3,T1,T2,T4,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.3.AND.IORD(2).EQ.1.AND.IORD(3).EQ.4.AND.
     &  IORD(4).EQ.2)THEN
*.3142
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T3,T1,T4,T2,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.3.AND.IORD(2).EQ.2.AND.IORD(3).EQ.1.AND.
     &  IORD(4).EQ.4)THEN
*.3214
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T3,T2,T1,T4,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.3.AND.IORD(2).EQ.2.AND.IORD(3).EQ.4.AND.
     &  IORD(4).EQ.1)THEN
*.3241
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T3,T2,T4,T1,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.3.AND.IORD(2).EQ.4.AND.IORD(3).EQ.1.AND.
     &  IORD(4).EQ.2)THEN
*.3412
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T3,T4,T1,T2,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.3.AND.IORD(2).EQ.4.AND.IORD(3).EQ.2.AND.
     &  IORD(4).EQ.1)THEN
*.3421
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T3,T4,T2,T1,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.4.AND.IORD(2).EQ.1.AND.IORD(3).EQ.2.AND.
     &  IORD(4).EQ.3)THEN
*. 4123
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T4,T1,T2,T3,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.4.AND.IORD(2).EQ.1.AND.IORD(3).EQ.3.AND.
     &  IORD(4).EQ.2)THEN
*. 4132
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T4,T1,T3,T2,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.4.AND.IORD(2).EQ.2.AND.IORD(3).EQ.1.AND.
     &  IORD(4).EQ.3)THEN
*. 4213
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T4,T2,T1,T3,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.4.AND.IORD(2).EQ.2.AND.IORD(3).EQ.3.AND.
     &  IORD(4).EQ.1)THEN
*. 4231
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T4,T2,T3,T1,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.4.AND.IORD(2).EQ.3.AND.IORD(3).EQ.1.AND.
     &  IORD(4).EQ.2)THEN
*. 4312
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T4,T3,T1,T2,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE IF (IORD(1).EQ.4.AND.IORD(2).EQ.3.AND.IORD(3).EQ.2.AND.
     &  IORD(4).EQ.1)THEN
*. 4321
          CALL OPCT1234M(IEX,ID1,ID2,ID3,ID4,
     &                   IT1_REO,IT2_REO,IT3_REO,IT4_REO,
     &                   T4,T3,T2,T1,HT1234,
     &                   IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE,FACX)
        ELSE 
          WRITE(6,*) ' Sorry, but jeppe has messed up again ... '
          WRITE(6,*) ' untested IORD sequence : ', (IORD(I),I=1,4)
          STOP       ' untested IORD sequence   '
        END IF
C     OPCT1234M(IOEX,IO1DX,IO2DX,IO3DX,IO4DX,
C    &                    IT1,IT2,IT3,IT4,T1,T2,T3,T4,OT1234,
C    &                    IT1SM,IT2SM,IT3SM,IT4SM,IOPSM,IOPINDX,
C    &                    ISIGN,FACX)
      END DO
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' HCT1234 : Output block of CC vector function '
         CALL WRTMAT(HT1234,1,LHT1234,1,LHT1234)
      END IF
*
      RETURN
      END 
      SUBROUTINE HTYPE_TO_ED(IHTYPE,IHINDEX,NHOP,ICONT,
     &           NTOP,ID1,ID2,ID3,ID4,IEX,IEXD1234,IORD,ISIGN_DE)
*    
*. A Hamilton operator (IHTYPE) is given as a set of C/A operators  
*. (in CA CB AA AB order), and a ordering array IHINDEX.
*. A contraction is given (ICONT), telling which T-operator 
*. each operator in H should be connected to.
*. Obtain 1 : The order in which the contraction should be done 
*.        2 : The individual excitations and deexcitations in CAAB form
*.        3 : Sign for bringing CA CB AA AB operator into EX D1 D2 D3 D4
*.            form
*
* Sign for bringing CA CB AA AB operator into EX D1 D2 D3 D4 form
*
* Jeppe Olsen, July 16, 2001, Finished Easter 2003 ( April)
*              Debugging initiated Nov. 2004 
*. Some phase modifications, Nov. 2004 ... 
*  ( I understand why somebody asks higher powers for signs,
*    they are difficult to get yourself ...)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IHTYPE(NGAS,4), IHINDEX(4),ICONT(4)
*. Output 
      INTEGER ID1(NGAS,4),ID2(NGAS,4),ID3(NGAS,4),ID4(NGAS,4)
      INTEGER IEX(NGAS,4),IEXD1234(4)
      INTEGER IORD(4)
* IORD(IOP) : New operator IOP is old operator IORD(IOP) ( new => old)
*. Local scratch 
      INTEGER IHOP(4,3)
      INTEGER IDOP(4,3,4),IEXOP(4,3), LD(4) 
C     INTEGER IT_AR(4)
      INTEGER IEX_INDEX(4),ID_INDEX(4,4)
*
      IZERO = 0
*
*. Reform IHTYPE from CAAB form to a string of NHOP operators
C     REFORM_HTYP(IHCAAB,IHOP,NOP,IWAY,NGAS,IDOREO,IREO)
C     CALL REFORM_HTYP(IHTYPE,IHOP,NHOP,1,NGAS,1,IHINDEX)
C     CALL REFORM_HTYP(IHTYPE,IHOP,NHOP,1,NGAS,1,IHINDEX)
*. Keep CAAB ordering ..
      CALL REFORM_HTYP(IHTYPE,IHOP,NHOP,1,NGAS,0,IHINDEX)
*. Obtain the deexcitations/excitations D1,D2,D3,D4,EX (they will now
*. be in CAAB order)
      DO I = 1, 4
        LD(I) = 0
      END DO
      NEX = 0
      DO IOP = 1, NHOP
        INDEX = ICONT(IOP)
*.      ^ Index tells to which T-operator an operator is connected, 0=>disc.
        IF(INDEX.EQ.0) THEN
*. Excitation part
          NEX = NEX + 1
          IEXOP(NEX,1) = IHOP(IOP,1)
          IEXOP(NEX,2) = IHOP(IOP,2)
          IEXOP(NEX,3) = IHOP(IOP,3)
C*. Original index for this operator
C         IEX_INDEX(NEX) = IHINDEX(IOP)
*. index for this operator in IHOP
          IEX_INDEX(NEX) = IOP
        ELSE 
*.Deexcitation part
          LD(INDEX) = LD(INDEX) + 1
          I = LD(INDEX)
          IDOP(I,1,INDEX) = IHOP(IOP,1)
          IDOP(I,2,INDEX) = IHOP(IOP,2)
          IDOP(I,3,INDEX) = IHOP(IOP,3)
C         ID_INDEX(I,INDEX) = IHINDEX(IOP)
          ID_INDEX(I,INDEX) = IOP
        END IF
      END DO
*. Obtain the execution order
C     CONTR_ORD4(NDOP,ID1,ID2,ID3,ID4,
C    &           LD1,LD2,LD3,LD4,IORD) 
       CALL CONTR_ORD4(NTOP,IDOP(1,1,1),IDOP(1,1,2),
     &                 IDOP(1,1,3),IDOP(1,1,4),
     &                 LD(1),LD(2),LD(3),LD(4),IORD)
*. I want IORD to contain 4 well-defined terms, also 
*. when there are less than 4 operators so
       IF(NTOP.LT.4) IORD(4) = 4
       IF(NTOP.LT.3) IORD(3) = 3
       IF(NTOP.LT.2) IORD(2) = 2
       IF(NTOP.LT.1) IORD(1) = 1
*...   (^I guess I do not need LT.0)
*. On return IORD : New => Old order
*. Set up IEXD1234 : Index in operator EXD1D2D3D4 to original order
      DO JEX = 1, NEX
        IEXD1234(JEX) = IEX_INDEX(JEX)
      END DO
      LEN = NEX
*. Loop over deexcitations in NEW order
      DO JDOP = 1, NTOP
        JDOP_OLD = IORD(JDOP)
C?      WRITE(6,*) 'JDOP, JDOP_OLD', JDOP, JDOP_OLD
C?      WRITE(6,*) ' LD(JDOP_OLD) ',  LD(JDOP_OLD)
        DO JOP = 1, LD(JDOP_OLD)
          LEN = LEN + 1
          IEXD1234(LEN) = ID_INDEX(JOP,JDOP_OLD)
        END DO
      END DO
*. IEXD1234 refers to indeces in IHOP, modify so they 
*. correspond to indeces in the original operator
      DO JOP = 1, NHOP
*. Set up the deexcitations numbered according to execution order 
        IEXD1234(JOP) = IHINDEX(IEXD1234(JOP))
      END DO
*. Excitation in CAAB form
      IDUM = 0
      CALL REFORM_HTYP(IEX,IEXOP,NEX,2,NGAS,0,IDUM)
C          REFORM_HTYP(IHCAAB,IHOP,NOP,IWAY,NGAS,0,IDUM)
*. First operator to be contracted
      CALL REFORM_HTYP(ID1,IDOP(1,1,IORD(1)),LD(IORD(1)),2,NGAS,0,IDUM)
*. Second op
      IF(NTOP.GE.2) THEN 
       CALL REFORM_HTYP(ID2,IDOP(1,1,IORD(2)),LD(IORD(2)),2,NGAS,0,IDUM)
      ELSE
       CALL ISETVC(ID2,IZERO,4*NGAS)
      END IF
*. third  op
      IF(NTOP.GE.3) THEN 
       CALL REFORM_HTYP(ID3,IDOP(1,1,IORD(3)),LD(IORD(3)),2,NGAS,0,IDUM)
      ELSE
       CALL ISETVC(ID3,IZERO,4*NGAS)
      END IF
*. fourth  op
      IF(NTOP.GE.4) THEN
       CALL REFORM_HTYP(ID4,IDOP(1,1,IORD(4)),LD(IORD(4)),2,NGAS,0,IDUM)
      ELSE
       CALL ISETVC(ID4,IZERO,4*NGAS)
      END IF
*. Sign for bringing H operator from CAAB order to EXD1234 order
C IPERM_PARITY(IPERM,NELMNT)
      ISIGN_DE = IPERM_PARITY(IEXD1234,NHOP)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================'
        WRITE(6,*) ' Output from HTYPE_TO_ED ' 
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        WRITE(6,*) ' Contraction of ops. (ICONT) '
        CALL IWRTMA(ICONT,1,NHOP,1,4)
        WRITE(6,*) ' Suggested order of contraction : '
        CALL IWRTMA(IORD,1,NTOP,1,NTOP)
        WRITE(6,*) ' IEX,D1,D2,D3,D4 : '
        CALL WRT_SPOX_TP(IEX,1)
        CALL WRT_SPOX_TP(ID1,1)
        IF(NTOP.GE.2) 
     &  CALL WRT_SPOX_TP(ID2,1)
        IF(NTOP.GE.3) 
     &  CALL WRT_SPOX_TP(ID3,1)
        IF(NTOP.GE.4) 
     &  CALL WRT_SPOX_TP(ID4,1) 
        WRITE(6,*) ' IEXD1234 : new to old index array '
        CALL IWRTMA(IEXD1234,1,NHOP,1,NHOP)
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_HTYP(IHCAAB,IHOP,NOP,IWAY,NGAS,IDOREO,IREO)
*
* Reform different forms of Hamiltonian operator 
*
* IWAY = 1 : IH_CAAB,IHINDX => IHOP
*      = 2 : IH_CAAB <= IHOP
* 
* If IDOREO = 1, reordering of indeces is done according to IREO
*
*. Note : Operators are ordered in CA CB AA AB order
*
* Jeppe Olsen, July 16, 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IREO(*)
*. Input and otput
      INTEGER IHCAAB(NGAS,4), IHOP(4,3)
*. Local scratch
      INTEGER JHOP(4,3)
*
      IF(IWAY.EQ.1) THEN
* HCAAB => HOP
       LOP = 0
       ICAAB = 0
       DO ICA = 1, 2
        DO IAB = 1, 2
          ICAAB = (ICA-1)*2 + IAB
          DO IGAS = 1, NGAS
            L = IHCAAB(IGAS,ICAAB)
            DO IOP = 1, L
              LOP = LOP + 1
              JHOP(LOP,1) = IGAS
              JHOP(LOP,2) = ICA
              JHOP(LOP,3) = IAB
            END DO
          END DO
        END DO
       END DO
       NOP = LOP
*. Reordering or transfer
       DO IOP = 1, NOP
         IF(IDOREO.EQ.0) THEN
           IOPEFF = IOP
         ELSE
           IOPEFF = IREO(IOP)
         END IF
         DO IND = 1, 3
           IHOP(IOP,IND) = JHOP(IOPEFF,IND)
         END DO
       END DO
      ELSE 
*. HOP => HCAAB
*. Nop must here be input
       IZERO = 0
       CALL ISETVC(IHCAAB,IZERO,4*NGAS)
       DO IOP = 1, NOP
         IGAS = IHOP(IOP,1)
         ICA  = IHOP(IOP,2)
         IAB  = IHOP(IOP,3)
         ICAAB = (ICA-1)*2+IAB
         IHCAAB(IGAS,ICAAB) = IHCAAB(IGAS,ICAAB) + 1
       END DO
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reform of Htype '
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' CAAB => OP reform '
        ELSE
         WRITE(6,*) ' OP => CAAB reform'
        END IF
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) ' The reorder array '
          CALL IWRTMA(IREO,1,NOP,1,NOP)
        END IF 
        WRITE(6,*) ' CAAB form '
        CALL WRT_SPOX_TP(IHCAAB,1)
        WRITE(6,*) ' Operator form '
C            WRT_CNTR(ICONT,NCONT,LDIM)
        CALL WRT_CNTR3(IHOP,NOP,4)
      END IF
*
      RETURN
      END
      SUBROUTINE CONTR_ORD4(NDOP,ID1,ID2,ID3,ID4,
     &           LD1,LD2,LD3,LD4,IORD) 
*
* NDOP contraction operators ID1, ID2, ... are given. 
* Find which order these contraction operators should be 
* applied 
*
* Jeppe Olsen, July15, 2001, debugging initialized nov. 2004
*
* Initial version : Contraction containing most elements
* are contracted first
* 
* Output : IORD : New operator I is original operator IORD(I)

*
      INCLUDE 'implicit.inc'
*. Input : Contraction operators 
      DIMENSION ID1(4,3),ID2(4,3),ID3(4,3),ID4(4,3)
*. Output 
      INTEGER IORD(4)
*. Local scratch
      INTEGER LEN_AR(4),INO(4), ISCR(4)
*
*. Length of different contractions
                    CALL DIM_CNTR(ID1,LD1,4,LEN_AR(1))
      IF(NDOP.GE.2) CALL DIM_CNTR(ID2,LD2,4,LEN_AR(2))
      IF(NDOP.GE.3) CALL DIM_CNTR(ID3,LD3,4,LEN_AR(3))
      IF(NDOP.GE.4) CALL DIM_CNTR(ID4,LD4,4,LEN_AR(4))
*. Multiply with -1 so numerically largest numbers 
*. will come first in sort
      DO I = 1, NDOP
       LEN_AR(I) = - LEN_AR(I)
      END DO
C           ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
       CALL ORDINT(LEN_AR,ISCR,NDOP,INO,0)   
*. INO now gives new order => old order which is what we want so
      CALL ICOPVE(INO,IORD,NDOP)
*. Fusk, enforce order 1 2 3 4
C     WRITE(6,*) ' Order of contraction is set to 1234'
C     IORD(1) = 1
C     IORD(2) = 2
C     IORD(3) = 3
C     IORD(4) = 4
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Contraction operators  : '
C  WRT_CNTR(ICONT,NCONT,LDIM)
       CALL WRT_CNTR3(ID1,LD1,4)
       IF(NDOP.GE.2) CALL WRT_CNTR3(ID2,LD2,4)
       IF(NDOP.GE.3) CALL WRT_CNTR3(ID3,LD3,4)
       IF(NDOP.GE.4) CALL WRT_CNTR3(ID4,LD4,4)
       WRITE(6,*) ' Order of contraction for operators '
       CALL IWRTMA(IORD,1,NDOP,1,NDOP)
      END IF
*
      RETURN
      END
      SUBROUTINE CONTR_POS(IHOP,NHOP,NTOP,IT1,IT2,IT3,IT4,
     &                     IONLY_CONN,NCONN,ICONN)
*
* A general operator IHOP containing NHOP operators are given.
*.Find possible contractions with the NTOP excitation 
*.operators IT1, IT2, ..
* Only connections where all  deexcitation operators in IHOP are 
* contracted with excitation operators are included
*
*. IF IONLY_CONN = 1, then connected terms are included
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IHOP(4,3)
*. GAS : 1, ICA : 2, IAB : 3
      INTEGER IT1(4*NGAS),IT2(4*NGAS),IT3(4*NGAS),IT4(4*NGAS)
*. Output
      INTEGER ICONN(4,*)
*. ICONN(JOP,JCONN) : Which T-operator is operator JOP in IHOP
*.                    connected with in connection JCONN
*. ICONN(JOP,JCONN) = 0 => No connection of operator JOP in IHOP
*
*. Local scratch 
      INTEGER ISCR(4), IMAXVAL(4)
      INTEGER IT(4*MXPNGAS,4)
      INTEGER JCONN(4)
      INTEGER IDE(4)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to CONTR_POS' 
        WRITE(6,*) ' Form of general operator : '
        CALL WRT_CNTR3(IHOP,NHOP,4)
        IF(NTOP.GE.1) CALL WRT_SPOX_TP(IT1,1)
        IF(NTOP.GE.2) CALL WRT_SPOX_TP(IT2,1)
        IF(NTOP.GE.3) CALL WRT_SPOX_TP(IT3,1)               
        IF(NTOP.GE.4) CALL WRT_SPOX_TP(IT4,1)
      END IF
*. Check if some of the operators are equivalent, i.e. 
*. have identical CA, AB, IGAS. 
*. Currently I am only checking whether operators 1 and 2 are 
*. identical and whether operators 3 and 4 are identical 
      N12EQV = 1
      DO ICOMP = 1, 3
        IF(IHOP(1,ICOMP).NE.IHOP(2,ICOMP)) N12EQV = 0
      END DO
*
      N34EQV = 0
      IF(NHOP.EQ.4) THEN
        N34EQV = 1
        DO ICOMP = 1, 3
          IF(IHOP(3,ICOMP).NE.IHOP(4,ICOMP)) N34EQV = 0
        END DO
      END IF

      I_OLD_OR_NEW = 2
      IF(I_OLD_OR_NEW.EQ.1) THEN       
*. Use IPHGAS1 to classify operators as excitation or deexcitation operators
*. Should be modified for open shell reference 
        DO JOP = 1, NHOP
          IF((IPHGAS1(IHOP(JOP,1)).EQ.1.AND.IHOP(JOP,2).EQ.2).OR.
     &       (IPHGAS1(IHOP(JOP,1)).EQ.2.AND.IHOP(JOP,2).EQ.1) ) THEN
             IDE(JOP) = 1
          ELSE
             IDE(JOP) = 2
          END IF
        END DO
*.    ELSE
*. Modifications allowing the use of High Spin Open shell
        DO JOP = 1, NHOP
         JGAS= IHOP(JOP,1)
         JCA = IHOP(JOP,2)
         JAB = IHOP(JOP,3)
         IF((IHPVGAS_AB(JGAS,JAB).EQ.2.AND.JCA.EQ.2).OR.
     &      (IHPVGAS_AB(JGAS,JAB).EQ.1.AND.JCA.EQ.1)) THEN
             IDE(JOP) = 1
          ELSE
             IDE(JOP) = 2
          END IF
        END DO
       END IF

C?    WRITE(6,*) ' The IDE array ', (IDE(I),I=1,NHOP)
*. A given connection is defined by NHOP integers between 1 and 
*. NTOP. Loop over these integers
      CALL ISETVC(IMAXVAL,NTOP,NHOP)
      IFIRST = 1
      IZERO = 0
      NCONN = 0
 1000 CONTINUE
        IF(IFIRST.EQ.1) THEN
          CALL ISETVC(ISCR,IZERO,NHOP)
          NONEW = 0
          IFIRST = 0
        ELSE
C              NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
          CALL NXTNUM2(ISCR,NHOP,IZERO,IMAXVAL,NONEW)
        END IF
C?      WRITE(6,*) ' Next suggested  connection '
C?      CALL IWRTMA(ISCR,1,NHOP,1,NHOP)
        IF(NONEW.EQ.0) THEN
*. Check if this connection is nonvanishing
          CALL ICOPVE(IT1,IT(1,1),4*NGAS)
          IF(NTOP.GE.2) CALL ICOPVE(IT2,IT(1,2),4*NGAS)
          IF(NTOP.GE.3) CALL ICOPVE(IT3,IT(1,3),4*NGAS)
          IF(NTOP.GE.4) CALL ICOPVE(IT4,IT(1,4),4*NGAS)
*
          I_AM_OKAY = 1
          DO IOP = 1, NHOP
            JTOP = ISCR(IOP)
            JGAS = IHOP(IOP,1)
            JCA =  IHOP(IOP,2)
            JAB =  IHOP(IOP,3)
*. Well, a creation operator in IHOP should connect with 
*. an annihilation operator in IT and vice versa so
            IF(JCA.EQ.1) THEN
             JTCA = 2
            ELSE
             JTCA = 1
            END IF
*
            JCAAB = (JTCA-1)*2+JAB
            J_AD = JGAS + (JCAAB-1)*NGAS
            IF(JTOP.GT.0) THEN
              IT(J_AD,JTOP) = IT(J_AD,JTOP)-1
              IF(IT(J_AD,JTOP).LT.0)  I_AM_OKAY = 0
            END IF
          END DO
C?        WRITE(6,*) ' I_AM_OKAY for this connection ', I_AM_OKAY
*. If required, check whether term is connected
          IF(IONLY_CONN.EQ.1.AND.I_AM_OKAY.EQ.1) THEN 
            CALL ISETVC(JCONN,IZERO,NTOP)
            DO JHOP = 1, NHOP
              IF(ISCR(JHOP).NE.0) JCONN(ISCR(JHOP)) = 1
            END DO
            ICONNECTED = 1
            DO JTOP = 1, NTOP
              IF(JCONN(JTOP).EQ.0) ICONNECTED = 0
            END DO
          ELSE
            ICONNECTED = 1
          END IF
C?        WRITE(6,*) ' ICONNECTED for this connection ', ICONNECTED
*         ^ End of check for connectedness is required
*. Check to see ensure that all deexcitation operators are 
*. connected
          DO JHOP = 1, NHOP
            IF(IDE(JHOP).EQ.1.AND.ISCR(JHOP).EQ.0) 
     &      I_AM_OKAY = 0       
          END DO
C?        WRITE(6,*) ' Check , All deexcitations connected ', I_AM_OKAY
*. Check that for equivalent operators the contracted operators 
*. occur in ascending order. This is to ensure that 
*. contractions splitting equivalent operators only 
*. are included once
          IF(N12EQV.EQ.1.AND.ISCR(1).LT.ISCR(2)) I_AM_OKAY = 0
          IF(NHOP.EQ.4) THEN
            IF(N34EQV.EQ.1.AND.ISCR(3).LT.ISCR(4)) I_AM_OKAY = 0
          END IF
*
          IF(I_AM_OKAY.EQ.1.AND.ICONNECTED.EQ.1) THEN
*. Well a new connection has been born 
             NCONN = NCONN + 1
             CALL ICOPVE(ISCR,ICONN(1,NCONN),NHOP)
          END IF
      GOTO 1000
        END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CONTR_POS '
*
        WRITE(6,*) ' Number of obtained connections = ', NCONN
        WRITE(6,*) ' Obtained connections : '
        CALL IWRTMA(ICONN,NHOP,NCONN,4,NCONN)
      END IF
*
      RETURN
      END
      SUBROUTINE EXP_MT_H_EXP_TC(T,CC_VEC)
*
* Very new CC vector function 
* (after new and probably before extremely new )
*
* Master routine 
*
* Jeppe Olsen, Initiated Jan. 2001
*              Finished  July 2001
*              Tested and debugged March-April 2003 (sic)
*              and in November 2004 
*
c      INCLUDE 'implicit.inc'
*. General input 
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      DIMENSION T(*)
*. Output
      DIMENSION CC_VEC(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXP_MT')
*
      LCCBD12 = LCCB 
      WRITE(6,*) ' ISIMTRH = ', ISIMTRH
*
C?    WRITE(6,*) ' Input T - amplitudes '
C?    CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
*
*. Find types of orbital- and spinorbital-excitations
      CALL H_TYPES(0,N1TP,N2TP,IDUM,IDUM,IDUM,
     &             N1OBTP,N2OBTP,IDUM,IDUM,IDUM,IDUM,
     &             IDUM,IDUM,IDUM)
      N12TP = N1TP + N2TP
      LEN = 4*NGAS*(N1TP+N2TP)
      CALL MEMMAN(KLHTP,LEN,'ADDL  ',1,'H_TYP ')
      CALL MEMMAN(KLHINDX,LEN,'ADDL ',1,'HINDX ')
      CALL MEMMAN(KLHSIGN,LEN,'ADDL ',1,'HSIGN ')
      LEN = (N1OBTP + N2OBTP)*2*NGAS
      CALL MEMMAN(KLHOBTP,LEN,'ADDL  ',1,'HOBTP ')
      N12OBTP = N1OBTP + N2OBTP
      CALL MEMMAN(KLNSOX_FOR_OX_H,N12OBTP,'ADDL  ',1,'NSOX_H')
      CALL MEMMAN(KLISOX_TO_OX_H,N12TP,'ADDL  ',1,'ISOXTH') 
      CALL MEMMAN(KLISOX_FOR_OX_H,N12TP,'ADDL  ',1,'ISOXFH') 
      CALL MEMMAN(KLIBSOX_FOR_OX_H,N12OBTP,'ADDL  ',1,'IBSOXH') 
      CALL MEMMAN(KLSOX_SPFLIP,N12TP,'ADDL  ',1,'HSPFLP') 
      CALL MEMMAN(KLH_EXC2,N12TP,'ADDL  ',1,'H_EXC2')
*
      CALL H_TYPES(1,N1TP,N2TP,WORK(KLHTP),
     &             WORK(KLHINDX),WORK(KLHSIGN),
     &             N1OBTP,N2OBTP,WORK(KLHOBTP),
     &             WORK(KLNSOX_FOR_OX_H),WORK(KLISOX_TO_OX_H),
     &             WORK(KLISOX_FOR_OX_H),WORK(KLIBSOX_FOR_OX_H),
     &             WORK(KLSOX_SPFLIP),WORK(KLH_EXC2) )
*. Allocate memory 
      CALL EXP_MT_H_EXP_TC_MEM
*.
      IHSM = 1
      NHTP = N1TP + N2TP
      ITSM = 1
*
      CALL EMTHET_COM(T,ITSM,CC_VEC,NSPOBEX_TP,
     &     WORK(KLSOBEX),WORK(KLIBSOBEX),
     &     WORK(KLLSOBEX),
     &     NOBEX_TP,WORK(KIBSOX_FOR_OX),WORK(KNSOX_FOR_OX),
     &     WORK(KISOX_FOR_OX),
     &     NHTP,WORK(KLHTP),WORK(KLHINDX),WORK(KLHSIGN),IHSM,
     &     WORK(KLSPOBEX_AC) )
*.
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXP_MT')
      RETURN
      END

      SUBROUTINE EMTHET_COM(T,ITSM,CC_VEC_FNC,NSPOBEX,
     &           ISPOBEX_TP,IBSPOBEX,LSPOBEX,NOBEX,
     &           IBSOX_FOR_OX,NSOX_FOR_OX,ISOX_FOR_OX,
     &           NHTP,IHTP,IHINDEX,IHSIGN,IHSM,ISPOBEX_AC)
*
* Calculate CC vector function using only 
* connected terms and commutators
*
* Jeppe Olsen, Jan. 2001, completed april 2003 (I hope ...)
*                        Well, I am sitting debugging nov. 2004
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Specific input
*. Info on included spinorbital excitations
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX+1),LSPOBEX(NSPOBEX+1)
      INTEGER IBSPOBEX(NSPOBEX+1)
*. Are given type of spin-orbital excitation active ?
      INTEGER ISPOBEX_AC(NSPOBEX)
*. Orbital excitations to spin-orbital excitations
      INTEGER IBSOX_FOR_OX(*),NSOX_FOR_OX(*),ISOX_FOR_OX(*)
*. Hamilton operator types
      INTEGER IHTP(4*NGAS,NHTP), IHINDEX(4,NHTP), IHSIGN(NHTP)
*. Input CC vector
      DIMENSION T(*)
*. Output CC vector function 
      DIMENSION CC_VEC_FNC(*)
*
*. Local scratch 
*
CT    INTEGER ITTOCC(4*MXPNGAS),ITTTOCC(4*MXPNGAS),ITTTTOCC(4*MXPNGAS)
      INTEGER IHTOCC(4*MXPNGAS), IHTTOCC(4*MXPNGAS)
      INTEGER IHTTTOCC(4*MXPNGAS), IHTTTTOCC(4*MXPNGAS)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'EMTHET')
*
      NTEST = 1000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Information from EMTHET_COM '
        WRITE(6,*) ' ============================='
        WRITE(6,*) ' Input T-amplitudes '
        CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
      END IF
*
      IUNIOP = NSPOBEX+1
      IBUNIOP = IBSPOBEX(IUNIOP)
      
*
      IDUM = -1
      XDUM = -1.0D0
      ONE = 1.0D0
      IONE = 1
*
*. Start by zeroing the CC-vector function 
      ZERO = 0.0D0
      CALL SETVEC(CC_VEC_FNC,ZERO,N_CC_AMP+1)
*. Loop over blocks of Hamilton operator
      DO IHOP = 1, NHTP
        FACTORH = FLOAT(IHSIGN(IHOP))
*. About signs : signs are all generated in HTYPE_TO_ED, so p.t
        FACTORH = 1.0D0
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Info for Hamiltonian op number ', IHOP 
          WRITE(6,*) ' Operator in CAAB form : '
          CALL WRT_SPOX_TP(IHTP(1,IHOP),1)
          WRITE(6,*) ' FACTORH = ', FACTORH
        END IF
*. This operator in explicit form 
*. Number of elementary excitation and deexcitation operators in this op
        CALL RANK_FOR_CAAB(IHTP(1,IHOP),NEX,NDEEX)
*. If the operator contains deexcitations, H(IHOP) !ref> = 0, so
        IF(NDEEX.EQ.0) THEN
           NCOM_MIN = 0
        ELSE
           NCOM_MIN = 1
        END IF
*. The largest number of commutators equals the number of 
*. deexcitation operators so 
        NCOM_MAX = NDEEX
        IF(NTEST.GE.10) 
     &  WRITE(6,*) ' NCOM_MIN, NCOM_MAX = ', NCOM_MIN, NCOM_MAX
*. Loop over number of commutators
        DO NCOMMU = NCOM_MIN, NCOM_MAX
*
          IF(NCOMMU.EQ.0) CALL LUCIAQENTER('0_COM')
          IF(NCOMMU.EQ.1) CALL LUCIAQENTER('1_COM')
          IF(NCOMMU.EQ.2) CALL LUCIAQENTER('2_COM')
          IF(NCOMMU.EQ.3) CALL LUCIAQENTER('3_COM')
          IF(NCOMMU.EQ.4) CALL LUCIAQENTER('4_COM')
*
          IF(NTEST.GE.100) 
     &    WRITE(6,*) ' Number of commutators : ', NCOMMU
*
          IF(NCOMMU.EQ.0) THEN
*
*. 0 commutators : H(IHTP) |ref>
*
C           INUM_FOR_OCC(IOCC_HTF_AR,IHTFTP)
            CALL INUM_FOR_OCC(IHTP(1,IHOP),IHCTP)
            IF(IHCTP.GE.1) THEN
              LCCFBLK = LSPOBEX(IHCTP)
              IF(NTEST.GE.100) 
     /        WRITE(6,*) '0: HCT1234 will be called for IHCTP = ', IHCTP
                   CALL HCT1234(IHTP(1,IHOP),IHINDEX(1,IHOP),0,
     &             ISPOBEX_TP(1,IUNIOP),ISPOBEX_TP(1,IUNIOP),
     &             ISPOBEX_TP(1,IUNIOP),ISPOBEX_TP(1,IUNIOP),
     &             IONE,IONE,IONE,IONE,IHSM,
     &             ONE,ONE,ONE,ONE,
     &             CC_VEC_FNC(IBSPOBEX(IHCTP)),1,FACTORH,LCCFBLK)
            END IF
          END IF
*
          IF(NCOMMU.EQ.1) THEN
*. Loop over operators in T
*. Loop over combinations of NCOMMU operators
           DO IT1 = 1, NSPOBEX
            IF(ISPOBEX_AC(IT1).EQ.1) THEN
              IT1_B = IBSPOBEX(IT1)
*. Occupation of HT1
             CALL OP_T_OCC(IHTP(1,IHOP),ISPOBEX_TP(1,IT1),
     &            IHTOCC,IMZEROHT)
             IHCTP = 0
*. address of HT1
             IF(IMZEROHT.EQ.0)  CALL INUM_FOR_OCC(IHTOCC,IHCTP)
             IF(IHCTP.GE.1) THEN
                LCCFBLK = LSPOBEX(IHCTP)
                IF(NTEST.GE.10) 
     /          WRITE(6,*)
     &                 '1: HCT1234 will be called for IHCTP = ', IHCTP
                     CALL HCT1234(IHTP(1,IHOP),IHINDEX(1,IHOP),1,
     &               ISPOBEX_TP(1,IT1),ISPOBEX_TP(1,IUNIOP),
     &               ISPOBEX_TP(1,IUNIOP),ISPOBEX_TP(1,IUNIOP),
     &               ITSM,IONE,IONE,IONE,IHSM,
     &               T(IT1_B),ONE,ONE,ONE,
     &               CC_VEC_FNC(IBSPOBEX(IHCTP)),1,FACTORH,LCCFBLK)
             END IF
            END IF
*           ^ End if excitation was active
           END DO
          END IF
*         ^ End of NCOMMU .EQ. 1
*
          IF(NCOMMU.EQ.2) THEN
*. Loop over combinations of NCOMMU operators
           DO IT1 = 1, NSPOBEX
           IF(ISPOBEX_AC(IT1).EQ.1) THEN
            IT1_B = IBSPOBEX(IT1)
*. Occupation of HT1
            CALL OP_T_OCC(IHTP(1,IHOP),ISPOBEX_TP(1,IT1),
     &           IHTOCC,IMZEROHT)
*. address of HT1
COLD        IHCTP = 0
COLD        IF(IMZEROHT.EQ.0) CALL INUM_FOR_OCC(IHTOCC,IHCTP)
COLD        IF(IHCTP.NE.0) THEN
C            DO IT2 = 1, IT1
             DO IT2 = IT1, NSPOBEX
             IF(ISPOBEX_AC(IT2).EQ.1) THEN
              IT2_B = IBSPOBEX(IT2)
*. Occupation of HT1T2
              CALL OP_T_OCC(IHTOCC,ISPOBEX_TP(1,IT2),
     &            IHTTOCC,IMZEROHTT)
              CALL INUM_FOR_OCC(IHTTOCC,IHCTP)
              IF(IHCTP.GE.1) THEN
*. Overall factor
C     FAC_IN_BCH(NT,IT1,IT2,IT3,IT4)
                XBCH = FAC_IN_BCH(2,IT1,IT2,IDUMMY,IDUMMY)
COLD            FACTOR = FACTOR*FACTORH
                FACTOR = XBCH*FACTORH
                LCCFBLK = LSPOBEX(IHCTP)
                IF(NTEST.GE.10) 
     /          WRITE(6,*)
     &                 '2: HCT1234 will be called for IHCTP = ', IHCTP
                     CALL HCT1234(IHTP(1,IHOP),IHINDEX(1,IHOP),2,
     &               ISPOBEX_TP(1,IT1),ISPOBEX_TP(1,IT2),
     &               ISPOBEX_TP(1,IUNIOP),ISPOBEX_TP(1,IUNIOP),
     &               ITSM,ITSM,IONE,IONE,IHSM,
     &               T(IT1_B),T(IT2_B),ONE,ONE,
     &               CC_VEC_FNC(IBSPOBEX(IHCTP)),1,FACTOR,LCCFBLK)
              END IF
*             ^ End if HT1T2 is nonvanishing.
             END IF
*            ^ End if excitation IT2 was active
             END DO
*            ^ End of loop over IT2
COLD        END IF
*           ^ End of HT1 is nonvanishing
           END IF
*          ^ End if excitation IT1 was active
           END DO
*          ^ End of loop over IT1
          END IF
*         ^ End of NCOMMU .EQ. 2
*
          IF(NCOMMU.EQ.3) THEN
*. Loop over combinations of NCOMMU operators
           DO IT1 = 1, NSPOBEX
            IF(ISPOBEX_AC(IT1).EQ.1) THEN
            IT1_B = IBSPOBEX(IT1)
*. Occupation of HT1
            CALL OP_T_OCC(IHTP(1,IHOP),ISPOBEX_TP(1,IT1),
     &           IHTOCC,IMZEROHT)
*. address of HT1
COLD        IHCTP = 0
COLD        IF(IMZEROHT.EQ.0) CALL INUM_FOR_OCC(IHTOCC,IHCTP)
COLD        IF(IHCTP.NE.0) THEN
C            DO IT2 = 1, IT1
             DO IT2 = IT1, NSPOBEX
             IF(ISPOBEX_AC(IT2).EQ.1) THEN
              IT2_B = IBSPOBEX(IT2)
*. Occupation of HT1T2
              CALL OP_T_OCC(IHTOCC,ISPOBEX_TP(1,IT2),
     &            IHTTOCC,IMZEROHTT)
COLD          IHCTTP = 0
COLD          IF(IMZEROHTT.EQ.0) CALL INUM_FOR_OCC(IHTTOCC,IHCTTP)
COLD          IF(IHCTTP.GE.1) THEN
C              DO IT3 = 1, IT2
               DO IT3 = IT2, NSPOBEX
               IF(ISPOBEX_AC(IT3).EQ.1) THEN
                IT3_B = IBSPOBEX(IT3)
*. Occupation for HT1T2T3
                CALL OP_T_OCC(IHTTOCC,ISPOBEX_TP(1,IT3),
     &          IHTTTOCC,IMZEROHTTT)
                CALL INUM_FOR_OCC(IHTTTOCC,IHCTP)
                IF(IHCTP.GE.1) THEN
                 LCCFBLK = LSPOBEX(IHCTP)
*. Overall factor
                 XBCH = FAC_IN_BCH(3,IT1,IT2,IT3,IDUMMY)
COLD             FACTOR = FACTOR*FACTORH
                 FACTOR = XBCH*FACTORH
                 IF(NTEST.GE.10) 
     /           WRITE(6,*)
     &                  '3: HCT1234 will be called for IHCTP = ',IHCTP
                      CALL HCT1234(IHTP(1,IHOP),IHINDEX(1,IHOP),3,
     &                ISPOBEX_TP(1,IT1),ISPOBEX_TP(1,IT2),
     &                ISPOBEX_TP(1,IT3),ISPOBEX_TP(1,IUNIOP),
     &                ITSM,ITSM,ITSM,IONE,IHSM,
     &                T(IT1_B),T(IT2_B),T(IT3_B),ONE,
     &                CC_VEC_FNC(IBSPOBEX(IHCTP)),1,FACTOR,LCCFBLK)
                END IF
*               ^ End if IHCTP ne 0
               END IF
*              ^ End if excitation IT3 was active
               END DO
*              ^ End of loop over IT3
COLD          END IF
*             ^ End if HT1T2 is nonvanishing.
             END IF
*            ^ End if excitation IT2 was active
             END DO
*            ^ End of loop over IT2
COLD        END IF
*           ^ End of HT1 is nonvanishing
           END IF
*          ^ End if excitation IT1 was active
           END DO
*          ^ End of loop over IT1
          END IF
*         ^ End of NCOMMU .EQ. 3
*
          IF(NCOMMU.EQ.4) THEN
           DO IT1 = 1, NSPOBEX
            IF(ISPOBEX_AC(IT1).EQ.1) THEN
            IT1_B = IBSPOBEX(IT1)
*. Occupation of HT1
            CALL OP_T_OCC(IHTP(1,IHOP),ISPOBEX_TP(1,IT1),
     &           IHTOCC,IMZEROHT)
*. address of HT1
COLD        IHCTP = 0
COLD        IF(IMZEROHT.EQ.0) CALL INUM_FOR_OCC(IHTOCC,IHCTP)
COLD        IF(IHCTP.NE.0) THEN
C            DO IT2 = 1, IT1
             DO IT2 = IT1, NSPOBEX
             IF(ISPOBEX_AC(IT2).EQ.1) THEN
              IT2_B = IBSPOBEX(IT2)
*. Occupation of HT1T2
              CALL OP_T_OCC(IHTOCC,ISPOBEX_TP(1,IT2),
     &            IHTTOCC,IMZEROHTT)
COLD          IHCTTP = 0
COLD          IF(IMZEROHTT.EQ.0) CALL INUM_FOR_OCC(IHTTOCC,IHCTTP)
COLD          IF(IHCTTP.GE.1) THEN
C              DO IT3 = 1, IT2
               DO IT3 = IT2, NSPOBEX
               IF(ISPOBEX_AC(IT3).EQ.1) THEN
                IT3_B = IBSPOBEX(IT3)
*. Occupation for HT1T2T3
                CALL OP_T_OCC(IHTTOCC,ISPOBEX_TP(1,IT3),
     &          IHTTTOCC,IMZEROHTTT)
COLD            IHCTTTP = 0
COLD            IF(IMZEROHTTT.EQ.0) CALL INUM_FOR_OCC(IHTTTOCC,IHCTTTP)
COLD            IF(IHCTTTP.NE.0) THEN
C                DO IT4 = 1, IT3
                 DO IT4 = IT3, NSPOBEX
                 IF(ISPOBEX_AC(IT4).EQ.1) THEN
                  IT4_B = IBSPOBEX(IT4)
*. Occupation for HT1T2T3T4
                  CALL OP_T_OCC(IHTTTOCC,ISPOBEX_TP(1,IT4),
     &            IHTTTTOCC,IMZEROHTTTT)
*. Overall factor
                  CALL INUM_FOR_OCC(IHTTTTOCC,IHCTP)
                  IF(IHCTP.GE.1) THEN
                   XBCH = FAC_IN_BCH(4,IT1,IT2,IT3,IT4)
                   FACTOR = XBCH*FACTORH
COLD               FACTOR = FACTOR*FACTORH
                   IF(NTEST.GE.10) 
     &             WRITE(6,*)
     &                  '4: HCT1234 will be called for IHCTP=',IHCTP
                   CALL HCT1234(IHTP(1,IHOP),IHINDEX(1,IHOP),4,
     &                  ISPOBEX_TP(1,IT1),ISPOBEX_TP(1,IT2),
     &                  ISPOBEX_TP(1,IT3),ISPOBEX_TP(1,IT4),
     &                  ITSM,ITSM,ITSM,ITSM,IHSM,
     &                  T(IT1_B),T(IT2_B),T(IT3_B),T(IT4_B),
     &                  CC_VEC_FNC(IBSPOBEX(IHCTP)),1,FACTOR)
                  END IF
*                 ^ End of IHCTP ne. 0
                 END IF
*                ^ End if excitation IT4 was active
                 END DO
*                ^ End of loop over IT4
COLD            END IF
*                ^ End if IHCTTTP ne 0
               END IF
*              ^ End if excitation IT3 was active
               END DO
*              ^ End of loop over IT3
COLD          END IF
*             ^ End if HT1T2 is nonvanishing.
             END IF
*            ^ End if excitation IT2 was active
             END DO
*            ^ End of loop over IT2
COLD        END IF
*           ^ End of HT1 is nonvanishing
           END IF
*          ^ End if excitation IT1 was active
           END DO
*          ^ End of loop over IT1
          END IF
*         ^ End of NCOMMU .EQ. 4
*
          IF(NCOMMU.EQ.0) CALL LUCIAQEXIT('0_COM')
          IF(NCOMMU.EQ.1) CALL LUCIAQEXIT('1_COM')
          IF(NCOMMU.EQ.2) CALL LUCIAQEXIT('2_COM')
          IF(NCOMMU.EQ.3) CALL LUCIAQEXIT('3_COM')
          IF(NCOMMU.EQ.4) CALL LUCIAQEXIT('4_COM')
        END DO
*       ^ End of loop over commutators
      END DO
*     ^ End of loop over H-operators
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'EMTHET')
      RETURN
      END 
      SUBROUTINE OPCT1234M(IOEX,IO1DX,IO2DX,IO3DX,IO4DX,
     &                    IT1,IT2,IT3,IT4,T1,T2,T3,T4,OT1234,
     &                    IT1SM,IT2SM,IT3SM,IT4SM,IOPSM,IOPINDX,
     &                    ISIGN,FACX)
*
* Master routine for contracting operator with upto 4 operators
*
* Jeppe Olsen, April 2003 ( I hope it will be working before May ..)
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
C and the local scratch 
C     COMMON/CC_SCR2/KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
C    &               KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTOCC7,
C    &               KLSTOCC8,
C    &               KLSTREO,
C    &               KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
C    &               KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
C    &               KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
C    &               KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
C    &               KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
C    &               KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
C    &               KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
C    &               KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
C    &               KIX5_CA,KSX5_CA,KIX5_CB,KSX5_CB,
C    &               KIX5_AA,KSX5_AA,KIX5_AB,KSX5_AB,
C    &               KIX6_CA,KSX6_CA,KIX6_CB,KSX6_CB,
C    &               KIX6_AA,KSX6_AA,KIX6_AB,KSX6_AB,
C    &               KIX7_CA,KSX7_CA,KIX7_CB,KSX7_CB,
C    &               KIX7_AA,KSX7_AA,KIX7_AB,KSX7_AB,
C    &               KIX8_CA,KSX8_CA,KIX8_CB,KSX8_CB,
C    &               KIX8_AA,KSX8_AA,KIX8_AB,KSX8_AB,
C    &               KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
C    &               KLTSCR5,KLTSCR6,
C    &               KLOPSCR,
C    &               KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
C    &               KLSMD1,KLSMD2,KLSMD3,KLSMD4,
C    &               KLSMEX,KLSMK1,KLSMK2,KLSML1,
C    &               KLNMD1,KLNMD2,KLNMD3,KLNMD4,
C    &               KLNMEX,KLNMK1,KLNMK2,KLNML1,
C    &               KLOCK1, KLOCK2, KLOCK3,KLOCK4, 
C    &               KLOCL1, KLOCL2, KLOCL3,KLOCL4, KL_IBF, 
C    &               KLEXEORD
C from MEM
      COMMON/CC_SCR2/KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &               KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTOCC7,
     &               KLSTOCC8,KLSTREO,
     &               KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &               KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &               KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &               KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &               KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &               KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &               KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &               KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &               KIX5_CA,KSX5_CA,KIX5_CB,KSX5_CB,
     &               KIX5_AA,KSX5_AA,KIX5_AB,KSX5_AB,
     &               KIX6_CA,KSX6_CA,KIX6_CB,KSX6_CB,
     &               KIX6_AA,KSX6_AA,KIX6_AB,KSX6_AB,
     &               KIX7_CA,KSX7_CA,KIX7_CB,KSX7_CB,
     &               KIX7_AA,KSX7_AA,KIX7_AB,KSX7_AB,
     &               KIX8_CA,KSX8_CA,KIX8_CB,KSX8_CB,
     &               KIX8_AA,KSX8_AA,KIX8_AB,KSX8_AB,
     &               KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
     &               KLTSCR5,KLTSCR6,
     &               KLOPSCR,
     &               KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &               KLSMD1,KLSMD2,KLSMD3,KLSMD4,
     &               KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &               KLNMD1,KLNMD2,KLNMD3,KLNMD4,
     &               KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &               KLOCK1, KLOCK2, KLOCK3,KLOCK4, 
     &               KLOCL1, KLOCL2, KLOCL3,KLOCL4, KL_IBF, 
     &               KLEXEORD
      CALL OPCT1234(IOEX,IO1DX,IO2DX,IO3DX,IO4DX,
     &           IT1,IT2,IT3,IT4,T1,T2,T3,T4,OT1234,
     &           IT1SM,IT2SM,IT3SM,IT4SM,IOPSM,
     &           LCCBD12,LCCB,     
     &    WORK(KIX1_CA),WORK(KSX1_CA),WORK(KIX1_CB),WORK(KSX1_CB),
     &    WORK(KIX1_AA),WORK(KSX1_AA),WORK(KIX1_AB),WORK(KSX1_AB),
     &    WORK(KIX2_CA),WORK(KSX2_CA),WORK(KIX2_CB),WORK(KSX2_CB),
     &    WORK(KIX2_AA),WORK(KSX2_AA),WORK(KIX2_AB),WORK(KSX2_AB),
     &    WORK(KIX3_CA),WORK(KSX3_CA),WORK(KIX3_CB),WORK(KSX3_CB),
     &    WORK(KIX3_AA),WORK(KSX3_AA),WORK(KIX3_AB),WORK(KSX3_AB),
     &    WORK(KIX4_CA),WORK(KSX4_CA),WORK(KIX4_CB),WORK(KSX4_CB),
     &    WORK(KIX4_AA),WORK(KSX4_AA),WORK(KIX4_AB),WORK(KSX4_AB),
     &    WORK(KIX5_CA),WORK(KSX5_CA),WORK(KIX5_CB),WORK(KSX5_CB),
     &    WORK(KIX5_AA),WORK(KSX5_AA),WORK(KIX5_AB),WORK(KSX5_AB),
     &    WORK(KIX6_CA),WORK(KSX6_CA),WORK(KIX6_CB),WORK(KSX6_CB),
     &    WORK(KIX6_AA),WORK(KSX6_AA),WORK(KIX6_AB),WORK(KSX6_AB),
     &    WORK(KIX7_CA),WORK(KSX7_CA),WORK(KIX7_CB),WORK(KSX7_CB),
     &    WORK(KIX7_AA),WORK(KSX7_AA),WORK(KIX7_AB),WORK(KSX7_AB),
     &    WORK(KIX8_CA),WORK(KSX8_CA),WORK(KIX8_CB),WORK(KSX8_CB),
     &    WORK(KIX8_AA),WORK(KSX8_AA),WORK(KIX8_AB),WORK(KSX8_AB),
     &    WORK(KLOCK1),WORK(KLOCK2),WORK(KLOCK3),
     &    WORK(KLOCK4),WORK(KLOCL1),WORK(KLOCL2),WORK(KLOCL3),
     &    WORK(KLOCL4),
     &    WORK(KLSTOCC1),WORK(KLSTOCC2),WORK(KLSTOCC3),
     &    WORK(KLSTOCC4),WORK(KLSTOCC5),WORK(KLSTOCC6),
     &    WORK(KLSTOCC7),
     &    WORK(KLTSCR1),WORK(KLTSCR2),WORK(KLTSCR3),WORK(KLTSCR4),
     &    WORK(KLTSCR5),WORK(KLTSCR6),WORK(KLOPSCR),
     &    WORK(KLSMD1),WORK(KLSMD2),WORK(KLSMD3),WORK(KLSMD4),
     &    WORK(KLSMEX),WORK(KLSMK1),WORK(KLSML1),
     &    WORK(KLNMD1),WORK(KLNMD2),WORK(KLNMD3),WORK(KLNMD4),
     &    WORK(KLNMEX),WORK(KLNMK1),WORK(KLNML1),
     &    IOPINDX,
     &    WORK(KLZ),WORK(KLZSCR),WORK(KLSTREO),
     &    ISIGN,FACX,N_TDL_MAX)
*
      RETURN
      END
*
      SUBROUTINE OPCT1234(IOEX,IO1DX,IO2DX,IO3DX,IO4DX,
     &           IT1,IT2,IT3,IT4,T1,T2,T3,T4,OT1234,
     &           IT1SM,IT2SM,IT3SM,IT4SM,IOPSM,
     &           LDB,LB,
     &           IX1_CA,SX1_CA,IX1_CB,SX1_CB,
     &           IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &           IX2_CA,SX2_CA,IX2_CB,SX2_CB,
     &           IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &           IX3_CA,SX3_CA,IX3_CB,SX3_CB,
     &           IX3_AA,SX3_AA,IX3_AB,SX3_AB,
     &           IX4_CA,SX4_CA,IX4_CB,SX4_CB,
     &           IX4_AA,SX4_AA,IX4_AB,SX4_AB,
     &           IX5_CA,SX5_CA,IX5_CB,SX5_CB,
     &           IX5_AA,SX5_AA,IX5_AB,SX5_AB,
     &           IX6_CA,SX6_CA,IX6_CB,SX6_CB,
     &           IX6_AA,SX6_AA,IX6_AB,SX6_AB,
     &           IX7_CA,SX7_CA,IX7_CB,SX7_CB,
     &           IX7_AA,SX7_AA,IX7_AB,SX7_AB,
     &           IX8_CA,SX8_CA,IX8_CB,SX8_CB,
     &           IX8_AA,SX8_AA,IX8_AB,SX8_AB,
     &           IOC_K1,IOC_K2,IOC_K3,IOC_K4,
     &           IOC_L1,IOC_L2,IOC_L3,IOC_L4,
     &           ISTR_D1,ISTR_D2,ISTR_D3,ISTR_D4,ISTR_EX,ISTR_K,ISTR_L,
     &           TSCR1,TSCR2,TSCR3,
     &           TSCR4,TSCR5,TSCR6,OPSCR,
     &           ISM_CAAB_D1,ISM_CAAB_D2,ISM_CAAB_D3,ISM_CAAB_D4,
     &           ISM_CAAB_EX,ISM_CAAB_K,ISM_CAAB_L,
     &           INM_CAAB_D1,INM_CAAB_D2,INM_CAAB_D3,INM_CAAB_D4,
     &           INM_CAAB_EX,INM_CAAB_K,INM_CAAB_L,
     &           IEXD1234_INDX,
     &           IZ,IZSCR,ISTREO,
     &           ISIGNG,FACX,N_TDL_MAX)
*
* Contract indeces of Operator O with indeces of 
* excitation operator T1,T2, T3, T4
*
* Operator O is defined by an operator part IOEX, and 
* deexcitation parts IO1DX,IO2DX,IO3DX,IO4DX.
*
* LDB : Batch size for D1, D2, D3, D4
* LB2   : Batch size for for remaining expansions
*
* Jeppe Olsen, May 2000, Finished March 2003 ( Well it became July 2003)
*                        and even November 2004
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'ctcc.inc'
*. Input   
      INTEGER IOEX(NGAS,4),IO1DX(NGAS,4),IO2DX(NGAS,4)
      INTEGER IO3DX(NGAS,4),IO4DX(NGAS,4)
      INTEGER IT1(NGAS,4),IT2(NGAS,4),IT3(NGAS,4),IT4(NGAS,4)
      DIMENSION T1(*),T2(*),T3(*),T4(*)
*. Map of index EXD1234 operator to original index of Hamiltonian
      INTEGER IEX1234_INDX(4)
*. Output
      DIMENSION OT1234(*)
*
*. Local scratch
*
*. Occupation of conjugated operators
      INTEGER IO1DX_DAG(4*MXPNGAS),IO2DX_DAG(4*MXPNGAS)
      INTEGER IO3DX_DAG(4*MXPNGAS),IO4DX_DAG(4*MXPNGAS)
*. Occupation of OP T1 T2
C     INTEGER IOC_OT1T2(NGAS,4)
*
*. Occupation of gasspaces for various strings 
*. Number of strings per sym of the various string supergroups
      INTEGER NOEX(8,4),NO1DX(8,4),NO2DX(8,4),NO3DX(8,4),NO4DX(8,4)
      INTEGER NT1(8,4), NT2(8,4), NT3(8,4), NT4(8,4)
      INTEGER NK1(8,4),NK2(8,4),NK3(8,4),NK4(8,4)
      INTEGER NL1(8,4),NL2(8,4),NL3(8,4),NL4(8,4)
C     INTEGER NOT1T2(8,4)
*. Offsets of strings with given sym for the strings of various CC ops
      INTEGER IBOEX(8,4),IBO1DX(8,4),IBO2DX(8,4),IBO3DX(8,4)
      INTEGER IBO4DX(8,4)
*. Offset in operators to strings with given sym
      INTEGER IBT1_TCC(8,8,8),IBT2_TCC(8,8,8)
      INTEGER IBT3_TCC(8,8,8),IBT4_TCC(8,8,8)
      INTEGER IBL1_TCC(8,8,8),IBL2_TCC(8,8,8)
      INTEGER IBL3_TCC(8,8,8),IBL4_TCC(8,8,8)
C     INTEGER IBOT1T2_TCC(8,8,8)
*
C     INTEGER IBOT1T2_TCC(8,8,8)
*
      INTEGER IB_D1K1(8,8,4),IB_D2K2(8,8,4) 
      INTEGER IB_D3K3(8,8,4),IB_D4K4(8,8,4)
      INTEGER IB_EXK1(8,8,4),IB_L1K2(8,8,4)
      INTEGER IB_L2K3(8,8,4),IB_L3K4(8,8,4)
*
*
*. Scratch through parameter list. 
*
*. IX1_* : Number of operators in operator * Largest C or A string
      DIMENSION IX1_CA(*),SX1_CA(*),IX1_CB(*),SX1_CB(*)
      DIMENSION IX1_AA(*),SX1_AA(*),IX1_AB(*),SX1_AB(*)
*. IX2, SX2
      DIMENSION IX2_CA(*),SX2_CA(*),IX2_CB(*),SX2_CB(*)
      DIMENSION IX2_AA(*),SX2_AA(*),IX2_AB(*),SX2_AB(*)
*. IX3, SX3
      DIMENSION IX3_CA(*),SX3_CA(*),IX3_CB(*),SX3_CB(*)
      DIMENSION IX3_AA(*),SX3_AA(*),IX3_AB(*),SX3_AB(*)
*. IX4, SX4
      DIMENSION IX4_CA(*),SX4_CA(*),IX4_CB(*),SX4_CB(*)
      DIMENSION IX4_AA(*),SX4_AA(*),IX4_AB(*),SX4_AB(*)
*. IX5, SX5
      DIMENSION IX5_CA(*),SX5_CA(*),IX5_CB(*),SX5_CB(*)
      DIMENSION IX5_AA(*),SX5_AA(*),IX5_AB(*),SX5_AB(*)
*. IX6, SX6
      DIMENSION IX6_CA(*),SX6_CA(*),IX6_CB(*),SX6_CB(*)
      DIMENSION IX6_AA(*),SX6_AA(*),IX6_AB(*),SX6_AB(*)
*. IX7, SX7
      DIMENSION IX7_CA(*),SX7_CA(*),IX7_CB(*),SX7_CB(*)
      DIMENSION IX7_AA(*),SX7_AA(*),IX7_AB(*),SX7_AB(*)
*. IX8, SX8
      DIMENSION IX8_CA(*),SX8_CA(*),IX8_CB(*),SX8_CB(*)
      DIMENSION IX8_AA(*),SX8_AA(*),IX8_AB(*),SX8_AB(*)
*
*. ISTR_D1, ISTR_K, ISTR_L : for occupations of strings  
      INTEGER ISTR_D1(*), ISTR_D2(*), ISTR_D3(*), ISTR_D4(*)
      INTEGER ISTR_EX(*)
      INTEGER ISTR_K(*), ISTR_L(*) 
*. For occupation of intermediate strings
      INTEGER IOC_K1(NGAS,4),IOC_K2(NGAS,4)
      INTEGER IOC_K3(NGAS,4),IOC_K4(NGAS,4)
      INTEGER IOC_L1(NGAS,4),IOC_L2(NGAS,4)
      INTEGER IOC_L3(NGAS,4),IOC_L4(NGAS,4)
*. For intermediates with both strings batched 
      DIMENSION TSCR1(LDB*LDB),TSCR2(LDB*LDB)
      DIMENSION TSCR3(LDB*LDB)
*. For an intermediate with only one string batched 
      DIMENSION TSCR4(N_TDL_MAX), TSCR5(N_TDL_MAX),TSCR6(N_TDL_MAX)
*. For a batch of coefficients for Operator
*. Operator accessed as OP(EX,D2,D1) so
      DIMENSION OPSCR(LDB*LDB*LB)
*. For part of Hamiltonian
      INTEGER ISM_CAAB_D1(4,*), INM_CAAB_D1(4,*)
      INTEGER ISM_CAAB_D2(4,*), INM_CAAB_D2(4,*)
      INTEGER ISM_CAAB_D3(4,*), INM_CAAB_D3(4,*)
      INTEGER ISM_CAAB_D4(4,*), INM_CAAB_D4(4,*)
      INTEGER ISM_CAAB_EX(4,*), INM_CAAB_EX(4,*)
*. For CC operators
      INTEGER ISM_CAAB_K(4,*), INM_CAAB_K(4,*)
      INTEGER ISM_CAAB_L(4,*), INM_CAAB_L(4,*)
*
      INTEGER IZ(*), IZSCR(*), ISTREO(*)
*. Collecting info on largest T(D,L) Block
      COMMON/CNTDL_MAX/NTDL_MAX_ACT,ID2_MAX_ACT(MXPNGAS*4),
     &      IL1_MAX_ACT(MXPNGAS*4),IEX_MAX_ACT(MXPNGAS*4),
     &      IT1_MAX_ACT(MXPNGAS*4),IT2_MAX_ACT(MXPNGAS*4)
*
* The story goes as
* Loop over symmetry of D4
* Loop over batches  of D4
*  Loop over symmetry of D3
*  Loop over batches  of D3
*   Loop over symmetry of D2
*   Loop over batches  of D2
*
* Part 1 : O, T1 => OT1(d2,d3,d4,l1) in TSCR4
*
*    Loop over symmetry of D1
*    Loop over batches  of D1
*
*    Loop over batches of Ex
*      Fetch O as O(d1,d2,d3,d4,ex)
*      Loop over batches of K1
*        T1(I) => T1(d1,k1)
*        OT1(d2,d3,d4,ex,k1) = sum(d1) O(d1,d2,d3,d4,ex)T1(d1,k1)
*        OT1(d2,d3,d4,ex,k1) => OT1(d2,d3,d4,l1)
*      End of loop over batches of K1
*    End of loop over bathces of Ex
*    End of loop over batches of d1
*    End of loop over symmetry of D1
* ( We now have OT1(l1,d4,d3,d2), all l1 !, d2,d3,d4 in batch)    
*
* Part 2 : OT1 * T2 => OT12(d3,d4,l2) in TSCR5
*
*    Loop over batches of k2
*      T2(i) => T2(d2,k2)
*      Loop over batches of l1
*        OT12(d3,d4,l1,k2) = OT1(d2,d3,d4,l1)*T2(d2,k2) 
*        OT12(d3,d4,l1,k2) => OT12(d3,d4,l2)
*      End of loop over batches of l1
*    End of loop over batches of K2
*   End of loop over batches  of D2
*   End of loop over symmetry of D2
*
* Part 3 : OT12 * T3 => OT123(d4,l3) in TSCR6
*
*   Loop over batches of K3
*     T3(I) => T2(D3,K3)
*     Loop over batches of L2
*       OT123(d4,l2,k3) = sum(d3) OT12(d3,d4,l2) T2(d3,k3) 
*       OT123(d4,l2,k3) => OT123(d4,l3)
*     End of loop over batches of L2
*   End of loop over batches of K3
*  End of loop over batches  of D3
*  End of loop over symmetry of D3
*
* Part 4 : OT123 * T4 => OT1234(I)
*
*  Loop over batches of K4
*    T4(I) => T4(d4,k4)
*    Loop over batches of L3
*     OT1234(l3,k4) = OT123(d4,l3)T4(d4,k4)
*     OT1234(l3,k4) => OT1234(I)
*    End of loop over batches of L3
*  End of loop over batches of K4
* End of loop over batches of D4
* End of loop over symmetry of D4
*
* - That's all she wrote
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'OPCTT ')
*
      ZERO = 0.0D0
      ONE = 1.0D0
*
      NTEST = 00 
      IF(NTEST.GE.100 )THEN
        WRITE(6,*) ' OPCT1234 entered '
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Ex Deex1 Deex2 Deex3 Deex4 ' 
        WRITE(6,*)
        CALL WRT_SPOX_TP(IOEX,1)
        CALL WRT_SPOX_TP(IO1DX,1)
        CALL WRT_SPOX_TP(IO2DX,1)
        CALL WRT_SPOX_TP(IO3DX,1)
        CALL WRT_SPOX_TP(IO4DX,1)
        WRITE(6,*) 
        WRITE(6,*) ' Form of T1 T2 T3 T4'
        CALL WRT_SPOX_TP(IT1,1)
        CALL WRT_SPOX_TP(IT2,1)
        CALL WRT_SPOX_TP(IT3,1)
        CALL WRT_SPOX_TP(IT4,1)
*
      END IF
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' N_TDL_MAX in OPCT1T2 ', N_TDL_MAX
      SIGNG = DBLE(ISIGNG)
*. Symmetry of final type  
      IF(NTEST.GE.100) WRITE(6,*) ' IT1SM, IT2SM, IT3SM, IT4SM : ',
     &             IT1SM, IT2SM, IT3SM, IT4SM
      IT12SM = MULTD2H(IT1SM,IT2SM)
      IT123SM= MULTD2H(IT12SM,IT3SM)
      IT1234SM= MULTD2H(IT123SM,IT4SM)
      IL4SM= MULTD2H(IOPSM,IT1234SM)
      IF(NTEST.GE.100) WRITE(6,*) 'IT12SM,IT123SM,IT1234SM,IL4SM ',
     &            IT12SM,IT123SM,IT1234SM,IL4SM
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Symmetry of OT1234 = ', IL4SM
      END IF
*
* =================================================
* Occupation of the various gas spaces for strings 
* =================================================
*
*. K1 = O1DX  I1 
      CALL NEW_CAAB_OC(IOC_K1,IT1,IO1DX,1,1,NGAS)
*. L1 = OEX K1 
      CALL NEW_CAAB_OC(IOC_L1,IOC_K1,IOEX,2,1,NGAS)
*. K2 = O2DX I2
      CALL NEW_CAAB_OC(IOC_K2,IT2,IO2DX,1,1,NGAS)
*. L2 = K2 L1
      CALL NEW_CAAB_OC(IOC_L2,IOC_K2,IOC_L1,2,1,NGAS)
*. K3 = O3DX I3
      CALL NEW_CAAB_OC(IOC_K3,IT3,IO3DX,1,1,NGAS)
*. L3 = K3 L2
      CALL NEW_CAAB_OC(IOC_L3,IOC_K3,IOC_L2,2,1,NGAS)
*. K4 = O2DX I2
      CALL NEW_CAAB_OC(IOC_K4,IT4,IO4DX,1,1,NGAS)
*. L4 = K4 L3
      CALL NEW_CAAB_OC(IOC_L4,IOC_K4,IOC_L3,2,1,NGAS)
*. (L4 is final type of OT1234)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation of K1, K2, K3, K4 '
        CALL WRT_SPOX_TP(IOC_K1,1)
        CALL WRT_SPOX_TP(IOC_K2,1)
        CALL WRT_SPOX_TP(IOC_K3,1)
        CALL WRT_SPOX_TP(IOC_K4,1)
        WRITE(6,*) ' Occupation of L1, L2, L3, L4 '
        CALL WRT_SPOX_TP(IOC_L1,1)
        CALL WRT_SPOX_TP(IOC_L2,1)
        CALL WRT_SPOX_TP(IOC_L3,1)
        CALL WRT_SPOX_TP(IOC_L4,1)
      END IF
*
* =============================================================
*. Obtain symmetry-dimensions and -offsets for various strings
* =============================================================
*
      DO I_CAAB = 1, 4
*. OEX 
       CALL NST_SPGP(IOEX(1,I_CAAB),NOEX(1,I_CAAB))
       CALL ZBASE(NOEX(1,I_CAAB),IBOEX(1,I_CAAB),NSMST)
C      ZBASE(NVEC,IVEC,NCLASS)
*. D1
       CALL NST_SPGP(IO1DX(1,I_CAAB),NO1DX(1,I_CAAB))
       CALL ZBASE(NO1DX(1,I_CAAB),IBO1DX(1,I_CAAB),NSMST)
*. D2
       CALL NST_SPGP(IO2DX(1,I_CAAB),NO2DX(1,I_CAAB))
       CALL ZBASE(NO2DX(1,I_CAAB),IBO2DX(1,I_CAAB),NSMST)
*. D3
       CALL NST_SPGP(IO3DX(1,I_CAAB),NO3DX(1,I_CAAB))
       CALL ZBASE(NO3DX(1,I_CAAB),IBO3DX(1,I_CAAB),NSMST)
*. D4
       CALL NST_SPGP(IO4DX(1,I_CAAB),NO4DX(1,I_CAAB))
       CALL ZBASE(NO4DX(1,I_CAAB),IBO4DX(1,I_CAAB),NSMST)
*. T1
       CALL NST_SPGP(IT1(1,I_CAAB),NT1(1,I_CAAB))
*. T2
       CALL NST_SPGP(IT2(1,I_CAAB),NT2(1,I_CAAB))
*. T3
       CALL NST_SPGP(IT3(1,I_CAAB),NT3(1,I_CAAB))
*. T4
       CALL NST_SPGP(IT4(1,I_CAAB),NT4(1,I_CAAB))
*. K1
       CALL NST_SPGP(IOC_K1(1,I_CAAB),NK1(1,I_CAAB))
*. K2
       CALL NST_SPGP(IOC_K2(1,I_CAAB),NK2(1,I_CAAB))
*. K3
       CALL NST_SPGP(IOC_K3(1,I_CAAB),NK3(1,I_CAAB))
*. K4
       CALL NST_SPGP(IOC_K4(1,I_CAAB),NK4(1,I_CAAB))
*. L1   
       CALL NST_SPGP(IOC_L1(1,I_CAAB),NL1(1,I_CAAB))
*. L2   
       CALL NST_SPGP(IOC_L2(1,I_CAAB),NL2(1,I_CAAB))
*. L3   
       CALL NST_SPGP(IOC_L3(1,I_CAAB),NL3(1,I_CAAB))
*. L4   
       CALL NST_SPGP(IOC_L4(1,I_CAAB),NL4(1,I_CAAB))
*
      END DO
      CALL MEMCHK2('A_NST_')
*. We now have the various dimensions, so we can write T1,T2,T3,T4 if 
*. required
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input T1 block '
        CALL WRT_TCC_BLK(T1,IT1SM,NT1(1,1),NT1(1,2),NT1(1,3),NT1(1,4),
     &                   NSMST)
        WRITE(6,*) ' Input T2 block '
        CALL WRT_TCC_BLK(T2,IT2SM,NT2(1,1),NT2(1,2),NT2(1,3),NT2(1,4),
     &                   NSMST)
        WRITE(6,*) ' Input T3 block '
        CALL WRT_TCC_BLK(T3,IT3SM,NT3(1,1),NT3(1,2),NT3(1,3),NT3(1,4),
     &                   NSMST)
        WRITE(6,*) ' Input T4 block '
        CALL WRT_TCC_BLK(T4,IT4SM,NT4(1,1),NT4(1,2),NT4(1,3),NT4(1,4),
     &                   NSMST)
        WRITE(6,*) ' Input OT1234 block '
        CALL WRT_TCC_BLK(OT1234,IL4SM,NL4(1,1),NL4(1,2),NL4(1,3),
     &        NL4(1,4),NSMST)
      END IF
*
*
* ============================
*. Offsets to T1, T2, T3, T4
* ============================
*
*. T1
      CALL Z_TCC_OFF2(IBT1_TCC,LEN_T1,NT1(1,1),NT1(1,2),NT1(1,3),
     &                NT1(1,4),IT1SM,NSMST)
*. T2
      CALL Z_TCC_OFF2(IBT2_TCC,LEN_T2,NT2(1,1),NT2(1,2),NT2(1,3),
     &                NT2(1,4),IT2SM,NSMST)
*. T3
      CALL Z_TCC_OFF2(IBT3_TCC,LEN_T3,NT3(1,1),NT3(1,2),NT3(1,3),
     &                NT3(1,4),IT3SM,NSMST)
*. T4
      CALL Z_TCC_OFF2(IBT4_TCC,LEN_T4,NT4(1,1),NT4(1,2),NT4(1,3),
     &                NT4(1,4),IT4SM,NSMST)
*. The following has been moved as symmetry of L1 - L4 is not known here
*. L1
C     CALL Z_TCC_OFF2(IBL1_TCC,LEN_L1,NL1(1,1),
C    &     NL1(1,2),NL1(1,3),NL1(1,4),IL1,NSMST)
*. L2
C     CALL Z_TCC_OFF2(IBL2_TCC,LEN_L2,NL2(1,1),
C    &     NL2(1,2),NL2(1,3),NL2(1,4),IL2,NSMST)
*. L4
C     CALL Z_TCC_OFF2(IBL3_TCC,LEN_L3,NL3(1,1),
C    &     NL3(1,2),NL3(1,3),NL3(1,4),IL3,NSMST)
*. L4
C     CALL Z_TCC_OFF2(IBL4_TCC,LEN_L4,NL4(1,1),
C    &     NL4(1,2),NL4(1,3),NL4(1,4),IL4,NSMST)
*
*
* ========================
*.  D1{\dagger} K1  => T1 mapping
* =========================
*
*. Obtain D1{\dagger}
      CALL CONJ_CAAB(IO1DX,IO1DX_DAG,NGAS,SP_D1K1)
      CALL T1T2_TO_T12_MAP(IO1DX_DAG,IOC_K1,IT1,
     &     IX1_CA,SX1_CA,IX1_CB,SX1_CB,IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &     IB_D1K1,ISTR_D1,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D1K1)
*
* ======================
* Ex K1 => L1 mapping
* ======================
*
      CALL T1T2_TO_T12_MAP(IOEX,IOC_K1,IOC_L1,
     &     IX2_CA,SX2_CA,IX2_CB,SX2_CB,IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &     IB_EXK1,ISTR_K,ISTR_L,ISTREO,IZ,IZSCR,SIGN_EXK1)
*
* ========================
*.  D2{\dagger} K2  => T2 mapping
* =========================
*
*. Obtain D2{\dagger}
      CALL CONJ_CAAB(IO2DX,IO2DX_DAG,NGAS,SP_D2K2)
      CALL T1T2_TO_T12_MAP(IO2DX_DAG,IOC_K2,IT2,
     &     IX3_CA,SX3_CA,IX3_CB,SX3_CB,IX3_AA,SX3_AA,IX3_AB,SX3_AB,
     &     IB_D2K2,ISTR_D2,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D2K2)
*
* ==========
* L1 K2 => L2
* ==========
*
      CALL T1T2_TO_T12_MAP(IOC_L1,IOC_K2,IOC_L2,
     &     IX4_CA,SX4_CA,IX4_CB,SX4_CB,IX4_AA,SX4_AA,IX4_AB,SX4_AB,
     &     IB_L1K2,ISTR_L,ISTR_K,ISTREO,IZ,IZSCR,SIGN_L1K2)
C      T1T2_TO_T12_MAP(I1SPOBEX,I2SPOBEX,I12SPOBEX,
C    &                       ICA_MAP,XCA_MAP,ICB_MAP,XCB_MAP,
C    &                       IAA_MAP,XAA_MAP,IAB_MAP,XAB_MAP,
C    &                       IB,I1OCC,I2OCC,I1REO,IZ,IZSCR,SIGNP) 
*
* ========================
*.  D3{\dagger} K3  => T3 mapping
* =========================
*
*. Obtain D3{\dagger}
      CALL CONJ_CAAB(IO3DX,IO3DX_DAG,NGAS,SP_D3K3)
      CALL T1T2_TO_T12_MAP(IO3DX_DAG,IOC_K3,IT3,
     &     IX5_CA,SX5_CA,IX5_CB,SX5_CB,IX5_AA,SX5_AA,IX5_AB,SX5_AB,
     &     IB_D3K3,ISTR_D3,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D3K3)
*
* ==========
* L2 K3 => L3
* ==========
*
      CALL T1T2_TO_T12_MAP(IOC_L2,IOC_K3,IOC_L3,
     &     IX6_CA,SX6_CA,IX6_CB,SX6_CB,IX6_AA,SX6_AA,IX6_AB,SX6_AB,
     &     IB_L2K3,ISTR_L,ISTR_K,ISTREO,IZ,IZSCR,SIGN_L2K3)
* ========================
*.  D4{\dagger} K4  => T4 mapping
* =========================
*
*. Obtain D4{\dagger}
      CALL CONJ_CAAB(IO4DX,IO4DX_DAG,NGAS,SP_D4K4)
      CALL T1T2_TO_T12_MAP(IO4DX_DAG,IOC_K4,IT4,
     &     IX7_CA,SX7_CA,IX7_CB,SX7_CB,IX7_AA,SX7_AA,IX7_AB,SX7_AB,
     &     IB_D4K4,ISTR_D4,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D4K4)
*
* ==========
* L3 K4 => L4
* ==========
*
      CALL T1T2_TO_T12_MAP(IOC_L3,IOC_K4,IOC_L4,
     &     IX8_CA,SX8_CA,IX8_CB,SX8_CB,IX8_AA,SX8_AA,IX8_AB,SX8_AB,
     &     IB_L3K4,ISTR_L,ISTR_K,ISTREO,IZ,IZSCR,SIGN_L3K4)
      CALL MEMCHK2('C     ')
*

*. And the individual strings : all symmetries constructed
*. D1 strings 
      CALL STR_CAAB(IO1DX,ISTR_D1)
*. D2 strings 
      CALL STR_CAAB(IO2DX,ISTR_D2)
*. D3 strings 
      CALL STR_CAAB(IO3DX,ISTR_D3)
*. D4 strings 
      CALL STR_CAAB(IO4DX,ISTR_D4)
*. Ex strings 
      CALL STR_CAAB(IOEX,ISTR_EX)
*. Batch length : Should allow to vary 
      LD1_BAT = LDB
      LD2_BAT = LDB
      LD3_BAT = LDB
      LD4_BAT = LDB
      LEX_BAT = LDB
      LKL_BAT = LDB 
*. Symmetry of Op T1, Op T1 T2, Op T1 T2 T3, Op T1 T2 T3 T4
      IOPT1SM    = MULTD2H(IOPSM,IT1SM)
      IOPT12SM   = MULTD2H(IOPT1SM,IT2SM)
      IOPT123SM  = MULTD2H(IOPT12SM,IT3SM)
      IOPT1234SM = MULTD2H(IOPT123SM,IT4SM)
*
      NL4_TOT = LEN_TCCBLK(NL4(1,1),NL4(1,2),NL4(1,3),
     &                 NL4(1,4),IOPT1234SM,NSMST)
*
      DO ID4SM = 1, NSMST
       K4SM = MULTD2H(ID4SM,IT4SM)
*. By def L4 = final operator so 
       L4SM = IOPT1234SM
       NL4_TOT = LEN_TCCBLK(NL4(1,1),NL4(1,2),NL4(1,3),
     &                NL4(1,4),L4SM,NSMST)
       ID4SM_NEW = 1
*. Offsets for L4 with sym L4SM
       CALL Z_TCC_OFF2(IBL4_TCC,LEN_L4,NL4(1,1),
     &      NL4(1,2),NL4(1,3),NL4(1,4),L4SM,NSMST)
*. Number of D4 strings with actual symmetry
       ND4_TOT = LEN_TCCBLK(NO4DX(1,1),NO4DX(1,2),NO4DX(1,3),
     &               NO4DX(1,4),ID4SM,NSMST)
*. Number of K4 strings with actual symmetry
       NK4_TOT = LEN_TCCBLK(NK4(1,1),NK4(1,2),NK4(1,3),
     &                 NK4(1,4),K4SM,NSMST)
*. Number of D4 batches 
       ND4_BAT = ND4_TOT/LD4_BAT
       IF(ND4_BAT*LD4_BAT.LT.ND4_TOT) ND4_BAT=ND4_BAT+1
*. And loop over batches of D4
       DO ID4_BAT = 1, ND4_BAT
        ID4_START = (ID4_BAT-1)*LD4_BAT + 1
        ID4_STOP  = MIN(ND4_TOT,ID4_START+LD4_BAT-1)
        ID4_BATLEN = ID4_STOP-ID4_START+1
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' ID4SM, ID4_BAT, ID4_BATLEN = ',
     &                 ID4SM, ID4_BAT, ID4_BATLEN
        END IF
*. Generate D4 strings for given sym and batch
        CALL ISMNM_FOR_TCC_BAT(NO4DX,ISM_CAAB_D4,INM_CAAB_D4,ID4SM,
     &       ID4_BATLEN,ID4SM_NEW,
     &       ISM_C1_D4,ISM_CA1_D4,ISM_AA1_D4,
     &       INM_AB1_D4,INM_AA1_D4,INM_CA1_D4,INM_CB1_D4,
     &       ISM_CINI_D4,ISM_CAINI_D4,ISM_AAINI_D4,      
     &       INM_ABINI_D4,INM_AAINI_D4,INM_CAINI_D4,INM_CBINI_D4,0)
        ID4SM_NEW = 0
* The results of this loop is the generation of OT123(d4,l3) in TSCR6
        L3SM = MULTD2H(IOPT123SM,ID4SM)
*. Number of L3 strings with actual symmetry
        NL3_TOT = LEN_TCCBLK(NL3(1,1),NL3(1,2),NL3(1,3),
     &                 NL3(1,4),L3SM,NSMST)
*. Offsets for L3 with sym L3SM
        CALL Z_TCC_OFF2(IBL3_TCC,LEN_L3,NL3(1,1),
     &       NL3(1,2),NL3(1,3),NL3(1,4),L3SM,NSMST)
        LENGTH = ID4_BATLEN*NL3_TOT
        IF(LENGTH.GT.N_TDL_MAX) THEN
           WRITE(6,*) ' N_TDL_MAX too small '
           WRITE(6,*) ' Current, required = ', N_TDL_MAX, LENGTH
           WRITE(6,*) 'ID4_BATLEN,NL3_TOT = ', ID4_BATLEN, NL3_TOT
           STOP ' N_TDL_MAX too small '
        END IF
        CALL SETVEC(TSCR6,ZERO,LENGTH)
*. Loop over D3
        DO ID3SM = 1, NSMST
         K3SM = MULTD2H(ID3SM,IT3SM)
         IK3SM_NEW = 1
         ID3SM_NEW = 1
         IF(NTEST.GE.100)  WRITE(6,*) ' ID3SM,K3SM = ',ID3SM, K3SM
*. Number of D3 strings with actual symmetry
         ND3_TOT = LEN_TCCBLK(NO3DX(1,1),NO3DX(1,2),NO3DX(1,3),
     &                 NO3DX(1,4),ID3SM,NSMST)
*. Number of K3 strings with actual symmetry
         NK3_TOT = LEN_TCCBLK(NK3(1,1),NK3(1,2),NK3(1,3),
     &                 NK3(1,4),K3SM,NSMST)
*. Number of D3 batches 
         IF(NK3_TOT.GT.0 .AND. NL3_TOT.GT.0) THEN
         ND3_BAT = ND3_TOT/LD3_BAT
         IF(ND3_BAT*LD3_BAT.LT.ND3_TOT) ND3_BAT=ND3_BAT+1
*. And loop over batches of D3
       IF(NTEST.GE.100) WRITE(6,*) 'ND3_TOT, ND3_BAT ',ND3_TOT,ND3_BAT
         DO ID3_BAT = 1, ND3_BAT
          IF(NTEST.GE.100) WRITE(6,*) ' ID3_BAT = ', ID3_BAT
          ID3_START = (ID3_BAT-1)*LD3_BAT + 1
          ID3_STOP  = MIN(ND3_TOT,ID3_START+LD3_BAT-1)
          ID3_BATLEN = ID3_STOP-ID3_START+1
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' ID3SM, ID3_BAT, ID3_BATLEN = ',
     &                   ID3SM, ID3_BAT, ID3_BATLEN
          END IF
*. Generate D3 strings for given sym and batch
          CALL ISMNM_FOR_TCC_BAT(NO3DX,ISM_CAAB_D3,INM_CAAB_D3,ID3SM,
     &         ID3_BATLEN,ID3SM_NEW,
     &         ISM_C1_D3,ISM_CA1_D3,ISM_AA1_D3,
     &         INM_AB1_D3,INM_AA1_D3,INM_CA1_D3,INM_CB1_D3,
     &         ISM_CINI_D3,ISM_CAINI_D3,ISM_AAINI_D3,      
     &         INM_ABINI_D3,INM_AAINI_D3,INM_CAINI_D3,INM_CBINI_D3,0)
          ID3SM_NEW = 0
*. In the loops inside, we will construct OT1T2(d3,d4,l2) in TSCR5 
*  by looping over symmetries and batches of D2 and D1. 
          ID34SM = MULTD2H(ID3SM,ID4SM)
          L2SM   = MULTD2H(IOPT12SM,ID34SM)
*. Number of L2 strings with sym L2SM
          NL2_TOT = LEN_TCCBLK(NL2(1,1),NL2(1,2),NL2(1,3),
     &                  NL2(1,4),L2SM,NSMST)
*. Offsets for L2 with sym L2SM
          CALL Z_TCC_OFF2(IBL2_TCC,LEN_L2,NL2(1,1),
     &         NL2(1,2),NL2(1,3),NL2(1,4),L2SM,NSMST)
          LENNY = ID3_BATLEN*ID4_BATLEN*NL2_TOT
          IF(NTEST.GE.100) WRITE(6,*) ' L2SM, NL2_TOT ',
     &                                  L2SM, NL2_TOT
          IF(LENNY .GT. N_TDL_MAX ) THEN
*. Problem, too small length for TSCR5 
             WRITE(6,*) ' Dimension of TSCR5 too small '
             WRITE(6,*) ' Required and allocated : ',LENNY, N_TDL_MAX 
             WRITE(6,*) 'NL2_TOT, ID3_BATLEN, ID4_BATLEN = ',
     &                   NL2_TOT, ID3_BATLEN, ID4_BATLEN
             STOP ' Dimension of TSCR5 too small '
          END IF
          CALL SETVEC(TSCR5,ZERO,ID3_BATLEN*ID4_BATLEN*NL2_TOT)
          DO ID2SM = 1, NSMST
           K2SM = MULTD2H(ID2SM,IT2SM)
*. The symmetry of L1 obtained from OT1(l1,d4,d3,d2) has symmetry OPT1SM
           ID234SM = MULTD2H(ID2SM,ID34SM)
           L1SM = MULTD2H(IOPT1SM,ID234SM)
*. Number of L1 strings with sym L1SM
           NL1_TOT = LEN_TCCBLK(NL1(1,1),NL1(1,2),NL1(1,3),
     &                   NL1(1,4),L1SM,NSMST)
*. Offsets for L1 with sym L1SM
           CALL Z_TCC_OFF2(IBL1_TCC,LEN_L1,NL1(1,1),
     &          NL1(1,2),NL1(1,3),NL1(1,4),L1SM,NSMST)
*
           IK2SM_NEW = 1
           ID2SM_NEW = 1
           IF(NTEST.GE.100)  WRITE(6,*) ' ID2SM,K2SM = ',ID2SM, K2SM
*. Number of D2 strings with actual symmetry
           ND2_TOT = LEN_TCCBLK(NO2DX(1,1),NO2DX(1,2),NO2DX(1,3),
     &                   NO2DX(1,4),ID2SM,NSMST)
           IF(NTEST.GE.100) WRITE(6,*) ' ND2_TOT = ', ND2_TOT 
*. Number of K2 strings with actual symmetry
           NK2_TOT = LEN_TCCBLK(NK2(1,1),NK2(1,2),NK2(1,3),
     &                   NK2(1,4),K2SM,NSMST)
           IF(NK2_TOT.GT.0.AND.ND2_TOT.NE.0) THEN
*. Number of D2 batches 
           ND2_BAT = ND2_TOT/LD2_BAT
           IF(ND2_BAT*LD2_BAT.LT.ND2_TOT) ND2_BAT=ND2_BAT+1
*. And loop over batches of D2
       IF(NTEST.GE.100) WRITE(6,*) 'ND2_TOT, ND2_BAT ', ND2_TOT, ND2_BAT
           DO ID2_BAT = 1, ND2_BAT
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' ID2SM, ID2_BAT, ID2_BATLEN = ',
     &                     ID2SM, ID2_BAT, ID2_BATLEN
            END IF
*. Generate D3 strings for given sym and batch
            ID2_START = (ID2_BAT-1)*LD2_BAT + 1
            ID2_STOP  = MIN(ND2_TOT,ID2_START+LD2_BAT-1)
            ID2_BATLEN = ID2_STOP-ID2_START+1
*. Generate D2 strings for given sym and batch
            CALL ISMNM_FOR_TCC_BAT(NO2DX,ISM_CAAB_D2,INM_CAAB_D2,ID2SM,
     &           ID2_BATLEN,ID2SM_NEW,
     &           ISM_C1_D2,ISM_CA1_D2,ISM_AA1_D2,
     &           INM_AB1_D2,INM_AA1_D2,INM_CA1_D2,INM_CB1_D2,
     &           ISM_CINI_D2,ISM_CAINI_D2,ISM_AAINI_D2,      
     &           INM_ABINI_D2,INM_AAINI_D2,INM_CAINI_D2,INM_CBINI_D2,0)
            ID2SM_NEW = 0
*
* ==============================================================
* Part 1 O1(d1,d2,d3,d4,ex) T1(i) => OT1(d2,d3,d4,l1) in TSCR4
* ==============================================================
*
* Part 1 : Obtain in TSCR4 OT1(d4,d3,d2,l1) for all l1 of correct sym
*
* Length of OT1(d4,d3,d2,l1)
            LEN_OT1 = NL1_TOT*ID2_BATLEN*ID3_BATLEN*ID4_BATLEN
            IF(LEN_OT1 .GT. N_TDL_MAX ) THEN
*. Problem, too small length for TSCR4 
             WRITE(6,*) ' Dimension of TSCR4 too small '
             WRITE(6,*) ' Required and allocated : ',LEN_OT1, N_TDL_MAX 
             WRITE(6,*) 'NL1_TOT,ID2_BATLEN, ID3_BATLEN, ID4_BATLEN = ',
     &                   NL1_TOT,ID2_BATLEN, ID3_BATLEN, ID4_BATLEN
             WRITE(6,*) ' Ex, D2,D3,D4 and L1 : '
             CALL WRT_SPOX_TP(IOEX,1)
             CALL WRT_SPOX_TP(IO2DX,1)
             CALL WRT_SPOX_TP(IO3DX,1)
             CALL WRT_SPOX_TP(IO4DX,1)
             CALL WRT_SPOX_TP(IOC_L1,1)
             WRITE(6,*) ' T1, T2, T3, T4 '
             CALL WRT_SPOX_TP(IT1,1)
             CALL WRT_SPOX_TP(IT2,1)
             CALL WRT_SPOX_TP(IT3,1)
             CALL WRT_SPOX_TP(IT4,1)
             STOP ' Dimension of TSCR4 too small '
            END IF
*
            ZERO = 0.0D0
            CALL SETVEC(TSCR4,ZERO,LEN_OT1)
*. Loop 1 can be realized in a number of different ways. 
*. The present approacj minimizes the number of times 
*  integrals must be fetched 
*
            DO ID1SM = 1, NSMST
*. Obtain exsm from OP(ex,d4,d3,d2,d1) have symmetry IOPSM
             IF(NTEST.GE.100) WRITE(6,*) ' ID1SM = ', ID1SM
             ID1234SM   = MULTD2H(ID1SM,ID234SM)
             IEXSM = MULTD2H(IOPSM,ID1234SM)
             K1SM = MULTD2H(IT1SM,ID1SM)
*. Number of K1 strings with sym K1SM
             NK1_TOT = LEN_TCCBLK(NK1(1,1),NK1(1,2),NK1(1,3),
     &                 NK1(1,4),K1SM,NSMST)
             IF(NTEST.GE.100) WRITE(6,*) 'K1SM, NK1_TOT ,IEXSM = ', 
     &                                    K1SM, NK1_TOT ,IEXSM
*. Number of excitations strings in O with this symmetry
             NEX_TOT = LEN_TCCBLK(NOEX(1,1),NOEX(1,2),NOEX(1,3),
     &                 NOEX(1,4),IEXSM,NSMST)
*. Number of D1 strings with given symmetry
             ND1_TOT = LEN_TCCBLK(NO1DX(1,1),NO1DX(1,2),NO1DX(1,3),
     &              NO1DX(1,4),ID1SM,NSMST)
             IF(ND1_TOT.GT.0.AND.NEX_TOT.GT.0.AND.NK1_TOT.GT.0) THEN
*. Number of batches of excitation part of O
             NEX_BAT = NEX_TOT/LEX_BAT
             IF(NEX_BAT*LEX_BAT.LT.NEX_TOT) NEX_BAT = NEX_BAT + 1
*. Number of D1 batches
             ND1_BAT = ND1_TOT/LD1_BAT
             IF(ND1_BAT*LD1_BAT.LT.ND1_TOT) ND1_BAT=ND1_BAT+1
*. Loop over D1 batches
             ID1SM_NEW = 1
       IF(NTEST.GE.100) WRITE(6,*) 'ND1_TOT, ND1_BAT ', ND1_TOT, ND1_BAT
             DO ID1_BAT = 1, ND1_BAT
              IF(NTEST.GE.100) WRITE(6,*) ' ID1_BAT = ', ID1_BAT
              ID1_START = (ID1_BAT-1)*LD1_BAT + 1
              ID1_STOP  = MIN(ND1_TOT,ID1_START+LD1_BAT-1)
              ID1_BATLEN = ID1_STOP-ID1_START+1
              IF(NTEST.GE.100) 
     &        WRITE(6,*) ' ID1_STOP, ID1_START, ID1_BATLEN ',
     &                     ID1_STOP, ID1_START, ID1_BATLEN
*. Generate ID1 strings for given sym and batch
              CALL ISMNM_FOR_TCC_BAT(NO1DX,ISM_CAAB_D1,INM_CAAB_D1,
     &             ID1SM,ID1_BATLEN,ID1SM_NEW,
     &             ISM_C1_D1,ISM_CA1_D1,ISM_AA1_D1,
     &             INM_AB1_D1,INM_AA1_D1,INM_CA1_D1,INM_CB1_D1,
     &             ISM_CINI_D1,ISM_CAINI_D1,ISM_AAINI_D1,      
     &             INM_ABINI_D1,INM_AAINI_D1,INM_CAINI_D1,INM_CBINI_D1,
     &             0)
C?            REWIND(8)
C?            WRITE(8,*)
C?   &        ' INM_ABINI_D1,INM_AAINI_D1,INM_CAINI_D1,INM_CBINI_D1',
C?   &          INM_ABINI_D1,INM_AAINI_D1,INM_CAINI_D1,INM_CBINI_D1
C?            REWIND(8)
              ID1SM_NEW = 0
*
* all excitation operators in one shot. 
* As we only have limited ISM, INM, the integrals (d1,d2,d3,d4,ex)
* are obtained in batches 
              IEXSM_NEW = 1
              DO IEX_BAT = 1, NEX_BAT
               IF(NTEST.GE.100) WRITE(6,*) ' IEX_BAT = ', IEX_BAT
               IEX_START = (IEX_BAT-1)*LB + 1
               IEX_STOP  = MIN(NEX_TOT,IEX_START+LB-1)
               IEX_BATLEN = IEX_STOP - IEX_START + 1
*. Generate IEX strings for given sym and batch
               CALL ISMNM_FOR_TCC_BAT(NOEX,ISM_CAAB_EX,INM_CAAB_EX,
     &              IEXSM,IEX_BATLEN,IEXSM_NEW,
     &              ISM_C1_EX,ISM_CA1_EX,ISM_AA1_EX,
     &              INM_AB1_EX,INM_AA1_EX,INM_CA1_EX,INM_CB1_EX,
     &              ISM_CINI_EX,ISM_CAINI_EX,ISM_AAINI_EX,      
     &              INM_ABINI_EX,INM_AAINI_EX,INM_CAINI_EX,INM_CBINI_EX,
     &              0)
               IEXSM_NEW = 0
*. Length of EXD2D3D4
               LEN_EXD234 = IEX_BATLEN*ID2_BATLEN*ID3_BATLEN*ID4_BATLEN
*. Obtain integrals OP(D1,D2,D3,D4,EX)
                 IOFF = 1 +
     &           (IEX_START-1)*ID1_BATLEN*ID2_BATLEN*ID3_BATLEN*
     &                         ID4_BATLEN
C?               WRITE(6,*) ' IOFF = ', IOFF
                 CALL GET_OPINT4(OPSCR(IOFF),
     &           IO4DX,ID4_BATLEN,INM_CAAB_D4(1,1),ISM_CAAB_D4(1,1),
     &           ISTR_D4,IBO4DX,
     &           IO3DX,ID3_BATLEN,INM_CAAB_D3(1,1),ISM_CAAB_D3(1,1),
     &           ISTR_D3,IBO3DX,
     &           IO2DX,ID2_BATLEN,INM_CAAB_D2(1,1),ISM_CAAB_D2(1,1),
     &           ISTR_D2,IBO2DX,
     &           IO1DX,ID1_BATLEN,INM_CAAB_D1(1,1),ISM_CAAB_D1(1,1),
     &           ISTR_D1,IBO1DX,
     &           IOEX,IEX_BATLEN,INM_CAAB_EX(1,1),ISM_CAAB_EX(1,1),
     &           ISTR_EX,IBOEX,IEXD1234_INDX,FACX)
C     GET_OPINT4(OPSCR,
C    &IOD4X,LD4,INM_CAAB_D4,ISM_CAAB_D4,ISTR_D4,IBO4DX,
C    &IOD3X,LD3,INM_CAAB_D3,ISM_CAAB_D3,ISTR_D3,IBO3DX,
C    &IOD2X,LD2,INM_CAAB_D2,ISM_CAAB_D2,ISTR_D2,IBO2DX,
C    &IOD1X,LD1,INM_CAAB_D1,ISM_CAAB_D1,ISTR_D1,IBO1DX,
C    &IOEX ,LEX,INM_CAAB_EX,ISM_CAAB_EX,ISTR_EX,IBOEX,
C    &IEXD1234_INDX,FACX)
*
              END DO
*            ^ End of loop over batches of excitation operators OEX
*. We now have OP(d1,d2,d3,d4,ex) for all ex 
*. OP(d1,d2,d3,d4,ex)*T1(I) => OPT(d2,d3,d4,l1) in TSCR4
*
              SIGNX = SP_D1K1*SIGN_D1K1
              SIGN_LK = SIGN_EXK1
              SIGNG = FLOAT(ISIGNG)
              IF(NTEST.GE.100) WRITE(6,*) ' SIGNX, SIGNG, SIGN_LK = ', 
     /                     SIGNX, SIGNG, SIGN_LK
              NDUM = ID2_BATLEN*ID3_BATLEN*ID4_BATLEN
              IF(NTEST.GE.100) 
     &        WRITE(6,*) ' OP(d1,d2,d3,d4,ex)*T1(I) => OPT(d2,d3,d4,l1)'
              CALL OT_T(OPSCR,T1,TSCR4,NDUM,K1SM,IEXSM,ID1SM,
     &           NK1,NOEX,NT1,NO1DX,NL1,LB,NK1_TOT,NEX_TOT,LEN_T1,
     &           ID1_BATLEN,
     &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,I_NM_CAAB_L,
     &           ISM_CINI_D1,ISM_CAINI_D1,ISM_AAINI_D1,      
     &           INM_ABINI_D1,INM_AAINI_D1,INM_CBINI_D1,INM_CAINI_D1,
     &           SIGNX,SIGNG,SIGN_LK, 
     &           IB_D1K1,IB_EXK1,
     &           IX1_CA,SX1_CA,IX1_CB,SX1_CB,
     &           IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &           IX2_CA,SX2_CA,IX2_CB,SX2_CB,
     &           IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &           IBT1_TCC, IBL1_TCC,
     &           TSCR1,TSCR2,TSCR3)
C     OT_T(OT1,T2,OT1T2,NDUM,
C    &           ISM_K,ISM_L,ISM_D,
C    &           NK,NL,NI,ND,NLK,LEN_KLBAT,NK_TOT,NL_TOT,NI_TOT,ND_TOT,
C    &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
C    &           ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,
C    &           INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
C    &           SIGNX,SIGNG,SIGNLK,
C    &           IB_DK,IB_LK,
C    &           IX_DK_CA,SX_DK_CA,IX_DK_CB,SX_DK_CB,
C    &           IX_DK_AA,SX_DK_AA,IX_DK_AB,SX_DK_AB,
C    &           IX_LK_CA,SX_LK_CA,IX_LK_CB,SX_LK_CB,
C    &           IX_LK_AA,SX_LK_AA,IX_LK_AB,SX_LK_AB,
C    &           IBT2_TCC,IBI_TCC,TSCR1,TSCR2,TSCR3)


*
             END DO
*             ^ End of loop over D1 batch
              IF(NTEST.GE.100) THEN 
                WRITE(6,*) ' OT1(D4D3D2,L1) matrix after ID1SM=  ',ID1SM
                ND432 = ID2_BATLEN*ID3_BATLEN*ID4_BATLEN 
                CALL WRTMAT(TSCR4,ND432,NL1_TOT,ND432,NL1_TOT)
              END IF
*
            END IF
*           ^ End if there were nonvanishing dimensions
            END DO
*           ^ End of loop over symmetries of D1
*
            IF(NTEST.GE.100) THEN 
              WRITE(6,*) ' OT1(D4D3D2,L1) matrix '
              ND432 = ID2_BATLEN*ID3_BATLEN*ID4_BATLEN 
              CALL WRTMAT(TSCR4,ND432,NL1_TOT,ND432,NL1_TOT)
            END IF
*
* =========================================================
* Part 2 OT1(d2,d3,d4,l1) T2(i) => OT12(d3,d4,l2) in TSCR5
* =========================================================
*

*
*  We have OT1(d4,d3,d2,L1) for d4,d3,d2 in batch and all L1 of given symmetry
*  Construct OT12(d4,d3,l2) for d4,d3 in batch and all l2
*
            LEND34 = ID3_BATLEN*ID4_BATLEN
            SIGNX = SIGN_D2K2*SP_D2K2
            SIGNY = 1.0D0
            SIGNLK = SIGN_L1K2
CM          CALL SETVEC(TSCR5,ZERO,LEND34*NL2_TOT)
            IF(NTEST.GE.100) 
     &      WRITE(6,*) '  OT1(d2,d3,d4,l1) T2(i) => OT12(d3,d4,l2) '
            CALL OT_T(TSCR4,T2,TSCR5,LEND34,
     &       K2SM,L1SM,ID2SM,
     &       NK2,NL1,NT2,NO2DX,NL2,LB,NK2_TOT,NL1_TOT,LEN_T2,
     &       ID2_BATLEN,
     &       ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
     &       ISM_CINI_D2,ISM_CAINI_D2,ISM_AAINI_D2,
     &       INM_ABINI_D2,INM_AAINI_D2,INM_CBINI_D2,INM_CAINI_D2,
     &       SIGNX,SIGNY,SIGNLK,
     &       IB_D2K2,IB_L1K2,
     &       IX3_CA,SX3_CA,IX3_CB,SX3_CB,
     &       IX3_AA,SX3_AA,IX3_AB,SX3_AB,
     &       IX4_CA,SX4_CA,IX4_CB,SX4_CB,
     &       IX4_AA,SX4_AA,IX4_AB,SX4_AB,
     &       IBT2_TCC,IBL2_TCC,TSCR1,TSCR2,TSCR3)  
C     OT_T(OT1,T2,OT1T2,NDUM,
C    &           ISM_K,ISM_L,ISM_D,
C    &           NK,NL,NI,ND,NLK,LEN_KLBAT,NK_TOT,NL_TOT,NI_TOT,ND_TOT,
C    &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
C    &           ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,
C    &           INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
C    &           SIGNX,SIGNG,SIGNLK,
C    &           IB_DK,IB_LK,
C    &           IX_DK_CA,SX_DK_CA,IX_DK_CB,SX_DK_CB,
C    &           IX_DK_AA,SX_DK_AA,IX_DK_AB,SX_DK_AB,
C    &           IX_LK_CA,SX_LK_CA,IX_LK_CB,SX_LK_CB,
C    &           IX_LK_AA,SX_LK_AA,IX_LK_AB,SX_LK_AB,
C    &           IBT2_TCC,IBI_TCC,TSCR1,TSCR2,TSCR3)
           END DO
*          ^ End of loop over D2 batches
          END IF 
*         ^ End of NK2_TOT*ND2_TOT .ne. 0
          END DO
*         ^ End of loop over symmetries of D2 batches
          IF(NTEST.GE.100) THEN 
            WRITE(6,*) ' OT12(D4D3,L2) matrix '
            ND43 = ID3_BATLEN*ID4_BATLEN 
            CALL WRTMAT(TSCR5,ND43,NL2_TOT,ND43,NL2_TOT)
          END IF
* We now have OPT12(d3,d4,l2) in TSCR5
*
* Part 3 : OPT123(d4,l3) in TSCR6
*
          SIGNX = SIGN_D3K3*SP_D3K3
          SIGNY = 1.0D0
          SIGNLK = SIGN_L2K3
C         CALL SETVEC(TSCR6,ZERO,ID4_BATLEN*NL3_TOT)
          IF(NTEST.GE.100) 
     &    WRITE(6,*) '  OT12(d3,d4,l2) T3(i) => OT123(d4,l3) '
          CALL OT_T(TSCR5,T3,TSCR6,ID4_BATLEN,
     &       K3SM,L2SM,ID3SM,
     &       NK3,NL2,NT3,NO3DX,NL3,LB,NK3_TOT,NL2_TOT,LEN_T3,
     &       ID3_BATLEN,
     &       ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
     &       ISM_CINI_D3,ISM_CAINI_D3,ISM_AAINI_D3,
     &       INM_ABINI_D3,INM_AAINI_D3,INM_CBINI_D3,INM_CAINI_D3,
     &       SIGNX,SIGNY,SIGNLK,
     &       IB_D3K3,IB_L2K3,
     &       IX5_CA,SX5_CA,IX5_CB,SX5_CB,
     &       IX5_AA,SX5_AA,IX5_AB,SX5_AB,
     &       IX6_CA,SX6_CA,IX6_CB,SX6_CB,
     &       IX6_AA,SX6_AA,IX6_AB,SX6_AB,
     &       IBT3_TCC,IBL3_TCC,TSCR1,TSCR2,TSCR3)  
C     OT_T(OT1,T2,OT1T2,NDUM,
C    &           ISM_K,ISM_L,ISM_D,
C    &           NK,NL,NI,ND,NLK,LEN_KLBAT,NK_TOT,NL_TOT,NI_TOT,ND_TOT,
C    &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
C    &           ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,
C    &           INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
C    &           SIGNX,SIGNG,SIGNLK,
C    &           IB_DK,IB_LK,
C    &           IX_DK_CA,SX_DK_CA,IX_DK_CB,SX_DK_CB,
C    &           IX_DK_AA,SX_DK_AA,IX_DK_AB,SX_DK_AB,
C    &           IX_LK_CA,SX_LK_CA,IX_LK_CB,SX_LK_CB,
C    &           IX_LK_AA,SX_LK_AA,IX_LK_AB,SX_LK_AB,
C    &           IBT2_TCC,IBI_TCC,TSCR1,TSCR2,TSCR3)
         END DO
*        ^ End of loop over D3 batches 
        END IF
*.      ^ End if there are K3 and L3 strings
        END DO
*       ^ End of loop over D3 symmetry
        IF(NTEST.GE.100) THEN 
           WRITE(6,*) ' OT123(D4,L3) matrix '
           ND4 = ID4_BATLEN 
           CALL WRTMAT(TSCR6,ND4,NL3_TOT,ND4,NL3_TOT)
         END IF
*
* Part 4 : OT123(d4,l3)*T4 => OT1234(I)
        SIGNX = SIGN_D4K4*SP_D4K4
        SIGNY = 1.0D0
        SIGNLK = SIGN_L3K4
        IF(NTEST.GE.100) 
     &  WRITE(6,*) '  OT123(d4,l1) T4(i) => OT1234(l4) '
        CALL OT_T(TSCR6,T4,OT1234,1,
     &       K4SM,L3SM,ID4SM,
     &       NK4,NL3,NT4,NO4DX,NL4,LB,NK4_TOT,NL3_TOT,LEN_T4,
     &       ID4_BATLEN,
     &       ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
     &       ISM_CINI_D4,ISM_CAINI_D4,ISM_AAINI_D4,
     &       INM_ABINI_D4,INM_AAINI_D4,INM_CBINI_D4,INM_CAINI_D4,
     &       SIGNX,SIGNY,SIGNLK,
     &       IB_D4K4,IB_L3K4,
     &       IX7_CA,SX7_CA,IX7_CB,SX7_CB,
     &       IX7_AA,SX7_AA,IX7_AB,SX7_AB,
     &       IX8_CA,SX8_CA,IX8_CB,SX8_CB,
     &       IX8_AA,SX8_AA,IX8_AB,SX8_AB,
     &       IBT4_TCC,IBL4_TCC,TSCR1,TSCR2,TSCR3)  
C     OT_T(OT1,T2,OT1T2,NDUM,
C    &           ISM_K,ISM_L,ISM_D,
C    &           NK,NL,NI,ND,NLK,LEN_KLBAT,NK_TOT,NL_TOT,NI_TOT,ND_TOT,
C    &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
C    &           ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,
C    &           INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
C    &           SIGNX,SIGNG,SIGNLK,
C    &           IB_DK,IB_LK,
C    &           IX_DK_CA,SX_DK_CA,IX_DK_CB,SX_DK_CB,
C    &           IX_DK_AA,SX_DK_AA,IX_DK_AB,SX_DK_AB,
C    &           IX_LK_CA,SX_LK_CA,IX_LK_CB,SX_LK_CB,
C    &           IX_LK_AA,SX_LK_AA,IX_LK_AB,SX_LK_AB,
C    &           IBT2_TCC,IBI_TCC,TSCR1,TSCR2,TSCR3)
       END DO
*      ^ End of D4 batches 
      END DO
*     ^ End of D4 symmetries
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated OT1234 block '
        CALL WRT_TCC_BLK(OT1234,IL4SM,NL4(1,1),NL4(1,2),NL4(1,3),
     &       NL4(1,4),NSMST)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'OPCTT ')
      RETURN
      END
      SUBROUTINE GET_OPINT4(OPSCR,
     &IOD4X,LD4,INM_CAAB_D4,ISM_CAAB_D4,ISTR_D4,IBO4DX,
     &IOD3X,LD3,INM_CAAB_D3,ISM_CAAB_D3,ISTR_D3,IBO3DX,
     &IOD2X,LD2,INM_CAAB_D2,ISM_CAAB_D2,ISTR_D2,IBO2DX,
     &IOD1X,LD1,INM_CAAB_D1,ISM_CAAB_D1,ISTR_D1,IBO1DX,
     &IOEX ,LEX,INM_CAAB_EX,ISM_CAAB_EX,ISTR_EX,IBOEX,
     &IEXD1234_INDX,FACX)
*
* Fetch batch of operator integrals ordered as O(D1,D2,D3,D4,EX) 
*
* Connected with standard integral input
*
* Jeppe Olsen, May of 2000
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'glbbas.inc'
*. Input 
*.. Occupation of gas spaces         
      INTEGER IOD4X(NGAS,4),IOD3X(NGAS,4)
      INTEGER IOD2X(NGAS,4),IOD1X(NGAS,4),IOEX(NGAS,4)
*. Index in EXD1D2D3D4 to original index in H
      INTEGER IEXD1234_INDX(4)
*.. The sym and number of the included strings 
      INTEGER ISM_CAAB_EX(4,*),INM_CAAB_EX(4,*)
      INTEGER ISM_CAAB_D1(4,*),INM_CAAB_D1(4,*)
      INTEGER ISM_CAAB_D2(4,*),INM_CAAB_D2(4,*)
      INTEGER ISM_CAAB_D3(4,*),INM_CAAB_D3(4,*)
      INTEGER ISM_CAAB_D4(4,*),INM_CAAB_D4(4,*)
*.. The actual occupation of the various strings
      INTEGER ISTR_D1(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D2(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D3(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D4(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_EX(MX_ST_TSOSO_BLK_MX*NSMST,4)
*.. Start of strings with given symmetry
      INTEGER IBO1DX(8,4),IBO2DX(8,4),IBO3DX(8,4), IBO4DX(8,4)
      INTEGER IBOEX(8,4)
* 
*. Local scratch
*
C     INTEGER IOCC(4),ICREA(4),IANNI(4)
      INTEGER JD1STR(4),JD2STR(4),JD3STR(4),JD4STR(4),JEXSTR(4)
      INTEGER IJKL_EDD(4), IJKL_ORIG(4)
      INTEGER NOP_D1_CAAB(4),NOP_D2_CAAB(4),NOP_EX_CAAB(4)
      INTEGER NOP_D3_CAAB(4),NOP_D4_CAAB(4)
*. Output
      DIMENSION OPSCR(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GET_OPINT speaking : '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' IOEX, IOD4X, IOD3X, IOD2X, IOD1X  : '
        CALL WRT_SPOX_TP(IOEX,1)
        CALL WRT_SPOX_TP(IOD4X,1)
        CALL WRT_SPOX_TP(IOD3X,1)
        CALL WRT_SPOX_TP(IOD2X,1)
        CALL WRT_SPOX_TP(IOD1X,1)
        WRITE(6,*)
      END IF
*
*. Number of C/A A/B indeces 
*. D1
      NOP_D1_CA = IELSUM(IOD1X(1,1),NGAS)
      NOP_D1_CB = IELSUM(IOD1X(1,2),NGAS)
      NOP_D1_AA = IELSUM(IOD1X(1,3),NGAS)
      NOP_D1_AB = IELSUM(IOD1X(1,4),NGAS)
      NOP_D1 = NOP_D1_CA+NOP_D1_CB+NOP_D1_AA+NOP_D1_AB
      NOP_D1_CAAB(1) = NOP_D1_CA
      NOP_D1_CAAB(2) = NOP_D1_CB
      NOP_D1_CAAB(3) = NOP_D1_AA
      NOP_D1_CAAB(4) = NOP_D1_AB
*. D2
      NOP_D2_CA = IELSUM(IOD2X(1,1),NGAS)
      NOP_D2_CB = IELSUM(IOD2X(1,2),NGAS)
      NOP_D2_AA = IELSUM(IOD2X(1,3),NGAS)
      NOP_D2_AB = IELSUM(IOD2X(1,4),NGAS)
      NOP_D2 = NOP_D2_CA+NOP_D2_CB+NOP_D2_AA+NOP_D2_AB
      NOP_D2_CAAB(1) = NOP_D2_CA
      NOP_D2_CAAB(2) = NOP_D2_CB
      NOP_D2_CAAB(3) = NOP_D2_AA
      NOP_D2_CAAB(4) = NOP_D2_AB
*. D3
      NOP_D3_CA = IELSUM(IOD3X(1,1),NGAS)
      NOP_D3_CB = IELSUM(IOD3X(1,2),NGAS)
      NOP_D3_AA = IELSUM(IOD3X(1,3),NGAS)
      NOP_D3_AB = IELSUM(IOD3X(1,4),NGAS)
      NOP_D3 = NOP_D3_CA+NOP_D3_CB+NOP_D3_AA+NOP_D3_AB
      NOP_D3_CAAB(1) = NOP_D3_CA
      NOP_D3_CAAB(2) = NOP_D3_CB
      NOP_D3_CAAB(3) = NOP_D3_AA
      NOP_D3_CAAB(4) = NOP_D3_AB
*. D4
      NOP_D4_CA = IELSUM(IOD4X(1,1),NGAS)
      NOP_D4_CB = IELSUM(IOD4X(1,2),NGAS)
      NOP_D4_AA = IELSUM(IOD4X(1,3),NGAS)
      NOP_D4_AB = IELSUM(IOD4X(1,4),NGAS)
      NOP_D4 = NOP_D4_CA+NOP_D4_CB+NOP_D4_AA+NOP_D4_AB
*
      NOP_D4_CAAB(1) = NOP_D4_CA
      NOP_D4_CAAB(2) = NOP_D4_CB
      NOP_D4_CAAB(3) = NOP_D4_AA
      NOP_D4_CAAB(4) = NOP_D4_AB
*
      NOP_EX_CA = IELSUM(IOEX(1,1),NGAS)
      NOP_EX_CB = IELSUM(IOEX(1,2),NGAS)
      NOP_EX_AA = IELSUM(IOEX(1,3),NGAS)
      NOP_EX_AB = IELSUM(IOEX(1,4),NGAS)
      NOP_EX = NOP_EX_CA+NOP_EX_CB+NOP_EX_AA+NOP_EX_AB
      NOP_EX_CAAB(1) = NOP_EX_CA
      NOP_EX_CAAB(2) = NOP_EX_CB
      NOP_EX_CAAB(3) = NOP_EX_AA
      NOP_EX_CAAB(4) = NOP_EX_AB
*
      NCREA_ALPHA  = NOP_D1_CA + NOP_D2_CA + NOP_EX_CA
     &             + NOP_D3_CA + NOP_D4_CA
      NCREA_BETA   = NOP_D1_CB + NOP_D2_CB + NOP_EX_CB
     &             + NOP_D3_CB + NOP_D4_CB
      NANNI_ALPHA  = NOP_D1_AA + NOP_D2_AA + NOP_EX_AA
     &             + NOP_D3_AA + NOP_D4_AA
      NANNI_BETA   = NOP_D1_AB + NOP_D2_AB + NOP_EX_AB
     &             + NOP_D3_AB + NOP_D4_AB
*
C?    WRITE(6,*) ' NANNI_ALPHA, NANNI_BETA = ', 
C?   &             NANNI_ALPHA, NANNI_BETA
*
      NCREA = NCREA_ALPHA + NCREA_BETA
      NANNI = NANNI_ALPHA + NANNI_BETA
*
      NALPHA = NANNI_ALPHA + NCREA_ALPHA
      NBETA  = NANNI_BETA  + NCREA_BETA
*
      NOP = NCREA + NANNI
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original index '
        CALL IWRTMA(IEXD1234_INDX,1,NOP,1,NOP)
      END IF
*
      IF(.NOT.((NCREA.EQ.1.AND.NANNI.EQ.1 ) .OR.
     &         (NCREA.EQ.2.AND.NANNI.EQ.2 )     )) THEN
        WRITE(6,*) ' Unknown operator in GET_OPINT '
        WRITE(6,*) ' NCREA, NANNI = ', NCREA, NANNI
        STOP       ' Unknown operator in GET_OPINT '
      END IF
*.
      IF(.NOT.(NCREA_ALPHA.EQ.NANNI_ALPHA.AND.
     &         NCREA_BETA.EQ.NANNI_BETA       ) ) THEN
        WRITE(6,*) ' Unknown operator in GET_OPINT '
        WRITE(6,*) 
     &  ' NCREA_ALPHA, NANNI_ALPHA, NCREA_BETA, NANNI_BETA  ', 
     &    NCREA_ALPHA, NANNI_ALPHA, NCREA_BETA, NANNI_BETA 
        STOP       ' Unknown operator in GET_OPINT '
      END IF
*. Now we know that a usual MS = 0, one- or two-electron 
*. operator is in operation
*. More complicated operators can be concocted
*
      INT = 0
      DO JEX = 1, LEX
C?    WRITE(6,*) ' EX string : '
      CALL GET_CCSTR_FROM_LIST 
     &     (JEXSTR,INM_CAAB_EX(1,JEX),ISM_CAAB_EX(1,JEX),
     &      NOP_EX_CAAB,ISTR_EX(1,1),ISTR_EX(1,2),ISTR_EX(1,3),
     &      ISTR_EX(1,4),IBOEX(1,1),IBOEX(1,2),IBOEX(1,3),
     &      IBOEX(1,4) )
      DO JD4 = 1, LD4
      CALL GET_CCSTR_FROM_LIST 
     &     (JD4STR,INM_CAAB_D4(1,JD4),ISM_CAAB_D4(1,JD4),
     &      NOP_D4_CAAB,ISTR_D4(1,1),ISTR_D4(1,2),ISTR_D4(1,3),
     &      ISTR_D4(1,4),IBO4DX(1,1),IBO4DX(1,2),IBO4DX(1,3),
     &      IBO4DX(1,4) )
      DO JD3 = 1, LD3
      CALL GET_CCSTR_FROM_LIST 
     &     (JD3STR,INM_CAAB_D3(1,JD3),ISM_CAAB_D3(1,JD3),
     &      NOP_D3_CAAB,ISTR_D3(1,1),ISTR_D3(1,2),ISTR_D3(1,3),
     &      ISTR_D3(1,4),IBO3DX(1,1),IBO3DX(1,2),IBO3DX(1,3),
     &      IBO3DX(1,4) )
      DO JD2 = 1, LD2
C?    WRITE(6,*) ' D2 string : '
      CALL GET_CCSTR_FROM_LIST 
     &     (JD2STR,INM_CAAB_D2(1,JD2),ISM_CAAB_D2(1,JD2),
     &      NOP_D2_CAAB,ISTR_D2(1,1),ISTR_D2(1,2),ISTR_D2(1,3),
     &      ISTR_D2(1,4),IBO2DX(1,1),IBO2DX(1,2),IBO2DX(1,3),
     &      IBO2DX(1,4) )
      DO JD1 = 1, LD1
*. Occupation of J1 string in JD1STR
C?    WRITE(6,*) ' D1 string : '
      CALL GET_CCSTR_FROM_LIST 
     &     (JD1STR,INM_CAAB_D1(1,JD1),ISM_CAAB_D1(1,JD1),
     &      NOP_D1_CAAB,ISTR_D1(1,1),ISTR_D1(1,2),ISTR_D1(1,3),
     &      ISTR_D1(1,4),IBO1DX(1,1),IBO1DX(1,2),IBO1DX(1,3),
     &      IBO1DX(1,4) )
*. Indeces of operator EXD1D2 
      DO IOP_EX = 1, NOP_EX
        IJKL_EDD(IOP_EX) = JEXSTR(IOP_EX)
      END DO
      DO IOP_D1 = 1, NOP_D1
        IJKL_EDD(NOP_EX+IOP_D1) = JD1STR(IOP_D1)
      END DO
      DO IOP_D2 = 1, NOP_D2
        IJKL_EDD(NOP_EX+NOP_D1+IOP_D2) = JD2STR(IOP_D2)
      END DO
      DO IOP_D3 = 1, NOP_D3
        IOP_EFF = NOP_EX+NOP_D1+NOP_D2 + IOP_D3
        IJKL_EDD(IOP_EFF) = JD3STR(IOP_D3)
      END DO
      DO IOP_D4 = 1, NOP_D4
        IOP_EFF = NOP_EX+NOP_D1+NOP_D2 + NOP_D3 + IOP_D4
        IJKL_EDD(IOP_EFF) = JD4STR(IOP_D4)
      END DO
*. Original order
      DO IOP = 1, NOP   
        IJKL_ORIG(IEXD1234_INDX(IOP)) = IJKL_EDD(IOP)
      END DO
*
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' IJKL_EDD and IJKL_ORIG '
        CALL IWRTMA(IJKL_EDD,1,NOP,1,NOP)
        CALL IWRTMA(IJKL_ORIG,1,NOP,1,NOP)    
      END IF
*
      INT = INT + 1
      IF(NOP.EQ.2) THEN
*. Normal one-electron integral 
        IF(IREFTYP.NE.2) THEN
*. Similarity transformed integrals in orbital basis 
          OPSCR(INT) =
     &    GETH1_B(IJKL_ORIG(1),IJKL_ORIG(2))*FACX
        ELSE
*. Similarity transformed integrals in spinorbital basis 
          IF(NALPHA.EQ.2) THEN
           OPSCR(INT) =
     &     GETH1_B2(IJKL_ORIG(1),IJKL_ORIG(2),WORK(KFI_AL))*FACX
          ELSE IF (NBETA.EQ.2) THEN
           OPSCR(INT) =
     &     GETH1_B2(IJKL_ORIG(1),IJKL_ORIG(2),WORK(KFI_BE))*FACX
          END IF
        END IF
      ELSE IF(NOP.EQ.4.AND.NALPHA.EQ.2.AND.NBETA.EQ.2) THEN
*. Normal two-electron integral 
        IF(IREFTYP.NE.2) THEN
        OPSCR(INT) = FACX*GTIJKL(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                      IJKL_ORIG(3))
        ELSE 
            OPSCR(INT)
     &    = GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                   IJKL_ORIG(3),2,2)*FACX
        END IF
      ELSE IF (NOP.EQ.4.AND.NALPHA.EQ.4) THEN
*. coulomb - exchange  
*. There can be a factor: when several indeces belongs to the same 
*. class, the same integral is included several times. Counter act this
        IF(IREFTYP.NE.2) THEN
          OPSCR(INT) = GTIJKL(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                        IJKL_ORIG(3)) 
     &               - GTIJKL(IJKL_ORIG(1),IJKL_ORIG(3),IJKL_ORIG(2),
     &                      IJKL_ORIG(4))
          OPSCR(INT) = FACX*OPSCR(INT)
        ELSE
          OPSCR(INT)
     &    = GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                   IJKL_ORIG(3),4,0)
     &    - GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(3),IJKL_ORIG(2),
     &                   IJKL_ORIG(4),4,0)
          OPSCR(INT) = FACX*OPSCR(INT)
        END IF
      ELSE IF (NOP.EQ.4.AND.NBETA.EQ.4) THEN
        IF(IREFTYP.NE.2) THEN
          OPSCR(INT) = GTIJKL(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                        IJKL_ORIG(3)) 
     &               - GTIJKL(IJKL_ORIG(1),IJKL_ORIG(3),IJKL_ORIG(2),
     &                      IJKL_ORIG(4))
          OPSCR(INT) = FACX*OPSCR(INT)
        ELSE
          OPSCR(INT)
     &    = GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                   IJKL_ORIG(3),0,4)
     &    - GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(3),IJKL_ORIG(2),
     &                   IJKL_ORIG(4),0,4)
          OPSCR(INT) = FACX*OPSCR(INT)
        END IF
      END IF
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' INT, OPSCR(INT) = ', INT, OPSCR(INT)
      END IF
*.    ^ End of switch between different integral types
      END DO
      END DO
      END DO
      END DO
      END DO
*
      IF(NTEST.GE.100) THEN 
         WRITE(6,*)
         WRITE(6,*) ' Output matrix from GET_OPINT as X(D4D3D2EX,D1)' 
         WRITE(6,*) ' ============================================= '
         WRITE(6,*)
         CALL WRTMAT(OPSCR,LD4*LD3*LD2*LEX,LD1,LD4*LD3*LD2*LEX,LD1)
      END IF
*
      RETURN
      END 
      SUBROUTINE OT_T(OT1,T2,OT1T2,NDUM,
     &           ISM_K,ISM_L,ISM_D,
     &           NK,NL,NI,ND,NLK,LEN_KLBAT,NK_TOT,NL_TOT,NI_TOT,ND_TOT,
     &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
     &           ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,
     &           INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
     &           SIGNX,SIGNG,SIGNLK,
     &           IB_DK,IB_LK,
     &           IX_DK_CA,SX_DK_CA,IX_DK_CB,SX_DK_CB,
     &           IX_DK_AA,SX_DK_AA,IX_DK_AB,SX_DK_AB,
     &           IX_LK_CA,SX_LK_CA,IX_LK_CB,SX_LK_CB,
     &           IX_LK_AA,SX_LK_AA,IX_LK_AB,SX_LK_AB,
     &           IBT2_TCC,IBI_TCC,TSCR1,TSCR2,TSCR3)
* Input / Argument list 
* =====================
*
* OT1 : OT1 coefficients
*  T2 : T2  coefficients
* OT1T2 (OUTPUT) : OT1T2 coefficients
* NDUM : Number of dummy indeces 
* ISM_K : Symmetry of K in T2(D,K)
* ISM_L : Symmetry of L-operator in OT1(IDUM,D,L)
* ISM_D : Symmetry of  D in T2(D,K)
* NK : Number of CA,CB,AA,AB operators per sym for K-strings in T2(D,K)
* NL : Number of CA,CB,AA,AB operators per sym for L-strings in  OT1(IDUM,D,L)
* NI : Number of CA,CB,AA,AB operators per sym for I-strings in  T2(I)
* ND : Number of CA,CB,AA,AB operators per sym for D-strings in OT1(IDUM,D,L)
* NLK : Number of CA,CB,AA,AB operators per sym for LK-strings in OT1T2(IDUM,LK)
* LEN_KLBAT : Length of batches over L and K
* NK_TOT : Total number of K-strings in T2(D,K)
* NI_TOT : Total number of I-strings in T2(I)
* ND_TOT : Total number of D-strings in T2(D,K) ( number of D's in batch)
* ISM_CAAB_K, INM_CAAB_K, ISM_CAAB_L, INM_CAAB_L (SCRATCH) : scratch 
* arrays for CAAB of K and L
* ISM_CINI_L,ISM_CAINI_L,ISM_AAINI_L,INM_ABINI_L,INM_AAINI_L,INM_CBINI_L,
* INM_CAINI_L : Initial value
* SIGNX, SIGNG, SIGNLK : General sign factors 
* IB_DK(IDSM,IKSM,ICAAB) Offset to mappings for D and K ops with given sym 
*                        for ICAAB part of operator
* IB_LK(ILSM,IKSM,ICAAB) Offset to mappings for L and K ops with given sym 
*                        for ICAAB part of operator
* IX_DK_CA, SX_DK_CA, : (D,K) = > I mapping for CA
* IX_DK_CB, SX_DK_CB, : (D,K) = > I mapping for CB
* IX_DK_AA, SX_DK_AA, : (D,K) = > I mapping for AA
* IX_DK_AB, SX_DK_AB, : (D,K) = > I mapping for AB
* IX_LK_CA, SX_LK_CA  : (L,K) => LK mapping for CA
* IX_LK_CB, SX_LK_CB  : (L,K) => LK mapping for CB
* IX_LK_AA, SX_LK_AA  : (L,K) => LK mapping for AA
* IX_LK_AB, SX_LK_AB  : (L,K) => LK mapping for AB
* IBT2_TCC : Offsets to symmetry blocks in T2
* IBI_TCC : Offset to symmetry blocks in T1
* TSCR1, TSCR2,TSCR3 (SCRATCH)
* 

*
* Include next operator in collected operator 
*
* OT1(D,IDUM,L)T2(I) => OT12(IDUM,LK) for all LK,I,L, D in batch
*
* Loop over batches of K
*  T2(I) => T2(D,K)
*  Loop over batches of L
*   OT1T2(IDUM,L,K) = sum(D) OT1(D,IDUM,L)*T2(D,K)  
*   OT1T2(IDUM,L,K) => OT1T2(IDUM,LK)
*  End of loop over L batches 
* End of loop over K batches
*
* Jeppe Olsen, July 14, 2001, completed and debugged July 2003
*                             (Well it became Nov. 2004)
*
      INCLUDE 'implicit.inc'
*. Number of strings per symmetry
      INTEGER NK(8,4),NL(8,4),NI(8,4),ND(8,4),NLK(8,4)
*. Maps 
      INTEGER IB_DK(8,8,4)
      DIMENSION IX_DK_CA(*),SX_DK_CA(*),IX_DK_CB(*),SX_DK_CB(*)
      DIMENSION IX_DK_AA(*),SX_DK_AA(*),IX_DK_AB(*),SX_DK_AB(*)
      INTEGER IB_LK(8,8,4)
      DIMENSION IX_LK_CA(*),SX_LK_CA(*),IX_LK_CB(*),SX_LK_CB(*)
      DIMENSION IX_LK_AA(*),SX_LK_AA(*),IX_LK_AB(*),SX_LK_AB(*)
*. Offset 
      INTEGER IBT2_TCC(8,8,8),IBI_TCC(8,8,8)
*. Input
      DIMENSION OT1(*), T2(*)
*. Output 
      DIMENSION OT1T2(*)
*. Scratch through input 
      INTEGER ISM_CAAB_K(*),INM_CAAB_K(*)
      INTEGER ISM_CAAB_L(*),INM_CAAB_L(*)
      DIMENSION TSCR1(*),TSCR2(*),TSCR3(*)
*
      INCLUDE 'csm.inc'
*
*
*
      NTEST = 00
      ZERO = 0.0D0
      ONE  = 1.0D0
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' OT_T entered '
        WRITE(6,*) ' NK_TOT, NL_TOT, ND_TOT, NI_TOT = ',
     &               NK_TOT, NL_TOT, ND_TOT, NI_TOT
      END IF
*. Number of K-batches
C?    WRITE(6,*) 'NK_TOT, LEN_KLBAT = ', NK_TOT, LEN_KLBAT
      NK_BAT = NK_TOT/LEN_KLBAT
      IF(NK_BAT*LEN_KLBAT.LT.NK_TOT) NK_BAT = NK_BAT + 1
*. Number of batches of L
      NL_BAT = NL_TOT/LEN_KLBAT
      IF(NL_BAT*LEN_KLBAT.LT.NL_TOT) NL_BAT = NL_BAT + 1
*. And loop over batches of K 
      IKSM_NEW = 1
      DO IK_BAT = 1, NK_BAT
        IK_START = (IK_BAT-1)*LEN_KLBAT + 1
        IK_STOP  = MIN(NK_TOT,IK_START+LEN_KLBAT-1)
        IK_BATLEN = IK_STOP-IK_START+1 
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' IK_BAT, IK_START, IK_STOP = ',
     &                 IK_BAT, IK_START, IK_STOP 
        END IF
*. Generate K strings for given sym and batch
        CALL ISMNM_FOR_TCC_BAT(NK,ISM_CAAB_K,INM_CAAB_K,ISM_K,
     &       IK_BATLEN,IKSM_NEW,
     &       ISM_C1_K,ISM_CA1_K,ISM_AA1_K,
     &       INM_AB1_K,INM_AA1_K,INM_CA1_K,INM_CB1_K,     
     &       ISM_CINI_K,ISM_CAINI_K,ISM_AAINI_K,      
     &       INM_ABINI_K,INM_AAINI_K,INM_CAINI_K,INM_CBINI_K,1)
        IKSM_NEW = 0
*. Obtain T2(D,K) in TSCR1
        IF(NTEST.GE.100) WRITE(6,*) '  T2(I) => T2(D,K) reordering'
C       SIGNX = SIGN_D2K2*SP_D2K2
        CALL TI_TO_TOKBN(NSMST,
     &       ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,ISM_D,
     &       INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
     &       ND_TOT,ND,
     &       ISM_CINI_K,ISM_CAINI_K,ISM_AAINI_K,ISM_K,
     &       INM_ABINI_K,INM_AAINI_K,INM_CBINI_K,INM_CAINI_K,
     &       IK_BATLEN,NK,
     &       IX_DK_CA,SX_DK_CA,IB_DK(1,1,1),
     &       IX_DK_CB,SX_DK_CB,IB_DK(1,1,2),
     &       IX_DK_AA,SX_DK_AA,IB_DK(1,1,3),
     &       IX_DK_AB,SX_DK_AB,IB_DK(1,1,4),
     &       NI(1,1),NI(1,2),NI(1,3),NI(1,4),
     &       TSCR1,T2,IBT2_TCC,1,1,SIGNX,1)
*. And loop over batches of L
        ILSM_NEW = 1
        DO IL_BAT = 1, NL_BAT
          IL_START = (IL_BAT-1)*LEN_KLBAT + 1
          IL_STOP  = MIN(NL_TOT,IL_START+LEN_KLBAT-1)
          IF(NTEST.GE.100) WRITE(6,*) ' IL_BAT, IL_START, IL_STOP ',
     &                 IL_BAT, IL_START, IL_STOP
          IL_BATLEN = IL_STOP-IL_START+1 
*. Generate L strings for given sym and batch
          CALL ISMNM_FOR_TCC_BAT(NL,ISM_CAAB_L,INM_CAAB_L,ISM_L,
     &         IL_BATLEN,ILSM_NEW,
     &         ISM_C1_L,ISM_CA1_L,ISM_AA1_L,
     &         INM_AB1_L,INM_AA1_L,INM_CA1_L,INM_CB1_L,     
     &         ISM_CINI_L,ISM_CAINI_L,ISM_AAINI_L,      
     &         INM_ABINI_L,INM_AAINI_L,INM_CAINI_L,INM_CBINI_L,1)
          ILSM_NEW = 0
*. OT12(IDUM,L,K) = sum(D) OT1(D,IDUM,L)*T2(D,K) in TSCR3
*. Signg is introduced here !!
          NR_OT12 = NDUM*IL_BATLEN
          NC_OT12 = IK_BATLEN
          NR_OT1 = ND_TOT 
          NC_OT1 = NR_OT12
          NR_T2  = ND_TOT
          NC_T2 = IK_BATLEN
CE        IOFF = (IL_START - 1)*ND_TOT*NL_TOT + 1
          IOFF = (IL_START - 1)*ND_TOT*NDUM + 1
C?        WRITE(6,*) ' IOFF in OT_T ', IOFF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) 'OT1(IOFF) block '
            CALL WRTMAT(OT1(IOFF),NR_OT1,NC_OT1,NR_OT1,NC_OT1)
            WRITE(6,*) ' TSCR1 '
            CALL WRTMAT(TSCR1,NR_T2,NC_T2,NR_T2,NC_T2)
          END IF
          ZERO = 0.0D0
          CALL MATML7(TSCR3,OT1(IOFF),TSCR1, NR_OT12, NC_OT12 ,
     &         NR_OT1,NC_OT1,NR_T2,NC_T2,ZERO,SIGNG,1)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' TSCR fresh from MATML7'
            CALL WRTMAT(TSCR3,NR_OT12,NC_OT12,NR_OT12,NC_OT12)
          END IF
*. Expand OT1T2(IDUM,L,K) to OT1T2(IDUM,I) 
          IF(NTEST.GE.100) 
     &    WRITE(6,*) ' OT1T2(L,K) =>  OT1T2(I) reordering '
C?        WRITE(6,*) ' SIGNLK before call to TI_TO_TOKBN(2) ', SIGNLK
          CALL TI_TO_TOKBN(NSMST,
     &         ISM_CINI_L,ISM_CAINI_L,ISM_AAINI_L,ISM_L,
     &         INM_ABINI_L,INM_AAINI_L,INM_CBINI_L,INM_CAINI_L,
     &         IL_BATLEN,NL,
     &         ISM_CINI_K,ISM_CAINI_K,ISM_AAINI_K,ISM_K,
     &         INM_ABINI_K,INM_AAINI_K,INM_CBINI_K,INM_CAINI_K,
     &         IK_BATLEN,NK,
     &         IX_LK_CA,SX_LK_CA,IB_LK(1,1,1),
     &         IX_LK_CB,SX_LK_CB,IB_LK(1,1,2),
     &         IX_LK_AA,SX_LK_AA,IB_LK(1,1,3),
     &         IX_LK_AB,SX_LK_AB,IB_LK(1,1,4),
     &         NLK(1,1),NLK(1,2),NLK(1,3),NLK(1,4),
     &         TSCR3,OT1T2,IBI_TCC,2,NDUM,SIGNLK,0)
        END DO
*       ^ End of loop over batches of L
      END DO
*.    ^ End of loop over batches of K
*
      RETURN
      END
C       RANK_FOR_CAAB(IHTP(1,IHOP),NEX,NDEEX)
      SUBROUTINE RANK_FOR_CAAB(IHTP,NEX,NDEEX)
*
* An operator IHTP is given in CAAB form. 
* Obtain number of excitation and deexcitation operators
*
* IHPVGAS_AB is used to determine the hp nature of spinorbitals
*
* Jeppe Olsen, March 2003
*
c      INCLUDE 'implicit.inc'
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER IHTP(NGAS,4)
*
      NEX = 0
      NDEEX = 0
*
      DO ICA = 1, 2
        DO IAB = 1, 2
          ICAAB = (ICA-1)*2+ IAB
          DO IGAS = 1, NGAS
            IF((ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.2).OR.
     &         (ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.1)    ) THEN
                NEX = NEX  + IHTP(IGAS,ICAAB)
            END IF
            IF((ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.2).OR.
     &         (ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.1)    ) THEN
                NDEEX = NDEEX  + IHTP(IGAS,ICAAB)
            END IF
          END DO
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CAAB operator : '
        CALL WRT_SPOX_TP(IHTP,1)
C             WRT_SPOX_TP(IEX_TP,NEX_TP)
        WRITE(6,*) ' Number of excitation and deexcitation ops : ',
     &  NEX,NDEEX
      END IF
*
      RETURN
      END
      SUBROUTINE OPCT1T2_M(IHEXOCC,IHD1OCC,IHD2OCC,IT1OCC,IT2OCC,
     &           T1,T2,HT1T2,IT1SM,IT2SM,IHSM,ISIGN1,FACX)
*
* Outer routine for contraction of Operator like Hamiltonian 
* with atmost two operators - hiding reference to scratch arrays
* and a few parameters
*
* Jeppe Olsen, March 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER IHEXOCC(NGAS,4),IHD1OCC(NGAS,4),IHD2OCC(NGAS,4) 
      INTEGER IT1OCC(NGAS,4),IT2OCC(NGAS,4)
      DIMENSION T1(*),T2(*)
C from MEM
      COMMON/CC_SCR2/KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &               KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTOCC7,
     &               KLSTOCC8,KLSTREO,
     &               KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &               KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &               KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &               KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &               KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &               KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &               KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &               KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &               KIX5_CA,KSX5_CA,KIX5_CB,KSX5_CB,
     &               KIX5_AA,KSX5_AA,KIX5_AB,KSX5_AB,
     &               KIX6_CA,KSX6_CA,KIX6_CB,KSX6_CB,
     &               KIX6_AA,KSX6_AA,KIX6_AB,KSX6_AB,
     &               KIX7_CA,KSX7_CA,KIX7_CB,KSX7_CB,
     &               KIX7_AA,KSX7_AA,KIX7_AB,KSX7_AB,
     &               KIX8_CA,KSX8_CA,KIX8_CB,KSX8_CB,
     &               KIX8_AA,KSX8_AA,KIX8_AB,KSX8_AB,
     &               KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
     &               KLTSCR5,KLTSCR6,
     &               KLOPSCR,
     &               KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &               KLSMD1,KLSMD2,KLSMD3,KLSMD4,
     &               KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &               KLNMD1,KLNMD2,KLNMD3,KLNMD4,
     &               KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &               KLOCK1, KLOCK2, KLOCK3,KLOCK4, 
     &               KLOCL1, KLOCL2, KLOCL3,KLOCL4, KL_IBF, 
     &               KLEXEORD
*. input and output
      DIMENSION HT1T2(*)
*. Call OPCT1T2 with all parameters ...
      CALL OPCT1T2(IHEXOCC,IHD1OCC,IHD2OCC,IT1OCC,IT2OCC,
     &     T1,T2,HT1T2,IT1SM,IT2SM,IHSM,LCCBD12,LCCB,
     &     WORK(KLIOD2_ST),WORK(KLIOD1_ST),WORK(KLIOEX_ST),
     &     WORK(KIX1_CA),WORK(KSX1_CA),WORK(KIX1_CB),WORK(KSX1_CB),
     &     WORK(KIX1_AA),WORK(KSX1_AA),WORK(KIX1_AB),WORK(KSX1_AB),
     &     WORK(KIX2_CA),WORK(KSX2_CA),WORK(KIX2_CB),WORK(KSX2_CB),
     &     WORK(KIX2_AA),WORK(KSX2_AA),WORK(KIX2_AB),WORK(KSX2_AB),
     &     WORK(KIX3_CA),WORK(KSX3_CA),WORK(KIX3_CB),WORK(KSX3_CB),
     &     WORK(KIX3_AA),WORK(KSX3_AA),WORK(KIX3_AB),WORK(KSX3_AB),
     &     WORK(KIX4_CA),WORK(KSX4_CA),WORK(KIX4_CB),WORK(KSX4_CB),
     &     WORK(KIX4_AA),WORK(KSX4_AA),WORK(KIX4_AB),WORK(KSX4_AB),
     &     WORK(KLSTOCC1),WORK(KLSTOCC2),WORK(KLSTOCC3),
     &     WORK(KLSTOCC4),WORK(KLSTOCC5),WORK(KLSTOCC6),
     &     WORK(KLTSCR1),WORK(KLTSCR2),WORK(KLTSCR3),
     &     WORK(KLTSCR4),WORK(KLOPSCR),
     &     WORK(KLSMD1),WORK(KLSMD2),WORK(KLSMK1),WORK(KLSMK2),
     &     WORK(KLSMEX),WORK(KLSML1),
     &     WORK(KLNMD1),WORK(KLNMD2),WORK(KLNMK1),WORK(KLNMK2),
     &     WORK(KLNMEX),WORK(KLNML1),IEXD1D2_INDX,
     &     WORK(KLOCK1),WORK(KLOCK2),WORK(KLOCL1),
     &     WORK(KLZ),WORK(KLZSCR),WORK(KLSTREO),ISIGN1,FACX,
     &     N_TDL_MAX)
*
      RETURN
      END
      SUBROUTINE EXP_MT_H_EXP_TC_MEM
*
* set up scratch space for direct calculation of H EXP T
* for commutator/connected approach 
*
* Pointers to scratch are returned in CC_SCR2
*
*. Jeppe Olsen, April 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
*
      COMMON/CC_SCR2/KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &               KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTOCC7,
     &               KLSTOCC8,KLSTREO,
     &               KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &               KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &               KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &               KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &               KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &               KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &               KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &               KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &               KIX5_CA,KSX5_CA,KIX5_CB,KSX5_CB,
     &               KIX5_AA,KSX5_AA,KIX5_AB,KSX5_AB,
     &               KIX6_CA,KSX6_CA,KIX6_CB,KSX6_CB,
     &               KIX6_AA,KSX6_AA,KIX6_AB,KSX6_AB,
     &               KIX7_CA,KSX7_CA,KIX7_CB,KSX7_CB,
     &               KIX7_AA,KSX7_AA,KIX7_AB,KSX7_AB,
     &               KIX8_CA,KSX8_CA,KIX8_CB,KSX8_CB,
     &               KIX8_AA,KSX8_AA,KIX8_AB,KSX8_AB,
     &               KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
     &               KLTSCR5,KLTSCR6,
     &               KLOPSCR,
     &               KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &               KLSMD1,KLSMD2,KLSMD3,KLSMD4,
     &               KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &               KLNMD1,KLNMD2,KLNMD3,KLNMD4,
     &               KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &               KLOCK1, KLOCK2, KLOCK3,KLOCK4, 
     &               KLOCL1, KLOCL2, KLOCL3,KLOCL4, KL_IBF, 
     &               KLEXEORD
C     OPCT1234(IOEX,IO1DX,IO2DX,IO3DX,IO4DX,
C    &           IT1,IT2,IT3,IT4,T1,T2,T3,T4,OT1234,
C    &           IT1SM,IT2SM,IT3SM,IT4SM,IOPSM,
C    &           LDB,LB,
C    &           IX1_CA,SX1_CA,IX1_CB,SX1_CB,  +
C    &           IX1_AA,SX1_AA,IX1_AB,SX1_AB,  +
C    &           IX2_CA,SX2_CA,IX2_CB,SX2_CB,  +
C    &           IX2_AA,SX2_AA,IX2_AB,SX2_AB,  +
C    &           IX3_CA,SX3_CA,IX3_CB,SX3_CB,  +
C    &           IX3_AA,SX3_AA,IX3_AB,SX3_AB,  +
C    &           IX4_CA,SX4_CA,IX4_CB,SX4_CB,  +
C    &           IX4_AA,SX4_AA,IX4_AB,SX4_AB,  +
C    &           IX5_CA,SX5_CA,IX5_CB,SX5_CB,  +
C    &           IX5_AA,SX5_AA,IX5_AB,SX5_AB,  +
C    &           IX6_CA,SX6_CA,IX6_CB,SX6_CB,  +
C    &           IX6_AA,SX6_AA,IX6_AB,SX6_AB,  +
C    &           IX7_CA,SX7_CA,IX7_CB,SX7_CB,  +
C    &           IX7_AA,SX7_AA,IX7_AB,SX7_AB,  +
C    &           IX8_CA,SX8_CA,IX8_CB,SX8_CB,  +
C    &           IX8_AA,SX8_AA,IX8_AB,SX8_AB,  +
C    &           IOC_K1,IOC_K2,IOC_K3,IOC_K4,  +
C    &           IOC_L1,IOC_L2,IOC_L3,IOC_L4,  +
C    &           ISTR_D1,ISTR_D2,ISTR_D3,ISTR_D4,ISTR_EX,ISTR_K,ISTR_L, +
C    &           TSCR1,TSCR2,TSCR3, +
C    &           TSCR4,TSCR5,TSCR6,OPSCR, +
C    &           ISM_CAAB_D1,ISM_CAAB_D2,ISM_CAAB_D3,ISM_CAAB_D4, +
C    &           ISM_CAAB_EX,ISM_CAAB_K,ISM_CAAB_L, +
C    &           INM_CAAB_D1,INM_CAAB_D2,INM_CAAB_D3,INM_CAAB_D4, +
C    &           INM_CAAB_EX,INM_CAAB_K,INM_CAAB_L,+
C    &           IEXD1234_INDX,
C    &           IOC_OT1T2,IZ,IZSCR,ISTREO,
C    &           ISIGNG,FACX,N_TDL_MAX)

*.1 : Maps for T1T2 => T12 for individual strings types 
      CALL DIM_T1T2_TO_T12_MAP(LEN_T1T2_STRING,LENT_T1T2_TCC)
      LEN = LEN_T1T2_STRING
*
      CALL MEMMAN(KIX1_CA,LEN,'ADDL  ',1,'IX1_CA')
      CALL MEMMAN(KSX1_CA,LEN,'ADDL  ',2,'SX1_CA')
*
      CALL MEMMAN(KIX1_CB,LEN,'ADDL  ',1,'IX1_CB')
      CALL MEMMAN(KSX1_CB,LEN,'ADDL  ',2,'SX1_CB')
*
      CALL MEMMAN(KIX1_AA,LEN,'ADDL  ',1,'IX1_AA')
      CALL MEMMAN(KSX1_AA,LEN,'ADDL  ',2,'SX1_AA')
*
      CALL MEMMAN(KIX1_AB,LEN,'ADDL  ',1,'IX1_AB')
      CALL MEMMAN(KSX1_AB,LEN,'ADDL  ',2,'SX1_AB')
*
      CALL MEMMAN(KIX2_CA,LEN,'ADDL  ',1,'IX2_CA')
      CALL MEMMAN(KSX2_CA,LEN,'ADDL  ',2,'SX2_CA')
*
      CALL MEMMAN(KIX2_CB,LEN,'ADDL  ',1,'IX2_CB')
      CALL MEMMAN(KSX2_CB,LEN,'ADDL  ',2,'SX2_CB')
*
      CALL MEMMAN(KIX2_AA,LEN,'ADDL  ',1,'IX2_AA')
      CALL MEMMAN(KSX2_AA,LEN,'ADDL  ',2,'SX2_AA')
*
      CALL MEMMAN(KIX2_AB,LEN,'ADDL  ',1,'IX2_AB')
      CALL MEMMAN(KSX2_AB,LEN,'ADDL  ',2,'SX2_AB')
*
      CALL MEMMAN(KIX3_CA,LEN,'ADDL  ',1,'IX3_CA')
      CALL MEMMAN(KSX3_CA,LEN,'ADDL  ',2,'SX3_CA')
*
      CALL MEMMAN(KIX3_CB,LEN,'ADDL  ',1,'IX3_CB')
      CALL MEMMAN(KSX3_CB,LEN,'ADDL  ',2,'SX3_CB')
*
      CALL MEMMAN(KIX3_AA,LEN,'ADDL  ',1,'IX3_AA')
      CALL MEMMAN(KSX3_AA,LEN,'ADDL  ',2,'SX3_AA')
*
      CALL MEMMAN(KIX3_AB,LEN,'ADDL  ',1,'IX3_AB')
      CALL MEMMAN(KSX3_AB,LEN,'ADDL  ',2,'SX3_AB')
*
      CALL MEMMAN(KIX4_CA,LEN,'ADDL  ',1,'IX4_CA')
      CALL MEMMAN(KSX4_CA,LEN,'ADDL  ',2,'SX4_CA')
*
      CALL MEMMAN(KIX4_CB,LEN,'ADDL  ',1,'IX4_CB')
      CALL MEMMAN(KSX4_CB,LEN,'ADDL  ',2,'SX4_CB')
*
      CALL MEMMAN(KIX4_AA,LEN,'ADDL  ',1,'IX4_AA')
      CALL MEMMAN(KSX4_AA,LEN,'ADDL  ',2,'SX4_AA')
*
      CALL MEMMAN(KIX4_AB,LEN,'ADDL  ',1,'IX4_AB')
      CALL MEMMAN(KSX4_AB,LEN,'ADDL  ',2,'SX4_AB')
*
      CALL MEMMAN(KIX5_CA,LEN,'ADDL  ',1,'IX5_CA')
      CALL MEMMAN(KSX5_CA,LEN,'ADDL  ',2,'SX5_CA')
*
      CALL MEMMAN(KIX5_CB,LEN,'ADDL  ',1,'IX5_CB')
      CALL MEMMAN(KSX5_CB,LEN,'ADDL  ',2,'SX5_CB')
*
      CALL MEMMAN(KIX5_AA,LEN,'ADDL  ',1,'IX5_AA')
      CALL MEMMAN(KSX5_AA,LEN,'ADDL  ',2,'SX5_AA')
*
      CALL MEMMAN(KIX5_AB,LEN,'ADDL  ',1,'IX5_AB')
      CALL MEMMAN(KSX5_AB,LEN,'ADDL  ',2,'SX5_AB')
*
      CALL MEMMAN(KIX6_CA,LEN,'ADDL  ',1,'IX6_CA')
      CALL MEMMAN(KSX6_CA,LEN,'ADDL  ',2,'SX6_CA')
*
      CALL MEMMAN(KIX6_CB,LEN,'ADDL  ',1,'IX6_CB')
      CALL MEMMAN(KSX6_CB,LEN,'ADDL  ',2,'SX6_CB')
*
      CALL MEMMAN(KIX6_AA,LEN,'ADDL  ',1,'IX6_AA')
      CALL MEMMAN(KSX6_AA,LEN,'ADDL  ',2,'SX6_AA')
*
      CALL MEMMAN(KIX6_AB,LEN,'ADDL  ',1,'IX6_AB')
      CALL MEMMAN(KSX6_AB,LEN,'ADDL  ',2,'SX6_AB')
*
      CALL MEMMAN(KIX7_CA,LEN,'ADDL  ',1,'IX7_CA')
      CALL MEMMAN(KSX7_CA,LEN,'ADDL  ',2,'SX7_CA')
*
      CALL MEMMAN(KIX7_CB,LEN,'ADDL  ',1,'IX7_CB')
      CALL MEMMAN(KSX7_CB,LEN,'ADDL  ',2,'SX7_CB')
*
      CALL MEMMAN(KIX7_AA,LEN,'ADDL  ',1,'IX7_AA')
      CALL MEMMAN(KSX7_AA,LEN,'ADDL  ',2,'SX7_AA')
*
      CALL MEMMAN(KIX7_AB,LEN,'ADDL  ',1,'IX7_AB')
      CALL MEMMAN(KSX7_AB,LEN,'ADDL  ',2,'SX7_AB')
*
      CALL MEMMAN(KIX8_CA,LEN,'ADDL  ',1,'IX8_CA')
      CALL MEMMAN(KSX8_CA,LEN,'ADDL  ',2,'SX8_CA')
*
      CALL MEMMAN(KIX8_CB,LEN,'ADDL  ',1,'IX8_CB')
      CALL MEMMAN(KSX8_CB,LEN,'ADDL  ',2,'SX8_CB')
*
      CALL MEMMAN(KIX8_AA,LEN,'ADDL  ',1,'IX8_AA')
      CALL MEMMAN(KSX8_AA,LEN,'ADDL  ',2,'SX8_AA')
*
      CALL MEMMAN(KIX8_AB,LEN,'ADDL  ',1,'IX8_AB')
      CALL MEMMAN(KSX8_AB,LEN,'ADDL  ',2,'SX8_AB')
*
*.4 : KLZ,KLZSCR : Memory for a Z matrix and scratch for 
*     constructing Z
*     
      IATP = 1
      IBTP = 2
      NAEL = NELFTP(IATP)
      NBEL = NELFTP(IBTP)
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      CALL MEMMAN(KLZ,LZ,'ADDL  ',1,'Z     ')
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'ZSCR  ')
*. String occupations for given CAAB, all symmetris
      LEN = MX_ST_TSOSO_BLK_MX*NSMST*4
      CALL MEMMAN(KLSTOCC1,LEN,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTOCC2,LEN,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTOCC3,LEN,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTOCC4,LEN,'ADDL  ',1,'STOCC4')
      CALL MEMMAN(KLSTOCC5,LEN,'ADDL  ',1,'STOCC5')
      CALL MEMMAN(KLSTOCC6,LEN,'ADDL  ',1,'STOCC6')
      CALL MEMMAN(KLSTOCC7,LEN,'ADDL  ',1,'STOCC7')
      CALL MEMMAN(KLSTOCC8,LEN,'ADDL  ',1,'STOCC8')

*. Reorder array for given CAAB, all symmetries 
      CALL MEMMAN(KLSTREO,MX_ST_TSOSO_MX*NSMST,'ADDL  ',1,'STREO ')
*. Intermediate blocks, (LD12B,LD12B)
*. For the moment 
      LD12B = LCCBD12
      LB = LCCB
      LEN = LD12B**2
      LEN3 = LD12B**3
      WRITE(6,*) ' LCCB, LCCBD12 = ', LCCB, LCCBD12
*. Should TSCR1-TSCT6 all have dim LEN3 ?
      CALL MEMMAN(KLTSCR1,LEN3,'ADDL  ',2,'TSCR1 ')
      CALL MEMMAN(KLTSCR2,LEN3,'ADDL  ',2,'TSCR2 ')
      CALL MEMMAN(KLTSCR3,LEN3,'ADDL  ',2,'TSCR3 ')
*
*
C          K_RES_DIM(ISPOBEX,NSPOBEX,MAXOP,NELMNT_MAX)
      CALL K_RES_DIM(WORK(KLSOBEX),NSPOBEX_TP,4,MX_KBLK)    
COLD  CALL T_DL_DIM(WORK(KLSOBEX),NSPOBEX_TP,4,N_TDL_MAX,LD12B)
C     WRITE(6,*) ' MX_KBLK in H_EXP ... = ', MX_KBLK
*. A big one for storing T(D2,L)
      N_TDL_MAX = 6000000
      LEN = N_TDL_MAX
      WRITE(6,*) ' N_TDL_MAX increased to  ', N_TDL_MAX
      CALL MEMMAN(KLTSCR4,LEN,'ADDL  ',2,'TSCR4 ')
      CALL MEMMAN(KLTSCR5,LEN,'ADDL  ',2,'TSCR5 ')
      CALL MEMMAN(KLTSCR6,LEN,'ADDL  ',2,'TSCR6 ')
*. For batches of indeces of Hamiltionian operator
      LEN = 4*LD12B
      CALL MEMMAN(KLIOD1_ST,LEN,'ADDL  ',2,'IOD1_S')
      CALL MEMMAN(KLIOD2_ST,LEN,'ADDL  ',2,'IOD2_S')
      CALL MEMMAN(KLIOEX_ST,LEN,'ADDL  ',2,'IOEX_S')
*and for a batch of Hamiltonian
      LOPSCR = LD12B*LD12B*LB
*. The above is not correct as the excitation part of the 
*. integrals are obtained as a single block. 
*. The largest block of integrals (hp!hp) should therefore
*. be contained in OPSCR. I will just pt just use an 
*. a block that is too large ...
      LOPSCR = MXTOB**4
      CALL MEMMAN(KLOPSCR,LOPSCR,'ADDL  ',2,'OPSCR ')
*. Number and symmetries of each substring for 6 complete T-blocks 
*. Largest number of strings in intermediate arrays
*
*. CA,CB,AA,AB for batches of strings 
      LEN = 4*LB
      CALL MEMMAN(KLSMD1,LEN,'ADDL  ',1,'SM_D1 ')
      CALL MEMMAN(KLNMD1,LEN,'ADDL  ',1,'NM_D1 ')
*
      CALL MEMMAN(KLSMD2,LEN,'ADDL  ',1,'SM_D2 ')
      CALL MEMMAN(KLNMD2,LEN,'ADDL  ',1,'NM_D2 ') 
*
      CALL MEMMAN(KLSMD3,LEN,'ADDL  ',1,'SM_D3 ')
      CALL MEMMAN(KLNMD3,LEN,'ADDL  ',1,'NM_D3 ') 
*
      CALL MEMMAN(KLSMD4,LEN,'ADDL  ',1,'SM_D4 ')
      CALL MEMMAN(KLNMD4,LEN,'ADDL  ',1,'NM_D4 ') 
*
      CALL MEMMAN(KLSMEX,LEN,'ADDL  ',1,'SM_EX ')
      CALL MEMMAN(KLNMEX,LEN,'ADDL  ',1,'NM_EX ') 
*. K1, K2, L1 holds general strings
C     LEN = 4*MX_KBLK
      LEN = 4*LB
      CALL MEMMAN(KLSMK1,LEN,'ADDL  ',1,'SM_K1 ')
      CALL MEMMAN(KLNMK1,LEN,'ADDL  ',1,'NM_K1 ') 
*
*
      CALL MEMMAN(KLSMK2,LEN,'ADDL  ',1,'SM_K2 ')
      CALL MEMMAN(KLNMK2,LEN,'ADDL  ',1,'NM_K2 ') 
*
      CALL MEMMAN(KLSML1,LEN,'ADDL  ',1,'SM_L1 ')
      CALL MEMMAN(KLNML1,LEN,'ADDL  ',1,'NM_L1 ') 
*     Occupation of intermediate strings 
      LEN = 4*NGAS
      CALL MEMMAN(KLOCK1,LEN,'ADDL  ',1,'IOC_K1')
      CALL MEMMAN(KLOCK2,LEN,'ADDL  ',1,'IOC_K2')
      CALL MEMMAN(KLOCK3,LEN,'ADDL  ',1,'IOC_K3')
      CALL MEMMAN(KLOCK4,LEN,'ADDL  ',1,'IOC_K4')
      CALL MEMMAN(KLOCL1,LEN,'ADDL  ',1,'IOC_L1')
      CALL MEMMAN(KLOCL2,LEN,'ADDL  ',1,'IOC_L2')
      CALL MEMMAN(KLOCL3,LEN,'ADDL  ',1,'IOC_L3')
      CALL MEMMAN(KLOCL4,LEN,'ADDL  ',1,'IOC_L4')
*. Array for offsets in F
      CALL MEMMAN(KL_IBF,NSPOBEX_TP,'ADDL  ',1,'IB_F  ')
*. Array for execution order of spinorbital types in Exp T
      CALL MEMMAN(KLEXEORD,NSPOBEX_TP,'ADDL  ',1,'EXEORD')
*
      RETURN
      END
      SUBROUTINE WRT_CNTR3(ICONT,NCONT,LDIM)
*
* Write contraction operator ICONT in the form used in VNEWCCV
*
      INCLUDE 'implicit.inc'
      INTEGER ICONT(LDIM,3)
*
      WRITE(6,*)
     &' Information about operators to be contracted'
      WRITE(6,*) ' Gasspace  Cr/An   Spin  '
      WRITE(6,*) ' ========================'
      DO JCONT = 1, NCONT
        WRITE(6,'(I4,6X,I2,7X,I2,7X,I2)')
     &  ICONT(JCONT,1),ICONT(JCONT,2),ICONT(JCONT,3)
      END DO
*
      RETURN
      END
      SUBROUTINE WRT_3SM_ARRAY(ISYM, NSYM, MXPSYM)
*
* Write array containing three symmetry labels 
*
      INCLUDE 'implicit.inc'
*
      INTEGER ISYM(MXPSYM,MXPSYM,MXPSYM)
*
      DO I3 = 1, NSYM
        WRITE(6,*) ' Value of third index = ', I3
        CALL IWRTMA(ISYM(1,1,I3),NSYM,NSYM,MXPSYM,MXPSYM)
      END DO
*
      RETURN
      END
      SUBROUTINE WRITE_NORM_BLOCKS_IN_VEC(VEC,NBLOCK,LBLOCK)
*
* A blocked vector VEC is given.
*.Find and print norm of each block
*
*. Jeppe Olsen ( For debugging another CC code)
*. Jan 2005
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Input
      INTEGER LBLOCK(NBLOCK)
      DIMENSION VEC(*)
*
      WRITE(6,*) ' Norm of various blocks of vector : '
      WRITE(6,*) ' ==================================='
      IOFF = 1
      DO JBLOCK = 1, NBLOCK
       LEN = LBLOCK(JBLOCK)
       X = INPROD(VEC(IOFF),VEC(IOFF),LEN)
       WRITE(6,'(A,I5,E22.15)')
     & ' Block and norm ', JBLOCK, SQRT(X)
       IOFF = IOFF + LEN
      END DO
*
      RETURN
      END

c $Id$
