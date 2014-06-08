      SUBROUTINE DIM_SPII(IASPGRP,IBSPGRP,IOBTP,IAB,IAC,NSPII)
*
* Obtain Dimension of S(P,Iac,Iob) array 
*
* Jeppe Olsen, december 1999
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
*. Local scratch
      INTEGER NPA_A(MXPNSMST),NPA_B(MXPNSMST), NPA_AB(MXPNSMST)
      INTEGER NAC_A(MXPNSMST),NAC_B(MXPNSMST)
      INTEGER KGRP(MXPNGAS)
     
*.
*. Obtain supergroup obtained by annihilating/creating electron
*  of type IOBTP to IASPGRP
      IF(IAB.EQ.1) THEN
        CALL ACOP_SPGRP(NGAS,ISPGPFTP(1,IASPGRP),IAC,IOBTP,
     &                  KGRP,ICODE)
C            ACOP_SPGRP(NNGRP,IGRP_IN,IAC,IGAS,IGRP_OUT,ICODE)
      ELSE
        CALL ACOP_SPGRP(NGAS,ISPGPFTP(1,IBSPGRP),IAC,IOBTP,
     &                  KGRP,ICODE)
      END IF
*
      IF(ICODE.EQ.-1) THEN
        WRITE(6,*) ' DIM_SPII : ACOP cannot generate group '
        STOP       ' DIM_SPII : ACOP cannot generate group '
      ELSE IF (ICODE.EQ.0) THEN
*. Group is trivially vanishing
        NSPII = 0
      ELSE   
*. Output group is nontrivial 
*. Dimension for active and passive parts of A and B
        IF(IAB.EQ.1) THEN
          CALL NST_PA(ISPGPFTP(1,IASPGRP),NGAS,1,IOBTP,NAC_A,NPA_A)
          CALL NST_PA(ISPGPFTP(1,IBSPGRP),NGAS,0,IOBTP,NAC_B,NPA_B)
C              NST_PA(ISPGP,NIGRP,NACTE,IACTE,NST_ACT,NST_PAS)
        ELSE
          CALL NST_PA(ISPGPFTP(1,IBSPGRP),NGAS,1,IOBTP,NAC_B,NPA_B)
          CALL NST_PA(ISPGPFTP(1,IASPGRP),NGAS,0,IOBTP,NAC_A,NPA_A)
        END IF
*. Total Number of passive strings per sym
        CALL NST_TWO_TYPES(NPA_A,NPA_B,NPA_AB)
*. and then : the largest block
        L_SPII_MX = 0
        DO ISSM = 1, NSMST
         DO KASM = 1, NSMST
           L_SPII = 0
           DO ISM = 1, NSMST
             NI = NOBPTS(IOBTP,ISM)
             DO IB_AC_SM = 1, NSMST
               IA_AC_SM = MULTD2H(KASM,ISM)
               IAB_AC_SM = MULTD2H(IA_AC_SM,IB_AC_SM)
               IAB_PA_SM = MULTD2H(ISSM,IAB_AC_SM)
               NP = NPA_AB(IAB_PA_SM)
               IF(IAB.EQ.1) THEN
                 NAC = NAC_B(IB_AC_SM)
               ELSE
                 NAC = NAC_A(IA_AC_SM)
               END IF
               L_SPII = L_SPII + NI*NP*NAC
             END DO
           END DO
           L_SPII_MX = MAX(L_SPII_MX,L_SPII)
         END DO
        END DO
        NSPII = L_SPII_MX
      END IF
*     ^ End of ICODE switches
*
*. Temporarily set to zero
      NSPII = 0
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' IASPGRP, IBSPGRP, IOBTP IAB IAC ', 
     &               IASPGRP, IBSPGRP, IOBTP,IAB,IAC
        WRITE(6,*) ' Number of SPII strings ', NSPII
      END IF
*
      RETURN
      END
      SUBROUTINE NST_TWO_TYPES(NST1,NST2,NST12)
*
* Number of strings per sym in string types 1,2 =>
* Number of strings per type in combined string
*
* Jeppe Olsen, December 1999
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
*. Input
      INTEGER NST1(NSMST),NST2(NSMST)
*. Output
      INTEGER NST12(NSMST)
* 
      DO ISM_12 = 1, NSMST
        NST = 0
        DO ISM_1 = 1, NSMST
          ISM_2 = MULTD2H(ISM_1, ISM_12)
          NST = NST + NST1(ISM_1)*NST2(ISM_2)
        END DO
        NST12(ISM_12) = NST
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings in combined type '
        CALL IWRTMA(NST12,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
      SUBROUTINE NST_PA(ISPGP,NIGRP,NACTE,IACTE,NST_ACT,NST_PAS)
*
* Number of active and passive strings 
*
* Jeppe Olsen, December 1999
*
*
* ======
*. Input
* ======
*   ISPGP   : Groups defining supergroup
*   NISPGP  : Number of groups in supergroup
*   NACTE    : Number of active orbital spaces(needs not all be distinct)
*   IACTE    : The active gas orbital spaces (needs not all be distinct)
*
* ======
* Output
* ======
*
*    NST_ACT: Number of active strings per symmetry
*    NST_PAS: Number of passive strings per symmetry
*
      INCLUDE 'wrkspc.inc'
*. Specific Input
      DIMENSION IACTE(NACTE)
      DIMENSION ISPGP(NIGRP)
*. General Input
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Output
*. Number of active strings per symmetry
      DIMENSION NST_ACT(NSMST)
*. Number of passive strings per symmetry 
      DIMENSION NST_PAS(NSMST)
*
*  ==============
*. Local scratch
*  ==============
*
      INTEGER IPAS(MXPNGAS),IACT(MXPNGAS)
      INTEGER IACTGP(MXPNGAS),IPASGP(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' DIM_PA in service '
        WRITE(6,*) ' ================= '
        WRITE(6,*)
        WRITE(6,*) '   Groups '
        WRITE(6,*) ' =============================================='
        WRITE(6,*) 
C       WRITE(6,*)
        WRITE(6,'(A,16(1X,I2))') 
     & '     ',(ISPGP(II),II=1,NIGRP)
        WRITE(6,*) ' Number of operators ', NACTE
        WRITE(6,*) ' Gaspaces of operators ', (IACTE(I),I=1,NACTE)
      END IF
*
*. Info on Active operators : There can be several operators in
*  a active space. Obtain distinct active spaces
      IZERO = 0
      CALL ISETVC(IACTGP,IZERO,MXPNGAS)
      DO JACT = 1, NACTE
        IACTGP(IACTE(JACT)) = 1
      END DO
      NACT = 0
*. Active groups 
      DO JGRP = 1, NIGRP
        IF(IACTGP(IGSFGP(ISPGP(JGRP))).EQ.1) THEN
          NACT = NACT + 1
          IACT(NACT) = JGRP
        END IF
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of active groups ', NACT
        WRITE(6,*) 
     &  ' Active gasspaces(index in IGRP) : ', (IACT(JACT),JACT=1,NACT)
      END IF
*. passive groups (we know these types)
      NPAS  = NIGRP - NACT
      JACT = 1
      JPAS = 0
      DO JGRP = 1, NIGRP
        IF(JGRP.EQ.IACT(JACT)) THEN
          IF (JACT.LT.NACT) JACT = JACT + 1
        ELSE
          JPAS = JPAS + 1
          IPAS(JPAS) = JGRP
        END IF
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Number of passive types ',NPAS
        WRITE(6,*) ' Passive indeces in IGRP'
        CALL IWRTMA(IPAS,1,NPAS,1,NPAS)
      END IF
*. Actual groups constituting active and passive parts
      JACT = 1
      JPAS = 1
      DO IIGRP = 1, NIGRP
        IF(IGSFGP(ISPGP(IIGRP)).EQ.IACT(JACT)) THEN
          IACTGP(JACT) = ISPGP(IIGRP)
          IF(JACT.LT.NACT) JACT = JACT + 1
        ELSE IF ( IGSFGP(ISPGP(IIGRP)).EQ.IPAS(JPAS) ) THEN
          IPASGP(JPAS) = ISPGP(IIGRP)
          IF(JPAS.LT.NPAS) JPAS = JPAS + 1
        END IF
      END DO
*. Number of passive and active groups in I strings (unmodified strings)
*. per sym
      DO ISTSYM = 1, NSMST
*
        CALL NST_SPGRP2(NACT,IACTGP,ISTSYM,NSMST,NSTRIN,NDIST)
        NST_ACT(ISTSYM) = NSTRIN
*
        CALL NST_SPGRP2(NPAS,IPASGP,ISTSYM,NSMST,NSTRIN,NDIST)
        NST_PAS(ISTSYM) = NSTRIN
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from DIM_PA : '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) ' Number of Active strings per sym '
        CALL IWRTMA(NST_ACT,1,NSMST,1,NSMST)
        WRITE(6,*) ' Number of Passive strings per sym '
        CALL IWRTMA(NST_PAS,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
c $Id$
