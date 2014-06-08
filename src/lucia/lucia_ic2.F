      SUBROUTINE GEN_IC_ORBOP3(IWAY,NIC_ORBOP,IC_ORBOP,MX_OP_RANK,
     &                     MN_OP_RANK,MX_CREA,MN_CREA,NCMA,
     &                     IONLY_EXCOP,IREFSPC,
     &                     ITREFSPC,IADD_UNI)
*
* Generate general orbital excitation types to be used in 
* various types of internal  contraction 
* The orbital excitations working on IREFSPC should contain 
* an component in space ITREFSPC.
*
* Generalization of GEN_IC_ORBOP to general number of elementary 
* operators ( instead of single and double excitations)
*
* If IADD_UNI = 1, the unit operator ( containing zero operators)
* is added at the end
*
* Jeppe Olsen, November 2005
*
*
* IWAY = 1 : Number of orbital excitations for internal contraction
* IWAY = 2 : Generate also the actual orbital excitations 
*
* MX_OP_RANK : Largest allowed orbital rank
* MN_OP_RANK : Smallest allowed orbital rank
* MN_CREA : Smallest allowed number of creation operators
* MX_CREA : Largest allowed number of creation operators
* NCMA    : Difference between number of creation and annihilation operators

* IONLY_EXCOP = 1 => only excitation operators ( no annihilation in particle 
*                    space, no creation in inactive space )
*
*. Rank is defined as # crea of particles + # anni of holes 
*                    -# crea of holes     - # anni of particles

      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Local scratch
      INTEGER ITREFOCC(MXPNGAS,2)
*. Output ( if IWAY .ne. 1 ) 
      INTEGER IC_ORBOP(2*NGAS,*)
*. Local scratch
      INTEGER IOP(2*MXPNGAS)
*
      NTEST =   05
      IZERO = 0
*
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' IREFSPC, ITREFSPC = ', IREFSPC, ITREFSPC 
      NIC_ORBOP =  0
C     DO NCREA = MN_CREA, MX_CREA
        NANNI = NCREA - NCMA
        IF(NANNI.GT.0) THEN
*. Loop over ways of generating 
*. 
*. Include single excitations
*. Single excitations a+i a j
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          CALL ISETVC(IOP,IZERO,2*NGAS)
          IOP(IGAS) = 1
          IOP(NGAS+JGAS) = 1
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Next Orbital excitation '
            CALL IWRTMA(IOP,NGAS,2,NGAS,2)
          END IF
C              IRANK_ORBOP(IOP,NEX,NDEEX)
C              COMPARE_OPDIM_ORBDIM(IOP,IOKAY)
          CALL COMPARE_OPDIM_ORBDIM(IOP,IOKAY)
          IF(NTEST.GE.100) WRITE(6,*) ' IOKAY from COMPARE..', IOKAY
*. Is the action of this operator on IREFSPC included in ITREFSPC
      CALL ORBOP_ACCOCC(IOP,IGSOCCX(1,1,IREFSPC),ITREFOCC,NGAS,MXPNGAS)
      CALL OVLAP_ACC_MINMAX(ITREFOCC,IGSOCCX(1,1,ITREFSPC),NGAS,MXPNGAS,
     &         IOVERLAP)
      IF(NTEST.GE.100) WRITE(6,*) ' IOVERLAP from OVLAP..',IOVERLAP
      IF(IOVERLAP.EQ.0) IOKAY = 0
C     ORBOP_ACCOCC(IORBOP,IACC_IN,IACC_OUT,NGAS,MXPNGAS)
C     OVLAP_ACC_MINMAX(IACC1,IACC2,NGAS,MXPNGAS,IOVERLAP)
*. is there any operators in spaces that are frozen or deleted in ITREFSPC
C     CHECK_EXC_FR_OR_DE(IOP,IOCC,NGAS,IOKAY)
          CALL CHECK_EXC_FR_OR_DE(IOP,IGSOCCX(1,1,ITREFSPC),NGAS,IOKAY2)
          IF(NTEST.GE.100) WRITE(6,*) ' IOKAY2 from CHECK ... ', IOKAY2
          IF(IOKAY2.EQ.0) IOKAY = 0
          IF(IOKAY.EQ.1) THEN
            CALL IRANK_ORBOP(IOP,NEX,NDEEX)
            IOKAY2 = 1
            IF(IONLY_EXCOP.EQ.1.AND.NDEEX.NE.0) IOKAY2 = 0
            IRANK = NEX - NDEEX
            IF(NTEST.GE.100) WRITE(6,*) ' IRANK = ', IRANK
            IF(MN_OP_RANK.LE.IRANK.AND.IRANK.LE.MX_OP_RANK
     &      .AND.IOKAY2.EQ.1)THEN
              NIC_ORBOP  = NIC_ORBOP + 1
              IF(NTEST.GE.100) WRITE(6,*) ' Operator included '
              IF(IWAY.NE.1) 
     &        CALL ICOPVE(IOP,IC_ORBOP(1,NIC_ORBOP),2*NGAS)
            END IF
          END IF
        END DO
      END DO
      END IF
*. Double excitations a+i a+j a k a l
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              CALL ISETVC(IOP,IZERO,2*NGAS)
              IOP(IGAS) = 1
              IOP(JGAS) = IOP(JGAS) + 1
              IOP(NGAS+KGAS) = 1
              IOP(NGAS+LGAS) = IOP(NGAS+LGAS) + 1
              CALL COMPARE_OPDIM_ORBDIM(IOP,IOKAY)
*. Is the action of this operator on IREFSPC included in ITREFSPC
      CALL ORBOP_ACCOCC(IOP,IGSOCCX(1,1,IREFSPC),ITREFOCC,NGAS,MXPNGAS)
      CALL OVLAP_ACC_MINMAX(ITREFOCC,IGSOCCX(1,1,ITREFSPC),NGAS,
     &         MXPNGAS,IOVERLAP)
      IF(IOVERLAP.EQ.0) IOKAY = 0
          CALL CHECK_EXC_FR_OR_DE(IOP,IGSOCCX(1,1,ITREFSPC),NGAS,IOKAY2)
              IF(IOKAY2.EQ.0) IOKAY = 0
              IF(IOKAY.EQ.1) THEN
                CALL IRANK_ORBOP(IOP,NEX,NDEEX)
                IRANK = NEX - NDEEX
                IOKAY2 = 1
                IF(IONLY_EXCOP.EQ.1.AND.NDEEX.NE.0) IOKAY2 = 0
                IF(MN_OP_RANK.LE.IRANK.AND.IRANK.LE.MX_OP_RANK.AND.
     &            IOKAY2.EQ.1) THEN
                  IF(NTEST.GE.100) WRITE(6,*) ' Operator included '
                  NIC_ORBOP  = NIC_ORBOP + 1
                  IF(IWAY.NE.1) 
     &            CALL ICOPVE(IOP,IC_ORBOP(1,NIC_ORBOP),2*NGAS)
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
      IF(IADD_UNI.EQ.1) THEN
        NIC_ORBOP = NIC_ORBOP + 1
        IF(IWAY.NE.1) THEN
           IZERO = 0
           CALL ISETVC(IC_ORBOP(1,NIC_ORBOP),IZERO,2*NGAS)
        END IF
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Number of orbitalexcitation types generated ',
     &               NIC_ORBOP
        IF(IWAY.NE.1) THEN
         WRITE(6,*) ' And the actual orbitalexcitation types : '
         DO JC = 1, NIC_ORBOP
           WRITE(6,*) ' Orbital excitation type ', JC
           CALL IWRTMA(IC_ORBOP(1,JC),NGAS,2,NGAS,2) 
         END DO
        END IF
      END IF
*
      RETURN
      END 
c $Id$
