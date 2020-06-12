*. A note on GICCI expansions
*
* |0> = C_0 |ref> + O_1|ref> +  O_2 O_1|ref> + ... +
*       O_N ....  O_1|ref> 
*
* The C_0 coefficient is saved as the last element in the 
* collected vector
* Each O operators consists of a excitation operator and a
* projection operator
* O_I = P_I T_I
* The projection operator is pt just projecting a single space out.
*.There is an indirect projection also as O_I ... I_1|ref> is/should be
* evaluated in CI-space I, but this is not pt included.
*
* When a given operator O_I is optimized, this corresponds to
* optimizing the linear expansion
*
* |O_new> = Delta_0(C_0 |ref> + O_1|ref> + ... O_{I-1} ... O_1|ref>
*         + sum_mu delta_{mu I} (O_{I+1} + .... O_N ... O_{I+1}) 
*                               tau_{mu I} O_{I-1} ... O_1 |ref>
*
* an optimization consists this of determining Delta_0 and delta_{mu I}.
*. Note that the vector to be multiplied by Delta_0 depends upon I.
*
* In practice: in the optimization of a given vector, Delta_0 is stored
* in the element corresponding to the unit-operator.
*
* The GICCI vector corresponding to a set of elements (delta_{mu I},
* delta_0) = (delta, delta_0) is obtained as
*
* I = 1:
* -----
* C_0(new) = delta_0 C_0
* T_1(new) = delta
* T_J(new) = T_J for J> 1
*
* I > 1:
* ------
* C_0(new) = delta_0 C_0
* T_1(new) = T_1*delta_0
* T_I(new) = delta/delta_0
* T_J(new) = T_J for J neq 1,I
*
* The optimization of a given GICCI operator, corresponds to a
* linear variational space spanned by the basisvectors
* (C_0 |ref> + O_1|ref> + ... O_{I-1} ... O_1|ref>
* and 
* (O_{I+1} + .... O_N ... O_{I+1}) tau_{mu I} O_{I-1} ... O_1 |ref>
*
* These vectors are fixed and do not depend on the expansion of O(I)
* 

      SUBROUTINE LUCIA_GIC(ICTYP,EREF,EFINAL,CONVER,VNFINAL)
*
*
* Master routine for General internally contracted CI calculations,
* Sprin 10 version 
*
*
* Jeppe Olsen, March 2010 looking into contracted CI with several
*              operators
*
* Assumed  spaces
*  Space 1: Reference HF or CAS
*  Space 2: Space where standard CI is performed
*  Space 3,4..: Spaces where internal contracted CI will be performed
* 
      INCLUDE 'wrkspc.inc'
      REAL*8 
     &INPROD
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
*. Transfer common block for communicating with H_EFF * vector routines
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
*. A bit of local scratch
      DIMENSION ICASCR(MXPNGAS)
      CHARACTER*6 ICTYP
      LOGICAL CONVER
*
      EXTERNAL MTV_FUSK, STV_FUSK
      EXTERNAL H_S_EFF_ICCI_TV,H_S_EXT_ICCI_TV
      EXTERNAL HOME_SD_INV_T_ICCI
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICCI  ')
*. I will play with spinadaptation in this routine so 
         I_SPIN_ADAPT = 1
*
      NTEST = 10
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' Generalized Internal contracted section entered '
        WRITE(6,*) ' =============================================== '
        WRITE(6,*)
        WRITE(6,'(A,A)') ' Form of calculation  ', ICTYP
        WRITE(6,*) '  Symmetri of reference vector ' , IREFSM 
        WRITE(6,*)
        WRITE(6,*) ' Number of external operators ', NTEXC_G
        WRITE(6,*) ' Parameters defining internal contraction '
*
        WRITE(6,*) ' Form of External operators: '
        WRITE(6,*) 
     &  ' Op.,  Min. and Max exc. rank, int-exc, Proj. and Final space'
        WRITE(6,*) 
     &  ' ------------------------------------------------------------'
        DO IEXC_G = 1, NTEXC_G
         IF(ICEXC_INT_G(IEXC_G).EQ.1) THEN
          WRITE(6,'(1H ,1X,I2,4X,I2,7X,I2,14X,A,3X,I2,8X,I2)')
     &    IEXC_G, ICEXC_RANK_MIN_G(IEXC_G),ICEXC_RANK_MAX_G(IEXC_G),
     &    '  +  ', IPTCSPC_G(IEXC_G),ITCSPC_G(IEXC_G)
         ELSE
          WRITE(6,'(1H ,1X,I2,4X,I2,7X,I2,14X,A,3X,I2,8X,I2)')
     &    IEXC_G, ICEXC_RANK_MIN_G(IEXC_G),ICEXC_RANK_MAX_G(IEXC_G),
     &    '  -  ', IPTCSPC_G(IEXC_G),ITCSPC_G(IEXC_G)
         END IF
        END DO
*  
C        IF(ICEXC_INT.EQ.1) THEN
C          WRITE(6,*) 
C    &   ' Internal (ina->ina, sec->sec) excitations allowed'
C        ELSE
C          WRITE(6,*) 
C    &   ' Internal (ina->ina, sec->sec) excitations not allowed'
C        END IF
        WRITE(6,*) 
     &  '  Largest number of vectors in iterative supspace ', MXCIV
        WRITE(6,*) 
     &  '  Largest initial number of vectors in iterative supspace ',
     &    MXVC_I
        IF(IRESTRT_IC.EQ.1) THEN
          WRITE(6,*) ' Restarted calculation : '
          WRITE(6,*) '      IC coefficients  read from LUSC54'
          WRITE(6,*) '      CI for reference read from LUSC54 '
        END IF
      END IF
*
      IDUM = 0
*. Divide orbital spaces into inactive, active, secondary using 
*. space 1
      CALL CC_AC_SPACES(1,IREFTYP)
*
      MX_ST_TSOSO_MX = 0
      MX_ST_TSOSO_BLK_MX = 0
      MX_TBLK_MX = 0
      MX_TBLK_AS_MX = 0
      MAXLEN_I1_MX = 0
*
* Generate information  about T-operators
*
      DO IEX_G = 1, NTEXC_G
        IF(NTEST.GE.10) WRITE(6,*) ' T-excitation type = ', IEX_G
*
        ICEXC_RANK_MIN = ICEXC_RANK_MIN_G(IEX_G)
        ICEXC_RANK_MAX = ICEXC_RANK_MAX_G(IEX_G)
        ICEXC_INT      = ICEXC_INT_G(IEX_G)
*. these are transferred through CRUN 
        IF(IEX_G.EQ.1) THEN
*. Initial reference space is first space by assumption
          IREFSPC = 1
        ELSE
          IREFSPC = ITREFSPC
        END IF
        ITREFSPC = ITCSPC_G(IEX_G)
C       GET_TEX_INFO(ICEXC_RANK_MIN,ICEXC_RANK_MAX,ICEXC_INT,
C                         IREFSPC,ITREFSPC,
C    &           MX_ST_TSOSO, MX_ST_TSOSO_BLK, MX_TBLK,  MX_TBLK_AS)
        CALL GET_TEX_INFO(IREFSPC,ITREFSPC,
     &       MX_ST_TSOSO, MX_ST_TSOSO_BLK, MX_TBLK,  MX_TBLK_AS)
*
        MX_ST_TSOSO_MX = MAX(MX_ST_TSOSO_MX,MX_ST_TSOSO)
        MX_ST_TSOSO_BLK_MX = MAX(MX_ST_TSOSO_BLK_MX,MX_ST_TSOSO_BLK)
        MX_TBLK_MX = MAX(MX_TBLK_MX,MX_TBLK)
        MX_TBLK_AS_MX = MAX(MX_TBLK_AS_MX,MX_TBLK_AS)
        MAXLEN_I1_MX = MAX(MAXLEN_I1_MX,MAXLEN_I1)
*
        I_FT_GLOBAL = 2
        CALL TRANSFER_T_OFFSETS(I_FT_GLOBAL,IEX_G)
      END DO
      MAXLEN_I1 = MAXLEN_I1_MX
*
      IF(I_SPIN_ADAPT.EQ.1) THEN
*. A bit of general info on prototype spin combinations
      CALL PROTO_SPIN_MAT
*. Set up information about partial spin adaptation
        DO IEX_G = 1, NTEXC_G
*. Put information about excitations in place
          I_FT_GLOBAL = 1
          CALL TRANSFER_T_OFFSETS(I_FT_GLOBAL,IEX_G)
*. Information about partial spin adaptation for this T excitation type
          CALL GET_SP_INFO
*. And save offsets and arrays
          I_FT_GLOBAL = 2
          CALL TRANSFER_SPIN_OFFSETS(I_FT_GLOBAL,IEX_G)
        END DO
      END IF
*. Prepare calculation with first T-operator
      I_FT_GLOBAL = 1
      CALL TRANSFER_T_OFFSETS(I_FT_GLOBAL,1)
      CALL TRANSFER_SPIN_OFFSETS(I_FT_GLOBAL,1)
*. Initial space is first space by assumption ( of Jeppe)
      IREFSPC = 1
      ITREFSPC = ITCSPC_G(1)
      WRITE(6,*) ' IREFSPC, ITREFSPC = ', IREFSPC, ITREFSPC
*
      IF(ICTYP(1:4).EQ.'ICCI') THEN
*
*                    ==============================
*                    Internal contracted CI section 
*                    ==============================
*
* Solve Internal contracted CI problem 
         CALL LUCIA_ICCI(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                 EFINAL,CONVER,VNFINAL)
*
      ELSE IF(ICTYP(1:5).EQ.'GICCI') THEN
*. Generalized intetnal contraction CI
         CALL LUCIA_GICCI(ICTYP,EREF,
     &                 EFINAL,CONVER,VNFINAL)

      ELSE IF(ICTYP(1:4).EQ.'ICPT') THEN
*
*                    ==========================================
*                    Internal contracted Perturbation expansion 
*                    ==========================================
*
        CALL LUCIA_ICPT(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                 EFINAL,CONVER,VNFINAL)
*
      ELSE IF(ICTYP(1:4).EQ.'ICCC') THEN
* Internal contracted coupled cluster 
*
*                    ======================================
*                    Internal contracted Coupled Cluster 
*                    =======================================
*
        CALL LUCIA_ICCC(IREFSPC,ITREFSPC,ICTYP,EREF,EFINAL,
     &                  CONVER,VNFINAL)
      END IF
*
*.
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICCI  ')
*
      RETURN
      END
      SUBROUTINE GEN_IC_IN_ORBSPC(IWAY,NIC_ORBOP,IC_ORBOP,MX_OP_NUM,
     &                               IORBSPC)
*
* Generate orbitalexcitations for a given  orbital space with 
* the restriction that the number of creation- or annihilationoperators
* is less or equal to MX_OP_NUM. No check are performed to see 
* whether operators are non-vanishing for given space.
*
* Jeppe Olsen, For generating cumulants in a given orbitalsubspace 
*
* IWAY = 1 : Number of orbital excitations for internal contraction
* IWAY = 2 : Generate also the actual orbital excitations 
*

      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Output ( if IWAY .ne. 1 ) 
      INTEGER IC_ORBOP(2*NGAS,*)
*. Local scratch
      INTEGER IOP(2*MXPNGAS)
*
      NTEST =   05
      IZERO = 0
*
      NIC_ORBOP =  0
      DO NOP = 1, MX_OP_NUM
        CALL ISETVC(IOP,IZERO,2*NGAS)
        IOP(IORBSPC) = NOP
        IOP(NGAS+IORBSPC) = NOP
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Next Orbital excitation '
          CALL IWRTMA(IOP,NGAS,2,NGAS,2)
        END IF
        NIC_ORBOP  = NIC_ORBOP + 1
        IF(IWAY.NE.1) CALL ICOPVE(IOP,IC_ORBOP(1,NIC_ORBOP),2*NGAS)
      END DO
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
      SUBROUTINE GEN_IC_ORBOP(IWAY,NIC_ORBOP,IC_ORBOP,MX_OP_RANK,
     &                     MN_OP_RANK,IONLY_EXCOP,IREFSPC,ITREFSPC,
     &                     IADD_UNI,IPRNT)
*
* Generate single and double 
* orbital excitation types corresponding to internal contraction  
* The orbital excitations working on IREFSPC should contain 
* an component in space ITREFSPC.
*
* If IADD_UNI = 1, the unit operator ( containing zero operators)
* is added at the end
*
* Jeppe Olsen, August 2002
*
*
* IWAY = 1 : Number of orbital excitations for internal contraction
* IWAY = 2 : Generate also the actual orbital excitations 
*
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
      NTEST =   0
      NTEST = MAX(NTEST,IPRNT)
      IZERO = 0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ------------------------------'
        WRITE(6,*) ' Information from GEN_IC_ORBOP '
        WRITE(6,*) ' ------------------------------'
        WRITE(6,*)
        WRITE(6,*) ' IREFSPC, ITREFSPC = ', IREFSPC, ITREFSPC 
      END IF
*
      NIC_ORBOP =  0
      I_INCLUDE_SX = 1
      IF(I_INCLUDE_SX.EQ.0) THEN
        DO I = 1, 200
          WRITE(6,*) ' Excitation operators are excluded '
       END DO
      ELSE
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
      IF(NTEST.GE.2) THEN
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
      SUBROUTINE IRANK_ORBOP(IOP,NEX,NDEEX)
*
*     An orbital operator is given in IOP 
*     Find RANK of the operator
*
*     Find number of excitation ops  (# crea of particles + # anni of holes )
*                  deexcitation ops  (# crea of holes     + # anni of particles)
*     IHPVGAS in CGAS is used to determine types of orbitals
*
* Jeppe Olsen, August 2002
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER IOP(NGAS,2)
*
      NEX = 0
      NDEEX = 0
*
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) THEN
            NDEEX = NDEEX + IOP(IGAS,1)
            NEX   = NEX   + IOP(IGAS,2)
         ELSE IF (IHPVGAS(IGAS).EQ.2) THEN
            NEX = NEX + IOP(IGAS,1)
            NDEEX = NDEEX + IOP(IGAS,2)
         END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
         WRITE(6,*) ' Orbital excitation operator '
         WRITE(6,*) ' =========================== '
         CALL IWRTMA(IOP,NGAS,2,NGAS,2)
         WRITE(6,*)
         WRITE(6,*) ' Number of excitation operators ', NEX
         WRITE(6,*) ' Number of deexcitation operators ', NDEEX
      END IF
*
      RETURN
      END
      SUBROUTINE COMPARE_OPDIM_ORBDIM(IOP,IOKAY)
*
* Compare dimensions of orbitaloperator in CA form and 
* orbitals, and check that number of crea- or anni-operators 
* is smaller than number of orbitals in each gas space
*
* Jeppe Olsen, August 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Integer 
      INTEGER IOP(NGAS,2)
*
      IOKAY = 1
      DO ICA = 1, 2
        DO IGAS = 1, NGAS
          IF(IOP(IGAS,ICA).GT.2*NOBPT(IGAS)) IOKAY = 0
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orbital operator '
        CALL IWRTMA(IOP,NGAS,2,NGAS,2)
        IF(IOKAY.EQ.1) THEN
           WRITE(6,*) ' Operator is nonvanishing '
        ELSE
           WRITE(6,*) ' Operator is vanishing '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE GET_NCA_FOR_ORBOP(NORBEX,IORBEX,NC_FOR_OBEX,
     &           NA_FOR_OBEX,NGAS)
*
* Find number of creation and annihilation operators for set 
* of orbital excitation operators
*
* Jeppe Olsen, September 2002
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IORBEX(NGAS,2,NORBEX)
*. Output
      INTEGER NC_FOR_OBEX(NORBEX),NA_FOR_OBEX(NORBEX)
*
      DO I = 1, NORBEX
        NC_FOR_OBEX(I) = IELSUM(IORBEX(1,1,I),NGAS)
        NA_FOR_OBEX(I) = IELSUM(IORBEX(1,2,I),NGAS)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of creations per orbital operator '
        CALL IWRTMA(NC_FOR_OBEX,1,NORBEX,1,NORBEX)
        WRITE(6,*) ' Number of annihilations per orbital operator '
        CALL IWRTMA(NA_FOR_OBEX,1,NORBEX,1,NORBEX)
      END IF
*
      RETURN
      END
      SUBROUTINE ORBOP_ACCOCC(IORBOP,IACC_IN,IACC_OUT,NGAS,MXPNGAS)
*
* An orbital excitation CA form and an CI space in the form of 
* an accumulated occupation are given. Find accumulated occupation 
* of product 
*
* Jeppe Olsen, September 2002
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IORBOP(NGAS,2), IACC_IN(MXPNGAS,2) 
*. Output
      INTEGER IACC_OUT(MXPNGAS,2)
*
      IDEL = 0
      DO IGAS = 1, NGAS
        IDEL = IDEL + IORBOP(IGAS,1) - IORBOP(IGAS,2)
        IACC_OUT(IGAS,1) = MAX(0,IACC_IN(IGAS,1) + IDEL)
        IACC_OUT(IGAS,2) = MAX(0,IACC_IN(IGAS,2) + IDEL)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input ORBOP in CA form '
        CALL IWRTMA(IORBOP,NGAS,2,NGAS,2)
        WRITE(6,*) ' Input OCC in acc min/max form '
        CALL IWRTMA(IACC_IN,NGAS,2,MXPNGAS,2)
        WRITE(6,*) ' Output OCC in acc min/max form '
        CALL IWRTMA(IACC_OUT,NGAS,2,MXPNGAS,2)
      END IF
*
      RETURN
      END
      SUBROUTINE OVLAP_ACC_MINMAX(IACC1,IACC2,NGAS,MXPNGAS,IOVERLAP)
*
* Two spaces are given in the form of accumulated MAX/MIN 
* occupations. Check if the two spaces overlap, ie. there is
* a nonvanishing space that is contained in both.
*
* Jeppe Olsen, Sept 2002
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IACC1(MXPNGAS,2), IACC2(MXPNGAS,2)
*
      IOVERLAP = 1
      DO IGAS = 1, NGAS
*. Find common Min  being the Max of the individual Mins
        IMIN_12 = MAX(IACC1(IGAS,1),IACC2(IGAS,1))
*. Find common Max  being the Min of the individual Maxs
        IMAX_12 = MIN(IACC1(IGAS,2),IACC2(IGAS,2))
        IF(IMIN_12.GT.IMAX_12) IOVERLAP = 0
CE      IF(.NOT.( (IACC2(IGAS,1).GE.IACC1(IGAS,1).AND.
CE   &      IACC2(IGAS,1).LE.IACC1(IGAS,2)     ) .OR.
CE   &     (IACC2(IGAS,2).GE.IACC1(IGAS,1).AND.
CE   &      IACC2(IGAS,2).LE.IACC1(IGAS,2))    )       ) THEN 
CE        IOVERLAP = 0
CE      END IF
      END DO
*
      NTEST = 00 
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Two accumulated min/max occupations '
        CALL IWRTMA(IACC1,NGAS,2,MXPNGAS,2)
        CALL IWRTMA(IACC2,NGAS,2,MXPNGAS,2)
        IF(IOVERLAP.EQ.1) THEN
          WRITE(6,*) ' The occupations overlap '
        ELSE 
          WRITE(6,*) ' The occupations do not overlap '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE CHECK_EXC_FR_OR_DE(IOP,IOCC,NGAS,IOKAY)
*
* An orbital operator IOP in CA form and and occupation space 
* IOCC in accumulated min/max form is given. Ensure that there
* are no operators in frozen, ie. completely occupied spaces 
* spaces and that no operators are in deleted orbspaces, 
*.that is spaces with zero electrons
*       IOKAY = 1 => No such operators
*             = 0 0>    such operators occurs in IOP
*
* Jeppe Olsen, Sept 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER IOCC(MXPNGAS,2),IOP(NGAS,2)
*
      IOKAYL = 1
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          NELMIN = IOCC(1,1)
        ELSE
          NELMIN = IOCC(IGAS,1)-IOCC(IGAS-1,2)
        END IF
        NOP = IOP(IGAS,1) + IOP(IGAS,2)
*. Check to see if orbital space is deleted, i.e. 
*. contains no electrons
        IDELETED = 0
        IF(IGAS.EQ.1) THEN 
          IF(IOCC(1,2).EQ.0) IDELETED = 1
        ELSE 
          IF(IOCC(IGAS,2).EQ.IOCC(IGAS-1,1)) IDELETED = 1
        END IF
  
        IF(NOP.NE.0.AND.IDELETED.EQ.1) IOKAYL = 0
        IF(NOP.NE.0.AND.NELMIN.EQ.2*NOBPT(IGAS)) IOKAYL = 0
      END DO
*
      IF(IOKAYL.EQ.1) THEN
        IOKAY = 1
      ELSE
        IOKAY = 0
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orbital operator in CA form '
        CALL IWRTMA(IOP,NGAS,2,NGAS,2)
        IF(IOKAY.EQ.1) THEN
          WRITE(6,*) ' No operators in frozen or deleted spaces '
        ELSE
          WRITE(6,*) ' Operators in frozen or deleted spaces '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE CHECK_EXC_FR(IOP,IOCC,NGAS,IOKAY)
*
* An orbital operator IOP in CA form and and occupation space 
* IOCC in accumulated min/max form is given. Ensure that there
* are no operators in frozen, ie. completely occupied spaces 
* spaces
*       IOKAY = 1 => No such operators
*             = 0 0>    such operators occurs in IOP
*
* Jeppe Olsen, Sept 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER IOCC(MXPNGAS,2),IOP(NGAS,2)
*
      IOKAYL = 1
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          NELMIN = IOCC(1,1)
        ELSE
          NELMIN = IOCC(IGAS,1)-IOCC(IGAS-1,2)
        END IF
        NOP = IOP(IGAS,1) + IOP(IGAS,2)
        IF(NOP.NE.0.AND.NELMIN.EQ.2*NOBPT(IGAS)) IOKAYL = 0
      END DO
*
      IF(IOKAYL.EQ.1) THEN
        IOKAY = 1
      ELSE
        IOKAY = 0
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orbital operator in CA form '
        CALL IWRTMA(IOP,NGAS,2,NGAS,2)
        IF(IOKAY.EQ.1) THEN
          WRITE(6,*) ' No operators in frozen or deleted spaces '
        ELSE
          WRITE(6,*) ' Operators in frozen or deleted spaces '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE ICCI_COMPLETE_MAT(IREFSPC,ITREFSPC,I_SPIN_ADAPT)
*
* Master routine for Internal contraction with complete incore 
* construction of all matrices
*
* Jeppe Olsen, Sept 2002
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
*. Scratch for CI 
*
      NTEST = 10
      WRITE(6,*) 
      WRITE(6,*) ' Complete H and S matrices will be constructed '
      WRITE(6,*) ' =============================================='
      WRITE(6,*)
      WRITE(6,*) ' Reference space is ', IREFSPC
      WRITE(6,*) ' Space of Operators times reference space ', ITREFSPC
      WRITE(6,*)
      WRITE(6,*) ' Number of parameters in spinuncoupled basis ', 
     &           N_CC_AMP
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'IC_CMP ')
*. Space for old fashioned CI behind the curtain
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
*
* Space for complete H and S matrices 
*
      LEN = N_CC_AMP ** 2
      CALL MEMMAN(KLSMAT,LEN,'ADDL  ',2,'SMAT  ')
      CALL MEMMAN(KLHMAT,LEN,'ADDL  ',2,'HMAT  ')
      CALL MEMMAN(KLSCR1,LEN,'ADDL  ',2,'SCR1_C')
      CALL MEMMAN(KLSCR2,LEN,'ADDL  ',2,'SCR2_C')
*. Add an extra matrix to allow for backtransformation to 
*. original basis as a test
      CALL MEMMAN(KLXORT,LEN,'ADDL  ', 2,'XORT  ')

*. And a few working vectors 
      CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,N_CC_AMP,'ADDL  ',2,'VCC3  ')
*. Identify the unit  operator i.e. the operator with 
*. zero creation and annihilation operators
      IDOPROJ = 1
      IF(IDOPROJ.EQ.1) THEN
        CALL GET_SPOBTP_FOR_EXC_LEVEL(0,WORK(KLCOBEX_TP),NSPOBEX_TP,
     &       NUNIOP,IUNITP,WORK(KLSOX_TO_OX))
*. And the position of the unitoperator in the list of SPOBEX operators
        WRITE(6,*) ' NUNIOP, IUNITP = ', NUNIOP,IUNITP
        IF(NUNIOP.EQ.0) THEN
          WRITE(6,*) ' Unitoperator not found in exc space '
          WRITE(6,*) ' I will proceed without projection '
          IDOPROJ = 0
        ELSE
C  IFRMR(WORK,IROFF,IELMNT)
          IUNIOP = IFRMR(WORK(KLIBSOBEX),1,IUNITP)
          WRITE(6,*) ' IUNIOP = ', IUNIOP
        END IF
      END IF
*
C     COM_SH(S,H,VCC1,VCC2,VCC3,VEC1,VEC2,
C    &                  N_CC_AMP,IREFSPC,ITREFSPC,
C    &                  LUC,LUHC,LUSCR,LUSCR2,IDOPROJ,IUNIOP)
      CALL COM_SH(WORK(KLSMAT),WORK(KLHMAT),WORK(KLVCC1),WORK(KLVCC2),
     &            WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &            N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &            IDOPROJ,IUNIOP,1,1,0,0,0,0,0)
*. Obtain singularities on S 
C     CHK_S_FOR_SING(S,NDIM,NSING,X,SCR)
      CALL CHK_S_FOR_SING(WORK(KLSMAT),N_CC_AMP,NSING,
     &                    WORK(KLSCR1),WORK(KLSCR2),WORK(KLVCC2))
*. On output the eigenvalues are residing in WORK(KLSCR2) and 
*. the corresponding eigenvectors in WORK(KLSCR1).
*. The singular subspace is defined by the first NSING eigenvectors
      NNONSING = N_CC_AMP - NSING
      WRITE(6,*) ' Number of nonsingular eigenvalues of S ', NNONSING
      KLNONSING = KLSCR1 + NSING*N_CC_AMP
*. For saving transformation matrix 
      CALL COPVEC(WORK(KLNONSING),WORK(KLXORT),NNONSING*N_CC_AMP)
*. Transform H to a nonsigular - and orthogonal basis
*. I use the transformation matrix 
*  X = U sigma^{-1/2}, where U are the nonsingular 
*. eigenvectors of S and sigma are the corresponding 
*. eigenvectors
*. This transformation matrix turns the nonsingular part of S into 
*. a unitmatrix
C?    WRITE(6,*) ' Unscaled transformation matrix '
C?    CALL WRTMAT(WORK(KLNONSING),N_CC_AMP,NNONSING,
C?   &                            N_CC_AMP,NNONSING)
      DO I = 1, NNONSING
        SCALE = 1/SQRT(WORK(KLSCR2-1+NSING+I))
        CALL SCALVE(WORK(KLNONSING+(I-1)*N_CC_AMP),SCALE,N_CC_AMP)
      END DO
C?    WRITE(6,*) ' Scaled transformation matrix '
C?    CALL WRTMAT(WORK(KLNONSING),N_CC_AMP,NNONSING,
C?   &                            N_CC_AMP,NNONSING)
*. Transform 
*. H Xin SCR2
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
C?    WRITE(6,*) ' H before transformation '
C?    CALL WRTMAT(WORK(KLHMAT),N_CC_AMP,N_CC_AMP,N_CC_AMP,N_CC_AMP)
      CALL MATML7(WORK(KLSCR2),WORK(KLHMAT),WORK(KLNONSING),
     &            N_CC_AMP,NNONSING,N_CC_AMP,N_CC_AMP,
     &            N_CC_AMP,NNONSING,FACTORC,FACTORAB,0)
C?    WRITE(6,*) ' H halftransformed '
C?    CALL WRTMAT(WORK(KLSCR2),N_CC_AMP,N_CC_AMP,N_CC_AMP,N_CC_AMP)
*. X(T) H X in HMAT
      CALL MATML7(WORK(KLHMAT),WORK(KLNONSING),WORK(KLSCR2),
     &            NNONSING,NNONSING,N_CC_AMP,NNONSING,
     &            N_CC_AMP,NNONSING,FACTORC,FACTORAB,1)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformed Hamiltonian matrix '
        CALL WRTMAT(WORK(KLHMAT),NNONSING,NNONSING,NNONSING,NNONSING)
      END IF
*
*. Diagonalize transformed Hamiltonian 
*
C DIAG_SYM_MAT(A,X,SCR,NDIM,ISYM)

      IOLD = 1
      IF(IOLD.EQ.0) THEN
      CALL DIAG_SYM_MAT(WORK(KLHMAT),WORK(KLSCR1),WORK(KLSCR2),
     &                  NNONSING,0)
      ELSE
        ZERO = 0.0D0
        ONE = 1.0D0
        CALL TRIPAK(WORK(KLHMAT),WORK(KLSCR1),1,NNONSING,NNONSING)
        CALL COPVEC(WORK(KLSCR1),WORK(KLHMAT),NNONSING*(NNONSING+1)/2)
        CALL SETVEC(WORK(KLSCR1),ZERO,NNONSING*NNONSING)
        CALL SETDIA(WORK(KLSCR1),ONE,NNONSING,0)
C            SETDIA(MATRIX,VALUE,NDIM,IPACK)
        CALL JACOBI(WORK(KLHMAT),WORK(KLSCR1),NNONSING,NNONSING)
C            JACOBI(F,V,NB,NMAX) 
        CALL COPDIA(WORK(KLHMAT),WORK(KLSCR2),NNONSING,1)
      END IF

*
      WRITE(6,*) ' Ecore in ICCI_COMPLETE.. ', ECORE
      DO I = 1, NNONSING
        WORK(KLSCR2-1+I) = WORK(KLSCR2-1+I) + ECORE
      END DO
*
      WRITE(6,*) ' Eigenvalues of H matrix in IC basis '
      WRITE(6,*) ' ===================================='
      CALL WRTMAT_EP(WORK(KLSCR2),1,NNONSING,1,NNONSING)
*
      IF(I_SPIN_ADAPT.EQ.1) THEN
*. First back transform first eigenvector to original basis 
        CALL MATML7(WORK(KLVCC2),WORK(KLXORT),WORK(KLSCR1),
     &              N_CC_AMP,1,N_CC_AMP,NNONSING,NNONSING,1,
     &              FACTORC,FACTORAB,0)
        WRITE(6,*) ' First eigenvector in CAAB basis '
        CALL WRTMAT(WORK(KLVCC2),1,N_CC_AMP,1,N_CC_AMP)
*. Reform to CSF basis 
        CALL REF_CCV_CAAB_SP(WORK(KLVCC2),WORK(KLVCC1),
     &                       WORK(KLVCC3),1)
*, And reform back to CAAB basis 
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLVCC2),ZERO,N_CC_AMP)
        CALL REF_CCV_CAAB_SP(WORK(KLVCC2),WORK(KLVCC1),
     &                       WORK(KLVCC3),2)
C     REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
*. Play a bit around with spin adaptation 
*. Reorder from CAAB to CONF order ICONF(I) = ICAAB(IREO(I))
*. corresponding to a gathering
C     (VECO,VECI,INDEX,NDIM)
C       CALL GATVEC(WORK(KLVCC1),WORK(KLVCC2),WORK(KLREORDER_CAAB),
C    &              N_CC_AMP)
C       WRITE(6,*) ' First eigenvector in conf order '
C       CALL WRTMAT(WORK(KLVCC1),1,N_CC_AMP,1,N_CC_AMP)
      END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'IC_CMP ')
      RETURN
      END 
      SUBROUTINE COM_SH(S,H,VCC1,VCC2,VCC3,VEC1,VEC2,
     &                  N_CC_AMP,IREFSPC,ITREFSPC,
     &                  LUC,LUHC,LUSCR,LUSCR2,IDOPROJ,IUNIOP,
     &                  IDO_S,IDO_H,IDO_SPA,I_DO_EI,NSPA,IDOSUB,
     &                  ISUB,NSUB)
*
* Construct complete S and M matrices for 
* Excitations defined in CC_TCC and 
* reference space on LUC
*
* If IDOPROJ = 1, then the reference space is projected out 
*                 for all operators except the unitoperator
*
* IF IDOSUB.NE.0, the matrix is constructed in the space 
* defined by the NSUB elements in ISUB
*
* IDO_S = 1 => S is constructed 
* IDO_H = 1 => H is constructed 
*
* If IDO_SPA = 1, the matrices are constructed in the spinadapted basis 
* If I_DO_EI = 1, the matrices are constructed in the orthonormal EI
* basis
*
* Jeppe Olsen, Sept 2002
*
* For IDOPROJ = 1 , we are interested in calculating the matrix 
*
*       ( <0!H!0>          <0!H!P Q_j!0>      )
*       ( <0!Q+(i)P!H!O>   <0!Q+(I)PH PQ(J)!0>)
*
* The projection operators in front of evrything but !0>
* induces some assymmetry that is organized by at the end calculating 
* explicitly 
*     <0!H!0> and <0!Q+(I)P!H0> and overwriting the corresponding column 
*      and row
*
      INCLUDE 'implicit.inc'
      REAL*8 INPRDD
* 
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
*. Input
      INTEGER ISUB(*)
*. Output
      DIMENSION S(*),H(*)
*. Scratch
      DIMENSION VCC1(*),VCC2(*),VCC3(*)
      DIMENSION VEC1(*),VEC2(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'COM_SH')
      IF(IDO_SPA.EQ.1.OR.I_DO_EI.EQ.1) THEN
         IUNIOP = NSPA
C?       WRITE(6,*) ' Unit operator is set to last operator '
      END IF
*
      NTEST = 1005
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' COM_SH speaking '
         WRITE(6,*) ' IDOPROJ, IUNIOP = ', IDOPROJ,IUNIOP
         WRITE(6,*) ' IDO_SPA, NSPA = ', IDO_SPA,NSPA
         WRITE(6,*) ' IDO_S, IDO_H, = ', IDO_S, IDO_H
         WRITE(6,*) ' IREFSPC, ITREFSPC = ', IREFSPC,ITREFSPC
         WRITE(6,*) '  LUC, LUHC, LUSCR = ', LUC, LUHC, LUSCR
      END IF
*. Number of excitations in calculation 
        NVAR = NSPA
*. Dimension of space in which S or H is constructed 
      IF(IDOSUB.EQ.0) THEN
        NSBVAR = NVAR
      ELSE
        NSBVAR = NSUB
      END IF
*
      IUNIOP_EFF = 0
      IF(IDOSUB.NE.0.AND.IUNIOP.NE.0) THEN
*. Check if unitoperator is included in list 
        CALL FIND_INTEGER_IN_VEC(IUNIOP,ISUB,NSUB,IUNIOP_EFF)
      ELSE IF(IUNIOP.NE.0) THEN
        IUNIOP_EFF = IUNIOP
      END IF
      WRITE(6,*) ' IUNIOP_EFF = ', IUNIOP_EFF
   
 
      LEN = NSBVAR**2

      ZERO = 0.0D0
      IF(IDO_S.EQ.1) CALL SETVEC(S,ZERO,LEN)
      IF(IDO_H.EQ.1) CALL SETVEC(H,ZERO,LEN)
*
*
*. Use new approach based on H,S times vector routines
*. It has not been checked with subspaces 
*.       
      WRITE(6,*) ' NEW route used to construct ICCI matrices '
      DO I = 1, NSBVAR
        IF(NTEST.GE.5) WRITE(6,*) 'Constructing row of S,H for I = ',I
        ZERO = 0.0D0
        CALL SETVEC(VCC1,ZERO,NVAR)
        IF(IDOSUB.EQ.0) THEN
          VCC1(I) = 1.0D0
        ELSE 
          VCC1(ISUB(I)) = 1.0D0
        END IF
*
*. Overlap terms 
*
        IF(IDO_S.EQ.1) THEN
          CALL H_S_EXT_ICCI_TV(VCC1,XDUM,VCC2,0,1)
          IF(IDOSUB.EQ.0) THEN
            CALL COPVEC(VCC2,S(1+(I-1)*NSBVAR),NSBVAR)
          ELSE
            CALL GATVEC(S(1+(I-1)*NSBVAR),VCC2,ISUB,NSBVAR)
          END IF
        END IF
*
*. Hamilton terms 
*
        IF(IDO_H.EQ.1) THEN
          CALL H_S_EXT_ICCI_TV(VCC1,VCC2,XDUM,1,0)
          IF(IDOSUB.EQ.0) THEN
            CALL COPVEC(VCC2,H(1+(I-1)*NSBVAR),NSBVAR)
          ELSE
            CALL GATVEC(H(1+(I-1)*NSBVAR),VCC2,ISUB,NSBVAR)
          END IF
        END IF
*
      END DO
*
      IF(NTEST.GE.100) THEN
         IF(IDO_S.EQ.1) THEN
           WRITE(6,*) ' Constructed S matrix '
           WRITE(6,*) ' ==================== '
           CALL WRTMAT(S,NSBVAR,NSBVAR,NSBVAR,NSBVAR)
         END IF
         IF(IDO_H.EQ.1) THEN
           WRITE(6,*) ' Constructed H matrix '
           WRITE(6,*) ' ======================'
           CALL WRTMAT(H,NSBVAR,NSBVAR,NSBVAR,NSBVAR)
         END IF
       END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'COM_SH')
*
      RETURN
      END 
      SUBROUTINE EXTR_CIV(ISM,ISPCIN,LUIN,
     &                  ISPCX,IEX_OR_DE,LUUT,LBLK,
     &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
* A vector of sym ISM and space ISPCIN is given in LUIN
* Extract(IEX_OR_DE=1) or delete (IEX_OR_DE = 2) the 
* parts of the CI vector that is in space ISPCX
*
* The output form is the same as the input form, only
* some blocks are zeroed.
*
* Jeppe Olsen, September 2002 from EXP_CIV
*
      INCLUDE 'wrkspc.inc'
C     IMPLICIT REAL*8(A-H,O-Z)
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'

*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXTR_C')
*
      NTESTL = 00
C     NTEST = MAX(NTESTG,NTESTL)
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' EXTR_CIV: Subspace to be modified ', ISPCX
      END IF
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
*
*. Allowed combinations of strings types for input and ISPCX
*. spaces
*
      CALL MEMMAN(KLIABI,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABI')
      CALL MEMMAN(KLIABX,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABU')
      CALL IAIBCM(ISPCIN,WORK(KLIABI))
      CALL IAIBCM(ISPCX,WORK(KLIABX))
*
* type of each symmetry block ( full, lower diagonal, absent )
*
      CALL MEMMAN(KLBLIN,NSMST,'ADDL  ',1,'KLBLIN')
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLIN),IDUMMY)
*. A scratch block 
      LENGTH = MXSOOB
      CALL MEMMAN(KLVEC,LENGTH,'ADDL  ',2,'LVEC  ')
*
      IF(NTEST.GE.1000) THEN
        CALL REWINO(LUIN)
        WRITE(6,*) ' Initial vectors in EXTR_CIV '
        DO IROOT = 1, NROOT
          WRITE(6,*) ' Root number ', IROOT 
          CALL WRTVCD(WORK(KLVEC),LUIN,0,-1)
        END DO
      END IF
*     ^ End of test
*
      CALL REWINO(LUIN)
      CALL REWINO(LUUT)
      DO IROOT = 1, NROOT
*. Input vector should be first vector on file so
        IF(IROOT.EQ.1) THEN
          LLUIN = LUIN
        ELSE
*. With the elegance of an elephant
          CALL REWINO(LUSCR)
          CALL REWINO(LUIN)
          DO JROOT = 1, IROOT
            CALL REWINO(LUSCR)
            CALL COPVCD(LUIN,LUSCR,WORK(KLVEC),0,-1)
          END DO
          CALL REWINO(LUSCR)
          LLUIN = LUSCR
        END IF
*. Expcivs may need the IAMPACK parameter ( in case it must write
*  a zero block before any blocks have been read in.
*  Use IDIAG to decide
        IF(IDIAG.EQ.1) THEN
          IAMPACK = 0
        ELSE
          IAMPACK = 1
        END IF
C       WRITE(6,*) ' IAMPACK in EXPCIV ', IAMPACK
*
        CALL EXTRCIVS(LLUIN,WORK(KLVEC),WORK(KLIABI),
     &       NOCTPA,NOCTPB,WORK(KLBLIN),
     &       LUUT,WORK(KLIABX),IEX_OR_DE,
     &       IDC,NSMST,LBLK,IAMPACK,ISMOST(1,ISM),
     &       WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)))
*
      END DO
*
      IF(ICOPY.NE.0) THEN
*. Copy expanded vectors to LUIN
        CALL REWINO(LUIN)
        CALL REWINO(LUUT)
        DO IROOT = 1, NROOT
          CALL COPVCD(LUUT,LUIN,WORK(KLVEC),0,-1)
        END DO
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output  vectors in EXTR_CIV '
*
        CALL REWINO(LUUT)
        DO IROOT = 1, NROOT
C?        WRITE(6,*) ' Root number ', IROOT 
            CALL WRTVCD(WORK(KLVEC),LUUT,0,-1)
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXTR_C')
*
      RETURN
      END 
C       CALL EXTRCIVS(LLUIN,WORK(KLBLI),WORK(KLIABIN),
C    &       NOCTPA,NOCTPB,WORK(KLBLIN),
C    &       LUUT,WORK(KLIABX),IEX_OR_DE,
C    &       IDC,NSMST,LBLK,IAMPACK,ISMOST(1,ISM),
C    &       WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)))
      SUBROUTINE EXTRCIVS(LUIN,VEC,IABIN,
     &                   NOCTPA,NOCTPB,IBLTPIN,
     &                   LUUT,IABX,IEX_OR_DE,
     &                   IDC,NSMST,LBLK,IAMPACKED_IN,
     &                   ISMOST,NSSOA,NSSOB)
*
* IEX_OR_DE = 1 : Copy those blocks of LUIN that are allowed according
*                 to IABX, set remaining blocks to 0
* IEX_OR_DE = 2 : Copy blocks of LUIN that are not allowed according
*                 to IABX, set remaining blocks to 0
*
* Input vector on LUIN, Output vector in LUUT
* Output vector is supposed on start of vector
*
* LUIN is assumed to be single vector file,
* so rewinding will place vector on start of vector
*
* Note that the form of the two files will be identical, 
* just that LUUT will contain some zero blocks 
*
* ALL ICISTR = 1 code has been removed
*
* Jeppe Olsen, September 2002 from EXPCIVS
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER IABIN(NOCTPA,NOCTPB),IABX(NOCTPA,NOCTPB)
      INTEGER IBLTPIN(NSMST)
*, Symmetry of other string, given total symmetry
      INTEGER ISMOST(NSMST)
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
*. Scratch 
      DIMENSION VEC(*)
*
*. Loop over TTS blocks of output vector
      IATP = 1
      IBTP = 1
      IASM = 0
 1000 CONTINUE
*. Next block 
        CALL NXTBLK(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSMST,
     &              IBLTPIN,IDC,NONEW,IABIN,ISMOST,
     &              NSSOA,NSSOB,LBLOCK,LBLOCKP)
        IF(IABX(IATP,IBTP).EQ.0) THEN
          IF(IEX_OR_DE.EQ.1) THEN
             ICOPY = 0
          ELSE 
             ICOPY = 1
          END IF
        ELSE 
          IF(IEX_OR_DE.EQ.1) THEN
             ICOPY = 1
          ELSE 
             ICOPY = 0
          END IF
        END IF
*
        IF(NONEW.EQ.0) THEN
          CALL IFRMDS(LENGTH,1,-1,LUIN)
          CALL FRMDSC(VEC,LENGTH,-1,LUIN,IMZERO,IAMPACK)
* 
          CALL ITODS(LENGTH,1,-1,LUUT)
          IF(ICOPY.EQ.0) THEN
            CALL ZERORC(-1,LUUT,IAMPACKED_IN)
          ELSE
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC(VEC,LENGTH,-1,LUUT)
            ELSE
              CALL TODSCP(VEC,LENGTH,-1,LUUT)
            END IF
          END IF
      GOTO 1000
        END IF
*. End of file on output vector 
      CALL ITODS(-1,1,-1,LUUT)
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' EXPTRCIVS Speaking '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' ============ '
        WRITE(6,*) ' Input Vector '
        WRITE(6,*) ' ============ '
        WRITE(6,*)
        CALL WRTVCD(VEC,LUIN,1,LBLK)
        WRITE(6,*)
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Output Vector '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        CALL WRTVCD(VEC,LUUT,1,LBLK)
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_CONF_FOR_ORBEX(NCOC_FSM,NAOC_FSM,ICOC,IAOC,
     &           NOP_C,NOP_A, IBCOC_FSM,IBAOC_FSM,NSMST,IOPSM,
     &           ICAOC)
*. Obtain the configurations for given C and A occupations
*
*. Jeppe Olsen, Sept. 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*
* ======
*. Input
* ======
*
*. Number of C and A occupations per symmetry
      INTEGER NCOC_FSM(NSMST), NAOC_FSM(NSMST)
*. Offset for C and A occupations of given sym
      INTEGER IBCOC_FSM(NSMST), IBAOC_FSM(NSMST)
*. And the actual C and A orbital configurations
      INTEGER ICOC(NOP_C,*), IAOC(NOP_A,*) 
*
* =======
*. Output
* =======
*
      INTEGER ICAOC(NOP_C+NOP_A,*)
* 
      NTEST = 10
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' C and A strings of sym 1 '
        CALL IWRTMA(ICOC,NOP_C,NCOC_FSM(1),NOP_C,NCOC_FSM(1))
        CALL IWRTMA(IAOC,NOP_A,NAOC_FSM(1),NOP_A,NAOC_FSM(1))
      END IF 
      JCONF = 0
      DO ICSM = 1, NSMST
        IASM = MULTD2H(IOPSM,ICSM)
        NC = NCOC_FSM(ICSM)
        NA = NAOC_FSM(IASM)
        DO IA = 1, NA
          DO IC = 1, NC
            IC_ABS = IBCOC_FSM(ICSM) - 1 + IC
            IA_ABS = IBAOC_FSM(IASM) - 1 + IA
            JCONF = JCONF + 1
            CALL ICOPVE(ICOC(1,IC_ABS),ICAOC(1,JCONF),NOP_C)
            CALL ICOPVE(IAOC(1,IA_ABS),ICAOC(1+NOP_C,JCONF),NOP_A)
          END DO
        END DO
*       ^ End of loop over C and A
      END DO
*     ^ End of loop over sym of C strings
      NCONF = JCONF  
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of operators in C and A ',NOP_C, NOP_A
        WRITE(6,*) ' List of CA configurations '
        WRITE(6,*) ' =========================='
        WRITE(6,*)
        WRITE(6,*) ' Creation part       Annihilation part '
        WRITE(6,*) ' ======================================'
        DO JCONF = 1, NCONF
          WRITE(6,'(1H , 20(1X,I3))') (ICAOC(I,JCONF),I=1,NOP_C+NOP_A)
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_CA_CONF_FOR_ORBEX(ICEX_TP,IAEX_TP,
     &           NCOC_FSM,NAOC_FSM,IBCOC_FSM,IBAOC_FSM,
     &           KCOC,KAOC,KZC,KZA,KCREO,KAREO)
*
* Obtain the occupations, Arc weights and reordering matrices 
* for a Creation and Annihilation types defined by ICEX_TP, IAEX_TP
*
*
* Jeppe Olsen, Sept 2002
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER ICEX_TP(NGAS),IAEX_TP(NGAS)
* 
*. Output 
*
*. Number of creation and annihilation occupations per symmetry
      INTEGER NCOC_FSM(MXPNSMST), NAOC_FSM(MXPNSMST)
*. Start of creation and annihilation occupations of given symmetry
      INTEGER IBCOC_FSM(MXPNSMST), IBAOC_FSM(MXPNSMST)
*
* A number of terms are delivered in arrays allocated in this 
* subroutine 
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' INFO from  GET_CA_CONF_FOR_ORBEX '
        WRITE(6,*) ' Creation excitation type '
        CALL IWRTMA(ICEX_TP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Annihilation excitation type '
        CALL IWRTMA(IAEX_TP,1,NGAS,1,NGAS)
      END IF
*
*  ================
*. Creation strings 
*  ================
*
*.Number of strings per symmetry

      IDUMMY = 0
      CALL GET_CONF_FOR_OCCLS(ICEX_TP,NCOC_FSM,IBCOC_FSM,IDUMMY,
     &                        NSMST,1)
*
*. the actual occupation
*
      NCOC_TOT = IELSUM(NCOC_FSM,NSMST)
      NELC = IELSUM(ICEX_TP,NGAS)
      CALL MEMMAN(KCOC,NELC*NCOC_TOT,'ADDL  ',2,'COC   ')
      CALL GET_CONF_FOR_OCCLS(ICEX_TP,NCOC_FSM,IBCOC_FSM,WORK(KCOC),
     &                        NSMST,2)
*
* Arc weights for addressing creation occupations
*
*. Memory for arc weights 
      CALL MEMMAN(KZC,2*NTOOB*NELC,'ADDL  ',2,'ZCconf')

*. Local scratch is needed for REO_CONFIGS, so
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'REO_C1')
*
      CALL MEMMAN(KLSCR,(NTOOB+1)*(NELC+1),'ADDL  ',2,'LSCR  ')
      CALL MEMMAN(KLOCMIN,NTOOB,'ADDL  ',2,'LOCMIN')
      CALL MEMMAN(KLOCMAX,NTOOB,'ADDL  ',2,'LOCMAX')
*. Min/Max occupation
      CALL MXMNOC_OCCLS(WORK(KLOCMIN),WORK(KLOCMAX),NGAS,NOBPT,
     &                  ICEX_TP,0,0)
*. and the arc weights
      CALL CONF_GRAPH(WORK(KLOCMIN),WORK(KLOCMAX),NTOOB,NELC,
     &     WORK(KZC),NCONFT,WORK(KLSCR))
*. And remove the local memory 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'REO_C1')
* Reorder array : Lexical to actual numbers 
      CALL MEMMAN(KCREO,NCOC_TOT,'ADDL  ',2,'COC_RE')
      CALL REO_CONFIGS(WORK(KCOC),NCOC_TOT,NELC,WORK(KZC),
     &                 NTOOB,WORK(KCREO),IBCOC_FSM)
*
*  ======================
*. Annihilation  strings 
*  ======================
*
*
*. Number per symmetry
      CALL GET_CONF_FOR_OCCLS(IAEX_TP,NAOC_FSM,IBAOC_FSM,IAOC,NSMST,
     &     1)
*. The actual occupations
      NAOC_TOT = IELSUM(NAOC_FSM,NSMST)
      NELA = IELSUM(IAEX_TP,NGAS)
      CALL MEMMAN(KAOC,NELA*NAOC_TOT,'ADDL  ',2,'AOC   ')
      CALL GET_CONF_FOR_OCCLS(IAEX_TP,NAOC_FSM,IBAOC_FSM,WORK(KAOC),
     &                        NSMST,2)
*
* Arc weights for addressing occupations
*
*. Memory for arc weights 
      CALL MEMMAN(KZA,2*NTOOB*NELA,'ADDL  ',2,'ZCconf')

*. Local scratch is needed for REO_CONFIGS, so
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'REO_C1')
*
      CALL MEMMAN(KLSCR,(NTOOB+1)*(NELA+1),'ADDL  ',2,'LSCR  ')
      CALL MEMMAN(KLOCMIN,NTOOB,'ADDL  ',1,'LOCMIN')
      CALL MEMMAN(KLOCMAX,NTOOB,'ADDL  ',1,'LOCMAX')
*. Min/Max occupation
      CALL MXMNOC_OCCLS(WORK(KLOCMIN),WORK(KLOCMAX),NGAS,NOBPT,
     &                  IAEX_TP,0,0)
*. and the arc weights
      CALL CONF_GRAPH(WORK(KLOCMIN),WORK(KLOCMAX),NTOOB,NELA,
     &     WORK(KZA),NCONFT,WORK(KLSCR))
*. And remove the local memory 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'REO_C1')
* Reorder array : Lexical to actual numbers 
      CALL MEMMAN(KAREO,NAOC_TOT,'ADDL  ',2,'COC_RE')
      CALL REO_CONFIGS(WORK(KAOC),NAOC_TOT,NELA,WORK(KZA),
     &                 NTOOB,WORK(KAREO),IBAOC_FSM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  Number of C occupations per symmetry '
        CALL IWRTMA(NCOC_FSM,1,NSMST,1,NSMST)
        WRITE(6,*) '  Number of A occupations per symmetry '
        CALL IWRTMA(NAOC_FSM,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
* 
      SUBROUTINE GET_CONF_FOR_OCCLS(IOC_TP,NOC_FSM,IBOC_FSM,IOC,
     &           NSMST,IWAY)
*
* Obtain the number of occupations and the actual occupations ( IWAY = 2)
* for given occupation type (IOC_TP)
*
* Jeppe Olsen, September 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Input : Number of electrons per GAS space
      INTEGER IOC_TP(NGAS)
*. Input  if IWAY = 2 , else output
*  Offset for occupations of given sym  
      INTEGER IBOC_FSM(NSMST)
*. Output : Number of occupations per symmetru
      INTEGER NOC_FSM(NSMST)
*. Output if IWAY = 2 : The actual occupations ordered by symmetry
      INTEGER IOC(*)
*. Scratch space 
      INTEGER ICONF(MXPNEL)
*
      NEL = IELSUM(IOC_TP,NGAS)
*
      IZERO = 0
      CALL ISETVC(NOC_FSM,IZERO,NSMST)
*. Loop over configurations
      INI = 1
      NONEW = 0
      NCONF_TEST = 0
 1000 CONTINUE
*. Next configuration 
C            NEXT_CONF_FOR_OCCLS(ICONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
        CALL NEXT_CONF_FOR_OCCLS(ICONF,IOC_TP,NGAS,NOBPT,INI,NONEW)
        INI = 0
        NCONF_TEST = NCONF_TEST + 1
C?      WRITE(6,*) ' Nonew = ', NONEW
C?      WRITE(6,*) ' Conf from NEXT_CONF = ' 
C?      CALL IWRTMA(ICONF,1,NEL,1,NEL)
*
C?      IF(NCONF_TEST.GE.100) THEN
C?         WRITE(6,*) ' Enforced stop in GET_CONF '
C?         STOP        ' Enforced stop in GET_CONF '
C?      END IF
*
        IF(NONEW.EQ.0) THEN
*. Another configuration has been delivered  
*. Find symmetry
          ISYM = ISYMST(ICONF,NEL)
          NOC_FSM(ISYM) = NOC_FSM(ISYM) + 1
          IF(IWAY.EQ.2) THEN
            NOC_TOT = IBOC_FSM(ISYM)-1 + NOC_FSM(ISYM)
            CALL ICOPVE(ICONF,IOC(1+(NOC_TOT-1)*NEL),NEL)
          END IF

      GOTO 1000
        END IF
*. Total number of configurations 
      NCONF_TOT = IELSUM(NOC_FSM,NSMST)
*. Offsets 
C          ZBASE(NVEC,IVEC,NCLASS)
      CALL ZBASE(NOC_FSM,IBOC_FSM,NSMST)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Occupation over gas spaces : '
         CALL IWRTMA(IOC_TP,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of configurations per symmetry '
         CALL IWRTMA(NOC_FSM,1,NSMST,1,NSMST)
*
         IF(IWAY.EQ.2) THEN
            WRITE(6,*) ' The actual configurations '
            CALL IWRTMA(IOC,NEL,NCONF_TOT,NEL,NCONF_TOT)
         END IF
      END IF
*
      RETURN
      END
      SUBROUTINE REO_CONFIGS(ICONF,NCONF,NEL,IZ,NORBT,IREO,IB_FSM)
*
* Obtain reorder array lexical order => actual order 
* for a set of configurations
*
* Offsets are defined with respect to start of symmetry
*
* Jeppe Olsen, Sept. 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*
*. Input
* =======
*
*. The occupation of configurations
      INTEGER ICONF(NEL,NCONF)
*. Arcweights
      INTEGER IZ(NORBT,NEL,2)
*. Offset for strings with given symmetry
      INTEGER IB_FSM(*)
*
*. Output
* =======
*
*. Reorder array lexical => actual order 
      DIMENSION IREO(*)
*. Local scratch : for configuration in truncated form 
      DIMENSION ICONF2(MXPORB)
*
C?    WRITE(6,*) ' In REO .. NORBT, NEL = ', NORBT, NEL
C?    WRITE(6,*) ' In REO, Number of configurations=', NCONF
      DO I = 1, NCONF
*. Obtain configuration in compact form -using negative numbers 
*. to flag double occupied orbitals 
C            REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY)
C?      WRITE(6,*) ' Config to be reordered ', 
C?   &  (ICONF(J,I),J=1,NEL)
        CALL REFORM_CONF_OCC(ICONF(1,I),ICONF2,NEL,NOCOB,1)
C               ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
C?      WRITE(6,*) ' NOCOB = ', NOCOB
*. Symmetry of this configuration 
        ISM = ISYMST(ICONF(1,I),NEL)
C?      WRITE(6,*) ' ISM = ', ISM
        ILEX =  ILEX_FOR_CONF(ICONF2,NOCOB,NORBT,NEL,IZ,0,IREO)
        IREO(ILEX) = I - IB_FSM(ISM) + 1
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reorder array, lexical => actual address '
        WRITE(6,*) ' Actual address is w.r.t. to start of block'
        CALL IWRTMA(IREO,1,NCONF,1,NCONF)
      END IF
*
      RETURN
      END
      SUBROUTINE IABS_TO_REL(IARRAY,NBLOCK,LBLOCK)
*
* An array IARRAY is given. Reform IARRAY, so each index 
* refers to start of block
*
      INCLUDE 'implicit.inc'
      INTEGER IARRAY(*), LBLOCK(NBLOCK)
*
      IOFF = 1
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFF = 1
        ELSE 
          IOFF = IOFF + LBLOCK(IBLOCK-1)
        END IF
        DO I = IOFF, IOFF + LBLOCK(IBLOCK-1)-1
           IARRAY(I) = IARRAY(I) - IOFF + 1
        END DO
      END DO
      NELMNT = IOFF + LBLOCK(NBLOCK)-1
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Array with relative indexing '
        WRITE(6,*) ' ============================ '
        CALL IWRTMA(IARRAY,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
      SUBROUTINE CAAB_TO_CA_OC(ISM,ISPOBEX_TP,IOBEX_TP,IOBEX_NUM,
     &           ISOX_FOR_OX,IBSOX_FOR_OX,NSOX_FOR_OX,
     &           IBSPOBEX,
     &           MX_ST_TSOSO_BLK_MX,NOP_CA,
     &           IZC, IZA, ICREO,IAREO,ICAOC,
     &           IBCA,NCOC_FSM,
     &           IBCAAB_FOR_CA,ICAAB_FOR_CA_OP,ICAAB_FOR_CA_NUM,
     &           LCAAB_FOR_CA,NCAAB_FOR_CA,
     &           NOBCONF,NSPOBOP,NCOMP_FOR_PROTO)


*
* Obtain the spinorbital excitations for each orbital excitation
*
*
* Jeppe Olsen, September 02
*. Modified to allow general prototypes, August 2004
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*
*  =====
*. Input
*  =====
*
*. The array of all spinorbital excitations
      INTEGER ISPOBEX_TP(4*NGAS,*)
*. All orbital orbital operators, the orbital excitation in action is 
*. IOBEX_NUM
      INTEGER IOBEX_TP(NGAS*2,*)
*. The arcweights for the C and A orbital occupations
      INTEGER IZC(*),IZA(*)
*. The reorder arrays for the C and A orbital occupations
      INTEGER ICREO(*), IAREO(*)
*. The occupation of the C and A orbital occupations
C     INTEGER ICOC(*), IAOC(*)
*. Offset to CA configurations with a given sym of C
      INTEGER IBCA(*)
*. Number of creation strings per symmetry
      INTEGER NCOC_FSM(*)
*. The list of orbital configurations 
       INTEGER ICAOC(NOP_CA,NOBCONF)
*. The spinorbital excitation types for a given orbital excitation type 
       INTEGER ISOX_FOR_OX(*)
*. The start of spinorbital excitations in ISOX_FOR_OX for 
*. a given orbital excitations
       INTEGER IBSOX_FOR_OX(*)
*. Number of spinorbital excitations for each orbital excitation 
       INTEGER NSOX_FOR_OX(*)
*.Base for coefficients for given spinorbital excitation type
       INTEGER IBSPOBEX(*)
*. Number of Components for the various prototype CA's
      INTEGER NCOMP_FOR_PROTO(*)

*
* =======
*. Output
* =======
*. The CAAB strings for a given CA configurations
*. ( LCAAB is the (max) number of elementary excitations in 
*    the CAAB operators)
      INTEGER ICAAB_FOR_CA_OP(NOP_CA,*)
*. The address in the spinorbital list for the CAABS belonging to a CAAB
      INTEGER ICAAB_FOR_CA_NUM(NSPOBOP)
*. The number of CAAB operators for each CA operators
      INTEGER NCAAB_FOR_CA(NOBCONF)
*. The number of operators in each of the CA CB AA AB operators 
      INTEGER LCAAB_FOR_CA(4,NSPOBOP)
*. The address of the first CAAB operator for a given CA operator
      INTEGER IBCAAB_FOR_CA(NOBCONF)
*. Offset in for the CAAB 
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GCC_FD')
*
      NTEST = 10
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' -------------------------------'
        WRITE(6,*) ' Information from CAAB_TO_CA_OC '
        WRITE(6,*) ' -------------------------------'
        WRITE(6,*)
        WRITE(6,*) ' CA => CAAB map for orbital excitation ', IOBEX_NUM
        WRITE(6,*) ' The corresponding CA operator '
        CALL IWRTMA(IOBEX_TP(1,IOBEX_NUM),NGAS,2,NGAS,2)
        WRITE(6,*) '  NOBCONF,NSPOBOP = ',  NOBCONF,NSPOBOP
        WRITE(6,*) ' NCOC_FSM(1) ', NCOC_FSM(1)
      END IF
*
*. Set up the the array IBCAAB_FOR_CA assuming that all 
*. spinorbital excitations belonging to a given orbital excitation 
* are given 
*. Number of operators in creation and annihilation part
      NOP_C = IELSUM(IOBEX_TP(1     ,IOBEX_NUM),NGAS)
      NOP_A = IELSUM(IOBEX_TP(1+NGAS,IOBEX_NUM),NGAS)
      NOP_CA = NOP_C + NOP_A
C?    WRITE(6,*) ' NOP_C, NOP_A = ', NOP_C, NOP_A
      IOFF = 1
      DO JOBEX = 1, NOBCONF
*. Obtain prototype for this CA ex
        IPROTO = IPROTO_TYPE_FOR_CA(ICAOC(1,JOBEX),IOBEX_NUM,
     &           NOP_C,NOP_A)
        NDET_FOR_CA = NCOMP_FOR_PROTO(IPROTO)
       IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Orbital excitation '
          CALL IWRTMA(ICAOC(1,JOBEX),1,NOP_CA,1,NOP_CA)
          WRITE(6,*) ' Prototype of orbexc ', IPROTO
          WRITE(6,*) ' Number of dets for conf ', NDET_FOR_CA
       END IF
       IBCAAB_FOR_CA(JOBEX) = IOFF
       IOFF = IOFF + NDET_FOR_CA
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IBCAAB_FOR_CA : '
        CALL IWRTMA(IBCAAB_FOR_CA,1, NOBCONF,1, NOBCONF)
      END IF 
*
      IZERO = 0
      CALL ISETVC(NCAAB_FOR_CA,IZERO,NOBCONF)
*. Four blocks of string occupations
      CALL MEMMAN(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
*
*. Loop over spinorbitaltypes for the given orbital excitations
      JSTART = IBSOX_FOR_OX(IOBEX_NUM)
      JSTOP  = JSTART + NSOX_FOR_OX(IOBEX_NUM) - 1 
      DO JJSPOBEX = JSTART, JSTOP
        JSPOBEX = ISOX_FOR_OX(JJSPOBEX)
C        WRITE(6,*) ' .. OCS will be called for JSPOBEX = ',
C    &   JSPOBEX
        JOFF = IBSPOBEX(JSPOBEX)
        CALL CAAB_TO_CA_OCS(ISPOBEX_TP(1,JSPOBEX),JOFF,1,NOP_CA,
     &     IZC,IZA,ICREO,IAREO,
     &     WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &     WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),IBCA,
     &     NCOC_FSM,IBCAAB_FOR_CA,ICAAB_FOR_CA_OP,ICAAB_FOR_CA_NUM,
     &     NCAAB_FOR_CA,LCAAB_FOR_CA)
      END DO
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Info on CAAB => CA relations '
         WRITE(6,*) ' ============================='
         WRITE(6,*)
         DO JOBCONF = 1, NOBCONF
            WRITE(6,*) ' CA conf ', JOBCONF, ' has ', 
     &      NCAAB_FOR_CA(JOBCONF), ' CAAB contributions '
            WRITE(6,*) ' Original order of the contributions '
            IOFF = IBCAAB_FOR_CA(JOBCONF)
            N = NCAAB_FOR_CA(JOBCONF)
            CALL IWRTMA(ICAAB_FOR_CA_NUM(IOFF),1,N,1,N)
         END DO
       END IF
          
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GCC_FD')
*
      RETURN
      END 
      SUBROUTINE CAAB_TO_CA_OCS(ITSS_TP,ITOFF,ISM,NOP_CA,
     &           IZC, IZA, ICREO,IAREO,
     &           IOCC_CA,IOCC_CB,IOCC_AA,IOCC_AB,IBCA,NCOC_FSM,
     &           IBCAAB_FOR_CA,ICAAB_FOR_CA_OP,ICAAB_FOR_CA_NUM,
     &           NCAAB_FOR_CA,LCAAB_FOR_CA)
*
* An  spin-orbital excitation type belonging to 
* a given orbital excitation type is given.
*
* ITOFF is offset for this type of spinorbital excitation 
*
* Obtain mapping Orbital excitation => spinorbital excitation 
*
* Jeppe Olsen, September 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER ITSS_TP(4*NGAS)
*. Arc weights for creation and annihilation occupations
      INTEGER IZC(*), IZA(*)
*. Reorder arrays for creation and annihilation occupations
      INTEGER ICREO(*),IAREO(*)
*, Number of creation occupations per symmetry
      INTEGER NCOC_FSM(*)
*. Offset of CA occupation with given symmetry of C string
      INTEGER IBCA(*)
*. First CAAB determinant for each CA operator
      INTEGER IBCAAB_FOR_CA(*)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*
      INTEGER IOCC_C(MXPNEL),IOCC_A(MXPNEL), IOCCX(MXPNEL)
      INTEGER IMS_C(MXPNEL),IMS_A(MXPNEL)
*. Output
*. Updated number of CAAB's for each CA 
      INTEGER NCAAB_FOR_CA(*)
*. Length of CA CB AA AB for each CAAB
      INTEGER LCAAB_FOR_CA(4,*)
*. The CA CB AA AB strings  
      INTEGER ICAAB_FOR_CA_OP(NOP_CA,*)
*. configuration => standard order of each SPOBEX 
      INTEGER ICAAB_FOR_CA_NUM(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ----------------------------'
        WRITE(6,*) ' Output from CAAB_TO_CA_OCS '
        WRITE(6,*) ' ----------------------------'
      END IF
C?     WRITE(6,*) ' The first two elements of IZC and IZA in CA_OCS'
C?     CALL IWRTMA(IZC,2,1,2,1)
C?     CALL IWRTMA(IZA,2,1,2,1)
      IT = ITOFF - 1
*. Transform from occupations to groups
      CALL OCC_TO_GRP(ITSS_TP(1+0*NGAS),IGRP_CA,1      )
      CALL OCC_TO_GRP(ITSS_TP(1+1*NGAS),IGRP_CB,1      )
      CALL OCC_TO_GRP(ITSS_TP(1+2*NGAS),IGRP_AA,1      )
      CALL OCC_TO_GRP(ITSS_TP(1+3*NGAS),IGRP_AB,1      )
*
      NEL_CA = IELSUM(ITSS_TP(1+0*NGAS),NGAS)
      NEL_CB = IELSUM(ITSS_TP(1+1*NGAS),NGAS)
      NEL_AA = IELSUM(ITSS_TP(1+2*NGAS),NGAS)
      NEL_AB = IELSUM(ITSS_TP(1+3*NGAS),NGAS)
      IF(NTEST.GE.100) THEN
      WRITE(6,'(A,4I4)') ' NEL_CA, NEL_CB, NEL_AA, NEL_AB = ',
     &                     NEL_CA, NEL_CB, NEL_AA, NEL_AB
      END IF
      DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        DO ISM_CA = 1, NSMST
          ISM_CB = MULTD2H(ISM_C,ISM_CA)
          DO ISM_AA = 1, NSMST
           ISM_AB =  MULTD2H(ISM_A,ISM_AA)
           IF(NTEST.GE.100) THEN
             WRITE(6,'(A,4I5)') ' ISM_CA, ISM_CB, ISM_AA, ISM_AB',
     &                            ISM_CA, ISM_CB, ISM_AA, ISM_AB
           END IF
*. obtain strings
           CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,NSTR_CA,
     &          IOCC_CA, NORBT,0,IDUM,IDUM)
           CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,NSTR_CB,
     &          IOCC_CB, NORBT,0,IDUM,IDUM)
           CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,NSTR_AA,
     &          IOCC_AA, NORBT,0,IDUM,IDUM)
           CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,NSTR_AB,
     &          IOCC_AB, NORBT,0,IDUM,IDUM)
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                             NORBT,IDOREO,IZ,IREO)
*. Loop over T elements as  matric T(I_CA, I_CB, IAA, I_AB)
            DO I_AB = 1, NSTR_AB
             DO I_AA = 1, NSTR_AA
              DO I_CB = 1, NSTR_CB
               DO I_CA = 1, NSTR_CA
                IT = IT + 1
                IF(NTEST.GE.100) THEN
                WRITE(6,*) ' CA CB  strings '
                  CALL IWRTMA(IOCC_CA(1+(I_CA-1)*NEL_CA),
     &                  1,NEL_CA,1,NEL_CA)
                  CALL IWRTMA(IOCC_CB(1+(I_CB-1)*NEL_CB),
     &                  1, NEL_CB,1,NEL_CB)
                END IF 
*
 
* Adress of Combined creation string in list of creation occupations
*
*. Obtain the AB occuption in IOCC_C
C               ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,IDET_OC,IDET_SP,ISIGN)
                CALL ABSTR_TO_ORDSTR(
     &          IOCC_CA(1+(I_CA-1)*NEL_CA),IOCC_CB(1+(I_CB-1)*NEL_CB),
     &          NEL_CA, NEL_CB, IOCC_C,IMS_C,ISIGN_C)
*. Reform Occupation to compressed form 
                NEL_C = NEL_CA + NEL_CB
C                    REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY)
                CALL REFORM_CONF_OCC(IOCC_C,IOCCX,NEL_C,NOCOBX,1)
*. Address of C string
C                        ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,
C                                      IDOREO,IREO)
C?              WRITE(6,*) ' Lexical adress for C '
                IC_NUM = ILEX_FOR_CONF(IOCCX,NOCOBX,NTOOB,NEL_C,IZC,
     &                   1, ICREO)
*
* Adress of Combined annihilation  string in list of creation occupations
*
*. Obtain the AB occuption in IOCC_A
C               ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,IDET_OC,IDET_SP,ISIGN)
                CALL ABSTR_TO_ORDSTR(
     &          IOCC_AA(1+(I_AA-1)*NEL_AA),IOCC_AB(1+(I_AB-1)*NEL_AB),
     &          NEL_AA, NEL_AB, IOCC_A,IMS_A,ISIGN_A)
*. Reform Occupation to compressed form 
                NEL_A = NEL_AA + NEL_AB
C                    REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY)
                CALL REFORM_CONF_OCC(IOCC_A,IOCCX,NEL_A,NOCOBX,1)
*. Address of A occupation 
C                        ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,
C                                      IDOREO,IREO)
C?              WRITE(6,*) ' Lexical adress for A '
                IA_NUM = ILEX_FOR_CONF(IOCCX,NOCOBX,NTOOB,NEL_A,IZA,
     &                   1, IAREO)
                IF(NTEST.GE.100) THEN
                  WRITE(6,'(A,4I4)') ' I_AB, I_AA, I_CB, I_CA', 
     &                                 I_AB, I_AA, I_CB, I_CA 
                END IF
*. And adress of the corresponding CA string 
                ICA_ADR = IBCA(ISM_C) - 1 
     &                  + (IA_NUM-1)*NCOC_FSM(ISM_C) + IC_NUM
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' IBCA(ISM_C) = ', IBCA(ISM_C)
                  WRITE(6,*) ' NCOC_FSM(ISM_C) = ',NCOC_FSM(ISM_C) 
                  WRITE(6,*) ' IA_NUM, IC_NUM, ISM_C, ICA_ADR = ',
     &                         IA_NUM, IC_NUM, ISM_C, ICA_ADR
                END IF
C       STOP ' Jeppe Stop '
*. And enroll this spinorbital excitation in the list for orbital 
*. excitation ICA_ADR
                NCAAB_FOR_CA(ICA_ADR) = NCAAB_FOR_CA(ICA_ADR) + 1
                ICAAB_ADR = IBCAAB_FOR_CA(ICA_ADR)-1
     &                    +  NCAAB_FOR_CA(ICA_ADR)
                IF(NTEST.GE.100) THEN
                WRITE(6,*) ' IBCAAB_FOR_CA(ICA_ADR) = ',
     &                       IBCAAB_FOR_CA(ICA_ADR)
                WRITE(6,*) ' NCAAB_FOR_CA(ICA_ADR) ',
     &                       NCAAB_FOR_CA(ICA_ADR)
                WRITE(6,*) ' ICAAB_ADR = ', ICAAB_ADR
                END IF
                ICAAB_FOR_CA_NUM(ICAAB_ADR) = IT
                IPLACE = 1
                CALL ICOPVE(IOCC_CA(1+(I_CA-1)*NEL_CA), 
     &                      ICAAB_FOR_CA_OP(IPLACE,ICAAB_ADR),NEL_CA)
                IPLACE = IPLACE + NEL_CA
                CALL ICOPVE(IOCC_CB(1+(I_CB-1)*NEL_CB), 
     &                      ICAAB_FOR_CA_OP(IPLACE,ICAAB_ADR),NEL_CB)
                IPLACE = IPLACE + NEL_CB
                CALL ICOPVE(IOCC_AA(1+(I_AA-1)*NEL_AA), 
     &                      ICAAB_FOR_CA_OP(IPLACE,ICAAB_ADR),NEL_AA)
                IPLACE = IPLACE + NEL_AA
                CALL ICOPVE(IOCC_AB(1+(I_AB-1)*NEL_AB), 
     &                      ICAAB_FOR_CA_OP(IPLACE,ICAAB_ADR),NEL_AB)
*
                LCAAB_FOR_CA(1,ICAAB_ADR) = NEL_CA
                LCAAB_FOR_CA(2,ICAAB_ADR) = NEL_CB
                LCAAB_FOR_CA(3,ICAAB_ADR) = NEL_AA
                LCAAB_FOR_CA(4,ICAAB_ADR) = NEL_AB
               END DO
              END DO
             END DO
            END DO
*           ^ End of loop over elements of block
           END DO
*          ^ End of loop over ISM_AA
        END DO
*       ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Number of elements ', IT-ITOFF + 1
      END IF
*
      RETURN
      END
      FUNCTION IGATSUM(IVEC,IGAT,IOFF,NELMNT)
*
* IGATSUM = SUM(I=IOFF,IOFF-1+NELMNT) IVEC(IGAT(I))
*
      INCLUDE 'implicit.inc'
*
      INTEGER IVEC(*),IGAT(*)
*
      ISUM = 0
      DO I = IOFF, IOFF-1+NELMNT
        ISUM = ISUM + IVEC(IGAT(I))
      END DO
*
      IGATSUM = ISUM
*
      RETURN
      END
      SUBROUTINE WRITE_CAAB_CONFM
*
* Print the spinorbital excitations as obtained from configurations 
* order
*
*
* Jeppe Olsen, September 2002
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'ctcc.inc'
*
*. Loop over the various types of orbital excitations
      DO IOBEX_TP = 1, NOBEX_TP
*. And let another routine do the work for a given 
*. orbital excitation type
        CALL WRITE_CAAB_CONF(
     &       NCAOC(IOBEX_TP),WORK(KIBCAAB_FOR_CA(IOBEX_TP)),
     &       WORK(KICAAB_FOR_CA_OP(IOBEX_TP)),
     &       WORK(KICAAB_FOR_CA_NUM(IOBEX_TP)),
     &       WORK(KLCAAB_FOR_CA(IOBEX_TP)),
     &       WORK(KNCAAB_FOR_CA(IOBEX_TP))                            )
      END DO
*
      RETURN
      END 
      SUBROUTINE WRITE_CAAB_CONF(NCAOC,
     &           IBCAAB_FOR_CA,ICAAB_FOR_CA_OP,ICAAB_FOR_CA_NUM,
     &           LCAAB_FOR_CA,NCAAB_FOR_CA)
*
* Print spinorbital excitations from configuration information 
*
*
* Jeppe Olsen, September 02
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*
*  =====
*. Input
*  =====
*
*. The spinorbital excitations (CAABS) belonging to a CA
      INTEGER ICAAB_FOR_CA_OP(*)
*. The address in the spinorbital list for the CAABS belonging to a CA
      INTEGER ICAAB_FOR_CA_NUM(*)
*. The number of CAAB operators for each CA operators
      INTEGER NCAAB_FOR_CA(*)
*. The number of operators in each of the CA CB AA AB operators 
      INTEGER LCAAB_FOR_CA(4,*)
*. The address of the first CAAB operator for a given CA operator
      INTEGER IBCAAB_FOR_CA(*)
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'WCAABC')
*
      WRITE(6,*) ' Number of configurations for orbex-type', NCAOC
*
      DO ICA = 1, NCAOC
        IBCA = IBCAAB_FOR_CA(ICA)
        NCAAB = NCAAB_FOR_CA(ICA)
*
        DO ICAAB = 1, NCAAB
          LCA = LCAAB_FOR_CA(1,IBCA-1+ICAAB)
          LCB = LCAAB_FOR_CA(2,IBCA-1+ICAAB)
          LAA = LCAAB_FOR_CA(3,IBCA-1+ICAAB)
          LAB = LCAAB_FOR_CA(4,IBCA-1+ICAAB)
          LCAAB = LCA + LCB + LAA + LAB
          ICAAB_ABS = IBCA-1+ICAAB
*
          WRITE(6,*) ' Info for CA configuration and component = ',
     &                 ICA, ICAAB
          WRITE(6,*) ' LCA LCB LAA LAB = ', LCA, LCB, LAA, LAB
          WRITE(6,*) ' The corresponding CA CB AA AB strings '
          CALL IWRTMA(ICAAB_FOR_CA_OP(1+(ICAAB_ABS-1)*LCAAB),
     &                1,LCAAB,1,LCAAB)
        END DO
      END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'WCAABC')
      RETURN
      END 
      SUBROUTINE GEN_REORDER_CAABM(ICAAB_REO)
*
* Outer routine for 
* Generating reorder array going from configuration order of 
* CAAB to standard CAAB order. The array is delivered in 
* ICAAB_REO, which should be located outside
*
* This routine exploits that ICAAB_FOR_CA_NUM gives reordering
* within a given CA type
*
* Jeppe Olsen, September 2002 ( 20 hours to take off to UTRECHT)
*
* This routine collects the informations stored 
* seperately for each occupation type in a single array
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'ctcc.inc' 
      INCLUDE 'crun.inc'
*. Output
      INTEGER ICAAB_REO(*)
*
      NTEST = 00
      IONEM = -1
      CALL ISETVC(ICAAB_REO,IONEM,N_CC_AMP)
*. Loop over the various types of orbital excitations
      IBCONF = 1
      DO IOBEX_TP = 1, NOBEX_TP
*. The number of CAABs for a given Orbital excitation
        NSOX = IFRMR(WORK(KNSOX_FOR_OX),1,IOBEX_TP)
        IBSOX = IFRMR(WORK(KIBSOX_FOR_OX),1,IOBEX_TP)
        NCAAB = IGATSUM(WORK(KLLSOBEX),WORK(KISOX_FOR_OX),
     &                     IBSOX,NSOX)
        IF(NTEST.GE.1000) THEN 
          WRITE(6,*) ' IOBEX_TP, NSOX, IBSOX, NCAAB, IBCONF = ',
     &                 IOBEX_TP, NSOX, IBSOX, NCAAB, IBCONF
        END IF
        CALL ICOPVE(WORK(KICAAB_FOR_CA_NUM(IOBEX_TP)),ICAAB_REO(IBCONF),
     &              NCAAB)
*
        IBCONF = IBCONF + NCAAB
      END DO
*
      I_DO_CHECK = 0
      IF( I_DO_CHECK.EQ.1) THEN
*. Check that the sum of all reorder elements = N_CC_AMP*(N_CC_AMP+1)/2
        ICHECKSUM = IELSUM(ICAAB_REO,N_CC_AMP)
        IF(ICHECKSUM.NE.N_CC_AMP*(N_CC_AMP+1)/2) THEN
          WRITE(6,*) '  CHECKSUM in REO failed ... '
          WRITE(6,*) ' Reorder array for CAAB, CONF => CAAB order '
          WRITE(6,*) ' =========================================== '
          CALL IWRTMA(ICAAB_REO,1,N_CC_AMP,1,N_CC_AMP)
          STOP ' CHECKSUM in REO failed ... '
        ELSE 
          WRITE(6,*) ' Check sum passed '
        END IF
      END IF
        
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Reorder array for CAAB, CONF => CAAB order '
         WRITE(6,*) ' =========================================== '
         CALL IWRTMA(ICAAB_REO,1,N_CC_AMP,1,N_CC_AMP)
      END IF
*
      RETURN
      END 
      SUBROUTINE PROTO_SPIN_MAT
*
* Set up matrices transforming between CAAB and spinadapted  operator 
* basis. Quick fix for results for the utrecht meeting
*
*. Jeppe Olsen, September 2002
*. 
*. Modified to include 4 det case, August 2004
*
      INCLUDE 'wrkspc.inc'
*. Output
      COMMON/PROTO_SP_MAT/NSPA_FOP(6),NCAAB_FOP(6),IB_FOP(6),XTRA(100),
     &                    NSPA_FOP_G(6,MXPCYC),NCAAB_FOP_G(6,MXPCYC),
     &                    IB_FOP_G(6,MXPCYC)
*
      FACTOR = 1.0D0/DSQRT(2.0D0)
*
* For one component  : Type 1
      NSPA_FOP(1) = 1
      NCAAB_FOP(1) = 1
      IB_FOP(1) = 1
      XTRA(1) = 1.0D0
*
*. For two components : type 2
*
      
      NSPA_FOP(2) = 1
      NCAAB_FOP(2) = 2
      IB_FOP(2) = 2
      XTRA(2) = FACTOR 
      XTRA(3) = FACTOR 
*
* for four components : type 4
*
      NSPA_FOP(4) = 2
      NCAAB_FOP(4) = 4
      IB_FOP(4) = 4
      ZERO = 0.0D0
      CALL SETVEC(XTRA(IB_FOP(4)),ZERO, NSPA_FOP(4)* NCAAB_FOP(4))
*. CAAB's related by time reversal are ( I hope ...)
*. 1 and 4
*. 2 and 3
*.. 1 : 1 + 4
      XTRA(IB_FOP(4)-1+1+(1-1)*4 ) = FACTOR
      XTRA(IB_FOP(4)-1+4+(1-1)*4 ) = FACTOR
*.. 2 : 2 + 3
      XTRA(IB_FOP(4)-1+2+(2-1)*4 ) = FACTOR
      XTRA(IB_FOP(4)-1+3+(2-1)*4 ) = FACTOR
*
* For six components : type 6
*
      NSPA_FOP(6) = 3
      NCAAB_FOP(6) = 6
      IB_FOP(6) = 12
      ZERO = 0.0D0
      CALL SETVEC(XTRA(IB_FOP(6)),ZERO, NSPA_FOP(6)* NCAAB_FOP(6))
*. CAAB's related by time reversal are ( I hope ...)
*. 1 and 4
*. 2 and 3
*. 5 and 6
*.. 1 : 1 + 4
      XTRA(IB_FOP(6)-1+1+(1-1)*6 ) = FACTOR
      XTRA(IB_FOP(6)-1+4+(1-1)*6 ) = FACTOR
*.. 2 : 2 + 3
      XTRA(IB_FOP(6)-1+2+(2-1)*6 ) = FACTOR
      XTRA(IB_FOP(6)-1+3+(2-1)*6 ) = FACTOR
*.. 3 : 5 + 6
      XTRA(IB_FOP(6)-1+5+(3-1)*6 ) = FACTOR
      XTRA(IB_FOP(6)-1+6+(3-1)*6 ) = FACTOR
*
      RETURN
      END 
      SUBROUTINE REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
*
* Transform vector between CAAB form and spinadapted form 
*
* IWAY = 1 : CAAB => Spin adapted form 
* IWAY = 2 : Spin adapted form => CAAB
*
* Jeppe Olsen, September 2002
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'crun.inc'
*. Input and output
      DIMENSION  VEC_CAAB(*),VEC_SP(*)
*. and a scratch vector 
      DIMENSION VEC_SCR(*)
*
      NTEST = 000
        IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' REF_CCV_CCAB speaking'
         WRITE(6,*) ' ---------------------'
         WRITE(6,*)
         IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' CAAB => spinadapted basis transformation '
         ELSE
          WRITE(6,*) ' spinadapted basis => CAABtransformation '
         END IF
        END IF
*
      IF(IWAY.EQ.1) THEN
* CAAB => Spin adapted : Reorder to conf and then transform
        CALL GATVEC(VEC_SCR,VEC_CAAB,WORK(KLREORDER_CAAB),
     &              N_CC_AMP)
*
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Result from GATVEC '
          CALL WRTMAT(VEC_SCR,1,N_CC_AMP,1,N_CC_AMP)
        END IF
*. Offsets for CAAB and Spin adapted form will be updated in the process
        IB_CAAB = 1
        IB_SP = 1
        DO JOBTP = 1, NOBEX_TP
           CALL CAAB_SP_FOR_OCTP(VEC_SCR(IB_CAAB),VEC_SP(IB_SP),
     &                     WORK(KNCAAB_FOR_CA(JOBTP)),NCAOC(JOBTP),
     &                     N_SP,N_CAAB,1)
           IB_CAAB = IB_CAAB + N_CAAB
           IB_SP   = IB_SP   + N_SP
        END DO 
      ELSE 
*. Spin-adapted => CAAB transformation 
        IB_CAAB = 1
        IB_SP = 1
        DO JOBTP = 1, NOBEX_TP
C?         WRITE(6,*) ' REF_CCV : JOBTP = ', JOBTP
           CALL CAAB_SP_FOR_OCTP(VEC_SCR(IB_CAAB),VEC_SP(IB_SP),
     &                     WORK(KNCAAB_FOR_CA(JOBTP)),NCAOC(JOBTP),
     &                     N_SP,N_CAAB,2)
           IB_CAAB = IB_CAAB + N_CAAB
           IB_SP   = IB_SP   + N_SP
        END DO 
C SCAVEC(VECO,VECI,INDEX,NDIM)
        CALL SCAVEC(VEC_CAAB,VEC_SCR,WORK(KLREORDER_CAAB),N_CC_AMP)
      END IF
      N_CAAB_TOT = IB_CAAB - 1
      N_SP_TOT   = IB_SP   - 1
*
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' Test, N_CAAB_TOT, N_SP_TOT = ', 
     &                   N_CAAB_TOT, N_SP_TOT
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector in spinadapted basis '
        CALL WRTMAT(VEC_SP,1,N_SP_TOT,1,N_SP_TOT)
        WRITE(6,*) ' Vector in CAAB basis '
        CALL WRTMAT(VEC_CAAB,1,N_CAAB_TOT,1,N_CAAB_TOT)
      END IF
*
      RETURN
      END 
      SUBROUTINE CAAB_SP_FOR_OCTP(VEC_CAAB,VEC_SP,NCAAB_FOR_CA,
     &                             NCONF,N_SP,N_CAAB,IWAY )
*
* Transforming between spinadapted  and CAAB form of 
* vector for given OCTP
*
* IWAY = 1 : CAAB => Spin
* IWAY = 2 : Spin => CAAB
*
* Jeppe Olsen, September 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'proto_sp_mat.inc'
*, Input or output
      DIMENSION VEC_CAAB(*),VEC_SP(*)
*. Number of dets per configuration
      INTEGER NCAAB_FOR_CA(NCONF)
      
*
      NTEST = 00
*
      IB_SP = 1
      IB_CAAB = 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCONF, IWAY  = ', NCONF, IWAY
      END IF
      DO ICONF = 1, NCONF
*. Use number of determinants is used to decide the type of open shells 
*. ( Yes dirty initial version)
        NDET = NCAAB_FOR_CA(ICONF)
        NCSF = NSPA_FOP(NDET)
        IB   = IB_FOP(NDET)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' ICONF, NDET ,NCSF = ',ICONF, NDET ,NCSF
        END IF
        IF(IWAY.EQ.1) THEN
*VEC_CSF(I) = SUM(J) XTRA(J,I) VEC_DET(J) 
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS
C?        WRITE(6,*) ' IB, IB_CAAB, IB_SP = ', IB,IB_CAAB,IB_SP
          CALL MATVCC(XTRA(IB),VEC_CAAB(IB_CAAB),VEC_SP(IB_SP),
     &                NDET,NCSF,1)
C?        WRITE(6,*) ' XTRA, VEC_CAAB, VEC_SP : '
C?        CALL WRTMAT(XTRA(IB),NDET,NCSF,NDET,NCSF)
C?        CALL WRTMAT(VEC_CAAB(IB_CAAB),1,NDET,1,NDET)
C?        CALL WRTMAT(VEC_SP(IB_SP),1,NCSF,1,NCSF)
        ELSE
* VEC_DET(J) = SUM(I) XTRA(J,I) VEC_CSF(I)
          CALL MATVCC(XTRA(IB),VEC_SP(IB_SP),VEC_CAAB(IB_CAAB),
     &                NDET,NCSF,0)
        END IF
        IB_SP = IB_SP + NCSF
        IB_CAAB = IB_CAAB + NDET
      END DO
*. Length of SP and CAAB expansions should be returned so 
         N_SP = IB_SP - 1
         N_CAAB =IB_CAAB-1
*
      RETURN
      END
      SUBROUTINE NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
*
* Number of CSF's in current expansion obtained by reading 
* number of CAABs in the CA expansion 
*
* Jeppe Olsen, Amsterdam airport Sept 20, 2002
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'cprnt.inc'
*. Local scratch
      INTEGER NCNF_FOP(MXPNEL), ISCR(MXPNEL)
*. Find number of configurations with the various number of open shells
*. at the moment I am here assuming atmost 4 open shells..
*. At the moment I assume only combinations so 3 csfs for 4 open shells..
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRCSF)
      MAXNDET = 6
      IZERO = 0
      CALL ISETVC(ISCR,IZERO,MAXNDET)
      CALL ISETVC(NCNF_FOP,IZERO,MAXNDET)
*
C?    WRITE(6,*) ' Number of orbitalexcitationtypes ', NOBEX_TP
      DO IOBEX_TP = 1, NOBEX_TP
*. Count the number of times the various number of dets for 
*. a given CA occurs
*  COUNT_OCCURENCE(IVEC,IOCC,NELMNT,MAXVAL)
        NCA = NCAOC(IOBEX_TP)
C?      WRITE(6,*) ' IOBEX_TP, NCA ', IOBEX_TP, NCA
C?      WRITE(6,*) ' And the types '
C?      CALL IWRTMA(WORK(KNCAAB_FOR_CA(IOBEX_TP)),1,NCA,1,NCA)
        CALL COUNT_OCCURENCE(WORK(KNCAAB_FOR_CA(IOBEX_TP)),ISCR,NCA,
     &                       MAXNDET)
        IONE = 1
        CALL IVCSUM(NCNF_FOP,NCNF_FOP,ISCR,IONE,IONE,MAXNDET)
      END DO
*
      NSPA = NCNF_FOP(1)*1 + NCNF_FOP(2)*1 + NCNF_FOP(4)*2 
     &     + NCNF_FOP(6)*3
      NCAAB= NCNF_FOP(1)*1 + NCNF_FOP(2)*2 + NCNF_FOP(4)*4
     &     + NCNF_FOP(6)*6
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Number of CA ops with 1 comp = ',
     &  NCNF_FOP(1)
        WRITE(6,*) ' Number of CA ops with two comps = ',
     &  NCNF_FOP(2)
        WRITE(6,*) ' Number of CA ops with four comps = ',
     &  NCNF_FOP(4)
        WRITE(6,*) ' Number of CA ops with six comps = ',
     &  NCNF_FOP(6)
        WRITE(6,*) ' Number of spinadapted operators = ', NSPA
        WRITE(6,*) ' Number of CAAB                  = ', NCAAB
      END IF
*
      RETURN
      END 
      SUBROUTINE ICCC_COMPLETE_MAT(
     &        IREFSPC,ITREFSPC,I_SPIN_ADAPT,
     &        IROOT,T_EXT,C_0,INI_IT,IFIN_IT,VEC1,VEC2,IDIIS)

*
* Master routine for Internal Contraction Coupled Cluster 
* with complete incore * construction of all matrices.
*
* It is assumed that the excitation manifold produces 
* states that are orthogonal to the reference so 
* no projection is carried out
*
* Routine is allowed to leave without turning the lights off,
* i.e. leave routine with all allocations and marks intact.
*: Thus : Allocations are only done if INI_IT = 1
*          Deallocations are only done if IFIN_IT = 1
*
* IF IDIIS.NE.0, DIIS is used to accelerate convergence 
*
* Jeppe Olsen, Aug. 2005 
*
*. for DIIS units LUSC35 and LUSC36 will be used for storing vectors
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
*
      REAL*8
     &INPROD
*. Output : Coefficients of external correlation 
      DIMENSION T_EXT(*)
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
      COMMON/CLOCAL/KVEC1,KVEC2,MXCJ,
     & KLVCC1,KLVCC2,KLVCC3,KLVCC4,KLVCC5,KLSMAT,KLXMAT,KLJMAT,KLU,KLL,
     & NSING,NNONSING,KLCDIIS,KLDIA
*. Scratch for CI behind the curtain 
       DIMENSION VEC1(*),VEC2(*)
       WRITE(6,*) ' Code has should be modified to new MRCC vecfnc '
       STOP ' Code has should be modified to new MRCC vecfnc '

*. Number of Spin adapted functions ( and NCAAB for a check)
      CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
*. We will not include the unit-operator so 
      NSPAM1 = NSPA - 1
*
      NTEST = 10
      WRITE(6,*) 
      WRITE(6,*) ' Complete J matrix will be used '
      WRITE(6,*) ' ==============================='
      WRITE(6,*)
      WRITE(6,*) ' Reference space is ', IREFSPC
      WRITE(6,*) ' Space of Operators times reference space ', ITREFSPC
      WRITE(6,*)
      WRITE(6,*) ' Number of parameters in spinuncoupled basis ', 
     &           N_CC_AMP
      WRITE(6,*) ' Number of parameters in spincoupled   basis ', 
     &           NSPA
      WRITE(6,*) ' INI_IT, IFIN_IT = ', INI_IT, IFIN_IT
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial T-amplitudes '
        CALL WRTMAT(T_EXT,1,N_CC_AMP,1,N_CC_AMP)
      END IF
*. Allowed number of iterations
      NNEW_MAX = 15
      MAXITL = NNEW_MAX
*
      IF(INI_IT.EQ.1) 
     &CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICC_CMP')
*
* Space for complete J matrices 
*
*. And a few working vectors 
      IF(INI_IT.EQ.1) THEN
*. Space for old fashioned CI behind the curtain
COLD    CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
COLD    KVEC1P = KVEC1
COLD    KVEC2P = KVEC2
        CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'VCC1  ')
        CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'VCC2  ')
        CALL MEMMAN(KLVCC3,N_CC_AMP,'ADDL  ',2,'VCC3  ')
        CALL MEMMAN(KLVCC4,N_CC_AMP,'ADDL  ',2,'VCC4  ')
        CALL MEMMAN(KLVCC5,N_CC_AMP,'ADDL  ',2,'VCC5  ')
        CALL MEMMAN(KLVCC6,N_CC_AMP,'ADDL  ',2,'VCC6  ')
*. For complete matrices, three used pt
        LEN = NSPA**2
        CALL MEMMAN(KLSMAT,LEN,'ADDL  ',2,'SMAT  ')
        CALL MEMMAN(KLXMAT,LEN,'ADDL  ',2,'XMAT  ')
        CALL MEMMAN(KLJMAT,LEN,'ADDL  ',2,'JMAT  ')
*. Storage for LU decomposition of J
        LEN = NSPA*(NSPA+1)/2
        CALL MEMMAN(KLL,LEN,'ADDL  ',2,'L     ')
        CALL MEMMAN(KLU,LEN,'ADDL  ',2,'U     ')
*. Space for DIIS
        CALL MEMMAN(KLCDIIS,MAXITL,'ADDL ',2,'CDIIS ')
      END IF
*
*. Identify the unit  operator i.e. the operator with 
*. zero creation and annihilation operators
      IDOPROJ = 0
*. Construct metric (once again ..)
*. Prepare the routines used in COM_SH
*. Not used here
      C_0X = 0.0D0
      KLTOPX = -1
*. Used
      NREFX = N_REF
      IREFSPCX = IREFSPC
      ITREFSPCX = ITREFSPC
      NCAABX = N_CC_AMP
      NSPAX = NSPA
      IPROJSPCX = IREFSPC
*. Unitoperator in SPA order ... Please check ..
      IUNIOPX = 0
*. Metric only evaluated in first macro-it
      IF(INI_IT.EQ.1) THEN
       CALL COM_SH(WORK(KLSMAT),WORK(KLSMAT),WORK(KLVCC1),WORK(KLVCC2),
     &             WORK(KLVCC3),VEC1,VEC2,
     &             N_CC_AMP,IREFSPC,ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &             IDOPROJ,IUNIOP,1,0,1,I_DO_EI,NSPA,0,0,0)
*. ELiminate part referring to unit operator
       CALL TRUNC_MAT(WORK(KLSMAT),NSPA,NSPA,NSPAM1,NSPAM1)
C      GET_ON_BASIS(S,NVEC,NSING,X,SCRVEC1,SCRVEC2)
       CALL GET_ON_BASIS(WORK(KLSMAT),NSPAM1,NSING,
     &                  WORK(KLXMAT),WORK(KLVCC1),WORK(KLVCC2) )
       WRITE(6,*) ' Number of singularities in S ', NSING
       NNONSING = NSPAM1 - NSING
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Transformation matrix to nonsingular basis '
         CALL WRTMAT(WORK(KLXMAT),NSPAM1,NNONSING,NSPAM1,
     &              NNONSING)
       END IF
      END IF
*     ^ End if it was initial iteration 
      IF(IDIIS.NE.0) THEN
        CALL REWINO(LUSC35)
        CALL REWINO(LUSC36)
      END IF
*. Loop over Newton iterations 
      DO IT = 1, NNEW_MAX
*. Construct CC vector function  in VCC5 
C?      WRITE(6,*) ' MRCC vector function at current point '
        CALL MRCC_VECFNC(WORK(KLVCC5),T_EXT,NCOMMU_V,I_APPROX_HCOM_V,
     &                   IREFSPC,ITREFSPC) 
*. The energy is returned as first element in CAAB basis, so
        E = WORK(KLVCC5)
*. And set energy term to zero
        WORK(KLVCC5) = 0.0D0
        VCFNORM = SQRT(INPROD(WORK(KLVCC5+1),WORK(KLVCC5+1),NCAAB-1))
        WRITE(6,'(A,1X,I4,2E22.15)')
     &  ' It, vecfnc : energy and norm ', IT, E, VCFNORM 
*
C       MRCC_VECFNC(CCVECFNC,T,NCOMMU,IREFSPC,ITREFSPC)
*. Vectors are stored in CAAB basis - not the smartest..
        IF(IDIIS.EQ.1) THEN
*. It is assumed that DIIS leaved the file at end of file 
*. T_ext on LUSC35, VECFNC on LUSC36
          CALL VEC_TO_DISC(T_EXT,NCAAB,0,-1,LUSC35)
          CALL VEC_TO_DISC(WORK(KLVCC5),NCAAB,0,-1,LUSC36)
*. We have now IT vectors in LUSC36, find combination with lowest 
*. Norm 
C DIIS_SIMPLE(LUEVEC,NVEC,NDIM,C)
          CALL DIIS_SIMPLE(LUSC36,IT,NCAAB,WORK(KLCDIIS))
*. Obtain combination as given in CDIIS
C  MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
          CALL MVCSMD(LUSC35,WORK(KLCDIIS),LUSC37,LUSC38,
     &                WORK(KLVCC1),WORK(KLVCC2),IT,1,-1)
          CALL VEC_FROM_DISC(T_EXT,NCAAB,1,-1,LUSC37)
*. Calculate new vectorfunction for T  or use sum
          I_NEW_OR_SUM = 1
          IF(I_NEW_OR_SUM.EQ.1) THEN
            WRITE(6,*) ' CC vector-function recalculated after DIIS '
            CALL MRCC_VECFNC(WORK(KLVCC5),T_EXT,NCOMMU_V,IREFSPC,
     &           ITREFSPC)
*. Note : I am not storing new vectors in DIIS queue - 
*         to have symmetry between case where vecfunc is 
*         obtained from sum.
            E = WORK(KLVCC5)
            VCFNORM = SQRT(INPROD(WORK(KLVCC5+1),WORK(KLVCC5+1),
     &                NCAAB-1))
            WRITE(6,'(A,I4,2E22.15)')
     &      ' From DIIS : It, vecfnc : energy and norm ',
     &        IT, E, VCFNORM 
          ELSE
            CALL MVCSMD(LUSC36,WORK(KLCDIIS),LUSC37,LUSC38,
     &                  WORK(KLVCC1),WORK(KLVCC2),IT,1,-1)
            CALL VEC_FROM_DISC(WORK(KLVCC5),NCAAB,1,-1,LUSC37)
            VCFNORM = SQRT(INPROD(WORK(KLVCC5+1),WORK(KLVCC5+1),
     &                NCAAB-1))
            WRITE(6,'(A,I4,2E22.15)')
     &      ' From DIIS : It, norm of approx vecfnc  ',
     &        IT,  VCFNORM 
          END IF
*.        ^ End if VECFNC should be recalculated or obtained as sum
        END IF
*. Transform to SPA basis
        CALL REF_CCV_CAAB_SP(WORK(KLVCC5),WORK(KLVCC1),WORK(KLVCC2),1)
C            REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
*. and to orthonormal basis, save in VCC5
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
        CALL MATVCC(WORK(KLXMAT),WORK(KLVCC1),WORK(KLVCC5),NSPAM1,
     &              NNONSING,1)
*. Transform to Nonsigular basis 
*. Construct Jacobian matrix in nonsingular basis
*. Here : Evaluate Jacobian in first IT, and use fewer commutators
*
* A further simplification is possible. If - As pt only one 
* commutator is used, one can restrict the space to be the MRSD space
* instead of the presently used MRSDTQ space. To accomplish this
* add the MRSD space as  third space after the refspc  and ITREFSPC
        IF(INI_IT.EQ.1.AND.IT.EQ.1) THEN
        IF(NCOMMU_J.EQ.1) THEN
*. I assume that the third space has been defined
         ITREFSPC_L = 3
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
         WRITE(6,*) ' NOTE : Space 3 is used for COM_JMRCC '
*. Jacobian independent of T, so use T = 0 for simplicity
         ZERO = 0.0D0
         CALL SETVEC(WORK(KLVCC6),ZERO,N_CC_AMP)
         CALL COM_JMRCC(WORK(KLVCC6),NCOMMU_J,WORK(KLJMAT),WORK(KLVCC1),
     &                  WORK(KLVCC2), WORK(KLVCC3), WORK(KLVCC4),
     &                  N_CC_AMP,NSPAM1,NNONSING,IREFSPC,ITREFSPC_L,
     &                  WORK(KLXMAT) )
         ELSE 
*. More than one commutator, so J depends on T
           CALL COM_JMRCC(T_EXT,NCOMMU_J,WORK(KLJMAT),WORK(KLVCC1),
     &                    WORK(KLVCC2), WORK(KLVCC3), WORK(KLVCC4),
     &                    N_CC_AMP,NSPAM1,NNONSING,IREFSPC,ITREFSPC_L,
     &                    WORK(KLXMAT) )
         END IF
*. Obtain LU-Decomposition of Jacobian 
         CALL LULU(WORK(KLJMAT),WORK(KLL),WORK(KLU),NNONSING)
        END IF
*. Solve Linear equations J Delta = - Vecfnc, store solution in VCC1
        ONEM = -1.0D0
        CALL SCALVE(WORK(KLVCC5),ONEM,NNONSING)
        CALL MEMCHK2('AFTSCA')
        CALL LINSOL_FROM_LUCOMP(WORK(KLL),WORK(KLU),WORK(KLVCC5),
     &       WORK(KLVCC1),NNONSING,WORK(KLVCC2))
C     LINSOL_FROM_LUCOMP(XL,XU,RHS,X,NDIM,SCR1)
*. Transform solution to SPA basis and store in VCC2
C  MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
        CALL MATVCC(WORK(KLXMAT),WORK(KLVCC1),WORK(KLVCC2),
     &              NSPAM1,NNONSING,0)
        CALL MEMCHK2('AFTVC2')
        WORK(KLVCC2-1+NSPA) = 0.0D0
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Solution in SPA basis '
          CALL WRTMAT(WORK(KLVCC2),1,NSPA,1,NSPA)
        END IF
*. And transform to CAAB basis  and save in VCC1
C   REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
        CALL REF_CCV_CAAB_SP(WORK(KLVCC1),WORK(KLVCC2),WORK(KLVCC3),2)
        CALL MEMCHK2('AFTRF2')
*. Norm of change
        XNORM = SQRT(INPROD(WORK(KLVCC1),WORK(KLVCC1),N_CC_AMP))
        WRITE(6,*) ' Norm of correction ', XNORM
*. And update the T-coefficients
        ONE = 1.0D0
        CALL VECSUM(T_EXT,T_EXT,WORK(KLVCC1),ONE,ONE,N_CC_AMP)
        CALL MEMCHK2('AFTSUM')
      END DO
*     ^ End of loop over Newton iterations
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from T optimization ', IROOT
        WRITE(6,*) ' Updated amplitudes '
        CALL WRTMAT(T_EXT,1,NCAAB,1,NCAAB)
      END IF
*
      IF(IFIN_IT.EQ.1) 
     &CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICC_CMP')
      RETURN
      END 
      SUBROUTINE ICCI_COMPLETE_MAT2(IREFSPC,ITREFSPC,I_SPIN_ADAPT,
     &        IROOT,T_EXT,C_0,E_IROOT)

*
* Master routine for Internal contraction with complete incore 
* construction of all matrices.
*
* Version using spin adapted basis functions  or EI basis functions
*
* Jeppe Olsen, Sept 2002
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cei.inc'
*. Output : Coefficients of external correlation 
      DIMENSION T_EXT(*)
*. Number of Spin adapted functions ( and NCAAB for a check)
      IF(I_DO_EI.EQ.0) THEN
        CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
      ELSE
        NSPA = N_ZERO_EI
        NCAAB = NDIM_EI
      END IF
      NTEST = 100
      WRITE(6,*) 
      WRITE(6,*) ' Complete H and S matrices will be constructed '
      WRITE(6,*) ' =============================================='
      WRITE(6,*)
      WRITE(6,*) ' Reference space is ', IREFSPC
      WRITE(6,*) ' Space of Operators times reference space ', ITREFSPC
      WRITE(6,*)
      WRITE(6,*) 
     &' Number of parameters in spinuncoupled/original basis ', 
     &           NCAAB
      WRITE(6,*) 
     &' Number of parameters in spincoupled/zero-order  basis ', 
     &           NSPA
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'IC_CMP ')
*. Space for old fashioned CI behind the curtain
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
*
* Space for complete H and S matrices 
*
*. And a few working vectors 
      CALL MEMMAN(KLVCC1,NCAAB,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,NCAAB,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,NCAAB,'ADDL  ',2,'VCC3  ')
      CALL MEMMAN(KLVCC4,NCAAB,'ADDL  ',2,'VCC4  ')
      LEN = NSPA**2
      CALL MEMMAN(KLSHMAT,LEN,'ADDL  ',2,'SHMAT ')
      CALL MEMMAN(KLXMAT,LEN,'ADDL  ',2,'XMAT  ')
*. Identify the unit  operator i.e. the operator with 
*. zero creation and annihilation operators
      IDOPROJ = 1
      IF(IDOPROJ.EQ.1) THEN
        CALL GET_SPOBTP_FOR_EXC_LEVEL(0,WORK(KLCOBEX_TP),NSPOBEX_TP+1,
     &       NUNIOP,IUNITP,WORK(KLSOX_TO_OX))
*. And the position of the unitoperator in the list of SPOBEX operators
        WRITE(6,*) ' NUNIOP, IUNITP = ', NUNIOP,IUNITP
        IF(NUNIOP.EQ.0) THEN
          WRITE(6,*) ' Unitoperator not found in exc space '
          WRITE(6,*) ' I will proceed without projection '
          IDOPROJ = 0
        ELSE
          IUNIOP = IFRMR(WORK(KLIBSOBEX),1,IUNITP)
          WRITE(6,*) ' IUNIOP = ', IUNIOP
        END IF
      END IF
*. Construct metric 
      CALL COM_SH(WORK(KLSHMAT),WORK(KLSHMAT),WORK(KLVCC1),WORK(KLVCC2),
     &            WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &            N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &            IDOPROJ,IUNIOP,1,0,1,I_DO_EI,NSPA,0,0,0)
*. Obtain singularities on S 
      CALL CHK_S_FOR_SING(WORK(KLSHMAT),NSPA,NSING,
     &                    WORK(KLXMAT),WORK(KLVCC1),WORK(KLVCC2))
*. On output the eigenvalues are residing in WORK(KLVCC1) and 
*. the corresponding eigenvectors in WORK(KLXMAT).
*. The singular subspace is defined by the first NSING eigenvectors
      NNONSING = NSPA - NSING
      WRITE(6,*) ' Number of nonsingular eigenvalues of S ', NNONSING
      KLNONSING = KLXMAT + NSING*NSPA
*
      I_ANALYZE_SUM_SING = 0
      IF(I_ANALYZE_SUM_SING.EQ.1) THEN
*. Analyze sum of singularities : Print out Sum(i:sing) C(j,i)**2, 
*. where C(J,I) is in the original basis 
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLVCC3),ZERO,N_CC_AMP)
        DO JSING = 1, NSING
*. Transform to Standard basis 
C     REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
          CALL REF_CCV_CAAB_SP(WORK(KLVCC4),
     &         WORK(KLXMAT-1+(JSING-1)*NSPA),WORK(KLVCC2),2)
*. Square Vector in CAAB basis and add to VCC3)
          CALL VVTOV(WORK(KLVCC4),WORK(KLVCC4),WORK(KLVCC2),N_CC_AMP)
          ONE = 1.0D0
          CALL VECSUM(WORK(KLVCC3),WORK(KLVCC3),WORK(KLVCC2),ONE,ONE,
     &                N_CC_AMP)
        END DO
*. Change so summed sqareed elements add up to one 
        FACTOR = 1.0D0/SQRT(DBLE(NSING))
        DO I = 1, N_CC_AMP
          WORK(KLVCC3-1+I) = SQRT(WORK(KLVCC3-1+I))*FACTOR
        END DO
*. And analyze vector 
        CALL ANA_GENCC(WORK(KLVCC3),1)
      END IF
*     ^ End if sum of singularities should be analyzed
*
*. Obtain transformation to orthonormal basis 
*  X = U sigma^{-1/2}, where U are the nonsingular 
*. eigenvectors of S and sigma are the corresponding 
*. eigenvectors
      DO I = 1, NNONSING
        SCALE = 1/SQRT(WORK(KLVCC1-1+NSING+I))
        CALL SCALVE(WORK(KLNONSING+(I-1)*NSPA),SCALE,NSPA)
      END DO
*. Construct H matrix 
      CALL COM_SH(WORK(KLSHMAT),WORK(KLSHMAT),WORK(KLVCC1),WORK(KLVCC2),
     &            WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &            N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &            IDOPROJ,IUNIOP,0,1,1,I_DO_EI,NSPA,0,0,0)
*. To save space we now need to play a bit around: First we 
*. write H and the needed part of X on disc -they will be 
*. destroyed during transformation 
      LUSCR = 36
      CALL REWINO(LUSCR)
C          TODSC(A,NDIM,MBLOCK,IFIL)
      CALL TODSC(WORK(KLNONSING),NSPA*NNONSING,-1,LUSCR)
      CALL ITODS(-1,1,-1,LUSCR)
      CALL TODSC(WORK(KLSHMAT),NSPA*NSPA,-1,LUSCR)
      CALL ITODS(-1,1,-1,LUSCR)
*. Use low memory routine overwriting the input matrices 
C  TRNMA_LM(XTAX,A,X,NRA,NCA,NRX,NCX,SCRVEC)
      CALL TRNMA_LM(WORK(KLNONSING),WORK(KLSHMAT),WORK(KLNONSING),
     &               NSPA,NSPA,NSPA,NNONSING,WORK(KLVCC1))
      CALL COPVEC(WORK(KLNONSING),WORK(KLSHMAT),NNONSING*NNONSING)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformed Hamiltonian matrix '
        CALL WRTMAT(WORK(KLSHMAT),NNONSING,NNONSING,NNONSING,NNONSING)
      END IF
C     STOP ' Enforced stop after TRANMA_LM'
*
*. Diagonalize transformed Hamiltonian 
*
*. using EISPACK TRED2-TQL2
      IOLD = 0
      IF(IOLD.EQ.0) THEN
        CALL DIAG_SYMMAT_EISPACK(WORK(KLSHMAT),WORK(KLVCC1),
     &                           WORK(KLVCC2),NNONSING,IEIG_RETURN)
      ELSE
        ZERO = 0.0D0
        ONE = 1.0D0
        CALL TRIPAK(WORK(KLSHMAT),WORK(KLXMAT),1,NNONSING,NNONSING)
        CALL COPVEC(WORK(KLXMAT),WORK(KLSHMAT),NNONSING*(NNONSING+1)/2)
        CALL SETVEC(WORK(KLXMAT),ZERO,NNONSING*NNONSING)
        CALL SETDIA(WORK(KLXMAT),ONE,NNONSING,0)
C            SETDIA(MATRIX,VALUE,NDIM,IPACK)
        CALL JACOBI(WORK(KLSHMAT),WORK(KLXMAT),NNONSING,NNONSING)
C            JACOBI(F,V,NB,NMAX) 
        CALL COPDIA(WORK(KLSHMAT),WORK(KLVCC1),NNONSING,1)
        WRITE(6,*) ' Diagonalize JACOBI was used '
        WRITE(6,*) ' This does not order eigenvalues so STOP '
        STOP ' Will not proceed after call to JACOBI '
      END IF
*
      WRITE(6,*) ' Ecore in ICCI_COMPLETE.. ', ECORE
      DO I = 1, NNONSING
        WORK(KLVCC1-1+I) = WORK(KLVCC1-1+I) + ECORE
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvalues of H matrix in IC basis '
        WRITE(6,*) ' ===================================='
        CALL WRTMAT_EP(WORK(KLVCC1),1,NNONSING,1,NNONSING)
      END IF
      E_IROOT = WORK(KLVCC1-1+IROOT) 
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Energy after reoptimization of external',E_IROOT
      END IF
*
      IF(IOLD.NE.0) THEN
       WRITE(6,*) ' Warning : Information for specific root '
       WRITE(6,*) ' can not be obtained as IOLD = 0 does not give '
       WRITE(6,*) ' ordered roots '
      END IF
*. Transform root IROOT to original spin-adapted basis
      CALL COPVEC(WORK(KLSHMAT+(IROOT-1)*NNONSING),WORK(KLVCC2),
     &            NNONSING)
      CALL REWINO(LUSCR)
C      FRMDSC(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
      CALL FRMDSC(WORK(KLNONSING),NSPA*NNONSING,-1,LUSCR,IMZERO,
     /            I_AM_PACKED)
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(WORK(KLNONSING),WORK(KLVCC2),WORK(KLVCC4),
     &            NSPA,NNONSING,0)
      C_0 = 0.0D0
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' NUNIOP, IUNIOP = ',  NUNIOP, IUNIOP
      IF(NUNIOP.NE.0) C_0 = WORK(KLVCC4-1+IUNIOP)
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' C_0 = ', C_0
*. And transform to CAAB basis 
      IF(I_DO_EI.EQ.0) THEN
        CALL REF_CCV_CAAB_SP(T_EXT,WORK(KLVCC4),WORK(KLVCC2),2)
C            REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
      ELSE
*. EI in VCC4 to CAAB in T_EXT
        CALL TRANS_CAAB_ORTN(T_EXT,WORK(KLVCC4),1,2,2,WORK(KLVCC2),2) 
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from IC root nr ', IROOT
        WRITE(6,*) ' Energy is ', WORK(KLVCC1-1+IROOT)
        WRITE(6,*) ' Coefficient of zero-order state ', C_0
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' eigenvector from ICCI eigenequations '
        CALL WRTMAT(T_EXT,1,NCAAB,1,NCAAB)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'IC_CMP ')
      RETURN
      END 
      SUBROUTINE JACOBI(F,V,NB,NMAX)                                    00015000
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(*),V(NMAX,NB)                                         00016000
C                                                                       00017000
C     PURPOSE: TO DIAGONALIZE AN NB*NB-SIZED SUBSPACE OF THE            00018000
C     MATRIX F, AND TO TRANSFORM THE NB VECTORS V OF LENGTH             00019000
C     NMAX BY THE SAME UNITARY MATRIX THAT DIAGONALIZED F.              00020000
C     (NORMAL USAGE: NB=NMAX, AND V IS A UNIT MATRIX WHEN CALLED,       00021000
C     SO THAT V CONTAINS THE EIGENVECTORS ON EXIT.)                     00022000
C     F IS STORED AS UNDER-TRIANGULAR ROWS, AND ON EXIT HAS             00023000
C     BEEN REPLACED BY A NEAR-DIAGONAL MATRIX. THE OUT-OF               00024000
C     DIAGONAL ELEMENTS ARE SMALLER IN SIZE THAN THE PARAMETER          00025000
C     EPS.                                                              00026000
C                                (MALMQUIST 85-02-05)                   00027000
      PARAMETER (EPS=1.E-12,EPS2=EPS*EPS)                               00028000
   1  FMAX=0.0                                                          00029000
      II0=1                                                             00030000
C --- SCAN ALL NON-DIAGONAL ELEMENTS. THIS IS JUST AS EFFICIENT AS      00031000
C --- TO ROTATE SELECTED PAIRS ONLY.                                    00032000
      DO 60 I=2,NB                                                      00033000
        II=II0+I                                                        00034000
        JJ0=0                                                           00035000
        DO 50 J=1,I-1                                                   00036000
          FII=F(II)                                                     00037000
C --- NOTE: FII CANNOT BE SET OUTSIDE THIS LOOP.                        00038000
          IJ=II0+J                                                      00039000
          FIJ=F(IJ)                                                     00040000
          JJ=JJ0+J                                                      00041000
          FJJ=F(JJ)                                                     00042000
          FSQ=FIJ**2                                                    00043000
          FMAX=MAX(FMAX,FSQ)                                            00044000
          IF(FSQ.LT.EPS2) GOTO 40                                       00045000
          DIFFR=FII-FJJ                                                 00046000
          SIGN=1.0                                                      00047000
          IF(DIFFR.LT.0) THEN                                           00048000
            DIFFR=-DIFFR                                                00049000
            SIGN=-SIGN                                                  00050000
          END IF                                                        00051000
          DUM=DIFFR+SQRT(DIFFR**2+4*FSQ)                                00052000
          T=2*SIGN*FIJ/DUM                                              00053000
          C=1.0/SQRT(1+T**2)                                            00054000
          S=C*T                                                         00055000
C --- T,C,S=TAN,COS AND SIN OF ROTATION ANGLE.                          00056000
C --- ROTATE VECTORS:                                                   00057000
          DO 10 K=1,NMAX                                                00058000
            DUM=C*V(K,J)-S*V(K,I)                                       00059000
            V(K,I)=S*V(K,J)+C*V(K,I)                                    00060000
            V(K,J)=DUM                                                  00061000
  10        CONTINUE                                                    00062000
C --- ROTATE F MATRIX COMPONENTS WITH ONE INDEX=I OR J:                 00063000
          DO 31 K=1,J-1                                                 00064000
            KI=II0+K                                                    00065000
            KJ=JJ0+K                                                    00066000
            DUM=C*F(KJ)-S*F(KI)                                         00067000
            F(KI)=S*F(KJ)+C*F(KI)                                       00068000
            F(KJ)=DUM                                                   00069000
  31        CONTINUE                                                    00070000
          KK0=JJ0+J                                                     00071000
          DO 32 K=J+1,I-1                                               00072000
            KI=II0+K                                                    00073000
            KJ=KK0+J                                                    00074000
            DUM=C*F(KJ)-S*F(KI)                                         00075000
            F(KI)=S*F(KJ)+C*F(KI)                                       00076000
            F(KJ)=DUM                                                   00077000
            KK0=KK0+K                                                   00078000
  32        CONTINUE                                                    00079000
          KK0=II0+I                                                     00080000
          DO 33 K=I+1,NB                                                00081000
            KI=KK0+I                                                    00082000
            KJ=KK0+J                                                    00083000
            DUM=C*F(KJ)-S*F(KI)                                         00084000
            F(KI)=S*F(KJ)+C*F(KI)                                       00085000
            F(KJ)=DUM                                                   00086000
            KK0=KK0+K                                                   00087000
  33        CONTINUE                                                    00088000
C--- ROTATE THE II,IJ, AND JJ COMPONENTS:                               00089000
          C2=C**2                                                       00090000
          S2=S**2                                                       00091000
          CIJ=2*C*S*FIJ                                                 00092000
          F(II)=C2*FII+S2*FJJ+CIJ                                       00093000
          F(JJ)=S2*FII+C2*FJJ-CIJ                                       00094000
          F(IJ)=0.0                                                     00095000
  40      JJ0=JJ0+J                                                     00096000
  50      CONTINUE                                                      00097000
        II0=II0+I                                                       00098000
  60    CONTINUE                                                        00099000
C --- CHECK IF CONVERGED:                                               00100000
      IF(FMAX.GT.EPS2) GOTO 1                                           00101000
      RETURN                                                            00102000
      END                                                               00103000
      SUBROUTINE DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
*
* Diagonalize symmetric matrix using eispack routines 
* TRED2 and TQL2
*
* Jeppe Olsen, September 2002  
*
*. Arguments
* ===========
*
* A  : On input :  The matrix in full form
*      On output:  The eigenvectors 
* EIGVAL : Contains eigenvalues on output
* SCRVEC : Scratch vector 
* NDIM   : Dimension of matrices 
* IRETURN : ne 0 => Diagonalization was not complete ...
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION A(NDIM*NDIM)
*. Output
      DIMENSION EIGVAL(*)
*. Scratch
      DIMENSION SCRVEC(*)
*
      CALL LUCIAQENTER('EIS_D')
*
* 1 : Bring matrix to tridiagonal form 
*
      CALL TRED2L(NDIM,NDIM,A,EIGVAL,SCRVEC,A)
* 
* 2 : Obtain eigenvalues from tridiagonal form
*
C     TQL2(NM,N,D,E,Z,IERR)
      CALL TQL2L(NDIM,NDIM,EIGVAL,SCRVEC,A,IRETURN)
*
      IF(IRETURN.NE.0) THEN
        WRITE(6,*) ' Problem in TQL2 diagonalization, IRETURN = ',
     &               IRETURN
        STOP       ' Problem in TQL2 diagonalization '
      END IF
*
      CALL LUCIAQEXIT('EIS_D')
      RETURN
      END 
      SUBROUTINE GET_SXLIKE_CAABM(NSXLIKE,ISXLIKE,IWAY,I_SPIN_ADAPT)
* 
* Obtain spinorbital excitations CAAB that may contain a part of 
* a single excitation
*
* a+ a a i
* a+ a a+x ax ai, where x refers to some orbital index 
*
*
* IWAY = 1 : Just the number of SXlike CAABS
* IWAY = 2 : Number and the actual SXLIKE CAABS
*
*
* Jeppe Olsen, September 2002, for understanding and isolating 
* singularities
*. Modified a bit to allow more general prototypes, aug. 2004
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
*. Output ( IF IWAY.NE.1) 
      INTEGER ISXLIKE(*)
*. Local scratch 
      INTEGER ICASCR(2*MXPNGAS)
*. Loop over the various types of orbital excitations
      IBSXLIKE = 1
      IBCOMP = 1
      DO IOBEX_TP = 1, NOBEX_TP
*. Integer arrays for creation and annihilation part 
          CALL ICOPVE2(WORK(KOBEX_TP),1+(IOBEX_TP-1)*2*NGAS,2*NGAS,
     &                  ICASCR)
          NOP_C = IELSUM(ICASCR,NGAS)
          NOP_A = IELSUM(ICASCR(1+NGAS),NGAS)
          NOP_CA = NOP_C + NOP_A

*. And let another routine do the work for a given 
*. orbital excitation type
*. Effective operator rank of this type of operator 
        CALL GET_SXLIKE_CAAB(IWAY,IBSXLIKE,ISXLIKE,
     &       NCAOC(IOBEX_TP),WORK(KCAOC(IOBEX_TP)),NOP_C,NOP_A,
     &       I_SPIN_ADAPT,IBCOMP,WORK(KNCAAB_FOR_CA(IOBEX_TP)) )
      END DO
*
      NSXLIKE = IBSXLIKE -1
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of SX like operators = ', NSXLIKE
        IF(IWAY.NE.1) THEN
          WRITE(6,*) ' The SX like operators '
          CALL IWRTMA(ISXLIKE,1,NSXLIKE,1,NSXLIKE)
        END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_SXLIKE_CAAB(IWAY,IBSXLIKE,ISXLIKE,
     &           NCA_FOR_TP,ICA_FOR_TP,NOP_C,NOP_A,
     &           I_SPIN_ADAPT,IBCOMP,NCOMP_FOR_CA)
* 
* Obtain -for a given occupation type - the configurations 
* that are effectively single excitations
*
* It is assumed that no operators are purely internal
*
* Jeppe Olsen, Sept. 2002
*. Modified a bit to allow more general prototypes, aug. 2004
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'proto_sp_mat.inc'
*
*. Input
*. The occupation of the configations 
      INTEGER ICA_FOR_TP(NOP_C+NOP_A,NCA_FOR_TP)
*. Number of components for each CA excs
      INTEGER NCOMP_FOR_CA(*)
*. Output (IWAY = 2)
      INTEGER ISXLIKE(*)
*
      NTEST = 10
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Info from  GET_SXLIKE_CAAB '
         WRITE(6,*) 'NOP_C, NOP_A, NCA_FOR_TP = ',
     &               NOP_C, NOP_A, NCA_FOR_TP
      END IF

      DO ICA = 1, NCA_FOR_TP
        IF(NTEST.GE.1000) THEN 
          WRITE(6,*) ' Next CA configuration '
          CALL IWRTMA(ICA_FOR_TP(1,ICA),1,NOP_C+NOP_A,1,NOP_C+NOP_A)
        END IF
        NCOMP_CAAB = NCOMP_FOR_CA(ICA)
        NCOMP_SPA = NSPA_FOP(NCOMP_CAAB)
        IF(I_SPIN_ADAPT.EQ.1) THEN
           NCOMP = NCOMP_SPA
        ELSE
           NCOMP = NCOMP_CAAB
        END IF
        LSX = 0
        IF(NOP_C.EQ.1) THEN 
*. Single excitation 
          LSX = 1
        ELSE IF (NOP_C.EQ.2) THEN
*. Twobody excitation a+ a+ a a, 
          IF(ICA_FOR_TP(1,ICA).EQ.ICA_FOR_TP(3,ICA).OR.
     &       ICA_FOR_TP(1,ICA).EQ.ICA_FOR_TP(4,ICA).OR.
     &       ICA_FOR_TP(2,ICA).EQ.ICA_FOR_TP(3,ICA).OR.
     &       ICA_FOR_TP(2,ICA).EQ.ICA_FOR_TP(4,ICA)) LSX = 1
*
        END IF
        IF(NTEST.GE.1000) THEN
          IF(LSX.EQ.1) THEN
            WRITE(6,*) ' Excitation is single like '
          ELSE
            WRITE(6,*) ' Excitation is not single-like'
          END IF
        END IF
        IF(NTEST.GE.1000) WRITE(6,*) '  NCOMP = ', NCOMP
*
        IF(LSX.EQ.1) THEN
          IF(IWAY.NE.1) THEN 
            DO J = 1, NCOMP
              ISXLIKE(IBSXLIKE-1+J) = IBCOMP-1+J
            END DO
            IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' Corresponding added operators '
               CALL IWRTMA(ISXLIKE(IBSXLIKE),1,NCOMP,1,NCOMP)
            END IF
          END IF
          IBSXLIKE = IBSXLIKE + NCOMP
        END IF
        IBCOMP = IBCOMP + NCOMP
      END DO
*
      RETURN
      END 
      SUBROUTINE SXLIKE_SING(IREFSPC,ITREFSPC,NSXLIKE,I_SPIN_ADAPT)
*
* Study the space of single-excitation like operators 
* and determine singularities in this space 
*
*
* Jeppe Olsen, Oct 1, 2002
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SXLIKE')
*. Dimension of the space of SXLIKE operators 
C     GET_SXLIKE_CAABM(NSXLIKE,ISXLIKE,IWAY,I_SPIN_ADAPT)
      CALL GET_SXLIKE_CAABM(NSXLIKE,IDUM,1,I_SPIN_ADAPT)
*. And the actual operators
      CALL MEMMAN(KLSXLIKE,NSXLIKE,'ADDL  ',2,'SXLIKE')
      CALL GET_SXLIKE_CAABM(NSXLIKE,WORK(KLSXLIKE),2,I_SPIN_ADAPT)
*. Construct the overlap over the SXLIKE operators 
      
      CALL MEMMAN(KLSMAT,NSXLIKE**2,'ADDL  ',2,'SMAT  ')
      CALL MEMMAN(KLX   ,NSXLIKE**2,'ADDL  ',2,'SMAT  ')
      CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,N_CC_AMP,'ADDL  ',2,'VCC3  ')
*. Space for old fashioned CI behind the curtain
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
      IDOPROJ = 1
      IUNIOP = 0
      IF(I_SPIN_ADAPT.EQ.0) THEN
         NSPA = 0
      ELSE 
         CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
      END IF
C     COM_SH(S,H,VCC1,VCC2,VCC3,VEC1,VEC2,
C    &                  N_CC_AMP,IREFSPC,ITREFSPC,
C    &                  LUC,LUHC,LUSCR,LUSCR2,IDOPROJ,IUNIOP,
C    &                  IDO_S,IDO_H,IDO_SPA,NSPA,IDOSUB,ISUB,NSUB)
      CALL COM_SH(WORK(KLSMAT),WORK(KLSMAT),WORK(KLVCC1),WORK(KLVCC2),
     &            WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &            N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &            IDOPROJ,IUNIOP,1,0,I_SPIN_ADAPT,I_DO_EI,NSPA,1,
     &            WORK(KLSXLIKE),NSXLIKE)
*
C?    WRITE(6,*) ' The first 5 rows of S '
C?    CALL WRTMAT(WORK(KLSMAT),5,NSXLIKE,NSXLIKE,NSXLIKE)
C?    WRITE(6,*) ' And the last column '
C?    CALL WRTMAT(WORK(KLSMAT+(NSXLIKE-1)*NSXLIKE),1,NSXLIKE,1,NSXLIKE)
*. Diagonalize  metric and count singularities 
C  CHK_S_FOR_SING(S,NDIM,NSING,X,SCR,SCR2)
      CALL CHK_S_FOR_SING(WORK(KLSMAT),NSXLIKE,NSXSING,WORK(KLX),
     &                  WORK(KLVCC2),WORK(KLVCC3))
      WRITE(6,*) ' Number of singularities in SX like space = ',
     &             NSXSING
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SXLIKE')
      RETURN
      END
      SUBROUTINE FIND_INTEGER_IN_VEC(IVAL,IVEC,NELMNT,IELMNT)
*
* A vector of NELMNT elements is given in IVEC.
* Find the element IELMNT in IVEC with value IVAL.
* If there are several elements with this value, the last element 
* with correct value is returned.
* If an element with the value IELMNT is not obtained, IELMNT is 
* returned as zero 
*
* Jeppe Olsen, Oct. 2002
*
      INCLUDE 'implicit.inc'
*
      INTEGER IVEC(NELMNT)
*
      IELMNT = 0
      DO JELMNT = 1, NELMNT
        IF(IVEC(JELMNT).EQ.IVAL) IELMNT = JELMNT
      END DO
*
      RETURN
      END
      SUBROUTINE GET_ADR_FOR_OCCLS(IOCCLS_SEL,NOCCLS_SEL,NOP,IOP)
*
* Find the number and addresses (in configuration order) of spinadapted 
* operators for the 
* NOCCLS_SEL occupation classes  given in IOCCLS_SEL
*
* The operators are returned in IOP
*
* Jeppe Olsen, Oct 2002, Milano Airport ( Malpensa to be more exact)
*
*. General Input
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'corbex.inc'
*. Specific Input 
      INTEGER IOCCLS_SEL(NOCCLS_SEL)
*. Output
      INTEGER IOP(*)
      IBOP = 1
      DO JJOCCLS = 1, NOCCLS_SEL
         JOCCLS = IOCCLS_SEL(JJOCCLS)
        DO JOP = 1, NSPA_FOR_OCCLS(JOCCLS)
          IOP(IBOP) =  IBSPA_FOR_OCCLS(JOCCLS)-1+JOP
          IBOP = IBOP + 1
        END DO
      END DO
      NOP = IBOP - 1
*
      NTEST = 10 
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' SPA operators for the Excitation types ',
     &  (IOCCLS_SEL(I),I=1, NOCCLS_SEL),' : '
        WRITE(6,*) ' Dimension = ', NOP 
        IF(NTEST.GE.100) CALL IWRTMA(IOP,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_FOR_OBEXTP
*
* Number of CSF's per ocupation class and number of 
* CAAB's per orbital excitation type. 
*
* At the moment the code is adapted to ICCI, so only single 
* and double excitations are considered ( giving atmost 6 dets for 
* a given CONF)
*
* The output is delivered in  NSPA_FOR_OCCLS,NCAAB_FOR_OCCLS 
* given in CORBEX
*
* Jeppe Olsen, Milano Airport, Oct 2002
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'cprnt.inc'
      COMMON/PROTO_SP_MAT/NSPA_FOP(6),NCAAB_FOP(6),IB_FOP(6),XTRA(100),
     &                    NSPA_FOP_G(6,MXPCYC),NCAAB_FOP_G(6,MXPCYC),
     &                    IB_FOP_G(6,MXPCYC)
*. Local scratch
      INTEGER ISCR(MXPNEL)
*. Output is given in CORBEX
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRCSF)
      MAXNDET = 6
      IZERO = 0
*
      DO IOBEX_TP = 1, NOBEX_TP
*. Count the number of times the various number of dets for 
*. a given CA occurs
        CALL ISETVC(ISCR,IZERO,MAXNDET)
        NCA = NCAOC(IOBEX_TP)
        CALL COUNT_OCCURENCE(WORK(KNCAAB_FOR_CA(IOBEX_TP)),ISCR,NCA,
     &                       MAXNDET)
*
        NSPA = ISCR(1)*NSPA_FOP(1) + ISCR(2)*NSPA_FOP(2) 
     &       + ISCR(4)*NSPA_FOP(4) + ISCR(6)*NSPA_FOP(6)
        NCAAB= ISCR(1)*NCAAB_FOP(1) + ISCR(2)*NCAAB_FOP(2) 
     &       + ISCR(4)*NCAAB_FOP(4) + ISCR(6)*NCAAB_FOP(6)
*
        NSPA_FOR_OCCLS(IOBEX_TP) = NSPA
        NCAAB_FOR_OCCLS(IOBEX_TP) = NCAAB
      END DO
*. Offsets for SPA operators belonging to a given occlass
C  ZBASE(NVEC,IVEC,NCLASS)
      CALL ZBASE(NSPA_FOR_OCCLS,IBSPA_FOR_OCCLS,NOBEX_TP)
C     IBSPA_FOR_OCCLS
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Information about operators per orb. exc. type '
        WRITE(6,*) '=================================================='
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Orb. exc. type  Configurations  Spin-adapted    CAAB '
        WRITE(6,*) 
     &  ' ====================================================='
        DO IOBEX_TP = 1, NOBEX_TP 
           WRITE(6,'(6X,I3,6X,I9,6X,I9,4X,I9)')
     &     IOBEX_TP, NCAOC(IOBEX_TP),NSPA_FOR_OCCLS(IOBEX_TP),
     &      NCAAB_FOR_OCCLS(IOBEX_TP)
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE SING_IN_OCCLS(IREFSPC,ITREFSPC,IOCCLS_SEL,NOCCLS_SEL)
*
* Analyze singularities in the space of the SPA operators of the 
* NOCCLS_SEL occupation classes given in IOCCLS_SEL'
*
*
* Jeppe Olsen, Oct 4, 2002
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc'
*. Specific input
      INTEGER IOCCLS_SEL(NOCCLS_SEL)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SING_I')
*. Allocate a vector that can contain addresses of all operators 
      NSPA_TOT = IELSUM(NSPA_FOR_OCCLS,NOBEX_TP)
C?    WRITE(6,*) ' NSPA_TOT = ', NSPA_TOT
      CALL MEMMAN(KLSPAOP,NSPA_TOT,'ADDL  ',1,'SPAOP ')
*. The operators of the specified occupation classes 
C     GET_ADR_FOR_OCCLS(IOCCLS_SEL,NOCCLS_SEL,NOP,IOP)
      CALL GET_ADR_FOR_OCCLS(IOCCLS_SEL,NOCCLS_SEL,NOP,
     &                        WORK(KLSPAOP))
*. Construct the overlap matrix over the these  operators 
      CALL MEMMAN(KLSMAT,NOP**2,'ADDL  ',2,'SMAT  ')
      CALL MEMMAN(KLX   ,NOP**2,'ADDL  ',2,'XMAT  ')
      CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,N_CC_AMP,'ADDL  ',2,'VCC3  ')
*. Space for old fashioned CI behind the curtain
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
      IDOPROJ = 1
      IUNIOP = 0
C     COM_SH(S,H,VCC1,VCC2,VCC3,VEC1,VEC2,
C    &                  N_CC_AMP,IREFSPC,ITREFSPC,
C    &                  LUC,LUHC,LUSCR,LUSCR2,IDOPROJ,IUNIOP,
C    &                  IDO_S,IDO_H,IDO_SPA,NSPA_TOT,IDOSUB,ISUB,NSUB)
      CALL COM_SH(WORK(KLSMAT),WORK(KLSMAT),WORK(KLVCC1),WORK(KLVCC2),
     &            WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &            N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &            IDOPROJ,IUNIOP,1,0,1,I_DO_EI,NSPA_TOT,1,WORK(KLSPAOP),
     &            NOP)
*
*. Diagonalize  metric and count singularities 
      CALL CHK_S_FOR_SING(WORK(KLSMAT),NOP,NSING,WORK(KLX),
     &                  WORK(KLVCC2),WORK(KLVCC3))
      WRITE(6,*) ' Number of singularities in choosen space ',
     &             NSING
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SING_I')
      RETURN
      END 
      SUBROUTINE TRNMA_LM(XTAX,A,X,NRA,NCA,NRX,NCX,SCRVEC)
*
* XTAX = X(T) * A * X
* Low memory version where XTAX may be identical to X, 
* and X and A are overwritten. This works only if the 
* number of columns in X is less than or equal to the 
* number of columns in A
*
* SCRVEC is a scratch vector of the max dimension as A and X
*
* Jeppe Olsen, October 4, 2002
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Input
      DIMENSION A(*),X(*)
*. Output - which may be identical to X
      DIMENSION XTAX(*)
*. Scratch vector 
      DIMENSION SCRVEC(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input matrices X and A to TRNMA_LM '
        CALL WRTMAT(X,NRX,NCX,NRX,NCX)
        CALL WRTMAT(A,NRA,NCA,NRA,NCA)
        WRITE(6,*) ' NRX, NCX, NRA, NCA = ',  NRX, NCX, NRA, NCA
      END IF
*
      IF(NCX.GT.NCA) THEN
        WRITE(6,*) ' TRNMA_LM: NCX gt  NCA: ', NCX,NCA
        STOP 'TRNMA_LM: NCX gt  NCA'
      END IF
*
        
*.1 :  X(T) A in A
      DO L = 1, NCA
*. To avoid compiler warnings
        IB_AKL = 0
        DO I = 1, NCX
          IB_XKI = (I-1)*NRX + 1
          IB_AKL = (L-1)*NRA + 1
          SCRVEC(I) = INPROD(X(IB_XKI),A(IB_AKL),NRA)
        END DO
*. Address of (1,L) in XTA
        IB_AKL = (L-1)*NCX + 1
        IF(NCX.NE.0) THEN
          CALL COPVEC(SCRVEC,A(IB_AKL),NCX)
        ELSE
          ZERO = 0.0D0
          CALL SETVEC(A(IB_AKL),ZERO,NCX)
        END IF
      END DO
* X(T) A X in XTAX
      DO J = 1, NCX
        ZERO = 0.0D0
        CALL SETVEC(SCRVEC,ZERO,NCX)
        DO L = 1, NRX
          XLJ = X((J-1)*NRX+L)
          IB_XTA_IL = (L-1)*NCX + 1
          ONE = 1.0D0
          CALL VECSUM(SCRVEC,SCRVEC,A(IB_XTA_IL),ONE,XLJ,NCX)
        END DO
        CALL COPVEC(SCRVEC,XTAX((J-1)*NCX+1),NCX)
      END DO 
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Outputmatrix from TRANMA_LM '
         CALL WRTMAT(XTAX,NCX,NCX,NCX,NCX)
      END IF
*
      RETURN
      END
      subroutine tranma_lm_test
*
* Test new low memory transformation of matrix
*
* Jeppe Olsen 
*
      INCLUDE 'implicit.inc'
      PARAMETER(MXPDIM = 100)
      DIMENSION A(MXPDIM*MXPDIM),X(MXPDIM*MXPDIM)
      DIMENSION VEC(MXPDIM)
*
      A(1) = 1.0D0
      A(2) = 2.0D0
      A(3) = 3.0D0
      A(4) = 4.0D0
      X(1) = 1.0D0
      X(2) = 2.0D0
      X(3) = 2.0D0
      X(4) = 1.0D0
C  TRNMA_LM(XTAX,A,X,NRA,NCA,NRX,NCX,SCRVEC)
      CALL TRNMA_LM(X,A,X,2,2,2,2,VEC)
*
      RETURN
      END 
      SUBROUTINE SXLIKE_SING2(IREFSPC,ITREFSPC,NSXLIKE,I_SPIN_ADAPT)
*
* 1 : Obtain the single like excitation by diagonalizing in the space 
*     single-like configurations 
*
* 2 : Diagonalize complete metric in space othogonal to SX like 
*     excitations and analyze remaining singulatities
*
* Jeppe Olsen, Oct 7, 2002 - a night session in Palermo 
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SXLIKA')
*. Dimension of the space of SXLIKE operators 
C     GET_SXLIKE_CAABM(NSXLIKE,ISXLIKE,IWAY,I_SPIN_ADAPT)
      CALL GET_SXLIKE_CAABM(NSXLIKE,IDUM,1,I_SPIN_ADAPT)
*. And the actual operators
      CALL MEMMAN(KLSXLIKE,NSXLIKE,'ADDL  ',2,'SXLIKE')
      CALL GET_SXLIKE_CAABM(NSXLIKE,WORK(KLSXLIKE),2,I_SPIN_ADAPT)
*. Construct the overlap over the SXLIKE operators 
      
      CALL MEMMAN(KLSMAT,NSXLIKE**2,'ADDL  ',2,'SMAT  ')
      CALL MEMMAN(KLX   ,NSXLIKE**2,'ADDL  ',2,'SMAT  ')
      CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,N_CC_AMP,'ADDL  ',2,'VCC3  ')
*. Space for old fashioned CI behind the curtain
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
      IDOPROJ = 1
      IUNIOP = 0
      IF(I_SPIN_ADAPT.EQ.0) THEN
         NSPA = 0
      ELSE 
         CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
      END IF
C     COM_SH(S,H,VCC1,VCC2,VCC3,VEC1,VEC2,
C    &                  N_CC_AMP,IREFSPC,ITREFSPC,
C    &                  LUC,LUHC,LUSCR,LUSCR2,IDOPROJ,IUNIOP,
C    &                  IDO_S,IDO_H,IDO_SPA,NSPA,IDOSUB,ISUB,NSUB)
      CALL COM_SH(WORK(KLSMAT),WORK(KLSMAT),WORK(KLVCC1),WORK(KLVCC2),
     &            WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &            N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &            IDOPROJ,IUNIOP,1,0,I_SPIN_ADAPT,I_DO_EI,NSPA,1,
     &            WORK(KLSXLIKE),NSXLIKE)
*
*. Diagonalize  metric and count singularities 
      CALL CHK_S_FOR_SING(WORK(KLSMAT),NSXLIKE,NSXSING,WORK(KLX),
     &                  WORK(KLVCC2),WORK(KLVCC3))
      WRITE(6,*) ' Number of singularities in SX like space = ',
     &             NSXSING
*. On output we have the singularities as the first NSXSING singularities
*. Write these to disc and remove current local allocation 
      CALL REWINO(LUSC1)
      CALL TODSC(WORK(KLX),NSXLIKE*NSXLIKE,-1,LUSC1)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SXLIKA')
*
* Part 2 : Construct complete metric and orthogonalize to 
*          SX like singularitues
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SXLIKB')
*. Memory for metrix and a eigenvector basis 
      NSPA_T = IELSUM(NSPA_FOR_OCCLS,NOBEX_TP) 
      CALL MEMMAN(KLSMAT,NSPA_T**2,'ADDL  ',2,'SMAT  ')
      CALL MEMMAN(KLX   ,NSPA_T**2,'ADDL  ',2,'XMAT  ')
      CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,N_CC_AMP,'ADDL  ',2,'VCC3  ')
*. Space for old fashioned CI behind the curtain
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
      IDOPROJ = 1
      IUNIOP = 0
C     COM_SH(S,H,VCC1,VCC2,VCC3,VEC1,VEC2,
C    &                  N_CC_AMP,IREFSPC,ITREFSPC,
C    &                  LUC,LUHC,LUSCR,LUSCR2,IDOPROJ,IUNIOP,
C    &                  IDO_S,IDO_H,IDO_SPA,NSPA,IDOSUB,ISUB,NSUB)
      CALL COM_SH(WORK(KLSMAT),WORK(KLSMAT),WORK(KLVCC1),WORK(KLVCC2),
     &            WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &            N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &            IDOPROJ,IUNIOP,1,0,I_SPIN_ADAPT,I_DO_EI,
     &            NSPA_T,0,IDUM,IDUM)
*. Recreate the SX like singularities (NSXLIKE is known)
      CALL MEMMAN(KLSXORD,NSPA_T,'ADDL  ',2,'SXORD ')
      CALL GET_SXLIKE_CAABM(NSXLIKE,WORK(KLSXORD),2,I_SPIN_ADAPT)
*. Add terms that are not SX at end of list
      CALL COMPL_LIST(WORK(KLSXORD),NSXLIKE,NSPA_T)
*. Find the configurations that not are single excitations 

*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SXLIKB')
      RETURN
      END
C     COMPL_LIST(WORK(KLSXORD),NSXLIKE,NSPA_T)
      SUBROUTINE COMPL_LIST(ILIST,NIN,NTOT)
* A list is given with NIN elements in 
* ascending order. Complete list so all integers 
* between 1 and NTOT occurs
*
* Jeppe Olsen, Palermo Oct 8 2002, a few hours before liftof
*
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER ILIST(NTOT)
*. Loop over intergers to be in list
      KPIN = 1
      KTOT = NIN 
      DO I = 1, NTOT
*. Is this integer next element included list ?
        IF(KPIN.GT.NIN.OR.I.NE.ILIST(KPIN)) THEN
* I is not in list
          KTOT = KTOT + 1
          ILIST(KTOT) = I
        ELSE
*. I is in list already
          KPIN = KPIN + 1
        END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' completed list from COMPL_LIST '
        WRITE(6,*) ' NIN, NTOT = ', NIN, NTOT
        CALL IWRTMA(ILIST,1,NTOT,1,NTOT)
      END IF
*
      RETURN
      END 
      SUBROUTINE ICCI_RELAX_REFCOEFS_COM(T_EXT,N_EXT,H_REF,S_REF,N_REF,
     &                               VEC1,VEC2,IDO_SPA,IREFSPC,ITREFSPC,
     &                               C_0,ECORE,C_REF_OUT,IREFROOT,NCAAB,
     &                               E_RELAX)
*. Relax internal coefficients in the presence of external 
*. correlation function
*
* Initial version generating complete matrices 
*
* NCAAB is number of operators including unitoperator, all in elementary
* form
*
*
* Redetermine coefficients in reference wavefunction for 
* a given Set of external coefficients given by T_EXT.
*
* The wave-function is given as 
*
* |ICCI > = (C_0 + P \sum_{\mu}T_EXT_{\mu} \hat 0_{\mu} |0 >
*
* where |0> is the reference wave function that we will 
* reoptimize
*
* |0> = \sum_i d_i |i>
*
* P is an projection operator projecting on the orthogonal
* complement space of the reference space
*
* T_EXT is required to be in the CAAB basis
*  
* The equations to be solved are 
*
* H_REF C = E S_REF C with
*
* H_REF_ij = <0_i!H!0_j>
* S_REF_ij = <0_i ! 0_j>
*
* |0_i> = (C_0 + P T) |i >
*. Jeppe Olsen, July 2004,  new way of calculating matrix added aug. 04
*
C     INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPRDD, INPROD
*
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
*. Transfer common block - all parameters have an X here - dirty 
*. and naughty ( in the boring way )
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
*. Input
      DIMENSION T_EXT(N_EXT)
*. Output
      DIMENSION H_REF(N_REF,N_REF),S_REF(N_REF,N_REF)
      DIMENSION C_REF_OUT(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'IC_REL')
*
      NTEST = 100
      WRITE(6,*) 'ICC_RELAX...: NCAAB= ', NCAAB
*
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
*
*
*. Scratch : 3 vectors that can hold T_EXT in expanded form 
*
*. Construct/copy T_EXT in CAAB form in VCC1
      CALL MEMMAN(KLVCC1,NCAAB,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,NCAAB,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLREF1,N_REF   ,'ADDL  ',2,'REF1  ')
*
      CALL COPVEC(T_EXT,WORK(KLVCC1),NCAAB)
*
*. Prepare the transfer common block 
C    &       C_0X,KLTOPX,NREFX,IREFSPX,ITREFSPCX,NCAABX
      C_0X = C_0
      KLTOPX = KLVCC1
      NREFX = N_REF
      IREFSPCX = IREFSPC
      ITREFSPCX = ITREFSPC
      NCAABX = NCAAB
*
      ZERO = 0.0D0
      ONE = 1.0D0
      DO I = 1, N_REF
        CALL SETVEC(WORK(KLREF1),ZERO,N_REF)
        WORK(KLREF1-1+I) = ONE
        CALL H_S_EFF_ICCI_TV(WORK(KLREF1),H_REF(1,I),S_REF(1,I),1,1)
C            H_S_EFF_ICCI_TV(VECIN,VECOUT_H,VECOUT_S)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Effective S-matrix in reference space '
        CALL WRTMAT(S_REF,N_REF,N_REF,N_REF,N_REF)
        WRITE(6,*) ' The Effective H-matrix in reference space '
        CALL WRTMAT(H_REF,N_REF,N_REF,N_REF,N_REF)
      END IF
*
** And diagonalize 
*
C     GENEIG_WITH_SING_CHECK(A,S,EIGVEC,EIGVAL,NVAR,NSING,
C    &                                  WORK)
      LWORK = 5*N_REF**2 + 2*N_REF
      CALL MEMMAN(KLSCR_FOR_GENEIG,LWORK,'ADDL  ',2,'SC_GEI')
      CALL MEMMAN(KLEIGVC,N_REF**2,'ADDL  ',2,'EIGVC ')
      CALL MEMMAN(KLEIGVA,N_REF   ,'ADDL  ',2,'EIGVA ')
      CALL GENEIG_WITH_SING_CHECK(H_REF,S_REF,WORK(KLEIGVC),
     &     WORK(KLEIGVA),N_REF,NSING,WORK(KLSCR_FOR_GENEIG),0)
*
      IF(NSING.NE.0) THEN
        WRITE(6,*) ' Warning : Singularities in Reference CI '
        WRITE(6,*) ' Warning : Singularities in Reference CI '
        WRITE(6,*) ' Warning : Singularities in Reference CI '
        WRITE(6,*) ' Number of singularities = ', NSING
      END IF
*
      NNONSING = N_REF - NSING
      DO I = 1, NNONSING
        WORK(KLEIGVA-1+I) = WORK(KLEIGVA-1+I) + ECORE
      END DO
*. Energy of root IREFROOT
      E_RELAX =  WORK(KLEIGVA-1+IREFROOT)
*. Copy the coefficients of root IROOT to C_REF_OUT
      CALL COPVEC(WORK(KLEIGVC+(IREFROOT-1)*N_REF),C_REF_OUT,N_REF)
*. The eigenvector is normalized with the general metric, 
*. but we want standard normalization so 
      XNORM = INPROD(C_REF_OUT,C_REF_OUT,N_REF)
      SCALE = 1.0D0/SQRT(XNORM)
      WRITE(6,*) ' NORM in ..RELAX.. ', XNORM
      CALL SCALVE(C_REF_OUT,SCALE,N_REF)
*
      WRITE(6,*) ' Eigenvalues of H_EFF matrix '
      WRITE(6,*) ' ============================'
      CALL WRTMAT_EP(WORK(KLEIGVA),1,NNONSING,1,NNONSING)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated coefficients of reference state'
        CALL WRTMAT(C_REF_OUT,1,N_REF,1,N_REF)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'IC_REL')
*
      RETURN
      END 
      SUBROUTINE GENEIG_WITH_SING_CHECK(A,S,EIGVEC,EIGVAL,NVAR,NSING,
     &                                  WORK,IASPACK)
*
* A generalized eigenvalue problem A X = Lambda S X is 
* given for S positive semidefinite. 
*
* Check for singularities, and find eigensolutions in nonsingular subspace
* Intended as subspace diagonalizer for iterative solver, therefore
* not extremely space conserving.
*
* If IASPACK = 1 the input matrices are packed in lower half form 
*            = 0 the input matrices are in complete quadratic form
*
* Jeppe Olsen, Palermo, Oct. 2002 
*
      INCLUDE 'implicit.inc'
*. Input - matrices are supposed to be given in symmetry  packed form
      DIMENSION A(*),S(*)
*. Output 
*. Eigenvectors in input basis 
      DIMENSION EIGVEC(*)
*. And the eigenvalues 
      DIMENSION EIGVAL(*)
*. Scratch : should atleast be 5*NVAR**2 + 2*NVAR
      DIMENSION WORK(*)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Wellcome to  GENEIG_WITH_SING_CHECK '
        WRITE(6,*) ' Dimension of problem = ', NVAR
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input A and S matrices '
        IF(IASPACK.EQ.0) THEN
          CALL WRTMAT(A,NVAR,NVAR,NVAR,NVAR)
          CALL WRTMAT(S,NVAR,NVAR,NVAR,NVAR)
        ELSE 
          CALL PRSYM(A,NVAR)
          CALL PRSYM(S,NVAR)
        END IF
      END IF
C     STOP ' Jeppe forced me to stop '
*. Partition WORK   
*
       KFREE = 1
*
       KSSUB = 1
       KFREE = KFREE + NVAR**2
*
       KMSUB = KFREE
       KFREE = KFREE + NVAR**2
*
       KXORTN = KFREE
       KFREE = KFREE + NVAR**2
*
       KSCRMAT = KFREE
       KFREE   = KFREE + NVAR**2
*
       KSCRMAT2 = KFREE
       KFREE   = KFREE + NVAR**2
*
       KVEC1 = KFREE
       KFREE = KFREE+ NVAR
*
       KVEC2 = KFREE
       KFREE = KFREE+ NVAR
*. Outpack S matrix to full form
       ONE = 1.0D0
C            TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       IF(IASPACK.EQ.1) THEN
         CALL TRIPK3(WORK(KSSUB),S,2,NVAR,NVAR,ONE)
       ELSE 
         CALL COPVEC(S,WORK(KSSUB),NVAR**2)
       END IF
C           GET_ON_BASIS(S,NVEC,NSING,X,SCRVEC1,SCRVEC2)
       CALL GET_ON_BASIS(WORK(KSSUB),NVAR,NSING,WORK(KXORTN),
     &                   WORK(KVEC1),WORK(KVEC2))
       NNONSING = NVAR - NSING
*. Transform A to orthonormal basis 
       IF(IASPACK.EQ.1) THEN
         CALL TRIPK3(WORK(KMSUB),A,2,NVAR,NVAR,ONE)
        ELSE 
          CALL COPVEC(A,WORK(KMSUB),NVAR**2)
        END IF
C       TRNMA_LM(XTAX,A,X,NRA,NCA,NRX,NCX,SCRVEC)
       CALL TRNMA_LM(WORK(KSCRMAT),WORK(KMSUB),WORK(KXORTN),
     &               NVAR,NVAR,NVAR,NNONSING,WORK(KVEC1))
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Matrix in orthonormal nonsingular basis '
         CALL WRTMAT(WORK(KSCRMAT),NNONSING,NNONSING,NNONSING,NNONSING)
        END IF
*. Transformed matrix is returved in KSCRMAT
*. Diagonalize transformed matrix
*
C      DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
        CALL DIAG_SYMMAT_EISPACK(WORK(KSCRMAT),WORK(KVEC1),
     &                           WORK(KVEC2),NNONSING,IRETURN)
        CALL COPVEC(WORK(KVEC1),EIGVAL,NNONSING)
*. Obtain the eigenvectors in the original basis 
       FACTORC = 0.0D0
       FACTORAB = 1.0D0
       CALL MATML7(EIGVEC,WORK(KXORTN),WORK(KSCRMAT),NVAR,NNONSING,
     &             NVAR,NNONSING,NNONSING,NNONSING,FACTORC,FACTORAB,0)
       IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvalues ' 
        CALL WRTMAT(WORK(KVEC1),1,NNONSING,1,NNONSING)
        WRITE(6,*) ' Lowest eigenvector '
        CALL WRTMAT(EIGVEC(1),1,NVAR,1,NVAR)
       END IF
       IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Eigenvectors in original basis '
        CALL WRTMAT(EIGVEC,NVAR,NNONSING,NVAR,NNONSING)
      END IF
     &             
*
      RETURN
      END
      SUBROUTINE GET_ON_BASIS(S,NVEC,NSING,X,SCRVEC1,SCRVEC2)
*
* NVEC vectors with overlap matrix S are given.
* Obtain transformation matrix to orthonormal basis
*
* NSING is the number of singularities obtained 
* If there are singularities, the nonsingular transformation 
* os obtained as a NVEC x (NVEC-NSING) matrix in X 
* First vectors. The eigenvectors corresponding to the 
* singular eigenvectors are lost. 
*
*
* Jeppe Olsen, Palermo, oct 2002
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION S(NVEC*NVEC)
*. Output
      DIMENSION X(NVEC*NVEC)
*. Local scratch
      DIMENSION SCRVEC1(*), SCRVEC2(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GET_ON_BASIS speaking '
        WRITE(6,*) ' Input overlap matrix '
        CALL WRTMAT(S,NVEC,NVEC,NVEC,NVEC)
      END IF
*1 : Diagonalize S and save eigenvalues in SCRVEC1
      CALL COPVEC(S,X,NVEC*NVEC)
C          DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
      CALL DIAG_SYMMAT_EISPACK(X,SCRVEC1,SCRVEC2,NVEC,IRETURN)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvalues of metric '
        CALL WRTMAT(SCRVEC1,1,NVEC,1,NVEC)
      END IF
*2 : Count number of nonsingularities
      NNONSING = 0
      THRES = 1.0D-14
      DO I = 1, NVEC
        IF(ABS(SCRVEC1(I)).GT.THRES) THEN
          NNONSING = NNONSING + 1
          IF(I.NE.NNONSING) THEN
            SCRVEC1(NNONSING) = SCRVEC1(I)
            CALL COPVEC(X((I-1)*NVEC+1), X((NNONSING-1)*NVEC+1),NVEC)
          END IF
        END IF
      END DO
      NSING = NVEC - NNONSING
*2 : Rearrange so the nonsingular
*    eigenvectors and eigenvalues are  the first parts of X and 
*    SCRVEC1
CE    ISING = 0
CE    INONSING = 0
CE    DO I = 1, NVEC
CE      IF(ABS(SCRVEC1(I)) .GT. THRES) THEN
*. A nonsingular eigenpair
CE        INONSING = INONSING + 1
CE        ITO = INONSING
CE      ELSE 
*. A singular eigenpair
CE        ISING = ISING + 1
CE        ITO = ISING + NNONSING
CE      END IF
CE      IF(ITO.NE.I) THEN
CE        SCRVEC1(ITO) = SCRVEC1(I)
CE        CALL COPVEC(X((I-1)*NVEC+1), X((ITO-1)*NVEC+1),NVEC)
CE      END IF
CE    END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Nonsingular eigenvalues of metric '
        CALL WRTMAT(SCRVEC1,1,NNONSING,1,NNONSING)
      END IF
*3 : Construct orthonormal basis using 
*  X = U sigma^{-1/2}, 
*  where U are the nonsingular 
*. eigenvectors of S and sigma are the corresponding eigenvalues
      DO I = 1, NNONSING
        SCALE = 1/SQRT(SCRVEC1(I))
        IBX = (I-1)*NVEC+1
        CALL SCALVE(X(IBX),SCALE,NVEC)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformation matrix to nonsingular basis '
        CALL WRTMAT(X,NVEC,NNONSING,NVEC,NNONSING)
      END IF
*
      RETURN
      END 
      SUBROUTINE INFO2_FOR_PROTO_CA(
     &           NOBEX_TP,IOBEX_TP,ISOX_FOR_OX,NSOX_FOR_OX,IBSOX_FOR_OX,
     &           ISPOBEX_TP,NGAS,
     &           IB_PROTO_CA, MX_DBL_C_CA, MX_DBL_A_CA,
     &           NCOMP_FOR_PROTO_CA,NPROTO_CA)
*
* Info on the number of CAAB excitations for a CA operator
* with due respect given to the number of double occupied orbitals
* in the CA operators
*
* To obtain the number of CAAB components belonging to a given 
* CA excitations two things must be taken into account
* 1) the types of spin-orbital excitations belonging to this type 
* 2) the number of doubly occuring indeces in the C and in the A
*    part of the orbextp
*    
*
* So a prototype spin-orbital excitation is defined by three 
* numbers
* 1) the orbital excitation type (JOBEX_TP)
* 2) the number of doubly occupied orbital in the C part  (NDBL_C)
* 3) the number of double occupied orbitals in the A part (NDLB_A)
*
*. A prototype CA will thus be given the number/adress
*  IB_PROTO_CA(JOBEX_TP) + NDBL_A*(MAX_DBL_C+1) + NDBL_C
*
*. Thus, presently a prototype does not distinguish
*. between CA operators having doubly occupied orbitals
*. in different orbital subspaces. THis may be a problem
*  when more than 2 e ex operators must be included.
* Jeppe Olsen, August 2004
*
      INCLUDE 'implicit.inc'
      INCLUDE 'cprnt.inc'
*. IPRCSF is printflag in charge
*.  Input
*. ======
*. The CA operators 
      INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
*. Number of spin-orbital excitations for each orbital excitations
      INTEGER NSOX_FOR_OX(NOBEX_TP)
*. And the number/address of the spinorbital excitations for each orbexc
*. the adress refers to ISPOBEX_TP
      INTEGER ISOX_FOR_OX(NOBEX_TP)
*. Start in ISOX_FOR_OC for spinorbital exc belonging to given orbexc
      INTEGER IBSOX_FOR_OX(NOBEX_TP)
*. and the actual spin-orbital excitations 
      INTEGER ISPOBEX_TP(4*NGAS,*)
*.========
*. Output
*.========
*
*. Offset for prototypes CA belonging to a given CA
      INTEGER IB_PROTO_CA(NOBEX_TP)
*. max number of double occupied orbital in C part for given CA type
      INTEGER MX_DBL_C_CA(NOBEX_TP)
*. max number of double occupied orbital in A part for given CA type
      INTEGER MX_DBL_A_CA(NOBEX_TP)
*. Number of CAAB components for given prototype of CA
      INTEGER NCOMP_FOR_PROTO_CA(NPROTO_CA)
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRCSF)
*
*. Number and offset for prototypes for given CA type
*
      IOFF = 1
      DO JOBEX_TP = 1, NOBEX_TP
         IB_PROTO_CA(JOBEX_TP) = IOFF
        DO ICA = 1, 2
          MXDBL = 0
          DO IGAS = 1, NGAS
            MXDBL = MXDBL + IOBEX_TP((ICA-1)*NGAS+IGAS,JOBEX_TP)/2
          END DO
          IF(ICA.EQ.1) THEN 
             MX_DBL_C_CA(JOBEX_TP) =  MXDBL
          ELSE
             MX_DBL_A_CA(JOBEX_TP) =  MXDBL
          END IF
        END DO
        IOFF = IOFF + 
     &  (MX_DBL_C_CA(JOBEX_TP)+1)*(MX_DBL_A_CA(JOBEX_TP)+1)
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Max number of double occ orbs in C part '
        CALL IWRTMA(MX_DBL_C_CA,1,NOBEX_TP,1,NOBEX_TP)
        WRITE(6,*) ' Max number of double occ orbs in A part '
        CALL IWRTMA(MX_DBL_A_CA,1,NOBEX_TP,1,NOBEX_TP)
        WRITE(6,*) ' Offset for proto CA types '
        CALL IWRTMA(IB_PROTO_CA,1,NOBEX_TP,1,NOBEX_TP)
      END IF
*
*. Number of CAAB components per prototype CA 
*
      DO JOBEX_TP  = 1, NOBEX_TP
        DO NDBL_C = 0, MX_DBL_C_CA(JOBEX_TP)
        DO NDBL_A = 0, MX_DBL_A_CA(JOBEX_TP)
          IPROTO = IB_PROTO_CA(JOBEX_TP) 
     &   + (MX_DBL_C_CA(JOBEX_TP)+1)*NDBL_A + NDBL_C 
C?        WRITE(6,*) ' Info for IPROTO = ', IPROTO
*. Loop over spin-components of this excitation 
          ISPOX_START = IBSOX_FOR_OX(JOBEX_TP)
          ISPOX_STOP  = ISPOX_START + NSOX_FOR_OX(JOBEX_TP)-1
C?        WRITE(6,*) ' JOBEX_TP, START, STOP ',
C?   &     JOBEX_TP , ISPOX_START,  ISPOX_STOP
          NCOMP_PROTO = 0
          DO JJSPOBEX_TP = ISPOX_START,ISPOX_STOP
            NDBL_C_LEFT = NDBL_C
            NDBL_A_LEFT = NDBL_A
            JSPOBEX_TP = ISOX_FOR_OX(JJSPOBEX_TP)
C?          WRITE(6,*) ' INFO2, JJSPOBEX_TP, JSPOBEX_TP ',
C?   &      JJSPOBEX_TP, JSPOBEX_TP
            NCOMP = 1
            DO JGAS = 1, NGAS
*.  Number CA,CB operators in this SPOX
C?            WRITE(6,*) ' in INFO2.. ', JGAS, JSPOBEX_TP,
C?   &        JGAS, JSPOBEX_TP
              NCA = ISPOBEX_TP(JGAS+0*NGAS,JSPOBEX_TP)
              NCB = ISPOBEX_TP(JGAS+1*NGAS,JSPOBEX_TP)
*. Put as many double occupied orbitals in this space
              ND_C = MIN(MIN(NCA,NCB),NDBL_C_LEFT)
              NDBL_C_LEFT = NDBL_C_LEFT - ND_C
              NCA_S = NCA - ND_C
              NCB_S = NCB - ND_C
C?            WRITE(6,*) ' NCA_S, NCB_S = ', NCA_S, NCB_S
              NC_COMP = IBION(NCA_S+NCB_S,NCB_S)
*.  Number AA,AB operators in this SPOX
              NAA = ISPOBEX_TP(JGAS+2*NGAS,JSPOBEX_TP)
              NAB = ISPOBEX_TP(JGAS+3*NGAS,JSPOBEX_TP)
C?            WRITE(6,*) ' NAA, NAB = ', NAA, NAB
*. Put as many double occupied orbitals in this space
              ND_A = MIN(MIN(NAA,NAB),NDBL_A_LEFT)
              NDBL_A_LEFT = NDBL_A_LEFT - ND_A
              NAA_S = NAA - ND_A
              NAB_S = NAB - ND_A
C?            WRITE(6,*) ' NAA_S, NAB_S = ', NAA_S, NAB_S
              NA_COMP = IBION(NAA_S+NAB_S,NAB_S)
              NCOMP = NCOMP*NC_COMP*NA_COMP
C?            WRITE(6,*) ' JGAS, NA_COMP,NC_COMP =',
C?   &        JGAS,NA_COMP,NC_COMP
            END DO
*           ^ End of loop over GAS spaces
C?          WRITE(6,*) ' Number of comps for this spox', NCOMP
            IF(NDBL_C_LEFT.EQ.0.AND.NDBL_A_LEFT.EQ.0)
     &      NCOMP_PROTO = NCOMP_PROTO + NCOMP
          END DO
*         ^ End of loop over spinorbitalexcitations
          NCOMP_FOR_PROTO_CA(IPROTO) = NCOMP_PROTO
        END DO
        END DO
*       ^ End of loop over number of doubly occ C and A operators
      END DO
*     ^ End of loop over orbital excitations
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of CAAB components per prototype '
        CALL IWRTMA(NCOMP_FOR_PROTO_CA,1,NPROTO_CA,1,NPROTO_CA)
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' Information about prototype CA excitations '
        WRITE(6,*) ' ==========================================='
        WRITE(6,*)
        WRITE(6,*) ' Number   Obextp   ndbl_c   ndbl_a   ncomp '
        WRITE(6,*) ' =========================================='
        IPROTO = 0
        DO JOBEX_TP = 1, NOBEX_TP
          DO NDBL_A = 0, MX_DBL_A_CA(JOBEX_TP)
          DO NDBL_C = 0, MX_DBL_C_CA(JOBEX_TP)
            IPROTO = IPROTO + 1
            NCOMP = NCOMP_FOR_PROTO_CA(IPROTO)
            WRITE(6,'(5(3X,I5))') 
     &      IPROTO, JOBEX_TP, NDBL_C, NDBL_A, NCOMP
          END DO
          END DO
        END DO
      END IF
*
      RETURN
      END 
      FUNCTION NPROTO_CA(NOBEX_TP,IOBEX_TP,NGAS)
*
* Find the number of prototype CA operators 
* A prototype CA is ( at least today, aug 5, 2004)
* defined by orbital excitation, and the number of 
* orbitals occuring twice in the C and A parts.
*. Thus, presently a prototype does not distinguish
*. between CA operators having doubly occupied orbitals
*. in different orbital subspaces. THis may be a problem
*  when more than 2 e ex operators must be included.
*
* Jeppe Olsen, Aug 2005
*
      INCLUDE 'implicit.inc'
*
*. Input
      INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to NPROTO_CA '
        WRITE(6,*) ' NOBEX_TP, NGAS = ', NOBEX_TP, NGAS
        WRITE(6,*) ' IOBEX:'
        CALL IWRTMA(IOBEX_TP,2*NGAS,NOBEX_TP,2*NGAS,NOBEX_TP)
      END IF
*
*. Compiler warnings ...
      MXDBL_C = -2810
      MXDBL_A = -2810
* 
      NPROTO = 0
      DO JOBEX_TP = 1, NOBEX_TP
        DO ICA = 1, 2
          MXDBL = 0
          DO IGAS = 1, NGAS
            MXDBL = MXDBL + IOBEX_TP((ICA-1)*NGAS+IGAS,JOBEX_TP)/2
          END DO
          IF(ICA.EQ.1) THEN 
             MXDBL_C = MXDBL
          ELSE
             MXDBL_A = MXDBL
          END IF
        END DO
        NPROTO = NPROTO + (MXDBL_C+1)*(MXDBL_A+1)
      END DO
*
      NPROTO_CA = NPROTO 
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of prototype CA''s ', NPROTO_CA
      END IF
*
      RETURN
      END
      FUNCTION IPROTO_TYPE_FOR_CA(ICAEX,IOBEX_TP,NOP_C,NOP_A)
*
*. Obtain prototype number for a given CAEX.
*
*. Jeppe Olsen, August 2004
*
C     INCLUDE 'implicit.inc'
*. General input
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*. Specific input
      DIMENSION ICAEX(*)
C K_MX_DLB_C,K_MX_DLB_A,K_IB_PROTO,K_NCOMP_FOR_PROTO
*. Number of double occupied orbital indeces in C and A part
       NCL_C = NCL_FOR_CONF(ICAEX(1),NOP_C)
       NCL_A = NCL_FOR_CONF(ICAEX(1+NOP_C),NOP_A)
*. Obtain MAX number of CL orbitals in C and A parts for this type
       MX_CL_C = IFRMR(WORK(K_MX_DLB_C),1,IOBEX_TP)
       MX_CL_A = IFRMR(WORK(K_MX_DLB_A),1,IOBEX_TP)
*. And offset to prototypes for this obextp
       IB = IFRMR(WORK(K_IB_PROTO),1,IOBEX_TP)
*
       IPROTO = IB + (MX_CL_C+1)*NCL_A + NCL_C
*
       IPROTO_TYPE_FOR_CA = IPROTO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' C and A parts of CA operator '
        CALL IWRTMA(ICAEX,1,NOP_C,1,NOP_C)
        CALL IWRTMA(ICAEX(1+NOP_C),1,NOP_A,1,NOP_A)
        WRITE(6,*) ' NCL_C and NCL_A ', NCL_C, NCL_A
        WRITE(6,*) ' CAex corresponds to protoype ', IPROTO
      END IF
*
      RETURN
      END
      SUBROUTINE LUCIA_ICPT(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                      EFINAL,CONVER,VFINAL)
*
* Master routine for Internal Contraction perturbation theory
*
* LUCIA_IC is assumed to have been called to do the 
* prepatory work for working with internal contraction
*
* It is assumed that spin-adaptation is used ( no flag anymore..)
*
* It is standard that the unitoperator is included in 
* the operator manifold, but in PT theory this should be 
* excluded. This is easily done as the unitoperator is the 
* last operator in CA order.
*
* Jeppe Olsen, August 2004
*
C     INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
      LOGICAL CONVER
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cei.inc'
*. Transfer common block for communicating with H_EFF * vector routines
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
*. Transfer block for communicating zero order energy to 
*. routien for performing H0-E0 * vector
      include 'cshift.inc'
*
      CHARACTER*6 ICTYP
      EXTERNAL H0ME0TV_EXT_IC
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICPT  ')
      NTEST = 1001
      WRITE(6,*)
      WRITE(6,*) ' ===================='
      WRITE(6,*) ' ICPT section entered '
      WRITE(6,*) ' ===================='
      WRITE(6,*)
*
*. Form of ICPT calculation 
*
      IF(ICTYP(1:5).EQ.'ICPT2') THEN
        WRITE(6,*) ' Second-order calculation '
      ELSE IF (ICTYP(1:5).EQ.'ICPT3') THEN
        WRITE(6,*) ' Third-order calculation '
      ELSE 
        WRITE(6,'(A,A)') ' Unknown ICPT form : ', ICTYP
        STOP ' Unknown ICPT form '
      END IF
*
      IF(I_DO_EI.EQ.1) THEN
       WRITE(6,*) ' EI approach in use'
      ELSE
       WRITE(6,*) ' Partial spin-adaptation in use'
      END IF
*

      WRITE(6,*) ' Energy of reference state ', EREF
*. Number of parameters with and without spinadaptation
      IF(I_DO_EI.EQ.0) THEN
        CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
      ELSE
*. zero-particle operator is included in N_ZERO_EI
        NSPA = N_ZERO_EI
*. Note: NCAAB includes unitop
        NCAAB = NDIM_EI
      END IF
      IF(I_DO_EI.EQ.0) THEN
          WRITE(6,*) ' Number of spin-adapted operators ', NSPA
      ELSE
          WRITE(6,*) ' Number of orthonormal zero-order states',
     &    N_ZERO_EI
      END IF
      WRITE(6,*) ' Number of CAAB operators         ', NCAAB
*. Number of spin adapted operators without the unitoperator
      I_DIR_OR_IT = 2
      IF(I_DIR_OR_IT.EQ.1) THEN
        WRITE(6,*) ' Explicit construction of all matrices'
      ELSE
        WRITE(6,*) ' Iterative solution of equations'
      END IF
*
      NSPAM1 = NSPA - 1
*
* ==================================================
* 1 : Set up zero-order Hamiltonian in WORK(KFIFA)
* ==================================================
*
*. It is assumed that one-body density over reference resides 
*  in WORK(KRHO1)
*. Calculate zero-order Hamiltonian: use either actual or Hartree-Fock density
      I_ACT_OR_HF = 1
*. Zero-offdiagonal elements ?
      I_ZERO_OFF = 0
      IF(I_ACT_OR_HF.EQ.1) THEN
        WRITE(6,*) ' Zero-order Hamiltonian with actual density '
*. Inactive Fock matrix and core-energy- with original def. of 
*  inactive terms 
        CALL COPVEC(WORK(KH),WORK(KHINA),NINT1)
        CALL FISM(WORK(KHINA),ECC)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' The (standard) inactive Fock matrix '
          CALL APRBLM2(WORK(KHINA),NTOOBS,NTOOBS,NSMOB,1)
        END IF
        CALL FAM(WORK(KFIFA))
*. and add active and inactive fock matrix
        ONE = 1.0D0
        CALL VECSUM(WORK(KFIFA),WORK(KFIFA),WORK(KHINA),
     &              ONE,ONE,NINT1)
      
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' FI + FA matrix '
          CALL APRBLM2(WORK(KFIFA),NTOOBS,NTOOBS,NSMOB,1)
        END IF
      ELSE
        WRITE(6,*) ' Zero-order Hamiltonian with zero-order density '
        STOP ' I doubt this route is working says Jeppe '
*. IPHGAS1 should be used to divide into H,P,V, but IPHGAS is used, so swap
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' IPHGAS1 : '
          CALL IWRTMA(IPHGAS1(1),1,NGAS,1,NGAS)
        END IF
        CALL ISWPVE(IPHGAS(1),IPHGAS1(1),NGAS)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' IHPGAS in use '
          CALL IWRTMA(IPHGAS(1),1,NGAS,1,NGAS)
        END IF
*
        CALL COPVEC(WORK(KINT1O),WORK(KFIFA),NINT1)
        CALL FI(WORK(KFIFA),ECC,1)
        IF(NTEST.GE.100)THEN
          WRITE(6,*) ' FI before zeroing : '
          CALL APRBLM2(WORK(KFIFA),NTOOBS,NTOOBS,NSMOB,1)
        END IF
*. And clean up
        CALL ISWPVE(IPHGAS,IPHGAS1,NGAS)
*. zero offdiagonal elements 
C            ZERO_OFFDIAG_BLM(A,NBLOCK,LBLOCK,IPACK)
        IF(I_ZERO_OFF.EQ.1) 
     &  CALL ZERO_OFFDIAG_BLM(WORK(KFIFA),NSMOB,NTOOBS,1)
      END IF
*     ^ End if we should use actual or Hartree-Fock density
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' One-body zero-order Hamiltonian '
        CALL APRBLM2(WORK(KFIFA),NTOOBS,NTOOBS,NSMOB,1)
      END IF
*. Obtain zero-order energy
      CALL COPVEC(WORK(KFIFA),WORK(KINT1),NINT1)
*. Contributions from inactive orbitals
      E0INA =  EXP_ONEEL_INACT(WORK(KFIFA),1)
*. Contributions from active orbitals
      CALL EN_FROM_DENS(E0ACT,1,0)
*. And the synthesis
      E0FIFA = ECORE_EXT + E0INA + E0ACT
      WRITE(6,'(A,4E15.8)') ' E0FIFA,ECORE_EXT,E0INA,E0ACT =',
     &                        E0FIFA,ECORE_EXT,E0INA,E0ACT
      E0 = E0FIFA
*. Scratch space for CI 
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
*
*
* =====================================================================
* Obtain metric matrix and nonsingular set of operators in WORK(KLXMAT)
* =====================================================================
*
*. Some additional scratch, dominated by two complete matrices !!
*. And a few working vectors 
      CALL MEMMAN(KLVCC1,NCAAB,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,NCAAB,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,NCAAB,'ADDL  ',2,'VCC3  ')
      CALL MEMMAN(KLVCC4,NCAAB,'ADDL  ',2,'VCC4  ')
      CALL MEMMAN(KLRHS ,NCAAB,'ADDL  ',2,'RHS   ')
      CALL MEMMAN(KLC1  ,NCAAB,'ADDL  ',2,'C1    ')
      CALL MEMMAN(KLC1O ,NCAAB,'ADDL  ',2,'C1    ')
*. Identify the unit  operator i.e. the operator with 
*. zero creation and annihilation operators
      IDOPROJ = 1
      IF(IDOPROJ.EQ.1) THEN
        CALL GET_SPOBTP_FOR_EXC_LEVEL(0,WORK(KLCOBEX_TP),NSPOBEX_TP,
     &       NUNIOP,IUNITP,WORK(KLSOX_TO_OX))
*. And the position of the unitoperator in the list of SPOBEX operators
*. that is, in the CAAB representation
        WRITE(6,*) ' NUNIOP, IUNITP = ', NUNIOP,IUNITP
        IF(NUNIOP.EQ.0) THEN
          WRITE(6,*) ' Unitoperator not found in exc space '
          WRITE(6,*) ' I will proceed without projection '
          IDOPROJ = 0
        ELSE
          IUNIOP = IFRMR(WORK(KLIBSOBEX),1,IUNITP)
          IF(NTEST.GE.100) WRITE(6,*) ' IUNIOP = ', IUNIOP
        END IF
      END IF
*. 
*. Prepare transfer common block used for H(ICCI) * v, S(ICCI) * v 
* ( also used for constructing H,S)
*. The First three entries below are not used 
      C_0X = 0.0D0
      KLTOPX = -1
      NREFX = -1
*. Used
      IREFSPCX = IREFSPC
      ITREFSPCX = ITREFSPC
      IPROJSPCX = IREFSPC
      NCAABX = N_CC_AMP
      NSPAX = NSPA
*. Unitoperator in SPA format
      IUNIOPX = NSPA
*
*
      IF(I_DIR_OR_IT.EQ.1) THEN
*
* Approach based on construction of all matrices.
* Matrices are constructed in the partial spin-adapted or 
* in the zero-order basis
*
*. Construct complete matrices in the SPA representation
        LEN = NSPA**2
        CALL MEMMAN(KLSHMAT,LEN,'ADDL  ',2,'SHMAT ')
        CALL MEMMAN(KLXMAT ,LEN,'ADDL  ',2,'XMAT  ')
        IF(I_DO_EI.EQ.1) THEN
          I_DO_SPA = 0
        ELSE
          I_DO_SPA = 1
        END IF
*. The metric
        CALL COM_SH(WORK(KLSHMAT),WORK(KLSHMAT),WORK(KLVCC1),
     &              WORK(KLVCC2),
     &              WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &              N_CC_AMP,IREFSPC,ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &              IDOPROJ,IUNIOP,1,0,I_DO_SPA,I_DO_EI,NSPA,0,0,0)
        IREFSPCX = IREFSPC
*. ELiminate part referring to unit operator 
       CALL TRUNC_MAT(WORK(KLSHMAT),NSPA,NSPA,NSPAM1,NSPAM1)
C     TRUNC_MAT(A,NRI,NCI,NRO,NCO)
*. Obtain orthonormal basis for nonsingular part of S
C     GET_ON_BASIS(S,NVEC,NSING,X,SCRVEC1,SCRVEC2)
        CALL GET_ON_BASIS(WORK(KLSHMAT),NSPAM1,NSING,
     &                  WORK(KLXMAT),WORK(KLVCC1),WORK(KLVCC2) )
        WRITE(6,*) ' Number of singularities in S ', NSING
        NNONSING = NSPAM1 - NSING 
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Transformation matrix to nonsingular basis '
          CALL WRTMAT(WORK(KLXMAT),NSPAM1,NNONSING,NSPAM1,
     &                NNONSING)
        END IF
*. Save transformation to orthonormal basis - WORK(KLXMAT) will be overwritten
        LU28 = IGETUNIT(28)
        CALL REWINO(LU28)
        CALL VEC_TO_DISC(WORK(KLXMAT),NSPAM1*NNONSING,1,-1,LU28)
* 
* =======================================================
* Set up RHS of first-order equations = <0!H P T_{\mu}!0>
* =======================================================
*
        I12 = 2
        CALL GET_ICPT_RHS1(WORK(KLRHS),IREFSPC,ITREFSPC,
     &                     NSPA,NCAAB,I_DO_EI,
     &                     WORK(KVEC1),WORK(KVEC2),
     &                     WORK(KLVCC1),WORK(KLVCC2)        )
C     GET_ICPT_RHS1(RHS,IREFSPC,ITREFSPC,
C    &                        NSPA,NCAAB,
C    &                        VEC1,VEC2,VIC1,VIC2)
*. Transform RHS to orthonormal basis 
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
        CALL MATVCC(WORK(KLXMAT),WORK(KLRHS),WORK(KLVCC1),NSPAM1,
     &              NNONSING,1)
        CALL COPVEC(WORK(KLVCC1),WORK(KLRHS),NNONSING)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' RHS in orthonormal basis '
          CALL WRTMAT(WORK(KLRHS),1,NNONSING,1,NNONSING)
        END IF
*
* 
* =======================================================
* Set up Zero-order Hamiltonian in WORK(KLSHMAT)
* =======================================================
*
*. Complete matrix including unitop
*
*. Make KINT1 the zero-order-hamiltonian
        CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
*. And tell CI only to work with one-electron operator
        I12 = 1
        CALL COM_SH(WORK(KLSHMAT),WORK(KLSHMAT),WORK(KLVCC1),
     &              WORK(KLVCC2),
     &              WORK(KLVCC3),WORK(KVEC1),WORK(KVEC2),
     &              N_CC_AMP,IREFSPC, ITREFSPC,LUC,LUHC,LUSC1,LUSC2,
     &              IDOPROJ,IUNIOP,0,1,I_DO_SPA,I_DO_EI,NSPA,0,0,0)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The zero-order Hamiltonian in SPA basis '
          CALL WRTMAT(WORK(KLSHMAT),NSPA,NSPA,NSPA,NSPA)
        END IF
*E0 is the last element of H so
        E0 = WORK(KLSHMAT-1+(NSPA-1)*NSPA+NSPA)
        WRITE(6,*) ' The zero-order energy ', E0
*. Eliminate the unit-operator from H0
         CALL TRUNC_MAT(WORK(KLSHMAT),NSPA,NSPA,NSPAM1,NSPAM1)
*. Transform H to orthonormal basis 
C       TRNMA_LM(XTAX,A,X,NRA,NCA,NRX,NCX,SCRVEC)
        CALL TRNMA_LM(WORK(KLXMAT),WORK(KLSHMAT),WORK(KLXMAT),
     &                NSPAM1,NSPAM1,NSPAM1,NNONSING,WORK(KLVCC1) )
        CALL COPVEC(WORK(KLXMAT),WORK(KLSHMAT),NNONSING*NNONSING)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The zero-order Hamiltonian in orthonormal basis '
          CALL WRTMAT(WORK(KLSHMAT),NNONSING,NNONSING,NNONSING,NNONSING)
        END IF
*
* 
* =====================================
* Obtain First order correction to wf 
* =====================================
*
*. H0 - E0*1
*
        FACTOR = - E0
        CALL ADDDIA(WORK(KLSHMAT),FACTOR,NNONSING,0)
*. Diagonalixe H0-E0 , eigenvectors are returned in WORK(KLSHMAT),
*. eigenvalues in WORK(KLVCC1)
C     DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
        CALL DIAG_SYMMAT_EISPACK(WORK(KLSHMAT),WORK(KLVCC1),
     &       WORK(KLVCC2),NNONSING,IRETURN)
C       IF(NTEST.GE.100) THEN 
          WRITE(6,*) ' Eigenvalues of H0 - E0*1 '
          CALL WRTMAT(WORK(KLVCC1),1,NNONSING,1,NNONSING)
C       END IF
*. Transform RHS to eigenvector basis and store in WORK(KLVCC2)
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
        CALL MATVCC(WORK(KLSHMAT),WORK(KLRHS),WORK(KLVCC2),NNONSING,
     &              NNONSING,1)
*. And divide with eigenvalues - with check for singularities 
        THRES = 1.0D-10
        NSING = 0
        DO I = 1, NNONSING
         IF(ABS(WORK(KLVCC1)).GT.THRES) THEN
          WORK(KLVCC3-1+I) = WORK(KLVCC2-1+I)/WORK(KLVCC1-1+I) 
         ELSE 
          NSING = NSING + 1
          WORK(KLVCC3-1+I) = 0.0D0
         END IF
        END DO
*. and remember the - : !1> = -(H0-E0)**-1 V |0>
        ONEM = -1.0D0
        CALL SCALVE(WORK(KLVCC3),ONEM,NNONSING)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' First order correction in eigenvector basis '
          CALL WRTMAT(WORK(KLVCC3),1,NNONSING,1,NNONSING)
        END IF
        WRITE(6,*) ' Number of encountered singularities ', NSING
*. And transform to orthonormal basis 
        CALL MATVCC(WORK(KLSHMAT),WORK(KLVCC3),WORK(KLC1),NNONSING,
     &              NNONSING,0)
        IF(NTEST.GE.100) THEN
         WRITE(6,*) ' First-order correction in orthonormal basis '
         CALL WRTMAT(WORK(KLC1),1,NNONSING,1,NNONSING)
        END IF
*. And obtain energy corrections
*. E2 = <0!V!1> = <0!H|1>
        E2 = INPROD(WORK(KLVCC2),WORK(KLVCC3),NNONSING)
        WRITE(6,*) ' Second order energy correction ', E2
        WRITE(6,*) ' Second order approximation to energy ',
     &             EREF+E2+ECORE
        E2TOT = EREF + E2 + ECORE
        EFINAL = E2TOT
*
        IF(ICTYP(1:5).EQ.'ICPT3') THEN
* Obtain also 3'rd order energy = <1!V-E1!1> = <1!H-F-E1!1>
*. transform first order correction to original SPA basis
          CALL VEC_FROM_DISC(WORK(KLXMAT),NSPAM1*NNONSING,1,-1,LU28)
          CALL MATVCC(WORK(KLXMAT),WORK(KLC1),WORK(KLVCC1),
     &                NSPAM1,NNONSING,0)
*. Insert a zero at the place of the unit-operator 
          WORK(KLVCC1-1+NSPA) = 0.0D0
*. And transform first order correction to CAAB basis 
          IF(I_DO_SPA.EQ.1) THEN
*. From SPA basis to CAAB basis
           CALL REF_CCV_CAAB_SP(WORK(KLC1O),WORK(KLVCC1),WORK(KLVCC3),2)
          ELSE
*. From zero-order to CAAB basis
C  TRANS_CAAB_ORTN(T_CAAB,T_ORTN,ITSYM,ICO,ILR,SCR,ICOCON)
           CALL TRANS_CAAB_ORTN(WORK(KLC1O),WORK(KLVCC1),1,2,2,
     &                          WORK(KLVCC3),2)
          END IF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' First order correction in CAAB basis '
            CALL WRTMAT(WORK(KLC1O),1,NCAAB,1,NCAAB)
          END IF
*. Modify one-electron integrals to h - f
*. (remember that f is in KINT1 and h is in KFIFA ...
          ONE = 1.0D0
          CALL VECSUM(WORK(KINT1),WORK(KFIFA),WORK(KINT1),
     &    ONE,ONEM,NINT1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' h - f 1-e operator '
            CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,1)
          END IF
          I12 = 2
*. And calculate <1|V|1> and <1!1> 
C     GET_IC_EXPECT(EXPVAL,IREFSPC,ITREFSPC,
C    &                        OP1,OP2,VEC1,VEC2)
*. 
          ECORE_SAVE = ECORE
          ECORE = 0.0D0
          CALL GET_IC_EXPECT(EXPVAL,OVLAP,IREFSPC,ITREFSPC,WORK(KLC1O),
     &                       WORK(KLC1O),WORK(KVEC1),WORK(KVEC2),
     &                       WORK(KLVCC1))
          ECORE = ECORE_SAVE
          E1 = EREF - E0
          E3 = EXPVAL - E1*OVLAP
          E3TOT = EREF + E2 + E3 + ECORE
          EFINAL = E3TOT
          WRITE(6,*) ' <1!V!1> = ', EXPVAL
          WRITE(6,*) ' <1|1>   = ', OVLAP
          WRITE(6,*) ' Third  order energy correction ', E3
          WRITE(6,*) ' Third order approximation to energy ',
     &               EREF+E2+E3+ECORE
        END IF
*. Report back to LUCIA
*. No iterative procedure, so
        CONVER = .TRUE.
        VFINAL = 0.0D0
      ELSE 
*
*. Use iterative method to solve first order equations
*
* 
* =======================================================
* Set up RHS of first-order equations = <0!H P T_{\mu}!0>
* =======================================================
*
        I12 = 2
        IPERTOP = 0
*. Use one-electron operator with inactive and ph contributions
        CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
        CALL GET_ICPT_RHS1(WORK(KLRHS),IREFSPC,ITREFSPC,
     &                     NSPA,NCAAB,I_DO_EI,
     &                     WORK(KVEC1),WORK(KVEC2),
     &                     WORK(KLVCC1),WORK(KLVCC2)        )
*. Make FIFA the one-body-hamiltonian
        CALL COPVEC(WORK(KFIFA),WORK(KINT1),NINT1)
*. And tell CI only to work with one-electron operator
        I12 = 1
*. Prepare for solution of first-order eqs by iterative techniques
*
*. The statement below is dirty, and I hope it will 
* not give me trouble in the future. The deal is that 
* the last operator ( in spinadapted order !!) is the 
* unit operator, and this is excluded from the 
* first order operator manifold, so...
        NVAR = NSPA - 1
*. Diagonal preconditioner, unit vector or diagonal of H0
        I_CALC_DIAG = 1
        IF(I_CALC_DIAG.EQ.1) THEN
          IF(I_DO_EI.EQ.1) THEN
C           GET_DIAG_H0_EI(DIAG,I_IN_TP)
            CALL GET_DIAG_H0_EI(WORK(KLVCC1))
*. The last element in KLDIA is the zero-order energy(without core)
            E0_FROMDIAG = WORK(KLVCC1-1+N_ZERO_EI)
            IF(NTEST.GE.10)
     &      WRITE(6,*) ' Zero-order energy from diag (with ecore)', 
     &      E0_FROMDIAG
            DO I = 1, NVAR
              WORK(KLVCC1-1+I) = WORK(KLVCC1-1+I) - E0_FROMDIAG
            END DO
          ELSE
            STOP ' Diagonal only programmed for EI-approach'
          END IF
        ELSE
          ONE = 1.0D0
          CALL SETVEC(WORK(KLVCC1),ONE,NVAR)
        END IF
        CALL VEC_TO_DISC(WORK(KLVCC1),NVAR,1,-1,LUSC53)
*. Initial guess - zero - to LUSC54
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLVCC1),ZERO,NVAR)
        CALL VEC_TO_DISC(WORK(KLVCC1),NVAR,1,-1,LUSC54)
*. And right hand side to LUSC37 
C       WRITE(6,*) ' RHS before written to DISC '
C       CALL WRTMAT(WORK(KLRHS),1,NVAR,1,NVAR)
        CALL VEC_TO_DISC(WORK(KLRHS),NVAR,1,-1,LUSC37)
*
        THRESH = 1.0D-8
        MAXITL = MAXIT
        MAXIT_MACRO = MAXITM
        WRITE(6,*) ' MAXITL, MAXIT_MACRO =', MAXITL, MAXIT_MACRO
*
        CALL MEMMAN(KLERROR,MAXITL+1,'ADDL  ',2,'ERROR ')
* The 0's in H0ME0TV are zeros ...
        NTESTL = 10
*. For communicating zero-order energy to routine for 
*. H0 - E0 * v
C       SHIFT = -E0
        SHIFT = -E0_FROMDIAG
*
        NTESTL =   3
        DO IMIC = 1, MAXIT_MACRO
*. Put RHS back on file
          IF(IMIC.NE.1) CALL VEC_TO_DISC(WORK(KLRHS),NVAR,1,-1,LUSC37)
          CALL MICGCG(H0ME0TV_EXT_IC,LUSC54,LUSC37,LUSC38,LUSC39,LUSC40,
     &                LUSC53,WORK(KLVCC1),WORK(KLVCC2),MAXITL,CONVER,
     &                THRESH,ZERO,WORK(KLERROR),NVAR,0,0,VFINAL,NTESTL)
C  MICGCG(MV8,LU1,LU2,LU3,LU4,LU5,LUDIA,VEC1,VEC2,
C    &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
C    &                  LUPROJ,LUPROJ2,VFINAL,IPRT)
          IF(CONVER) GOTO 1001
        END DO
 1001   CONTINUE
*. The solution to the first-order eqs, without a minus, is now 
*. on LUSC54
        CALL VEC_FROM_DISC(WORK(KLVCC1),NVAR,1,-1,LUSC54)
*. and add a minus to obtain the first-order corrections
        ONEM = -1.0D0
        CALL SCALVE(WORK(KLVCC1),ONEM,NVAR)
*. E2 = <0!V!1> = <0!H|1>
        E2 = INPROD(WORK(KLVCC1),WORK(KLRHS),NVAR)
        WRITE(6,*) ' Second order energy correction ', E2
        WRITE(6,*) ' Second order approximation to energy ',
     &             EREF+E2
        EFINAL = EREF+E2
        IF(ICTYP(1:5).EQ.'ICPT3') THEN
* Obtain also 3'rd order energy = <1!V-E1!1> = <1!H-F-E1!1>
*. Insert a zero at the place of the unit-operator 
          WORK(KLVCC1-1+NSPA) = 0.0D0
*. And transform first order correction to CAAB basis 
C   REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
          IF(I_DO_EI.EQ.0) THEN
          CALL REF_CCV_CAAB_SP(WORK(KLC1O),WORK(KLVCC1),WORK(KLVCC3),2)
          ELSE
            CALL TRANS_CAAB_ORTN(WORK(KLC1O),WORK(KLVCC1),1,2,2,
     &                           WORK(KLVCC3),2)
          END IF

          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' First order correction in CAAB basis '
            CALL WRTMAT(WORK(KLC1O),1,NCAAB,1,NCAAB)
          END IF
C         IF(NTEST.GE.2) THEN
C           WRITE(6,*) ' Analysis of first-order correction'
C           CALL ANA_GENCC(WORK(KLC1O),1)
C         END IF
*. Modify one-electron integrals to h - f
*. (remember that f is in KINT1 and h is in KFIFA ...
          ONE = 1.0D0
          CALL VECSUM(WORK(KINT1),WORK(KFIFA),WORK(KINT1),
     &    ONE,ONEM,NINT1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' h - f 1-e operator '
            CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,1)
          END IF
          I12 = 2
*. And calculate <1|V|1> and <1!1> 
C     GET_IC_EXPECT(EXPVAL,IREFSPC,ITREFSPC,
C    &                        OP1,OP2,VEC1,VEC2)
          ECORE_SAVE = ECORE
          ECORE = 0.0D0
          CALL GET_IC_EXPECT(EXPVAL,OVLAP,IREFSPC,ITREFSPC,WORK(KLC1O),
     &                       WORK(KLC1O),WORK(KVEC1),WORK(KVEC2),
     &                       WORK(KLVCC1))
          ECORE = ECORE_SAVE
          E1 = EREF -ECORE - E0
          E3 = EXPVAL - E1*OVLAP
          WRITE(6,*) ' <1!V!1> = ', EXPVAL
          WRITE(6,*) ' <1|1>   = ', OVLAP
          WRITE(6,*) ' Third  order energy correction ', E3
          WRITE(6,*) ' Third order approximation to energy ',
     &               EREF+E2+E3
          EFINAL = EREF+E2+E3
*
          IF(NTEST.GE.2) THEN
            WRITE(6,*) ' Analysis of first-order correction'
            CALL ANA_GENCC(WORK(KLC1O),1)
          END IF
        END IF
      END IF
*     ^ End of swith between direct and iterative method for solving
*     first order eqs.
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICPT  ')
      RETURN
      END 
      SUBROUTINE GET_ICPT_RHS1(RHS,IREFSPC,ITREFSPC,
     &                        NSPA,NCAAB,I_DO_EI,
     &                        VEC1,VEC2,VIC1,VIC2)
*
* Obtain RHS side vector for first order ICPT equations 
*
*. RHS_{\mu} = <0|T+_{\mu}PH|0> = <0!HP T_{\mu}|0>
*
* I_DO_EI = 1 => EI approach used, output vector is in zero-order basis
* I_DO_EI = 0 => SPA approach used, output vector is in SPA basis
*
*. Jeppe Olsen, August 2004
*               October 2009: I_DO_EI added
*

C     INCLUDE 'implicit.inc'
*. General input
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'clunit.inc'
*. Scratch for CI 
      DIMENSION VEC1(*),VEC2(*)
*. Scratch space for IC vectors 
      DIMENSION VIC1(*),VIC2(*)
*. Output 
      DIMENSION RHS(*)
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output form GET_ICPT_RHS1'
        WRITE(6,*) ' -------------------------'
        WRITE(6,*) ' I_DO_EI, NCAAB, NSPA =', I_DO_EI,NCAAB,NSPA
      END IF
*
* RHS will be calculated as density <L|T_{\mu}|0> 
* with |L> = P H|0>

*. Obtain H|0> on LUHC
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
C?    WRITE(6,*) ' Test : ICSPC, ISSPC = ', ICSPC,ISSPC
      CALL MV7(VEC1,VEC2,LUC,LUHC,0,0)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' H !Ref> as delivered in LUHC '
        CALL WRTVCD(VEC1,LUHC,1,-1)
      END IF
*. P H  !0> on LUHC 
      CALL REWINO(LUHC)
      CALL EXTR_CIV(IREFSM,ITREFSPC,LUHC,IREFSPC,2,
     &      LUSC1,-1,LUSC2,1,1,IDC,NTEST)
C          EXTR_CIV(ISM,ISPCIN,LUIN,
C    &               ISPCX,IEX_OR_DE,LUUT,LBLK,
C    &               LUSCR,NROOT,ICOPY,IDC,NTESTG)
      IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' P H !Ref> as delivered in LUHC '
           CALL WRTVCD(VEC1,LUHC,1,-1)
      END IF
*     <0!T+(I)P H  !0>  = <LUHC!T(I)!LUC> 
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
      ZERO = 0.0D0
      CALL SETVEC(VIC1,ZERO,NCAAB)
      CALL SIGDEN_CC(VEC1,VEC2,LUC,LUHC,VIC1,2)
      IF(I_DO_EI.EQ.0) THEN
        CALL REF_CCV_CAAB_SP(VIC1,RHS,VIC2,1)
      ELSE
C             TRANS_CAAB_ORTN(T_CAAB,T_ORTN,ITSYM,ICO,ILR,SCR,ICOCON)
         CALL TRANS_CAAB_ORTN(VIC1,RHS,1,1,2,VIC2,1)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' RHS for first order correction '
        CALL WRTMAT(RHS,1,NSPA,1,NSPA)
      END IF
*
      RETURN
      END 
      SUBROUTINE TRUNC_MAT(A,NRI,NCI,NRO,NCO)
*
* Truncate a matrix A by deleting some of the last rows and columns
*
*. Jeppe Olsen, Aug. 2004
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION A(*)
      IJO = 0
      DO ICO = 1, NCO
       DO  IRO = 1, NRO
         IJO = IJO + 1
         IJI = (ICO-1)*NRI+IRO
         A(IJO) = A(IJI)
       END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' truncated matrix '
        CALL WRTMAT(A,NRO,NCO,NRO,NCO)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_IC_EXPECT(EXPVAL,OVLAP,IREFSPC,ITREFSPC,
     &                        OP1,OP2,VEC1,VEC2,VIC1)
*. Obtain expectation value 
*   <0!O1+ H O2 |0>
* for two operators delivered in CAAB form 
* Jeppe Olsen, August 2004
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      REAL*8 INPRDD 
*. For communicating with routines below 
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'clunit.inc'
*. Input : Two operators in CAAB format 
      DIMENSION OP1(*),OP2(*)
*. Scratch for CI
      DIMENSION VEC1(*), VEC2(*)
*. and a vector of the size of the IC expansion
      DIMENSION VIC1(*)
*
*. 1 : Obtain Op2 |0> on LUSC1
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
      CALL SIGDEN_CC(VEC1,VEC2,LUC,LUSC1,OP2,1)
*. Obtain P Op2 !0> on LUSC1
      CALL REWINO(LUSC1)
      CALL EXTR_CIV(IREFSM,ITREFSPC,LUSC1,IREFSPC,2,
     &      LUSC2,-1,LUSC3,1,1,IDC,NTEST)
*. Obtain H P Op2 |0> on LUHC
      ICSPC = ITREFSPC
      ISSPC = ITREFSPC
      CALL REWINO(LUHC)
      CALL MV7(VEC1,VEC2,LUSC1,LUHC,0,0)
*
* Two ways to proceed.
*
      I_NEW_OR_OLD = 1
      IF(I_NEW_OR_OLD.EQ.2) THEN
*. Obtain Op1 |0> on LUSC2
        ICSPC = IREFSPC
        ISSPC = ITREFSPC
        CALL SIGDEN_CC(VEC1,VEC2,LUC,LUSC2,OP1,1)
*. Obtain P  Op1 |0> on LUSC2
        CALL REWINO(LUSC2)
        CALL EXTR_CIV(IREFSM,ITREFSPC,LUSC2,IREFSPC,2,
     &        LUSC3,-1,LUSC34,1,1,IDC,NTEST)
*. Obtain <O| Op1+P H P Op 2 |0> as inner product
        EXPVAL = INPRDD(VEC1,VEC2,LUHC,LUSC2,1,-1)
*. and the overlap  <O| Op1+P P Op 2 |0>
        OVLAP = INPRDD(VEC1,VEC2,LUSC1,LUSC2,1,-1)
      ELSE 
*. Op1 => Op1+ 
        CALL CONJ_CCAMP(OP1,1,VIC1)
        CALL CONJ_T
*. Op1+ P Op2 |0> on LUSC2
        ICSPC = ITREFSPC
        ISSPC = IREFSPC
        CALL REWINO(LUSC2)
        CALL REWINO(LUSC1)
        CALL SIGDEN_CC(VEC1,VEC2,LUSC1,LUSC2,VIC1,1)
        OVLAP = INPRDD(VEC1,VEC2,LUC,LUSC2,1,-1)
*. H P Op2 |0> => P H P Op 2|0>  on LUHC
        CALL REWINO(LUHC)
        CALL EXTR_CIV(IREFSM,ITREFSPC,LUHC,IREFSPC,2,
     &        LUSC3,-1,LUSC34,1,1,IDC,NTEST)
*  P H P Op 2|0> => Op1+  P H P Op 2|0>  on LUSC1
        CALL SIGDEN_CC(VEC1,VEC2,LUHC,LUSC1,VIC1,1)
        EXPVAL = INPRDD(VEC1,VEC2,LUSC1,LUC,1,-1)
*. And clean up
        CALL CONJ_T
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Expectation value <0|Op1+ P H P Op2 |0> = ', EXPVAL
        WRITE(6,*) ' Overlap           <0|Op1+ P   P Op2 |0> = ', OVLAP
      END IF

*
      RETURN
      END 
      SUBROUTINE H_S_EFF_ICCI_TV(VECIN,VECOUT_H,VECOUT_S,
     &           I_DO_H,I_DO_S)
*
* Obtain effective H and S- matrices (in reference space )
* times vector ( in reference space ) for given external CI 
* vector.
* if (I_DO_H.EQ.1) 
* vecout_h(i) = <i!(C_0 + T+ P) H (C_0 + P T)|in>, |in> = sum(j) vecin(j) |j>
*
* If (I_DO_S.EQ.1) 
* vecout_s(i) = <i!(C_0 + T+ P)   (C_0 + P T)|in>  
* 
* it is assumed that space for CI (Work(kvec1p) etc has been 
* defined .., and that common block COM_H_S_EFF_ICCI_TV
* has been initialized 
* Jeppe Olsen, Aug. 2004
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
*. Scratch units in use : LUHC, LUSC1, LUSC2, LUSC3, LUSC34,LUSC35
*. Transfer common 
      COMMON/COM_H_S_EFF_ICCI_TV/C_0,KLTOP,NREF,IREFSPC,ITREFSPC,NCAAB,
     &       IUNIOP,NSPA,IPROJSPC
     &
* C0 : Coefficient of reference function 
* KLTOP : Pointer to T vector in WORK
*         T is assumed to be in CAAB form 
* NREF  : Number of parameters in reference vector 
*. Input : Vector in refence space 
      DIMENSION VECIN(*)
*. And output, also a vector in reference space
      DIMENSION VECOUT_H(*)
      DIMENSION VECOUT_S(*)
*
      NTEST = 00
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HS_EFV')
      CALL MEMMAN(KL_REFV1,NREF,'ADDL  ',2,'REFV1 ')
      CALL MEMMAN(KL_ICV1,NCAAB,'ADDL  ',2,'ICV1  ')
*
C?    WRITE(6,*) ' Start of H_S .... '
*     
*
*. Transfer Vecin to discfile LUSC1 using the format of LUDIA
*
*. Use VECOUT_H to write integer list 1,2,3, ... NREF ( A bit unesthetic ..)
       CALL ISTVC2(VECOUT_H,0,1,NREF)
       CALL REWINO(LUSC1)
       CALL REWINO(LUDIA)
       CALL WRSVCD(LUSC1,-1,WORK(KVEC1P),VECOUT_H,VECIN,NREF,NREF,
     &             LUDIA,1)
C            WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
*. Obtain T !vecin> on LUSC2
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
      CALL REWINO(LUSC1)
      CALL REWINO(LUSC2)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUSC2,
     &               WORK(KLTOP),1)
*. T |vecin> => P T |vecin> on LUSC2
      CALL REWINO(LUSC2)
      CALL REWINO(LUSC3)
      CALL EXTR_CIV(IREFSM,ITREFSPC,LUSC2,IPROJSPC,2,
     &                    LUSC3,-1,LUSCR34,1,1,IDC,NTEST)
C     EXTR_CIV(ISM,ISPCIN,LUIN,
C    &                  ISPCX,IEX_OR_DE,LUUT,LBLK,
C    &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
*. Expand vecin from IREFSPC to ITREFSPC  on LUSC34
      CALL REWINO(LUSC1)
      CALL REWINO(LUSC34)
      CALL EXPCIV(IREFSM,IREFSPC,LUSC1,ITREFSPC,LUSC34,-1,
     /                 LUSC35,1,0,IDC,NTEST)
C               EXPCIV(ISM,ISPCIN,LUIN,
C     &                 ISPCUT,LUUT,LBLK,
C     &                 LUSCR,NROOT,ICOPY,IDC,NTESTG)
*. And add C_0 !vecin> to P T |Vecin>, save result on  LUSC1
      ONE = 1.0D0
C?    WRITE(6,*) ' The LUSC2 and LUSC34 files '
C?    CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
C?    CALL WRTVCD(WORK(KVEC1P),LUSC34,1,-1)
      CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),ONE,C_0,LUSC2,LUSC34,
     &            LUSC1,1,-1)
C              VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
*. Now we have ( C_0 + P T ) |vecin> on LUSC1
C?     WRITE(6,*) '(C_0 + P T) |Vecin> '
C?     CALL WRTVCD(WORK(KVEC1P),LUSC1,1,-1)
*
* ================
*. Overlap terms 
* ================
*
*. obtain ( C_0 + P T ) |vecin> in reference space on LUSC2, LUSC1 => LUSC2
C?    WRITE(6,*) ' Start of overlap terms '
      IF(I_DO_S.EQ.1) THEN
        CALL REWINO(LUSC1)
        CALL REWINO(LUSC2)
        CALL EXPCIV(IREFSM,ITREFSPC,LUSC1,IREFSPC,LUSC2,-1,
     /                  LUSC3,1,0,IDC,NTEST)
*.  ( C_0 + P T ) |vecin> => P  ( C_0 + P T ) |vecin>, LUSC1 => LUSC3
        CALL REWINO(LUSC1)
        CALL REWINO(LUSC3)
        CALL EXTR_CIV(IREFSM,ITREFSPC,LUSC1,IPROJSPC,2,
     &                      LUSC3,-1,LUSC34,1,0,IDC,NTEST)
*.  P  ( C_0 + P T ) |vecin> => T+  P  ( C_0 + P T ) |vecin>, LUSC3 => LUSC34
*. Conjugate T
        CALL CONJ_CCAMP(WORK(KLTOP),1,WORK(KL_ICV1))
        CALL CONJ_T
        CALL REWINO(LUSC3)
        CALL REWINO(LUSC34)
        ICSPC = ITREFSPC
        ISSPC = IREFSPC
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC34,
     &                 WORK(KL_ICV1),1)
*. C_0 ( C_0 + P T ) |vecin> + T+  P  ( C_0 + P T ) |vecin> on LUSC35
*. C_0 * LUSC2 + LUSC34 => LUSC35
        CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),C_0,ONE,LUSC2,LUSC34,
     &              LUSC35,1,-1)
*. And now read in form LUSC35
        CALL REWINO(LUSC35)
        CALL FRMDSCN(VECOUT_S,-1,-1,LUSC35)
      ELSE 
*. It was assumed that T => T+ was done in connection with overlap so
        CALL CONJ_CCAMP(WORK(KLTOP),1,WORK(KL_ICV1))
        CALL CONJ_T
      END IF
*
* ==========
*  H terms 
* ==========
C?    WRITE(6,*) ' Start of Hamilton terms '
*. (C_0 + P T ) |vecin> =>  H (C_0 + P T ) |vecin>, LUSC1 => LUHC
      IF(I_DO_H.EQ.1) THEN
        CALL REWINO(LUSC1)
        CALL REWINO(LUHC)
        ICSPC = ITREFSPC
        ISSPC = ITREFSPC
        CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUHC,0,0)
*. Obtain H!(C_0+PT)!vecin> in LUSC2, just in reference space 
*  (obtained by contracting from ITREFSPC to IREFSPC), LUHC => LUSC2
       CALL REWINO(LUHC)
       CALL REWINO(LUSC2)
       CALL EXPCIV(IREFSM,ITREFSPC,LUHC,IREFSPC,LUSC2,-1,
     /                  LUSC3,1,0,IDC,NTEST)
*. H (C_0 + P T) |vecin> => P H (C_0 + P T) |vecin>, LUHC  = > LUHC via LUSC1
       CALL REWINO(LUHC)
       CALL REWINO(LUSC1)
C?    WRITE(6,*) ' LUHC before call to EXTR_CIV '
C?    CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
       CALL EXTR_CIV(IREFSM,ITREFSPC,LUHC,IPROJSPC,2,
     &                     LUSC1,-1,LUSC3,1,1,IDC,NTEST)
C     EXTR_CIV(ISM,ISPCIN,LUIN,
C    &                  ISPCX,IEX_OR_DE,LUUT,LBLK,
C    &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
*. T => T+ operator have been done in overlap part 
*. P H (C_0 + P T) |vecin> => T+ P H (C_0 + P T) |vecin>, LUHC => LUSC1
       ICSPC = ITREFSPC
       ISSPC = IREFSPC
       CALL REWINO(LUHC)
       CALL REWINO(LUSC1)
       CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUHC,LUSC1,
     &                WORK(KL_ICV1),1)
*. Clean up, conjugate so we get the standard T operator back
       CALL CONJ_T
*. add C_O *  H!(C_0+PT)!vecin> and  T+ P H (C_0 + P T) |vecin>,
*   C_0 * LUSC2 + LUSC1 => LUSC3
       CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),C_0,ONE,LUSC2,LUSC1,
     &             LUSC3,1,-1)
C              VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
*. And we must now just read the result from LUSC3
C      FRMDSCN(VEC,NREC,LBLK,LU)
        CALL REWINO(LUSC3)
        CALL FRMDSCN(VECOUT_H,-1,-1,LUSC3)
      ELSE 
*. It is assumed that T-ops are conjugated back in the above so
        CALL CONJ_T
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vecin, Vecout_H, Vecout_S from  H_S_EFF_ICCI_...'
        CALL WRTMAT(VECIN,1,NREF,1,NREF)
        IF(I_DO_H.EQ.1) THEN
          WRITE(6,*) ' Vecout_H '
          CALL WRTMAT(VECOUT_H,1,NREF,1,NREF)
        END IF
        IF(I_DO_S.EQ.1) THEN
          WRITE(6,*) ' Vecout_S '
          CALL WRTMAT(VECOUT_S,1,NREF,1,NREF)
        END IF
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'HS_EFV')
*
      RETURN
      END 
      SUBROUTINE H_S_EXT_ICCI_TV(VECIN,VECOUT_H,VECOUT_S,
     &                           I_DO_H,I_DO_S)
*
* Obtain ICCI Hamiltonian and metric times vector, 
* external part
*
* If(I_DO_H.eq.1) vecout_h(i) :
*     <0!       H (V_0 + P sum_j vecin(j) O(j)) |0>, V_0 = vecin(iuniop)
*     <0!O+(i) PH (V_0 + P sum_j vecin(j) O(j))|0>
* if(I_DO_S.eq.1) vecout_s(i) : 
*     <0!         (V_0 + P sum_j vecin(j) O(j)) |0> = V_0 
*     <0!O+(i) P  (V_0 + P sum_j vecin(j) O(j))|0>
*
* <0!0> is assumed normalized
*
* Vecin is supposed to be delivered in SPA basis (if I_DO_EI = 0)
* or in the Zeroorder basis (if I_DO_EI = 1)
*
* Jeppe Olsen, August 2004
* I_DO_EI added, August 2009
*
      INCLUDE 'wrkspc.inc'
      REAL*8
     &INPRDD
      INCLUDE 'clunit.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc'
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT_H(*), VECOUT_S(*)
*. For transfer of data
      COMMON/COM_H_S_EFF_ICCI_TV/C_0,KLTOP,NREF,IREFSPC,ITREFSPC,NCAAB,
     &                           IUNIOP,NSPA,IPROJSPC
      NTEST = 5
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' H_S_EXT_ICCI_TV  entered '
      ELSE IF(NTEST.GE.10) THEN
        WRITE(6,*) '---------------------------------'
        WRITE(6,*) ' Reporting from  H_S_EXT_ICCI_TV '
        WRITE(6,*) '---------------------------------'
        WRITE(6,*)
        WRITE(6,*) ' NSPA, NCAAB = ', NSPA, NCAAB
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input vector '
        CALL WRTMAT(VECIN,1,NSPA,1,NSPA)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HSE_TV')
      CALL MEMMAN(KL_VIC1,NCAAB+1,'ADDL  ',2,'VIC1  ')
      CALL MEMMAN(KL_VIC2,NCAAB+1,'ADDL  ',2,'VIC2  ')
*
      IF(IUNIOP.NE.0) THEN 
        V_0 = VECIN(IUNIOP)
      ELSE 
        V_0 = 0.0D0
      END IF
C?    WRITE(6,*) ' IUNIOP = ', IUNIOP
*
* =======================================================
* 1 : Obtain  (V_0 + P sum_j vecin(j) O(j)) |0> on LUSC1
* =======================================================
* 
*. Reform VECIN to CAAB basis and store in WORK(KL_VIC1)
      
      IF(I_DO_EI.EQ.0) THEN
        CALL REF_CCV_CAAB_SP(WORK(KL_VIC1),VECIN,WORK(KL_VIC2),2) 
      ELSE
        CALL TRANS_CAAB_ORTN(WORK(KL_VIC1),VECIN,1,2,2,WORK(KL_VIC2),2)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input vector in CAAB basis '
        CALL WRTMAT(WORK(KL_VIC1),1,NCAAB,1,NCAAB)
      END IF
*. Obtain T !0> on LUSC2
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
      CALL REWINO(LUC)
      CALL REWINO(LUSC2)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUSC2,
     &               WORK(KL_VIC1),1)
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' T |0> '
        CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
      END IF
*. T |0> => P T |0> on LUSC2
      CALL REWINO(LUSC2)
      CALL REWINO(LUSC3)
C?    WRITE(6,*) ' IREFSM, ITREFSPC, IPROJSPC = ',
C?   &             IREFSM, ITREFSPC, IPROJSPC 
      CALL EXTR_CIV(IREFSM,ITREFSPC,LUSC2,IPROJSPC,2,
     &                    LUSC3,-1,LUSCR34,1,1,IDC,NTEST)
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' P T |0> '
        CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
      END IF
*. Expand |0>  from IREFSPC to ITREFSPC  on LUSC34
      CALL REWINO(LUC)
      CALL REWINO(LUSC34)
C?    WRITE(6,*) ' IREFSM, IREFSPC, ITREFSPC ',
C?   &             IREFSM, IREFSPC, ITREFSPC
      CALL EXPCIV(IREFSM,IREFSPC,LUC,ITREFSPC,LUSC34,-1,
     /                 LUSC35,1,0,IDC,NTEST)
*. And add V_0 !0> to P T |0>, save result on  LUSC1
      ONE = 1.0D0
      CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),ONE,V_0,LUSC2,LUSC34,
     &            LUSC1,1,-1)
*. We now we have ( V_0 + P T ) |0> on LUSC1
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) '(V_0 + P T) |0> '
       CALL WRTVCD(WORK(KVEC1P),LUSC1,1,-1)
      END IF
CM    CALL MEMCHK2('BEF_OVL')
*
* ================
*. Overlap terms 
* ================
*
*.  ( V_0 + P T ) |0> => P  ( V_0 + P T ) |0>, LUSC1 => LUSC3
      IF(I_DO_S.EQ.1) THEN
        CALL REWINO(LUSC1)
        CALL REWINO(LUSC3)
        CALL EXTR_CIV(IREFSM,ITREFSPC,LUSC1,IPROJSPC,2,
     &                      LUSC3,-1,LUSC34,1,0,IDC,NTEST)
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) 'P (V_0 + P T) |0> '
         CALL WRTVCD(WORK(KVEC1P),LUSC1,1,-1)
        END IF
CM    CALL MEMCHK2('AFT_EX')
*. Obtain density <0!O+(i)  P  ( V_0 + P T ) |0>
        ICSPC = IREFSPC
        ISSPC = ITREFSPC
        ZERO = 0.0D0
        CALL SETVEC(WORK(KL_VIC1),ZERO,NCAAB)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUSC3,
     &                 WORK(KL_VIC1),2)
*. Transfer to SPA or EI basis
        IF(I_DO_EI.EQ.1) THEN
C              TRANS_CAAB_ORTN(T_CAAB,T_ORTN,ITSYM,ICO,ILR,SCR,ICOCON)
          CALL TRANS_CAAB_ORTN(WORK(KL_VIC1),VECOUT_S,1,1,2,
     &                        WORK(KL_VIC2),1)
        ELSE
          CALL REF_CCV_CAAB_SP(WORK(KL_VIC1),VECOUT_S,WORK(KL_VIC2),1) 
        END IF
*. and the unit terms
      IF(IUNIOP.NE.0) VECOUT_S(IUNIOP) = V_0 
      END IF
CM    CALL MEMCHK2('AFT_OVL')
*
* ================
*. Hamilton  terms 
* ================
*
*. (V_0 + P T ) |0> =>  H (V_0 + P T ) |0>, LUSC1 => LUHC
      IF(I_DO_H.EQ.1) THEN
        CALL REWINO(LUSC1)
        CALL REWINO(LUHC)
        ICSPC = ITREFSPC
        ISSPC = ITREFSPC
        CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUHC,0,0)
*. Obtain H!(V_0+PT)!0> in LUSC2, just in reference space 
*  (obtained by contracting from ITREFSPC to IREFSPC), LUHC => LUSC2
        CALL REWINO(LUHC)
        CALL REWINO(LUSC2)
        CALL EXPCIV(IREFSM,ITREFSPC,LUHC,IREFSPC,LUSC2,-1,
     /                 LUSC3,1,0,IDC,NTEST)
*. Obtain    <0! H!(V_0+PT)!0>
        H_UNI = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUC,LUSC2,1,-1)
*. H (V_0 + P T) |0> => P H (V_0 + P T) |0>, LUHC  = > LUHC via LUSC1
        CALL REWINO(LUHC)
        CALL REWINO(LUSC1)
        CALL EXTR_CIV(IREFSM,ITREFSPC,LUHC,IPROJSPC,2,
     &                    LUSC1,-1,LUSC3,1,1,IDC,NTEST)
*. <LUHC!T(I)!LUC> 
        ICSPC = IREFSPC
        ISSPC = ITREFSPC
        ZERO = 0.0D0
        CALL SETVEC(WORK(KL_VIC1),ZERO,NCAAB)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,
     &                 WORK(KL_VIC1),2)
*. Transfer to SPA or EI basis
        IF(I_DO_EI.EQ.1) THEN
C             TRANS_CAAB_ORTN(T_CAAB,T_ORTN,ITSYM,ICO,ILR,SCR,ICOCON)
         CALL TRANS_CAAB_ORTN(WORK(KL_VIC1),VECOUT_H,1,1,2,
     &                        WORK(KL_VIC2),1)
        ELSE
         CALL REF_CCV_CAAB_SP(WORK(KL_VIC1),VECOUT_H,WORK(KL_VIC2),1) 
        END IF
        VECOUT_H(IUNIOP) = H_UNI
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Direct ICCI, external part '
        WRITE(6,*) ' Input vector '
        CALL WRTMAT(VECIN,1,NSPA,1,NSPA)
        IF(I_DO_H.EQ.1) THEN
          WRITE(6,*) ' H(ICCI) times input vector '
          CALL WRTMAT(VECOUT_H,1,NSPA,1,NSPA)
        END IF
        IF(I_DO_S.EQ.1) THEN
          WRITE(6,*) ' S(ICCI) times input vector '
          CALL WRTMAT(VECOUT_S,1,NSPA,1,NSPA)
        END IF
      END IF
* 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HSE_TV')
      RETURN
      END
      SUBROUTINE GET_HS_DIA(HDIA,SDIA,IDO_H,IDO_S,IFORM,
     &                      VCC1,VCC2,VEC1,VEC2,
     &                      IREFSPC,ITREFSPC,
     &                      IUNIOP,NSPA,IDOSUB,ISUB,NSUB)
*
* Obtain some form of Diagonal of H and S
*
*  IFORM = 1 : Obtain diagonal of Hamiltonian 
*  IFORM = 2 : Obtain diagonal of number-conserving part of H
*
* reference space on LUC
*
* If IDOPROJ = 1, then the reference space is projected out 
*                 for all operators except the unitoperator
*
* IF IDOSUB.NE.0, the matrix is constructed in the space 
* defined by the NSUB elements in ISUB
* NOTE : CODE HAS NOT BEEN TESTED FOR IDOSUB = 1 !!!!
*
* IDO_S = 1 => Diagonal of S is constructed 
* IDO_H = 1 => Diagonal of H is constructed 
*
* Jeppe Olsen, August 2004
*
*
      INCLUDE 'implicit.inc'
* 
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
*. Input
      INTEGER ISUB(*)
*. Output
      DIMENSION HDIA(*),SDIA(*) 
*. Scratch
      DIMENSION VCC1(*),VCC2(*)
      DIMENSION VEC1(*),VEC2(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HS_DIA')
      IUNIOP = NSPA
*     ^ Unit operator is assumed to be last operator
*       as it is in configuration ordered approach
*
      NTEST = 205
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' GET_HS_DIA speaking '
         WRITE(6,*) ' IDO_S, IDO_H, = ', IDO_S, IDO_H
         WRITE(6,*) ' IREFSPC, ITREFSPC = ', IREFSPC,ITREFSPC
      END IF
*. Number of excitations in calculation 
      NVAR = NSPA
*. Dimension of space in which S or H is constructed 
      IF(IDOSUB.EQ.0) THEN
        NSBVAR = NVAR
      ELSE
        NSBVAR = NSUB
      END IF
*
      IUNIOP_EFF = 0
      IF(IDOSUB.NE.0.AND.IUNIOP.NE.0) THEN
*. Check if unitoperator is included in list 
        CALL FIND_INTEGER_IN_VEC(IUNIOP,ISUB,NSUB,IUNIOP_EFF)
      ELSE IF(IUNIOP.NE.0) THEN
        IUNIOP_EFF = IUNIOP
      END IF
C?    WRITE(6,*) ' IUNIOP_EFF = ', IUNIOP_EFF
*.       
      IF(IFORM.EQ.1) THEN
*. Calculate Diagonal of H and S by calculating complete matrix ..
        WRITE(6,*) ' Complete matrix approach to obtaining diagonals'
        DO I = 1, NSBVAR
        IF(NTEST.GE.5) WRITE(6,*) 'Constructing row of S,H for I = ',I
        ZERO = 0.0D0
        CALL SETVEC(VCC1,ZERO,NVAR)
        IF(IDOSUB.EQ.0) THEN
          VCC1(I) = 1.0D0
        ELSE 
          VCC1(ISUB(I)) = 1.0D0
        END IF
*
*. Overlap terms 
*
        IF(IDO_S.EQ.1) THEN
          CALL H_S_EXT_ICCI_TV(VCC1,XDUM,VCC2,0,1)
          IF(IDOSUB.EQ.0) THEN
            SDIA(I) = VCC2(I)
          ELSE
            SDIA(I) = VCC2(ISUB(I))
          END IF
        END IF
*
*. Hamilton terms 
*
        IF(IDO_H.EQ.1) THEN
          CALL H_S_EXT_ICCI_TV(VCC1,VCC2,XDUM,1,0)
          IF(IDOSUB.EQ.0) THEN
            HDIA(I) = VCC2(I)
          ELSE
            HDIA(I) = VCC2(ISUB(I))
          END IF
        END IF
*
       END DO
      END IF
*     ^ Switch between various IFORMS
*
      IF(NTEST.GE.100) THEN
         IF(IDO_S.EQ.1) THEN
           WRITE(6,*) ' Diagonal of S '
           WRITE(6,*) ' ============== '
           CALL WRTMAT(SDIA,1,NSBVAR,1,NSBVAR)
         END IF
         IF(IDO_H.EQ.1) THEN
           WRITE(6,*) ' Diagonal of H '
           WRITE(6,*) ' =============='
           CALL WRTMAT(HDIA,1,NSBVAR,1,NSBVAR)
         END IF
       END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HS_DIA')
*
      RETURN
      END 
      SUBROUTINE HOME_SD_INV_T_ICCI(VECIN,VECOUT,E0,LUL1,LUL2)
*
* Obtain Inverted diagonal operator times ICCI vector
*
* VECOUT(I) = sum_j <0!O+i (sum_I |I><I|(H0-E0)!I><I| ) Oj |0> VECIN(J)
*
* Note that this does not correspond to the solution of the equations
*
* sum(j) <0!O+i(H0-E0)O j|0> Vecout(j) = Vecin(i)
*
* For getting better preconditioners ( without too much human labor)
*
*
* Vecin and Vecout are in (partial ) spinadapted basis 
* 
*
* Jeppe Olsen, Sept. 2004
      INCLUDE 'wrkspc.inc'
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      REAL*8
     &INPRDD, INPROD
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT(*)
*. Transfer block 
      COMMON/COM_H_S_EFF_ICCI_TV/C_0,KLTOP,NREF,IREFSPC,ITREFSPC,NCAAB,
     &                           IUNIOP,NSPA,IPROJSPC 
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HOMESD')
*. 2 vectors that should hold IC expansion in CAAB format
      CALL MEMMAN(KLVIC1,NCAAB,'ADDL  ',2,'VIC1  ')
      CALL MEMMAN(KLVIC2,NCAAB,'ADDL  ',2,'VIC2  ')
*
* Obtain VECIN in CAAB basis 
*
C      REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
       CALL REF_CCV_CAAB_SP(WORK(KLVIC1),VECIN,WORK(KLVIC2),2)
*
*. Obtain sum_j Vecin(j) O_j !0> in SD basis and save on LUL1
*
      ICSPC = IREFSPC 
      ISSPC = ITREFSPC
      CALL REWINO(LUC)
      CALL REWINO(LUL1)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUL1,
     &               WORK(KLVIC1),1)
*. Norm of assumed residual 
      X1NORM = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUL1,LUL1,1,-1)
*
*. And then Multiply LU1 with (H0-E0)**-1, save result on LUL2
*
      FACTOR = -1.0D0*E0
      CALL REWINO(LUL1)
      CALL REWINO(LUL2)
      CALL  DIA0TRM_GAS(2,LUL1,LUL2,WORK(KVEC1P),WORK(KVEC2P),FACTOR)
*. Norm of (H0-E0)**-1 * residual 
      X2NORM = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUL2,LUL2,1,-1)
      WRITE(6,*) ' Norm of residual and (H0-E0)**-1 * resid ',
     &            X1NORM, X2NORM

C          DIATRM(ITASK,LUIN,LUOUT,VECIN,VECOUT,FACTOR)
*. We are interested in <0!0+i (H0-E0)**-1(SD) O_j!0> Vecin(j) =
*.         <LUL2!O_i!LUC>
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
      CALL REWINO(LUC) 
      CALL REWINO(LUL2)
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLVIC1),ZERO,NCAAB)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUL2,WORK(KLVIC1),2)
      X4NORM = INPROD(WORK(KLVIC1),WORK(KLVIC1),NCAAB)
      WRITE(6,*) ' Norm (H0-E0)**-1 * resid in ICCI(CAAB) basis ',
     &            X4NORM
*. And reformat to SP basis
      CALL REF_CCV_CAAB_SP(WORK(KLVIC1),VECOUT,WORK(KLVIC2),1)
*. Norm of (H0-E0)**-1 * residual 
      X3NORM = INPROD(VECOUT,VECOUT,NSPA)
      WRITE(6,*) ' Norm (H0-E0)**-1 * resid in ICCI(SPA) basis ',
     &            X3NORM
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input and output vectors from  HOME_SD_INV_T_ICCI'
        CALL WRTMAT(VECIN ,1,NSPA,1,NSPA)
        CALL WRTMAT(VECOUT,1,NSPA,1,NSPA)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HOMESD')
*
      RETURN
      END 
      SUBROUTINE H0ME0TV_EXT_IC(VEC1,VEC2,LU1,LU2)
*
*. Obtain H0 - E0 * vector for external part in IC formalism
*
*. Jeppe Olsen, Sept. 2004
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
* Info from transfer arrays
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0,KLTOP,NREF,IREFSPC,ITREFSPC,NCAAB,
     &       IUNIOP,NSPA,IPROJSPC
      include 'cshift.inc'

*. 1 : Read input vector in from disc : Remember that
*      unit operator is excluded, so NVAR = NSPA - 1
*. Obtain H0 and S times vectors
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'H0ME0 ')
      CALL MEMMAN(KLVEC3,NSPA,'ADDL  ',2,'VEC3IC')
*
      CALL VEC_FROM_DISC(VEC1,NSPA-1,1,-1,LU1)
      VEC1(NSPA) = 0.0D0
*
      CALL  H_S_EXT_ICCI_TV(VEC1,VEC2,WORK(KLVEC3),1,1)
C           H_S_EXT_ICCI_TV(VECIN,VECOUT_H,VECOUT_S,
C    &                           I_DO_H,I_DO_S)
* H0 * v, S *v => (H0-E0S)*V
      ONE = 1.0D0
      CALL VECSUM(VEC2,VEC2,WORK(KLVEC3),ONE,SHIFT,NSPA-1)
      CALL VEC_TO_DISC(VEC2,NSPA-1,1,-1,LU2)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Input and output vectors from H0ME0_EXT_IC '
       CALL WRTMAT(VEC1,1,NSPA-1,1,NSPA-1)
       CALL WRTMAT(VEC2,1,NSPA-1,1,NSPA-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'H0ME0 ')
*
      RETURN
      END
      SUBROUTINE GET_SING_IN_SX_SPACE(IREFSPC)
*
* Analyze singularities in space of single-excitations
*
* Jeppe Olsen, Comfort Inn in Oak Ridge, Sept. 17 2004, 5 am (to be precise)
*.Continued Dec 2004 at Korshøjen before Warwick meeting( 30 hours to take-off)
*                    
*
*. It is assumed that spin-densities have been calculated 
*  for reference state - although spin-densities may 
* be recalculated here ...
*
      INCLUDE 'wrkspc.inc'
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc' 
      INCLUDE 'csm.inc'
*. Local list of single excitations, atmost 100 orbitals 
      INTEGER ISX(2,100*100)
*
      NTEST = 10
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_SI')
*
      I_CALC_DENS = 1
      IF(I_CALC_DENS.EQ.1) THEN
*
*. Space for CI behind the curtain 
           CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
           KVEC1P = KVEC1
           KVEC2P = KVEC2
*. Recalculate density matrices 
*. Should the densities be calculated with original CI-vector 
*. or projected density matrix ?
        I_DO_PROJ = 1
        IF(I_DO_PROJ.EQ.1) THEN
*. Project part of CI-vector belonging to IREFSPC - 1
          IPROJSPC = IREFSPC - 1
          WRITE(6,*) 
     &    ' Space to be projected out from reference ',IPROJSPC
          IF(IPROJSPC.EQ.0) THEN
            WRITE(6,*) ' No projection will be done '
            WRITE(6,*) ' As suggested projection space is undefined'
            LUPROJ = LUC
          ELSE
*. Project IPROJSPC out, save on LUHC
*. P T(I) !Ref> back on LUSCR
C                EXTR_CIV(ISM,ISPCIN,LUIN,
C    &                    ISPCX,IEX_OR_DE,LUUT,LBLK,
C    &                    LUSCR,NROOT,ICOPY,IDC,NTESTG)
            CALL EXTR_CIV(IREFSM,IREFSPC,LUC,IPROJSPC,2,
     &                    LUHC,-1,LUSC2,1,0,IDC,NTEST)
            LUPROJ = LUHC
          END IF
        ELSE
          LUPROJ = LUC
        END IF
*. And do the densities 
        ISPNDEN = 2
        CALL COPVCD(LUPROJ,LUSC2,WORK(KVEC1),1,-1)
        CALL DENSI2(IDENSI,WORK(KRHO1),WORK(KRHO2),
     &       WORK(KVEC1),WORK(KVEC2),LUPROJ,LUSC2,EXPS2,ISPNDEN,
     &       WORK(KSRHO1),WORK(KRHO2AA),WORK(KRHO2AB),WORK(KRHO2BB),1)
      END IF
*.    ^ End if densities should be recalculated ..
*
* ==========================================================
* Very simple first try, just diagonalize using no symmetry
* ==========================================================
*
*. Allocate space for two scratch matrices - each of length 2*NTOOB**2
*
      LEN = NTOOB**2
      LEN2 = NTOOB**4
      CALL MEMMAN(KLVEC1,2*LEN ,'ADDL  ',2,'LVEC1  ')
      CALL MEMMAN(KLVEC2,2*LEN ,'ADDL  ',2,'LVEC2  ')
      CALL MEMMAN(KLMAT1,4*LEN2,'ADDL  ',2,'LMAT1  ')
      CALL MEMMAN(KLMAT2,4*LEN2,'ADDL  ',2,'LMAT2  ')
      CALL MEMMAN(KLMAT3,4*LEN2,'ADDL  ',2,'LMAT3  ')
      CALL MEMMAN(KLISX, LEN   ,'ADDL  ',1,'ISX    ')
*
      I_DIAG_AAOP = 1
      IF(I_DIAG_AAOP.EQ.1) THEN 
*. Diagonalize space of double annihilations 
 
*
* 1 : Double annihilation operators 
*

*
*. Diagonalize RHO2AB
*
*. The form of RHO2AB is <0!a+ia a+kb alb aja!0> written as 
*. rho2(ik,lj) ik=(k-1)*NORB+i, lj=(j-1)*NORB+l.
*. with the addressing of ik and jl this is not 
*. an overlap matrix !! (It took me some hours to figure this out)
*. 
*. If we define operator lj ( with above def of lj) to be
*  alb aja!0>, then the conjugated operator (lj)+ is
*  <0!a+ja a+lb - which in rho2ab is given address jl.
*. so reorganize row indeces
*
      DO L = 1, NTOOB
      DO J = 1, NTOOB
       LJ_IN = (J-1)*NTOOB + L
       LJ_OUT = (L-1)*NTOOB + J
*. looping in the wrong direction, but this is not timedefining
       DO ICOL = 1, NTOOB**2
         WORK(KLMAT1-1+(ICOL-1)*LEN+LJ_OUT) =
     &   WORK(KRHO2AB-1+(ICOL-1)*LEN+LJ_IN)
       END DO
      END DO
      END DO
C     CALL COPVEC(WORK(KRHO2AB),WORK(KLMAT1),LEN2)
      WRITE(6,*) ' Info for diagonalization of RHO2AB '
      CALL CHK_S_FOR_SING(WORK(KLMAT1),LEN,NSING,WORK(KLMAT2),
     &                    WORK(KLVEC1),WORK(KLVEC2)           )
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The eigenvectors for zero-eigenvalues '
        CALL WRTMAT(WORK(KLMAT1),LEN,NSING,LEN,NSING)
      END IF
C     CHK_S_FOR_SING(S,NDIM,NSING,X,SCR,SCR2)
*
*. Diagonalize RHO2AA
*
      LENS = NTOOB*(NTOOB+1)/2
      LENS2 = LENS**2
      CALL COPVEC(WORK(KRHO2AA),WORK(KLMAT1),LENS2)
*. Actually RHO2SS are organized so they are  minus the overlap so
      ONEM = -1.0D0
      CALL SCALVE(WORK(KLMAT1),ONEM,LENS2)
      WRITE(6,*) ' Info for diagonalization of RHO2AA '
      CALL CHK_S_FOR_SING(WORK(KLMAT1),LENS,NSING,WORK(KLMAT2),
     &                    WORK(KLVEC1),WORK(KLVEC2)           )
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The eigenvectors for zero-eigenvalues '
        CALL WRTMAT(WORK(KLMAT1),LENS,NSING,LENS,NSING)
      END IF
*
*. Diagonalize RHO2BB
*
      LENS = NTOOB*(NTOOB+1)/2
      LENS2 = LENS**2
      CALL COPVEC(WORK(KRHO2BB),WORK(KLMAT1),LENS2)
*. Actually RHO2SS are organized so they are  minus the overlap so
      ONEM = -1.0D0
      CALL SCALVE(WORK(KLMAT1),ONEM,LENS2)
      WRITE(6,*) ' Info for diagonalization of RHO2BB '
      CALL CHK_S_FOR_SING(WORK(KLMAT1),LENS,NSING,WORK(KLMAT2),
     &                    WORK(KLVEC1),WORK(KLVEC2)           )
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The eigenvectors for zero-eigenvalues '
        CALL WRTMAT(WORK(KLMAT1),LENS,NSING,LENS,NSING)
      END IF
*
      END IF
*.    ^ End if double annihilations should be diagonalized
*
      I_DIAG_FULLSX = 0
*
* 2 : And the single excitation operators 
*
* a : MS = 1 operators : a+ia ajb
*
*. The overlap is S_ij,kl 
* = <0!(a+ia ajb)^+ (a+ka alb)!0> 
* = - <0!a+ka a+jb alb aia!0>  + delta(i,k)<0!a+jb alb!0>
* = -RHO2AB(kj,li) + delta(i,k)(RHO1(jl)-RHO1S(jl))/2
*
      DO I = 1, NTOOB
       DO J = 1, NTOOB
        DO K = 1, NTOOB
         DO L = 1, NTOOB
           KJ = (J-1)*NTOOB + K
           LI = (I-1)*NTOOB + L
           JL = (L-1)*NTOOB + J
           KJLI = (KJ-1)*NTOOB**2 + LI
           IJKL = (L-1)*NTOOB**3 + (K-1)*NTOOB**2 + (J-1)*NTOOB + I
           WORK(KLMAT1-1+IJKL) = -WORK(KRHO2AB-1+KJLI)
           IF(I.EQ.K) WORK(KLMAT1-1+IJKL) = WORK(KLMAT1-1+IJKL)
     &               +(WORK(KRHO1-1+JL)-WORK(KSRHO1-1+JL))/2
         END DO
        END DO
       END DO
      END DO
      CALL COPVEC(WORK(KLMAT1),WORK(KLMAT3),LEN*LEN)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The MS=1 SX metric '
        CALL WRTMAT(WORK(KLMAT1),LEN,LEN,LEN,LEN)
      END IF
      IF(I_DIAG_FULLSX.EQ.1) THEN
        WRITE(6,*) ' Info for diagonalization of metric of MS=1 SX '
        CALL CHK_S_FOR_SING(WORK(KLMAT1),LEN,NSING,WORK(KLMAT2),
     &                      WORK(KLVEC1),WORK(KLVEC2)           )
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 
     &    'The eigenvectors of zero-eigenvalues as NORB X NORB matrices'
            DO I = 1, NSING
              ILOFF = KLMAT1 + (I-1)*LEN
              CALL WRTMAT(WORK(ILOFF),NTOOB,NTOOB,NTOOB,NTOOB)
            END DO
        END IF
      END IF
*     ^ End if full space of SX should be diagonalized
*. Divide orbital excitations according to symmetry and 
*. diagonalize subblocks
      DO ISYM = 1, NSMST
C      DO IRANK = -1,1
       DO IRANK =  0,0
*. Obtain single excitations of this symmetry and rank 
C        GET_SX_FOR_SYM_AND_EXCRANK(ISYM_SX,IRANK2_SX,NSX,ISX)
         IRANK2 = 2*IRANK
         CALL GET_SX_FOR_SYM_AND_EXCRANK(ISYM,IRANK2,NSX,ISX)
*. Obtain matrix of excitations of this symmetry and rank
         DO IEX = 1, NSX
           DO JEX = 1, NSX
             IC = ISX(1,IEX)
             IA = ISX(2,IEX)
             JC = ISX(1,JEX)
             JA = ISX(2,JEX)
             IADR_IN = (JA-1)*NTOOB**3 + (JC-1)*NTOOB**2 
     /               + (IA-1)*NTOOB + IC
             IADR_OUT = (JEX-1)*NSX + IEX
             WORK(KLMAT1-1+IADR_OUT) = WORK(KLMAT3-1+IADR_IN)
           END DO
         END DO
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Metric for MS, SYM, RANK = ', 1,ISYM,IRANK2
           CALL WRTMAT(WORK(KLMAT1),NSX,NSX,NSX,NSX)
         END IF
         WRITE(6,*)  
     &   ' Info for diagonalization of metric of SX for MS,SYM,RANK ',
     &     1,ISYM,IRANK2
         CALL CHK_S_FOR_SING(WORK(KLMAT1),NSX,NSING,WORK(KLMAT2),
     &                       WORK(KLVEC1),WORK(KLVEC2)           )
         IF(NTEST.GE.10) THEN
           WRITE(6,*) 
     &     ' The eigenvectors for zero-eigenvalues'
           CALL WRTMAT(WORK(KLMAT1),NSX,NSING,NSX,NSING)
         END IF
       END DO
      END DO

*
* b : MS = -1 operators : a+ib aja
*
*. The overlap is S_ij,kl 
* = <0!(a+ib aja)^+ (a+kb ala)!0> 
* = - <0!a+ja a+kb aib ala!0>  + delta(i,k)<0!a+ja ala!0>
* = -RHO2AB(jk,il) + delta(i,k)(RHO1(jl)+RHO1S(jl))/2
*
      DO I = 1, NTOOB
       DO J = 1, NTOOB
        DO K = 1, NTOOB
         DO L = 1, NTOOB
           JK = (K-1)*NTOOB + J
           IL = (L-1)*NTOOB + I
           JL = (L-1)*NTOOB + J
           JKIL = (JK-1)*NTOOB**2 + IL
           IJKL = (L-1)*NTOOB**3 + (K-1)*NTOOB**2 + (J-1)*NTOOB + I
           WORK(KLMAT1-1+IJKL) = -WORK(KRHO2AB-1+JKIL)
           IF(I.EQ.K) WORK(KLMAT1-1+IJKL) = WORK(KLMAT1-1+IJKL)
     &               +(WORK(KRHO1-1+JL)+WORK(KSRHO1-1+JL))/2
         END DO
        END DO
       END DO
      END DO
      CALL COPVEC(WORK(KLMAT1),WORK(KLMAT3),LEN*LEN)
      IF(I_DIAG_FULLSX.EQ.1) THEN
       WRITE(6,*) ' Info for diagonalization of metric of MS=-1 SX '
       CALL CHK_S_FOR_SING(WORK(KLMAT1),LEN,NSING,WORK(KLMAT2),
     &                     WORK(KLVEC1),WORK(KLVEC2)           )
       IF(NTEST.GE.10) THEN
         WRITE(6,*) 
     &   ' Eigenvectors for zero-eigenvalues as NORB X NORB matrices'
           DO I = 1, NSING
             ILOFF = KLMAT1 + (I-1)*LEN
             CALL WRTMAT(WORK(ILOFF),NTOOB,NTOOB,NTOOB,NTOOB)
           END DO
       END IF
      END IF
*. Divide orbital excitations according to symmetry and 
*. diagonalize subblocks
      DO ISYM = 1, NSMST
C      DO IRANK = -1,1
       DO IRANK =  0,0
*. Obtain single excitations of this symmetry and rank 
C        GET_SX_FOR_SYM_AND_EXCRANK(ISYM_SX,IRANK2_SX,NSX,ISX)
         IRANK2 = 2*IRANK
         CALL GET_SX_FOR_SYM_AND_EXCRANK(ISYM,IRANK2,NSX,ISX)
*. Obtain matrix of excitations of this symmetry and rank
         DO IEX = 1, NSX
           DO JEX = 1, NSX
             IC = ISX(1,IEX)
             IA = ISX(2,IEX)
             JC = ISX(1,JEX)
             JA = ISX(2,JEX)
             IADR_IN = (JA-1)*NTOOB**3 + (JC-1)*NTOOB**2 
     /               + (IA-1)*NTOOB + IC
             IADR_OUT = (JEX-1)*NSX + IEX
             WORK(KLMAT1-1+IADR_OUT) = WORK(KLMAT3-1+IADR_IN)
           END DO
         END DO
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Metric for MS, SYM, RANK = ', -1,ISYM,IRANK2
           CALL WRTMAT(WORK(KLMAT1),NSX,NSX,NSX,NSX)
         END IF
         WRITE(6,*)  
     &   ' Info for diagonalization of metric of SX for MS,SYM,RANK ',
     &     -1,ISYM,IRANK2
         CALL CHK_S_FOR_SING(WORK(KLMAT1),NSX,NSING,WORK(KLMAT2),
     &                       WORK(KLVEC1),WORK(KLVEC2)           )
         IF(NTEST.GE.10) THEN
           WRITE(6,*) 
     &     ' The eigenvectors for zero-eigenvalues'
           CALL WRTMAT(WORK(KLMAT1),NSX,NSING,NSX,NSING)
         END IF
       END DO
      END DO
*
* MS = 0 
*
* There are two types of operators : a+ia aja and a+ib ajb
*
* This leads to a 2*NTOOB matrix 
* S_ij,kl = 
* (<0!(a+ia aja)^+ a+ka ala |0> | <0!(a+ia aja)^+ a+kb alb !0> )
* ( ----------------------------| -----------------------------)
* (<0!(a+ib ajb)^+ a+ka ala |0> | <0!(a+ib ajb)^ a+kb alb !0>  )
*
* The aaaa part 
*
*  <0!(a+ia aja)^+ a+ka ala |0>
*=-<0!a+ja a+ka aia ala!0> + delta(i,k) <0!a+ja ala!0>
*
      LEND = 2*NTOOB**2
      VALUE = -1234
      CALL SETVEC(WORK(KLMAT1),VALUE,LEND**2)
      DO I = 1, NTOOB
        DO J = 1, NTOOB
          DO K = 1, NTOOB
            DO L = 1, NTOOB
              IF(J.GT.K) THEN
                JK = J*(J-1)/2+K
                SIGN_JK =-1.0D0
              ELSE 
                JK = K*(K-1)/2 + J
                SIGN_JK = 1.0D0
              END IF
              IF(I.GT.L) THEN
                IL = I*(I-1)/2 + L
                SIGN_IL = -1.0D0
              ELSE
                IL = L*(L-1)/2 + I
                SIGN_IL =1.0D0
              END IF
              JKIL = (IL-1)*NTOOB*(NTOOB+1)/2 + JK
              IJKL = ((L-1)*NTOOB+K-1)*2*NTOOB**2 + (J-1)*NTOOB + I
              JL = (L-1)*NTOOB + J
              WORK(KLMAT1-1+IJKL) =-SIGN_JK*SIGN_IL*WORK(KRHO2AA-1+JKIL)
              IF(I.EQ.K)   WORK(KLMAT1-1+IJKL) =   WORK(KLMAT1-1+IJKL)
     &                   +(WORK(KRHO1-1+JL)+WORK(KSRHO1-1+JL))/2
            END DO
          END DO
        END DO
      END DO
*
* the bbbb part
*
*  <0!(a+ib ajb)^+ a+kb alb |0>
*=-<0!a+jb a+kb aib alb!0> + delta(i,k) <0!a+jb alb!0>
*
      DO I = 1, NTOOB
        DO J = 1, NTOOB
          DO K = 1, NTOOB
            DO L = 1, NTOOB
              IF(J.GT.K) THEN
                JK = J*(J-1)/2+K
                SIGN_JK = 1.0D0
              ELSE 
                JK = K*(K-1)/2 + J
                SIGN_JK = -1.0D0
              END IF
              IF(I.GT.L) THEN
                IL = I*(I-1)/2 + L
                SIGN_IL = 1.0D0
              ELSE
                IL = L*(L-1)/2 + I
                SIGN_IL =-1.0D0
              END IF
              JKIL = (IL-1)*NTOOB*(NTOOB+1)/2 + JK
              IJKL = ((L-1)*NTOOB+K-1+NTOOB**2 )*2*NTOOB**2 
     &             + (J-1)*NTOOB + I + NTOOB**2
              JL = (L-1)*NTOOB + J
              WORK(KLMAT1-1+IJKL) =-SIGN_JK*SIGN_IL*WORK(KRHO2BB-1+JKIL)
              IF(I.EQ.K)   WORK(KLMAT1-1+IJKL) =   WORK(KLMAT1-1+IJKL)
     &                   +(WORK(KRHO1-1+JL)-WORK(KSRHO1-1+JL))/2
            END DO
          END DO
        END DO
      END DO
*
* the aabb and bbaa part 
*
* S_ijkl(aabb) = <0!a+ja a+kb alb aia!0>
* S_ijkl(bbaa) = S_klij(aabb)
*
      DO I = 1, NTOOB
        DO J = 1, NTOOB
          DO K = 1, NTOOB
            DO L = 1, NTOOB
              JKLI = (I-1)*NTOOB**3 + (L-1)*NTOOB**2 + (K-1)*NTOOB + J
              IJKL = ((L-1)*NTOOB+K-1+NTOOB**2)*2*NTOOB**2
     &             + (J-1)*NTOOB + I
              WORK(KLMAT1-1+IJKL) = WORK(KRHO2AB-1+JKLI)
              KLIJ = ((J-1)*NTOOB + I-1)*2*NTOOB**2
     /             +  (L-1)*NTOOB + K + NTOOB**2
              WORK(KLMAT1-1+KLIJ) = WORK(KLMAT1-1+IJKL)
            END DO
          END DO
        END DO
      END DO
*
      LEND = 2*NTOOB**2
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The metric for MS = 0 '
        CALL WRTMAT(WORK(KLMAT1),LEND,LEND,LEND,LEND)
      END IF
*
      CALL COPVEC(WORK(KLMAT1),WORK(KLMAT3),LEND*LEND)
      IF(I_DIAG_FULLSX.EQ.1) THEN
       WRITE(6,*) ' Info for diagonalization of metric of MS = 0 SX '
       CALL CHK_S_FOR_SING(WORK(KLMAT1),LEND,NSING,WORK(KLMAT2),
     &                     WORK(KLVEC1),WORK(KLVEC2)           )
        IF(NTEST.GE.10) THEN
         WRITE(6,*) 
     &   ' Eigenvectors for zero-eigenvalues as 2 NORB X NORB matrices'
           DO I = 1, NSING
             ILOFF = KLMAT1 + (I-1)*LEND
             CALL WRTMAT(WORK(ILOFF),NTOOB,NTOOB,NTOOB,NTOOB)
             ILOFF = KLMAT1 + (I-1)*LEND + LEN
             CALL WRTMAT(WORK(ILOFF),NTOOB,NTOOB,NTOOB,NTOOB)
           END DO
        END IF
      END IF
*     ^ End if diag should be performed in full space
*. Divide orbital excitations according to symmetry and 
*. diagonalize subblocks
      DO ISYM = 1, NSMST
C      DO IRANK = -1,1
       DO IRANK =  0,0
*. Obtain single excitations of this symmetry and rank 
C        GET_SX_FOR_SYM_AND_EXCRANK(ISYM_SX,IRANK2_SX,NSX,ISX)
         IRANK2 = 2*IRANK
         CALL GET_SX_FOR_SYM_AND_EXCRANK(ISYM,IRANK2,NSX,ISX)
*. Obtain matrix of excitations of this symmetry and rank
         DO IEX = 1, NSX
           DO JEX = 1, NSX
             IC = ISX(1,IEX)
             IA = ISX(2,IEX)
             JC = ISX(1,JEX)
             JA = ISX(2,JEX)
*.aaaa
             IADR_IN = ((JA-1)*NTOOB+JC-1)*2*NTOOB**2 
     /               +  (IA-1)*NTOOB+IC
             IADR_OUT = (JEX-1)*2*NSX + IEX
             WORK(KLMAT1-1+IADR_OUT) = WORK(KLMAT3-1+IADR_IN)
*.aabb
             IADR_IN = ((JA-1)*NTOOB+JC+NTOOB**2-1)*2*NTOOB**2 
     /               +  (IA-1)*NTOOB+IC 
             IADR_OUT = (JEX+NSX-1)*2*NSX + IEX
             WORK(KLMAT1-1+IADR_OUT) = WORK(KLMAT3-1+IADR_IN)
*.bbaa
             IADR_IN = ((JA-1)*NTOOB+JC-1)*2*NTOOB**2 
     /               +  (IA-1)*NTOOB+IC + NTOOB**2
             IADR_OUT = (JEX-1)*2*NSX + IEX+NSX
             WORK(KLMAT1-1+IADR_OUT) = WORK(KLMAT3-1+IADR_IN)
*.bbbb
             IADR_IN = ((JA-1)*NTOOB+JC+NTOOB**2-1)*2*NTOOB**2 
     /               +  (IA-1)*NTOOB+IC + NTOOB**2
             IADR_OUT = (JEX+NSX-1)*2*NSX + IEX+NSX
             WORK(KLMAT1-1+IADR_OUT) = WORK(KLMAT3-1+IADR_IN)
           END DO
         END DO
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Metric for MS, SYM, RANK = ',  0,ISYM,IRANK2
           CALL WRTMAT(WORK(KLMAT1),2*NSX,2*NSX,2*NSX,2*NSX)
         END IF
         WRITE(6,*)  
     &   ' Info for diagonalization of metric of SX for MS,SYM,RANK ',
     &      0,ISYM,IRANK2
         CALL CHK_S_FOR_SING(WORK(KLMAT1),2*NSX,NSING,WORK(KLMAT2),
     &                         WORK(KLVEC1),WORK(KLVEC2)           )
         IF(NTEST.GE.10) THEN
           WRITE(6,*) 
     &   ' The eigenvectors for zero-eigenvalues'
           CALL WRTMAT(WORK(KLMAT1),2*NSX,NSING,2*NSX,NSING)
         END IF
       END DO
      END DO
* 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_SI')
*
      STOP ' Enforced stop in GET_SING_IN_SX_LIKE'
      RETURN
      END 
      SUBROUTINE MINGENEIG(MSTV,PRECTV,IPREC_FORM,THRES_E,THRES_R,
     &                  I_ER_CONV,
     &                  VEC1,VEC2,VEC3,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,
     &                  LU3,LU4,LU5,LUDIAM,LUDIAS,LUS,NROOT,MAXVEC,
     &                  NINVEC,
     &                  APROJ,AVEC,SPROJ,WORK,IPRT,EIGSHF,AVECP,
     &                  I_DO_PRECOND,CONVER,EFINAL,VFINAL)
*
* Iterative routine for generalized eigenvalue  problem
*
* M X = Lambda S X
*
* Version requiring 3 vectors in core 
*
* Jeppe Olsen Oct 2002 from MINDA4
*             Finished June 2004 at Korshoejen 53
*
* Input :
* =======
*        MSTV : Name of routine performing matrix*vector calculations
*        PRECTV : Name of precondtioner used if IPREC_FORM = 1
*        IPREC_FORM = 1 : use simple diagonal preconditioner 
*                   = 2 : Use external routine PRECTV to perform precond.
*        THRES_E: Convergence threshold for eigenvalue
*        THRES_R: Convergence threshold for residual norm
*        I_ER_CONV= 1 => Change in eigenvalue is used as conv. criterium
*                 = 2 => Norm or residual     is used as conv. criterium
*        LU1 : Initial set of vectors
*        VEC1,VEC2,VEC3 : Vectors,each must be dimensioned to hold
*                    complete vector
*        LU2,LU3   : Scatch files
*        LUDIAM    : File containing diagonal of matrix M
*        LUDIAS    : File containing diagonal of matrix S
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
* On input LU1 is supposed to hold initial guess to eigenvectors
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*), VEC3(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),SPROJ(*),WORK(*),AVECP(*)
*. Scratch files that may be used by matrix times vector
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3,
     &       LUCBIO_SAVE, LUHCBIO_SAVE, LUC_SAVE
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*MAXVEC
*      SPROJ  : MAXVEC*MAXVEC
*      AVEC   : MAXVEC*MAXVEC
*      WORK   : MAXVEC*MAXVEC
*      AVECP  : MAXVEC*MAXVEC
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(10)
* MSTV : Routine for matrix and metric times vector
* PRECTV : Routine for preconditioner times vector
       EXTERNAL MSTV, PRECTV
*
C?     WRITE(6,*) ' MINGENEIG: I_ER_CONV, THRES_E, THRES_R = ',
C?   &                         I_ER_CONV, THRES_E, THRES_R
       ONE = 1.0D0
       ZERO = 0.0D0
*. And the scratch files
       LUSCR1 = LU3
       LUSCR2 = LU4
       LUSCR3 = LU5
       LUCBIO_SAVE  = 0
       LUHCBIO_SAVE = 0
       LUC_SAVE = 0
*
*. Current code always reset to 2*NROOT so : 
       IF( MAXVEC .LT. 3 * NROOT ) THEN
         WRITE(6,*) ' SORRY MINGENEIG WOUNDED , MAXVEC .LT. 3*NROOT '
         STOP ' ENFORCED STOP IN MINGENEIG'
       END IF
*
       KFREE = 1
*
       KSSUB = 1
       KFREE = KFREE + MAXVEC*MAXVEC
*
       KMSUB = KFREE
       KFREE = KFREE + MAXVEC*MAXVEC
*
       KXORTN = KFREE
       KFREE = KFREE + MAXVEC*MAXVEC
*
       KSCRMAT = KFREE
       KFREE   = KFREE + MAXVEC*MAXVEC
*
       KSCRMAT2 = KFREE
       KFREE   = KFREE + MAXVEC*MAXVEC
*
       KVEC1 = KFREE
       KFREE = KFREE+ MAXVEC
*
       KVEC2 = KFREE
       KFREE = KFREE+ MAXVEC
       CONVER = .FALSE.
*
*.   INITAL ITERATION
*
       ITER = 1
*
       IPRT = 10000
       WRITE(6,*) 
     & ' MINGENEIG: IPRT, NVAR,MAXVEC  = ' , IPRT, NVAR, MAXVEC
       WRITE(6,'(A,I2,2(2X,E8.3))')
     & ' MINGENEIG: I_ER_CONV, THRES_E, THRES_R', 
     &              I_ER_CONV, THRES_E, THRES_R
       IF(IPRT.GE.200) THEN
        WRITE(6,*) ' Initial vectors in LU1 '
        CALL REWINO(LU1)
        DO IVEC = 1, NINVEC
         CALL WRTVCD(VEC1,LU1,1,-1)
        END DO
       END IF
       CALL GFLUSH(6)
*
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       CALL REWINO(LUS)
       WRITE(6,*) ' NVAR at start of MINGENEIG = ', NVAR
       DO  IVEC = 1,NINVEC
*. M and S times initial vector IVEC
         CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LU1)
         WRITE(6,*) ' Before MSTV '
         CALL MSTV(VEC1,VEC2,VEC3,1,1)
         WRITE(6,*) ' After MSTV ' 
         WRITE(6,*) ' NVAR, LU2, LUS = ', NVAR, LU2,LUS
*
         CALL VEC_TO_DISC(VEC2,NVAR,0,-1,LU2)
         CALL VEC_TO_DISC(VEC3,NVAR,0,-1,LUS)
* Update projected matrix 
         CALL REWINO(LU1)
         DO  JVEC = 1, IVEC
           CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LU1)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
           SPROJ(IJ) = INPROD(VEC1,VEC3,NVAR)
         END DO
       END DO
*
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' Initial matrix in subspace '
         CALL PRSYM(APROJ,NINVEC)
         WRITE(6,*) ' Initial metric in subspace '
         CALL PRSYM(SPROJ,NINVEC)
       END IF
*. Check for singularities in subspace matrix
C           TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       CALL TRIPAK(WORK(KSSUB),SPROJ,2,NINVEC,NINVEC)
C           GET_ON_BASIS(S,NVEC,NSING,X,SCRVEC1,SCRVEC2)
       CALL GET_ON_BASIS(WORK(KSSUB),NINVEC,NSING,WORK(KXORTN),
     &                   WORK(KVEC1),WORK(KVEC2))
       NNONSING = NINVEC - NSING
*. Transform Subspace M to orthonormal basis 
       CALL TRIPAK(WORK(KMSUB),APROJ,2,NINVEC,NINVEC)
       CALL COPVEC(WORK(KXORTN),WORK(KSCRMAT),
     &             NNONSING*NINVEC)
C           TRNMA_LM(XTAX,A,X,NRA,NCA,NRX,NCX,SCRVEC)
       CALL TRNMA_LM(WORK(KSCRMAT),WORK(KMSUB),WORK(KSCRMAT),
     &               NINVEC,NINVEC,NINVEC,NNONSING,WORK(KVEC1))
*. Transformed matrix is returved in KSCRMAT
        IF(IPRT.GE.20) THEN
          WRITE(6,*) ' NNONSING = ', NNONSING
          WRITE(6,*) ' Matrix in ON basis '
          CALL WRTMAT(WORK(KSCRMAT),NINVEC,NINVEC,NINVEC,NINVEC)
        END IF
*. Diagonalize transformed matrix
C            DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
       CALL DIAG_SYMMAT_EISPACK(WORK(KSCRMAT),WORK(KVEC1),
     &                          WORK(KVEC2),NNONSING,IRETURN)
C?     WRITE(6,*) ' Eigenvalues on return from DIAG_SYM .... '
C?     CALL WRTMAT(WORK(KVEC1),1,NNONSING,1,NNONSING)
*. Obtain the eigenvectors in the original basis 
       FACTORC = 0.0D0
       FACTORAB = 1.0D0
       CALL MATML7(AVEC,WORK(KXORTN),WORK(KSCRMAT),NINVEC,NNONSING,
     &             NINVEC,NNONSING,NNONSING,NNONSING,FACTORC,FACTORAB,0)
     &             
       DO IROOT = 1, NROOT 
         EIG(1,IROOT) = WORK(KVEC1-1+IROOT)
       END DO
*
       IF( IPRT  .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Initial set of eigenvalues '
         WRITE(6,'(5F22.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF),IROOT=1,NNONSING)
         WRITE(6,*) ' Initial subspace eigenvectors '
         CALL WRTMAT(AVEC,NINVEC,NROOT,NINVEC,NROOT)
       END IF
       NVEC = NINVEC
       NROOT_EFF = MIN(NROOT,NNONSING)
       IF(NNONSING.LT.NROOT) THEN
         WRITE(6,*) ' Linear dependencies in initial set of vectors '
         WRITE(6,*) ' NROOT, NNONSING = ', NROOT, NNONSING
         WRITE(6,*) ' Linear dependencies in initial set of vectors '
       END IF
*
      ITER_EFF = 1
      DO ITER = 2, MAXIT+1
        CALL GFLUSH(6)
*. In iteration MAXIT + 1, only the residuals are obtained ...
        IF(IPRT  .GE. 10 ) 
     &  WRITE(6,*) ' INFO FORM ITERATION .... ', ITER
*
** 1  New directions to be included
*
*   R = H*X - EIGAPR*S*X
        IADD = 0
        CONVER = .TRUE.
C?      WRITE(6,*) ' NROOT_EFF = ' , NROOT_EFF
        DO 100 IROOT = 1, NROOT_EFF
*. H*X in VEC3
C  MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),LU3,LU4,
     &                VEC1,VEC2,NVEC,1,-1)
          CALL VEC_FROM_DISC(VEC3,NVAR,1,-1,LU3)
          IF(IPRT.GE.600) THEN
            WRITE(6,*) ' MX '
            CALL WRTMAT(VEC3,1,NVAR,1,NVAR)
          END IF
*. S*X in VEC2
          CALL MVCSMD(LUS,AVEC((IROOT-1)*NVEC+1),LU3,LU4,
     &                VEC1,VEC2,NVEC,1,-1)
          CALL VEC_FROM_DISC(VEC2,NVAR,1,-1,LU3)
          IF(IPRT.GE.600) THEN
            WRITE(6,*) ' SX '
            CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
          END IF
*. MX - ESX in VEC1
          FACTOR = -EIG(ITER-1,IROOT)
          CALL VECSUM(VEC1,VEC3,VEC2,ONE,FACTOR,NVAR)
          IF ( IPRT  .GE.600 ) THEN
            WRITE(6,*) '  ( MX - ESX ) '
            CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
          END IF
          RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
          RNRM(ITER-1,IROOT) = RNORM
*  STRANGE PLACE TO TEST CONVERGENCE , BUT ....
          RTCNV(IROOT) = .FALSE.
          IF(I_ER_CONV.EQ.2) THEN
            IF(RNORM.LT. THRES_R) THEN
               RTCNV(IROOT) = .TRUE.
            ELSE
               RTCNV(IROOT) = .FALSE.
               CONVER = .FALSE.
            END IF
          ELSE 
           IF(ITER.EQ.2) THEN
              CONVER = . FALSE.
           ELSE 
            IF(ABS(EIG(ITER-1,IROOT)-EIG(ITER-2,IROOT)).LT.THRES_E)
     &      THEN
              RTCNV(IROOT) = .TRUE.
            ELSE
              RTCNV(IROOT) = .FALSE.
              CONVER = .FALSE.
            END IF
           END IF
          END IF
*
          IF(ITER.LE.MAXIT.AND. .NOT. RTCNV(IROOT) ) THEN
            IADD = IADD + 1
            IF(I_DO_PRECOND.EQ.1) THEN 
            IF(IPREC_FORM.EQ.1) THEN
*. Just use simple diagonal preconditioner
*.Multiply with diag(M-eig*S) to get new direction 
                CALL VEC_FROM_DISC(VEC2,NVAR,1,-1,LUDIAM)
                CALL VEC_FROM_DISC(VEC3,NVAR,1,-1,LUDIAS)
                FACTOR = -EIG(ITER-1,IROOT)
                CALL VECSUM(VEC2,VEC2,VEC3,ONE,FACTOR,NVAR)
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) ' Diagonal(M) - E*DIAG(S) '
                  CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
                  END IF
                CALL DIAVC2(VEC2,VEC1,VEC2,ZERO,NVAR)
C                    DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
                CALL COPVEC(VEC2,VEC1,NVAR)
                IF ( IPRT  .GE. 600) THEN
                  WRITE(6,*) '  (Diag(M)-E*Diag(S))-1 *( MX - ESX ) '
                  CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
                END IF
            ELSE 
*.  Perform more advanced preconditioning by using a 
*. external preconditionings routine 
               E = EIG(ITER-1,IROOT) + EIGSHF
               CALL PRECTV(VEC1,VEC2,E,LUDIAM,LUDIAS,VEC3)
               CALL COPVEC(VEC2,VEC1,NVAR)
            END IF
            END IF
*. VEC1 contains now new direction 
*. 1.3 ORTHOGONALIZE TO ALL PREVIOUS VECTORS
*. Should one use the S-metric or the standard metric?
*. I think one can argue for both. Therefore a swith here
*
            I_USE_1_OR_S = 2
            IF(I_USE_1_OR_S.EQ.1) THEN
              CALL COPVEC(VEC1,VEC2,NVAR)
            ELSE 
              WRITE(6,*) ' Before MSTV2'
              CALL MSTV(VEC1,VEC3,VEC2,0,1)
              WRITE(6,*) ' After MSTV2'
            END IF
            XNRMI = INPROD(VEC1,VEC2,NVAR)
            CALL REWINO( LU1 )
            DO IVEC = 1,NVEC+IADD-1
              CALL VEC_FROM_DISC(VEC3,NVAR,0,-1,LU1)
              OVLAP = INPROD(VEC3,VEC2,NVAR)
              CALL VECSUM(VEC1,VEC1,VEC3,1.0D0,-OVLAP,NVAR)
            END DO
*. 1.4 Normalize vector and check for linear dependency
            IF(I_USE_1_OR_S.EQ.1) THEN
              CALL COPVEC(VEC1,VEC2,NVAR)
            ELSE 
              WRITE(6,*) '  Before MSTV3'
              CALL MSTV(VEC1,VEC3,VEC2,0,1)
              WRITE(6,*) ' After MSTV3'
            END IF
            SCALE = INPROD(VEC1,VEC2,NVAR)
            IF(ABS(SCALE)/XNRMI .LT. 1.0D-10) THEN
*. Linear dependency
              IADD = IADD - 1
              IF ( IPRT  .GE. 10 ) 
     ¥        WRITE(6,*) '  Trial vector linear dependent so OUT !!'
            ELSE
              FACTOR = 1.0D0/SQRT(SCALE)
              CALL SCALVE(VEC1,FACTOR,NVAR)
              CALL VEC_TO_DISC(VEC1,NVAR,0,-1,LU1)

              IF ( IPRT  .GE.600 ) THEN
                WRITE(6,*) 
     &          ' Orthonormalized (Diag(M)-E*Diag(S))-1 *( MX - ESX ) '
                CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
              END IF
            END IF
*           ^ End if no singularity
          END IF
*         ^ End if this root was not converged
  100   CONTINUE
* 
        IF( CONVER ) GOTO  1001
*
**  2 : OPTIMAL COMBINATION OF NEW AND OLD DIRECTION
*
        IF(.NOT.CONVER.AND.ITER.LE.MAXIT) THEN 
          ITER_EFF = ITER_EFF + 1
*   Augment projected matrices
          CALL REWINO( LU1)
          CALL REWINO( LU2)
          CALL REWINO( LUS)
          DO IVEC = 1, NVEC
            CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LU1)
            CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LU2)
            CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LUS)
          END DO
*
          DO IVEC = 1, IADD
           CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LU1)
              WRITE(6,*) ' Before MSTV4'
           CALL MSTV(VEC1,VEC2,VEC3,1,1)
              WRITE(6,*) ' After MSTV4'
           CALL VEC_TO_DISC(VEC2,NVAR,0,-1,LU2)
           CALL VEC_TO_DISC(VEC3,NVAR,0,-1,LUS)
           CALL REWINO( LU1)
           DO JVEC = 1, NVEC+IVEC
             IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
             CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LU1)
             APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
             SPROJ(IJ) = INPROD(VEC1,VEC3,NVAR)
           END DO
          END DO
          IF(IPRT.GE.10) THEN
            WRITE(6,*) ' Subspace M and S matrices '
            CALL PRSYM(APROJ,NVEC+IADD)
            CALL PRSYM(SPROJ,NVEC+IADD)
          END IF
*
        I_DO_SYMTEST = 1
        IF(I_DO_SYMTEST.EQ.1) THEN
          WRITE(6,*) ' Symmetry of subspace matrices tested'
* Test: Construct complete subspace matrices without assuming
*       Hermiticity
          CALL REWINO(LU1)
          NVECA = NVEC + IADD
          DO IVEC = 1, NVECA
            CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LU1)
            CALL REWINO(LU2)
            CALL REWINO(LUS)
            DO JVEC = 1, NVECA
              IJ = (JVEC-1)*(NVECA) + IVEC
              CALL VEC_FROM_DISC(VEC2,NVAR,0,-1,LU2)
              WORK(KSCRMAT-1+IJ) = INPROD(VEC1,VEC2,NVAR)
              CALL VEC_FROM_DISC(VEC2,NVAR,0,-1,LUS)
              WORK(KSCRMAT2-1+IJ) = INPROD(VEC1,VEC2,NVAR)
            END DO
          END DO
          WRITE(6,*) ' Full A and S subspace matrices '
          CALL WRTMAT(WORK(KSCRMAT),NVECA,NVECA,NVECA,NVECA)
          WRITE(6,*)
          CALL WRTMAT(WORK(KSCRMAT2),NVECA,NVECA,NVECA,NVECA)
        END IF ! End if hermiticity of submatrices should be tested

       
    

*. Save the previous set of eigenvectors in AVECP
C              COPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
          CALL COPMT2(AVEC,AVECP,NVEC,NNONSING,NVEC+IADD,NNONSING,1)
*. We now have new subspace matrices, so diagonalize
          NVEC = NVEC + IADD
*. Check for singularities in subspace matrix
          ONE = 1.0D0
          CALL TRIPAK(WORK(KSSUB),SPROJ,2,NVEC,NVEC)
C?        WRITE(6,*) ' Projected S matrix in expanded form '
C?        CALL WRTMAT(WORK(KSSUB),NVEC,NVEC,NVEC,NVEC)
          CALL GET_ON_BASIS(WORK(KSSUB),NVEC,NSING,WORK(KXORTN),
     &                      WORK(KVEC1),WORK(KVEC2))
          NNONSING = NVEC - NSING
          IF(NNONSING.LT.NROOT) THEN
            WRITE(6,*) ' Number of roots in nonsing problem '
            WRITE(6,*) ' Is lower than the required number of roots'
            WRITE(6,*) NNONSING, NROOT
            WRITE(6,*) ' I expect trouble but will continue '
          END IF
*. Transform Subspace M to orthonormal basis 
          CALL TRIPAK(WORK(KMSUB),APROJ,2,NVEC,NVEC)
          CALL COPVEC(WORK(KXORTN),WORK(KSCRMAT),
     &                NNONSING*NVEC)
C              TRNMA_LM(XTAX,A,X,NRA,NCA,NRX,NCX,SCRVEC)
          CALL TRNMA_LM(WORK(KSCRMAT),WORK(KMSUB),WORK(KSCRMAT),
     &                  NVEC,NVEC,NVEC,NNONSING,WORK(KVEC1))
*. Transformed matrix is returved in KSCRMAT
*. Diagonalize transformed matrix
C              DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
          IF(IPRT.GE.20) THEN
            WRITE(6,*) ' Matrix in orthonormal basis '
            CALL WRTMAT(WORK(KSCRMAT),NNONSING,NNONSING,NNONSING,
     &                  NNONSING)
           END IF
          CALL DIAG_SYMMAT_EISPACK(WORK(KSCRMAT),WORK(KVEC1),
     &                          WORK(KVEC2),NNONSING,IRETURN)
*. Obtain the eigenvectors in the original basis 
          FACTORC = 0.0D0
          FACTORAB = 1.0D0
          CALL MATML7(AVEC,WORK(KXORTN),WORK(KSCRMAT),NVEC,NNONSING,
     &               NVEC,NNONSING,NNONSING,NNONSING,FACTORC,FACTORAB,0)
          DO IROOT = 1, NROOT
            EIG(ITER,IROOT) = WORK(KVEC1-1+IROOT)
          END DO
*
          IF(IPRT .GE. 3 ) THEN
            WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
            WRITE(6,'(5F22.13)')
     &      ( (EIG(ITER,IROOT)+EIGSHF) ,IROOT=1,NROOT)
          END IF
*
          IF( IPRT  .GE. 5 ) THEN
            WRITE(6,*) ' Projected M-and S-matrices'
            CALL PRSYM(APROJ,NVEC)
            CALL PRSYM(SPROJ,NVEC)
            WRITE(6,*) ' Subspace eigen-values and -vectors'
            WRITE(6,'(2X,E20.13)') 
     &      (EIG(ITER,IROOT)+EIGSHF,IROOT = 1, NROOT)
            CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
          END IF
        END IF
*       ^ End if not converged
*
**  Reset / Assemble current eigenvectors if 
*   space for another set of NROOT vectors is not possible 
        IF(NVEC+NROOT.GT.MAXVEC.AND..NOT.CONVER) THEN
*. Orthogonalize previous set of eigenvectors on current 
*. set using normal metric !
           CALL COPVEC(AVECP,AVEC(NROOT*NVEC+1),NROOT*NVEC)
           IF(IPRT.GE.20) THEN
             WRITE(6,*) ' Nonorthonormal basis for reset '
             CALL WRTMAT(AVEC,NVEC,2*NROOT,NVEC,2*NROOT)
           END IF
*. Overlap matrix of the 2*NROOT vectors : All vectors on file 
* are orthonormal, so overlap matrix is simple to obtain.
           CALL MATML7(WORK(KSCRMAT),AVEC,AVEC,2*NROOT,2*NROOT,
     &                 NVEC,2*NROOT,NVEC,2*NROOT,ZERO,ONE,1)
            IF(IPRT.GE.20) THEN
              WRITE(6,*) ' Overlap of nonorthonormal reset vecs '
              CALL WRTMAT(WORK(KSCRMAT),2*NROOT,2*NROOT,
     &                    2*NROOT,2*NROOT)
            END IF
*. Orthogonalize vectors by forward Gram-Schmidt diagonalization
           CALL MGS3(WORK(KSCRMAT2),WORK(KSCRMAT),2*NROOT,WORK(KVEC1))
           IF(IPRT.GE.20) THEN
             WRITE(6,*) ' Transformation matrix to orthonormal basis '
             CALL WRTMAT(WORK(KSCRMAT2),2*NROOT,2*NROOT,
     &                    2*NROOT,2*NROOT)
           END IF
*. In KSCRMAT2 we now have the expansion of the orthogonal 
*. eigenvectors in terms of the new and the previous eigenvectors.
*. Obtain the expansion of the orthogonal eigenvectors in terms of 
*. the vectors on disc
           CALL MATML7(AVECP,AVEC,WORK(KSCRMAT2),NVEC,2*NROOT,
     &                 NVEC,2*NROOT,2*NROOT,2*NROOT,ZERO,ONE,ZERO)
           CALL COPVEC(AVECP,AVEC,NVEC*2*NROOT)
           IF(IPRT.GE.20) THEN
             WRITE(6,*) ' Orthonormal basis for reset vectors '
             CALL WRTMAT(AVEC,NVEC,2*NROOT,NVEC,2*NROOT)
           END IF
*. Obtain the corresponding Vectors on Disc
*. The c-Vectors 
           CALL REWINO(LU3)
           DO IROOT = 1, 2*NROOT
             CALL MVECSUM(AVEC((IROOT-1)*NVEC+1),NVEC,NVAR,VEC1,VEC2,
     &                    LU1,1,1)
             CALL VEC_TO_DISC(VEC1,NVAR,0,-1,LU3)
           END DO
           CALL REWINO(LU3)
           CALL REWINO(LU1)
           DO IROOT = 1, 2*NROOT
             CALL COPVCD(LU3,LU1,VEC1,0,-1)
            END DO
*. and the sigma-vectors
           CALL REWINO(LU3)
           DO IROOT = 1, 2*NROOT
             CALL MVECSUM(AVEC((IROOT-1)*NVEC+1),NVEC,NVAR,VEC1,VEC2,
     &                    LU2,1,1)
             CALL VEC_TO_DISC(VEC1,NVAR,0,-1,LU3)
           END DO
           CALL REWINO(LU3)
           CALL REWINO(LU2)
           DO IROOT = 1, 2*NROOT
             CALL COPVCD(LU3,LU2,VEC1,0,-1)
            END DO
*. And the S-vectors 
           CALL REWINO(LU3)
           DO IROOT = 1, 2*NROOT
             CALL MVECSUM(AVEC((IROOT-1)*NVEC+1),NVEC,NVAR,VEC1,VEC2,
     &                    LUS,1,1)
             CALL VEC_TO_DISC(VEC1,NVAR,0,-1,LU3)
           END DO
           CALL REWINO(LU3)
           CALL REWINO(LUS)
           DO IROOT = 1, 2*NROOT
             CALL COPVCD(LU3,LUS,VEC1,0,-1)
            END DO
*
           IF(IPRT.GE.20) THEN
             WRITE(6,*) ' Reset set of 2*NROOT eigenvectors '
             CALL WRTMAT(AVEC,NVEC,2*NROOT,NVEC,2*NROOT)
           END IF
*. Subspace matrices for the new basis-vectors 
C     SUBSPC_MAT_FROM_VECTORS(LUV,LUAV,NVECP,NVEC,ASUB,
C    &           ISYM,VEC1,VEC2,NVAR)
           CALL SUBSPC_MAT_FROM_VECTORS(LU1,LU2,0,2*NROOT,APROJ,
     &          1,VEC1,VEC2,NVAR)
           CALL SUBSPC_MAT_FROM_VECTORS(LU1,LUS,0,2*NROOT,SPROJ,
     &          1,VEC1,VEC2,NVAR)
*
           NVEC = 2*NROOT
*. and reset the matrix defining the roots
           CALL SETVEC(AVEC,ZERO,NVEC**2)
           CALL SETDIA(AVEC,ONE,NVEC,0)
        END IF
*       ^ End if Reset was required
      END DO
*     ^ End of loop over iterations
 1001 CONTINUE
*     ^ Statement to which we skip if converged
*. Well, the last iteration was used to to construct the residual, 
*. and does therefore not really count so 
      ITER = ITER_EFF 
*
*. construct the first NROOT approximations to the 
*. eigenvectors on LU1 and the corresponding sigmavectors on LU2
*
*. The c-Vectors 
      CALL REWINO(LU3)
      DO IROOT = 1, NROOT
        CALL MVECSUM(AVEC((IROOT-1)*NVEC+1),NVEC,NVAR,VEC1,VEC2,
     &               LU1,1,1)
        CALL VEC_TO_DISC(VEC1,NVAR,0,-1,LU3)
      END DO
      CALL REWINO(LU3)
      CALL REWINO(LU1)
      DO IROOT = 1, NROOT
        CALL COPVCD(LU3,LU1,VEC1,0,-1)
      END DO
*. and the sigma-vectors
      CALL REWINO(LU3)
      DO IROOT = 1, NROOT
        CALL MVECSUM(AVEC((IROOT-1)*NVEC+1),NVEC,NVAR,VEC1,VEC2,
     &               LU2,1,1)
        CALL VEC_TO_DISC(VEC1,NVAR,0,-1,LU3)
      END DO
      CALL REWINO(LU3)
      CALL REWINO(LU2)
      DO IROOT = 1, NROOT
        CALL COPVCD(LU3,LU2,VEC1,0,-1)
      END DO
*. Obtain the Final C-vector in VEC1
      CALL VEC_FROM_DISC(VEC1,NVAR,1,-1,LU1)
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
C        ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
      END IF
*. Final eigenvalues
      DO IROOT = 1, NROOT
         FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
      END DO
*
      EFINAL = FINEIG(NROOT)
      VFINAL = RNRM(ITER,NROOT)
*
      IF ( IPRT .GT. 1 ) THEN
        DO IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,'(A,F18.10)')  
     &    ' The final approximation to eigenvalue ', FINEIG(IROOT)
          IF(IPRT.GE.1000) THEN
            WRITE(6,*) '  The final approximation to eigenvector'
            CALL WRTVCD(VEC1,LU1,1,-1)
          END IF
          WRITE(6,'(A)') ' Summary of iterations '
          WRITE(6,'(A)') ' ----------------------'
          WRITE(6,'(A)')
     &    ' Iteration point        Eigenvalue         Residual '
          DO I=1,ITER
            WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          END DO
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
        END DO
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),
     &                 RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
        END DO
      END IF
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE SUBSPC_MAT_FROM_VECTORS(LUV,LUAV,NVECP,NVEC,ASUB,
     &           ISYM,VEC1,VEC2,NVAR)
*
* Obtain subspace matrix from a set of vectors (on file LUV) and matrix times
* vectors ( on file LUAV)
*
*. Input 
*  LUV : file containing vectors 
*  LUAV: file containing matrix times vectors 
*  NVECP : Number of vectors for which subspace matrix already 
*          have been constructed 
*  NVEC   : Number of vectors
*  ISYM   : = 1 => matrix is symmetric, only lower half of ASUB
*           is calculated
*            =0 => matrix is not symmetric, complete SUB is obtained
*
*. Output
*  ASUB : Updated subspace matrix 
*
* Scratch
* ======
* VEC1, VEC2, Should be able to hold vectors 
*
* Jeppe Olsen, June 2004, trying to get back to work ....
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Output
      DIMENSION ASUB(*)
*. Scratch
      DIMENSION VEC1(NVAR),VEC2(NVAR)
*
      IF(ISYM.EQ.1) THEN
* Calculate A(i,j) = Vec(i)(T) A Vec(j) for i.le.j.
        CALL REWINO(LUV)
        DO I = 1, NVECP
         CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LUV)
        END DO
        DO I = NVECP+1,NVEC
          CALL VEC_FROM_DISC(VEC1,NVAR,0,-1,LUV)
*
          CALL REWINO(LUAV)
          DO J = 1, NVECP
            CALL VEC_FROM_DISC(VEC2,NVAR,0,-1,LUAV)
          END DO
          DO J = NVECP+1,I
            CALL VEC_FROM_DISC(VEC2,NVAR,0,-1,LUAV)
            IJ = I*(I-1)/2 + J
            ASUB(IJ) = INPROD(VEC1,VEC2,NVAR)
          END DO
        END DO
       ELSE 
          WRITE(6,*) ' Sorry ISYM = 0 option not yet implemented '
          STOP '  SUBSPC_MAT_FROM_VECTORS : ISYM = 0 not implemented '
       END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated subspace matrix '
        CALL PRSYM(ASUB,NVEC)
      END IF
*
      RETURN
      END  
      SUBROUTINE MTV_FUSK(VECIN,VECOUT)
*
* Fusk version of vector * matrix
*
      INCLUDE 'implicit.inc'
*
      PARAMETER(NDIM_FUSK = 4)
      DIMENSION A(NDIM_FUSK*NDIM_FUSK)
*
      DO I = 1, NDIM_FUSK ** 2 
       A(I) = 1.1D0
      END DO
      DO I = 1, NDIM_FUSK
        A((I-1)*NDIM_FUSK+I) = DBLE(I)
      END DO
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(A,VECIN,VECOUT,NDIM_FUSK,NDIM_FUSK,0)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input and output form MTV_FUSK '
        CALL WRTMAT(VECIN,1,NDIM_FUSK,1,NDIM_FUSK)
        CALL WRTMAT(VECOUT,1,NDIM_FUSK,1,NDIM_FUSK)
      END IF
*
      RETURN
      END 
      SUBROUTINE STV_FUSK(VECIN,VECOUT)
*
* Fusk version of Metric * vector 
*
      INCLUDE 'implicit.inc'
*
      PARAMETER(NDIM_FUSK = 4)
      DIMENSION S(NDIM_FUSK*NDIM_FUSK)
*
      DO I = 1, NDIM_FUSK ** 2 
       S(I) = 0.0D0
      END DO
      DO I = 1, NDIM_FUSK
        S((I-1)*NDIM_FUSK+I) = 1.0D0 + 0.1*FLOAT(I-1) 
      END DO
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(S,VECIN,VECOUT,NDIM_FUSK,NDIM_FUSK,0)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input and output form STV_FUSK '
        CALL WRTMAT(VECIN,1,NDIM_FUSK,1,NDIM_FUSK)
        CALL WRTMAT(VECOUT,1,NDIM_FUSK,1,NDIM_FUSK)
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_SX_FOR_SYM_AND_EXCRANK(ISYM_SX,IRANK2_SX,NSX,ISX)
*
* Obtain single excitations of given symmetry and excitation rank
* Orbital numbers are in TS order 
* IHPVGAS is used to decide excitation rank
*
*. Jeppe Olsen, Dec. 2004
*
      INCLUDE 'implicit.inc'
*. Input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*. Output : Creation and annihilation part of SX
      INTEGER ISX(2,*)
*
      NSX = 0
      DO ICOB = 1, NTOOB
        DO IAOB = 1, NTOOB
          ISYM = MULTD2H(ISMFTO(ICOB),ISMFTO(IAOB))
          IHPV_C = IHPVGAS(ITPFTO(ICOB))
          IHPV_A = IHPVGAS(ITPFTO(IAOB))
          IF(IHPV_C.EQ.1) THEN
*. Creation of hole, corresponds to deexcitaion
            IR_C = -1
          ELSE IF(IHPV_C.EQ.2) THEN
*. creation of particle, corresponds to excitation 
            IR_C = 1
          ELSE 
*. Valence 
            IR_C = 0
          END IF
          IF(IHPV_A.EQ.1) THEN
*. Annihilation of hole, corresponds to excitation
            IR_A = 1
          ELSE IF(IHPV_A.EQ.2) THEN
*. Annihilation of particle, corresponds to de-excitation 
            IR_A =-1
          ELSE 
*. Valence 
            IR_A = 0
          END IF
          IRANK2 = IR_C + IR_A
          IF(IRANK2.EQ.IRANK2_SX.AND.ISYM.EQ.ISYM_SX) THEN
            NSX = NSX + 1
            ISX(1,NSX) = ICOB
            ISX(2,NSX) = IAOB
          END IF
        END DO
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' SX for rank*2 and symmetry ', IRANK2_SX,ISYM_SX
        WRITE(6,*) ' Number of excitations obtained ', NSX
        CALL WRT_SXLIST(ISX,NSX)
      END IF
*
      RETURN
      END 
      SUBROUTINE WRT_SXLIST(ISX,NSX)
*
* Write list of single excitations 
*
*. Jeppe Olsen, Dec. 2004
*
      INCLUDE 'implicit.inc'
      INTEGER ISX(2,NSX)
*
      DO JSX = 1, NSX
       WRITE(6,'(A,I3,A,I3,A)') '(',ISX(1,JSX),',',ISX(2,JSX),')'
      END DO 
*
      RETURN
      END
      SUBROUTINE REFORM_RDM_TO_CUMULANTS(CUMULANTS,ISPOBEX_TP,LSOBEX_TP)
*
* Reform density matrices to cumulants
* 
* On input CUMULANTS is asumed to contain the RDM, on 
* output it will contain the cumulants
*
*. Jeppe Olsen
*
      INCLUDE 'wrkspc.inc'
*
      INCLUDE 'glbbas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cprnt.inc'
*. Type and length of the various spinorbitalexcitationtypes
      INTEGER ISPOBEX_TP(4*NGAS,*), LSOBEX_TP(*)
*
      NTEST = 100
*. Loop over types of spinorbital excitations
      DO IXTP = 1, NSPOBEX_TP
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Type of spin-orbital excitations : '
          CALL WRT_SPOX_TP(ISPOBEX_TP(1,IXTP),1)
        END IF
*. Rank of type (here : just number of creation operators )
        IRANK = IELSUM(ISPOBEX_TP(1,IXTP),2*NGAS)
        WRITE(6,*) ' Rank of operator ', IRANK
*
        IF(IRANK.EQ.1) THEN
*. Reduced density matrices are directly cumulants so no reforming 
        ELSE IF(IRANK.EQ.2) THEN
*. Two-particle cumulant, C(ic1,ic2,ia1,ia2) = D(ic1,ic2,ia1,ia2)
*                 -D(ic1,ia1)*D(ic2,ia2) + D(ic1,ia2)D(ic2,ia1)
*. spinsubtype : aa, ab,bb
           IAOP = IELSUM(ISPOBEX_TP(1,IXTP),NGAS)
           IF(IAOP.EQ.2) THEN
*. AA type
           ELSE IF(IAOP.EQ.1) THEN
*. AB type 
           ELSE IF(IAOP.EQ.0) THEN
*. AB type 
           END IF
        END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' And here comes : The cumulants '
        IPRNCIV_SAVE = IPRNCIV
        IPRNCIV = 1
        CALL ANA_GENCC(CUMULANTS,1)
        IPRNCIV = IPRNCIV_SAVE
      END IF
*
      RETURN
      END
*    |||||
*     '('  
*     \ /
* CLONE:
      SUBROUTINE GEN_IC_ORBOP2(IWAY,NIC_ORBOP,IC_ORBOP,
     &                     INC_SING, INC_DOUB,
     &                     IONLY_EXCOP,I_IGN_OVL,
     &                     IREFSPC,ITREFSPC,IADD_UNI)
*
* Generate single and double 
* orbital excitation types corresponding to internal contraction  
* The orbital excitations working on IREFSPC should contain 
* an component in space ITREFSPC.
*
* Operator-manifold is specified by the arrays
*
*  inc_sing = ( <+2> ,  <0>, <-2> )
*  inc_doub = ( <+4> , <+2>, <0>, <-2>, <-4> )
*
*  the indices can be calculated as
*          idx1 = 2 - rank/2   and   idx2 = 3 - rank/2
*
* where an entry of 1 means inclusion of operators of this rank
* and a zero means to skip this type of operators
*
* If IADD_UNI = 1, the unit operator ( containing zero operators)
* is added at the end
*
* Jeppe Olsen, August 2002
*
*
* IWAY = 1 : Number of orbital excitations for internal contraction
* IWAY = 2 : Generate also the actual orbital excitations 
*
* IONLY_EXCOP = 1 => only excitation operators ( no annihilation in particle 
*                    space, no creation in inactive space )
*
* I_IGN_OVL = 1   => we ignore the overlap criterion and include operators
*                    that in first order vanish, but which in higher order
*                    may contribute
*
*. Rank is defined as # crea of particles + # anni of holes 
*                    -# crea of holes     - # anni of particles

      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input array
      INTEGER INC_SING(3), INC_DOUB(5)
*. Local scratch
      INTEGER ITREFOCC(MXPNGAS,2)
*. Output ( if IWAY .ne. 1 ) 
      INTEGER IC_ORBOP(2*NGAS,*)
*. Local scratch
      INTEGER IOP(2*MXPNGAS)
*
      NTEST =  100 
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' IREFSPC, ITREFSPC = ', IREFSPC, ITREFSPC 
        WRITE(6,'(X,A,3I2)') ' INC_SING = ', INC_SING(1:3) 
        WRITE(6,'(X,A,5I2)') ' INC_DOUB = ', INC_DOUB(1:5) 
      END IF
      NIC_ORBOP =  0
      IF (NTEST.GE.100) WRITE(6,*) ' output for singles:'
*. Single excitations a+i a j
      DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          IZERO = 0
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
          IF (I_IGN_OVL.NE.1) THEN
      CALL ORBOP_ACCOCC(IOP,IGSOCCX(1,1,IREFSPC),ITREFOCC,NGAS,MXPNGAS)
      CALL OVLAP_ACC_MINMAX(ITREFOCC,IGSOCCX(1,1,ITREFSPC),NGAS,MXPNGAS,
     &         IOVERLAP)
      IF(NTEST.GE.100) WRITE(6,*) ' IOVERLAP from OVLAP..',IOVERLAP
      IF(IOVERLAP.EQ.0) IOKAY = 0
           ELSE
             IOKAY = 1
           END IF
C     ORBOP_ACCOCC(IORBOP,IACC_IN,IACC_OUT,NGAS,MXPNGAS)
C     OVLAP_ACC_MINMAX(IACC1,IACC2,NGAS,MXPNGAS,IOVERLAP)
*. is there any operators in spaces that are frozen or deleted in ITREFSPC
C     CHECK_EXC_FR_OR_DE(IOP,IOCC,NGAS,IOKAY)
          CALL CHECK_EXC_FR(IOP,IGSOCCX(1,1,ITREFSPC),NGAS,IOKAY2)
          IF(NTEST.GE.100) WRITE(6,*) ' IOKAY2 from CHECK ... ', IOKAY2
          IF(IOKAY2.EQ.0) IOKAY = 0
          IF(IOKAY.EQ.1) THEN
            CALL IRANK_ORBOP(IOP,NEX,NDEEX)
            IOKAY2 = 1
            IF(IONLY_EXCOP.EQ.1.AND.NDEEX.NE.0) IOKAY2 = 0
            IRANK = NEX - NDEEX
            IF(NTEST.GE.100) WRITE(6,*) ' IRANK = ', IRANK
            IF(INC_SING(2-IRANK/2).NE.0
c test
c            IF(INC_SING(2-IRANK).NE.0
     &      .AND.IOKAY2.EQ.1)THEN
              NIC_ORBOP  = NIC_ORBOP + 1
              IF(NTEST.GE.100) WRITE(6,*) ' Operator included '
              IF(IWAY.NE.1) 
     &        CALL ICOPVE(IOP,IC_ORBOP(1,NIC_ORBOP),2*NGAS)
            END IF
          END IF
        END DO
      END DO
*. Double excitations a+i a+j a k a l
      IF (NTEST.GE.100) WRITE(6,*) ' output for doubles:'
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS
          DO KGAS = 1, NGAS
            DO LGAS = 1, KGAS
              CALL ISETVC(IOP,IZERO,2*NGAS)
              IOP(IGAS) = 1
              IOP(JGAS) = IOP(JGAS) + 1
              IOP(NGAS+KGAS) = 1
              IOP(NGAS+LGAS) = IOP(NGAS+LGAS) + 1
              IF(NTEST.GE.200) THEN
                WRITE(6,*) ' Next Orbital excitation '
                CALL IWRTMA(IOP,NGAS,2,NGAS,2)
              END IF
              CALL COMPARE_OPDIM_ORBDIM(IOP,IOKAY)
              IF(NTEST.GE.200) WRITE(6,*) ' IOKAY from COMPARE..', IOKAY
*. Is the action of this operator on IREFSPC included in ITREFSPC
              IF (I_IGN_OVL.NE.1) THEN
      CALL ORBOP_ACCOCC(IOP,IGSOCCX(1,1,IREFSPC),ITREFOCC,NGAS,MXPNGAS)
      CALL OVLAP_ACC_MINMAX(ITREFOCC,IGSOCCX(1,1,ITREFSPC),NGAS,
     &         MXPNGAS,IOVERLAP)
      IF(NTEST.GE.200) WRITE(6,*) ' IOVERLAP from OVLAP..',IOVERLAP
      IF(IOVERLAP.EQ.0) IOKAY = 0
              ELSE
                IOKAY = 1
              END IF
              CALL CHECK_EXC_FR(IOP,IGSOCCX(1,1,ITREFSPC),NGAS,IOKAY2)
              IF(NTEST.GE.200)
     &             WRITE(6,*) ' IOKAY2 from CHECK ... ', IOKAY2
              IF(IOKAY2.EQ.0) IOKAY = 0
              IF(IOKAY.EQ.1) THEN
                CALL IRANK_ORBOP(IOP,NEX,NDEEX)
                IOKAY2 = 1
                IF(IONLY_EXCOP.EQ.1.AND.NDEEX.NE.0) IOKAY2 = 0
                IRANK = NEX - NDEEX
                IF(NTEST.GE.100) WRITE(6,*) ' IRANK = ', IRANK
                IF(INC_DOUB(3-IRANK/2).NE.0 .AND.
c test
c                IF(INC_DOUB(3-IRANK).NE.0 .AND.
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
* END OF CLONE
      SUBROUTINE PROJ_VEC_TO_ICSPC(LUREF,LUIN,LUOUT,VEC1_CI,VEC2_CI,
     &           VEC1_IC,VEC2_IC,VEC3_IC,RMAT_IC,
     &           IREFSPC,ITREFSPC,NSPA,N_IC_OP,N_NONSING,S_IC,
     &           X_IC_NONSING,LUSCR)
*
* A vector is given in uncontracted basis (Determinant basis)
* on LUIN. Project this vector to the space given by the 
* internal contracted operators O_i |ref> where |ref> is 
* the vector on LUREF
*
* Jeppe Olsen, May 2005 for settling whether the IC triples 
* correction is the exact second order MP triples correction
*
* The projected vector is 
*
* sum_ij O_i|ref> S_{ij}^-1 <ref|O+j|LUIN>
*
* So the procedure is 
* 1 : Calculate  <ref|O+j|LUIN> as density 
* 2 : Invert S and multiply on <ref|O+j|LUIN>
* 3 : Expand resulting vector in SD space
* 4 : And compare
*
      INCLUDE 'wrkspc.inc'
      REAL*8 INPRDD
      INCLUDE 'cands.inc'
* ========
*.  Input 
* ========
*. Metric in IC basis - unitoperator excluded  IS DESTROYED IN THIS ROUTINE !!!
      DIMENSION S_IC((NSPA-1)**2)
*.Transformation basis IC=> Non-sing basis  (minus unit operator)
      DIMENSION X_IC_NONSING(NSPA-1,N_NONSING)
* =========
*. Scratch 
* =========
*. Scratch for CI
      DIMENSION VEC1_CI(*), VEC2_CI(*)
*. For holding IC vectors
       DIMENSION VEC1_IC(N_IC_OP), VEC2_IC(N_IC_OP)
*. and an matrix in IC basis
       DIMENSION RMAT_IC(N_IC_OP,N_IC_OP)
*
      NTEST = 00
*
*     <REF!T+(I)P H  !0>  = <LUIN!T(I)!LUREF>
*
      IF(NTEST.GE.10) THEN
      WRITE(6,*) ' PROJ ..., LUIN, LUOUT, LUSCR = ', LUIN,LUOUT,LUSCR
      WRITE(6,*) ' PROJ ... N_IC_OP, NSPA, N_NONSING = ',
     &                       N_IC_OP, NSPA, N_NONSING 
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input vector in SD basis '
        CALL WRTVCD(VEC1_CI,LUIN,1,-1)
      END IF
*. Both sides are in the form of the ITREFSPC so :
      ICSPC = ITREFSPC
      ISSPC = ITREFSPC
      ZERO = 0.0D0
      CALL SETVEC(VEC1_IC,ZERO,N_IC_OP)
      CALL SIGDEN_CC(VEC1_CI,VEC2_CI,LUREF,LUIN,VEC1_IC,2)
      CALL REF_CCV_CAAB_SP(VEC1_IC,VEC2_IC,VEC3_IC,1)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transition density <ref|O+j|LUIN> in IC basis '
        CALL WRTMAT(VEC2_IC,1,NSPA,1,NSPA)
      END IF
*. and transform to nonsingular basis
      CALL MATVCC(X_IC_NONSING,VEC2_IC,VEC1_IC,NSPA-1,N_NONSING,1)
*. Transform the metric to the nonsingular space 
C     TRNMAD(A,X,SCR,NDIMI,NDIMO)
      CALL TRNMAD(S_IC,X_IC_NONSING,RMAT_IC,NSPA-1,N_NONSING)
* Obtain inverse metric  in S_IC
      CALL  INVMAT(S_IC,RMAT_IC,N_NONSING,N_NONSING,ISING)
*. Multiply  <ref|O+j|LUIN> with inverse metric
      CALL MATVCC(S_IC,VEC1_IC,VEC2_IC,N_NONSING,N_NONSING,0)
*. Transform back to SPA basis 
      CALL MATVCC(X_IC_NONSING,VEC2_IC,VEC1_IC,NSPA-1,N_NONSING,0)
*. We have left out the coefficient corresponding to the 
*. zero-order state. Set this to zero 
      VEC1_IC(NSPA) = 0.0D0
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Projected vector in IC basis '
        CALL WRTMAT(VEC1_IC,1,NSPA,1,NSPA)
      END IF
*. We now have projected vector in IC basis, expand in SD
*. basis to allow comparison
C     REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY) 
      CALL REF_CCV_CAAB_SP(VEC2_IC,VEC1_IC,VEC3_IC,2)
      CALL SIGDEN_CC(VEC1_CI,VEC2_CI,LUREF,LUOUT,VEC2_IC,1)
*. Obtain difference between the two vectors on LUSCR
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      FAC1 = 1.0D0
      FAC2 = -1.0D0
      CALL VECSMD(VEC1_CI,VEC2_CI,FAC1,FAC2,LUIN,LUOUT,LUSCR,1,-1)
*. Norm of LUIN and of LUIN-LUOUT
      XNORM_IN = INPRDD(VEC1_CI,VEC2_CI,LUIN,LUIN,1,-1)
      XNORM_OUT = INPRDD(VEC1_CI,VEC2_CI,LUOUT,LUOUT,1,-1)
      XNORM_DIFF = INPRDD(VEC1_CI,VEC2_CI,LUSCR,LUSCR,1,-1)
*. And compare individual elements
      WRITE(6,*) ' Comparison of LUIN and LUOUT '
      CALL CMP2VCD(VEC1_CI,VEC2_CI,LUIN,LUOUT,0.0D0,1,-1)
*
      WRITE(6,*) ' Comparing vector and vector projected to IC space '
      WRITE(6,*) ' Squared norm of input vector = ', XNORM_IN
      WRITE(6,*) ' Squared norm of output vector = ', XNORM_OUT
      WRITE(6,*) ' Squared norm of difference    = ', XNORM_DIFF
*
      RETURN
      END
      SUBROUTINE TRNMAD(A,X,SCR,NDIMI,NDIMO)
*
* Obtain X(T) A X and store it in A
* Allows different dimensions in input and output matrices 
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION A(*), X(NDIMI,NDIMO)
*. Scratch
      DIMENSION SCR(NDIMI*NDIMO)
      NTEST = 000
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Info from TRNMAD '
        WRITE(6,*) '   NDIMI, NDIMO = ', NDIMI,NDIMO
        WRITE(6,*) ' Input X matrix '
        CALL WRTMAT(X,NDIMI,NDIMO,NDIMI,NDIMO)
        WRITE(6,*) ' Input A matrix '
        CALL WRTMAT(A,NDIMI,NDIMI,NDIMI,NDIMI)
       END IF
     
*
*. 1 : X(T) A in SCR 
      ZERO = 0.D0
      CALL SETVEC(SCR,ZERO,NDIMI*NDIMO)
      CALL MATML7(SCR,X,A,NDIMO,NDIMI,NDIMI,NDIMO,NDIMI,NDIMI,
     &              0.0D0,1.0D0,1)
*. X(T) A X in A
      CALL MATML7(A,SCR,X,NDIMO,NDIMO,NDIMO,NDIMI,NDIMI,NDIMO,
     &            0.0D0,1.0D0,0)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformed matrix : '
        CALL WRTMAT(A,NDIMO,NDIMO,NDIMO,NDIMO)
      END IF
*
      RETURN
      END
      SUBROUTINE EXPND_T_TO_NOSYM(XIN,XOUT,ICAAB)
*
* A matrix XIN is given in symmetry packed form XIN(CA,CB,AA,AB)
* Expand to form without symmetry
*
* Jeppe Olsen
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
* Specific input
      INTEGER ICAAB(NGAS,4)
      DIMENSION XIN(*)
*. Output
      DIMENSION XOUT(*)
*
      IDUM = -1
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'EXPNOS')
*. 
      NOP_CA = IELSUM(ICAAB(1,1),NGAS)
      NOP_CB = IELSUM(ICAAB(1,2),NGAS)
      NOP_AA = IELSUM(ICAAB(1,3),NGAS)
      NOP_AB = IELSUM(ICAAB(1,4),NGAS)
*
      NOP_MX = MAX(NOP_CA,NOP_CB,NOP_AA,NOP_AB)
    
*. Set up arrays for indexing ICA, ICB, IAA, IAB without symmetry
      CALL MEMMAN(KLZ_CA,NOP_CA*NTOOB,'ADDL  ',2,'Z_CA  ')
      CALL MEMMAN(KLZ_CB,NOP_CB*NTOOB,'ADDL  ',2,'Z_CB  ')
      CALL MEMMAN(KLZ_AA,NOP_AA*NTOOB,'ADDL  ',2,'Z_AA  ')
      CALL MEMMAN(KLZ_AB,NOP_AB*NTOOB,'ADDL  ',2,'Z_AB  ')
      LSCR = 2*NTOOB + (NOP_MX+1)*(NTOOB+1)
      CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'ZLSCR')
C          WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,ISCR,NTEST)
      CALL WEIGHT_SPGP(WORK(KLZ_CA),NGAS,ICAAB(1,1),NOBPT,WORK(KLSCR),0)
      CALL WEIGHT_SPGP(WORK(KLZ_CB),NGAS,ICAAB(1,2),NOBPT,WORK(KLSCR),0)
      CALL WEIGHT_SPGP(WORK(KLZ_AA),NGAS,ICAAB(1,3),NOBPT,WORK(KLSCR),0)
      CALL WEIGHT_SPGP(WORK(KLZ_AB),NGAS,ICAAB(1,4),NOBPT,WORK(KLSCR),0)
*. Total number of strings per ICAAB ( is also given in last elements of Z's)
      NST_CA = NST_FOR_OCC(ICAAB(1,1),NOBPT,NGAS)
      NST_CB = NST_FOR_OCC(ICAAB(1,2),NOBPT,NGAS)
      NST_AA = NST_FOR_OCC(ICAAB(1,3),NOBPT,NGAS)
      NST_AB = NST_FOR_OCC(ICAAB(1,4),NOBPT,NGAS)
*. In the general form, a string is XOUT(ICA,ICB,IAA,IAB) will be adressed 
*. as a standard fortran array
*. We are now ready to do the reordering
      ZERO = 0.0D0
      NELMNT = NST_CA*NST_CB*NST_AA*NST_AB 
      CALL SETVEC(XOUT,ZERO,NELMNT)
*. Four scratch blocks for holding blocks of 

      CALL  EXPND_T_TO_NOSYMS(XIN,XOUT,ICAAB,ISM,
     &      WORK(KLZ_CA),WORK(KLZ_CB),WORK(KLZ_AA),WORK(KLZ_AB),
     &      IOCC_CA, IOCC_CB, IOCC_AA, IOCC_AB,NORB,MSCOMB_CC)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'EXPNOS')
      RETURN
      END
      SUBROUTINE EXPND_T_TO_NOSYMS(XIN,XOUT,ICAAB,ISM,
     &      IZ_CA,IZ_CB,IZ_AA,IZ_AB,
     &      IOCC_CA, IOCC_CB, IOCC_AA, IOCC_AB,NORB,MSCOMB_CC)
*
*. An array T(ICA,ICB,IAA,IAB) is given in symmetry-ordered form.
*. Unpack to form without symmetry
*. 
*. Jeppe Olsen, April 2005
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER ICAAB(NGAS,4)
      DIMENSION XIN(*)
      INTEGER IZ_CA(*),IZ_CB(*),IZ_AA(*),IZ_AB(*)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*. Output
      DIMENSION XOUT(*)
*. Total number of strings for the various groups
      NST_CA_TOT = NST_FOR_OCC(ICAAB(1,1),NOBPT,NGAS)
      NST_CB_TOT = NST_FOR_OCC(ICAAB(1,2),NOBPT,NGAS)
      NST_AA_TOT = NST_FOR_OCC(ICAAB(1,3),NOBPT,NGAS)
      NST_AB_TOT = NST_FOR_OCC(ICAAB(1,4),NOBPT,NGAS)

*
*. Transform from occupations to groups
      CALL OCC_TO_GRP(ICAAB(1,1),IGRP_CA,1)
      CALL OCC_TO_GRP(ICAAB(1,2),IGRP_CB,1)
      CALL OCC_TO_GRP(ICAAB(1,3),IGRP_AA,1)
      CALL OCC_TO_GRP(ICAAB(1,4),IGRP_AB,1)
*
      NEL_CA = IELSUM(ICAAB(1,1),NGAS)
      NEL_CB = IELSUM(ICAAB(1,2),NGAS)
      NEL_AA = IELSUM(ICAAB(1,3),NGAS)
      NEL_AB = IELSUM(ICAAB(1,4),NGAS)
*. It is assumed that no reduction due to spin symmetri is used.
      DO ISM_C = 1, NSMST
       ISM_A = MULTD2H(ISM,ISM_C) 
       DO ISM_CA = 1, NSMST
        ISM_CB = MULTD2H(ISM_C,ISM_CA)
        DO ISM_AA = 1, NSMST
         ISM_AB =  MULTD2H(ISM_A,ISM_AA)
         ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
         ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
*. obtain strings
         CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,NSTR_CA,
     &        IOCC_CA, NORB,0,IDUM,IDUM)
         CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,NSTR_CB,
     &        IOCC_CB, NORB,0,IDUM,IDUM)
         CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,NSTR_AA,
     &        IOCC_AA, NORB,0,IDUM,IDUM)
         CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,NSTR_AB,
     &        IOCC_AB, NORB,0,IDUM,IDUM)
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
         DO I_AB = 1, NSTR_AB
*. Number in nonsymmetric form
C  ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
          I_AB_EXP = ISTRNM(IOCC_AB(1+(I_AB-1)*NEL_AB),NORB,IZ_AB,
     &               IDUM,0)
          DO I_AA = 1, NSTR_AA
           I_AA_EXP = ISTRNM(IOCC_AA(1+(I_AA-1)*NEL_AA),NORB,IZ_AA,
     &                IDUM,0)
           DO I_CB = 1, NSTR_CB
            I_AB_EXP = ISTRNM(IOCC_CB(1+(I_CB-1)*NEL_CB),NORB,IZ_CB,
     &                 IDUM,0)
            DO I_CA = 1, NSTR_CA
             I_CA_EXP = ISTRNM(IOCC_CA(1+(I_CA-1)*NEL_CA),NORB,IZ_CA,
     &                  IDUM,0)
             IT = IT + 1
             IT_EXP = (IAB_EXP-1)*NST_CA_TOT*NST_CB_TOT*NST_AA_TOT 
     &              + (IAB_EXP-1)*NST_CA_TOT*NST_CB_TOT
     &              + (ICB_EXP-1)*NST_CA_TOT
     &              + ICA_EXP
             XOUT(IT_EXP) = XIN(IT)
            END DO
*           ^ End of loop over alpha creation strings
           END DO
*          ^ End of loop over beta creation strings
          END DO
*         ^ End of loop over alpha annihilation 
         END DO 
*        ^ End of loop over beta annihilation 
  777   CONTINUE
        END DO
       END DO
      END DO
*      ^ End of loop over symmetry blocks
      RETURN
      END
      SUBROUTINE LUCIA_ICCC(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                      EFINAL,CONVER,VNFINAL)
*
* Master routine for Internal Contraction multireference coupled cluster theory
*
* LUCIA_IC is assumed to have been called to do the 
* prepatory work for working with internal contraction
*
* It is assumed that spin-adaptation is used ( no flag anymore..)
*
* It is standard that the unitoperator is included in 
* the operator manifold, but in CC ( and PT)  theory this should be 
* excluded. This is easily done as the unitoperator is the 
* last operator in CA order.
*
* Jeppe Olsen, August 2005
*
C     INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
      LOGICAL CONVER,CONVERL
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cei.inc'
*. Transfer common block for communicating with H_EFF * vector routines
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
*. Transfer block for communicating zero order energy to 
*. routine for performing H0-E0 * vector
      INCLUDE 'cshift.inc'
*
      CHARACTER*6 ICTYP
      EXTERNAL H0ME0TV_EXT_IC
*. Number of commutators used in approach
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICCC  ')
      NTEST = 5
*
*. a bit of dirty work before print:
*  I will add unitoperator to the spin-orbital excitations-
* evrything prepared, I just need to increase number of types
* Is already done in old non-IE-route, add for IE-route
      IF(I_DO_EI.EQ.1) THEN
        NSPOBEX_TP = NSPOBEX_TP + 1
      END IF
      WRITE(6,*)
      WRITE(6,*) ' ===================='
      WRITE(6,*) ' ICCC section entered '
      WRITE(6,*) ' ===================='
      WRITE(6,*)
*
*. Form of ICPT calculation 
*
      WRITE(6,'(A,A)') ' Type of ICCC calculation : ', ICTYP
      WRITE(6,*) ' Energy of reference state ', EREF
      WRITE(6,*) ' Reference space ', IREFSPC
      WRITE(6,*) ' Extended space (ITREFSPC) ', ITREFSPC
      WRITE(6,*) ' Number of commutators employed : '
      WRITE(6,*) '    In energy evaluation     ', NCOMMU_E
      WRITE(6,*) '    In approximate Jacobian  ', NCOMMU_J
      WRITE(6,*) '    In vector function       ', NCOMMU_V
*
      IF(I_FIX_INTERNAL.EQ.0) THEN
        WRITE(6,*) ' Internal (reference) wave-function reoptimized'
      ELSE
        WRITE(6,*) ' Internal (reference) wave-function frozen'
      END IF
      IF(I_INT_HAM.EQ.1) THEN
        WRITE(6,*) ' One-body H0 used for internal zero-order states'
      ELSE
        WRITE(6,*) ' One-body H used for internal zero-order states'
      END IF
*
*. Approximate highest commutator 
      N_APPROX_HCOM = I_APPROX_HCOM_E + I_APPROX_HCOM_V 
     &              + I_APPROX_HCOM_J
      IF(N_APPROX_HCOM.NE.0) THEN
        WRITE(6,*) ' Highest commutator approximated for '
        IF(I_APPROX_HCOM_E.EQ.1) WRITE(6,*) '    energy-function'
        IF(I_APPROX_HCOM_V.EQ.1) WRITE(6,*) '    vector-function'
        IF(I_APPROX_HCOM_J.EQ.1) WRITE(6,*) '    approximate Jacobian'
      END IF
      IF(I_DO_EI.EQ.1) THEN
       WRITE(6,*) ' EI approach in use'
      ELSE
       WRITE(6,*) ' Partial spin-adaptation in use'
      END IF
*
      WRITE(6,*) ' LUCIA_ICCC: IREFSPC, ITREFSPC =', IREFSPC, ITREFSPC
      WRITE(6,*) ' Number of spinorbitalexctypes (inc. unit)'
     &           , NSPOBEX_TP
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' The list of spinorbitalexcitations'
        CALL WRT_SPOX_TP_JEPPE(WORK(KLSOBEX),NSPOBEX_TP)
      END IF
*. Number of parameters with and without spinadaptation
      IF(I_DO_EI.EQ.0) THEN
        CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
      ELSE
*. zero-particle operator is included in N_ZERO_EI
        NSPA = N_ZERO_EI 
*. Note: NCAAB and N_CC_AMP below now both includes unitop
        NCAAB = NDIM_EI
        N_CC_AMP = NCAAB 
      END IF
      IF(NTEST.GE.10) THEN
        IF(I_DO_EI.EQ.0) THEN
          WRITE(6,*) ' Number of spin-adapted operators ', NSPA
        ELSE
          WRITE(6,*) ' Number of orthonormal zero-order states',
     &                 N_ZERO_EI
        END IF
        WRITE(6,*) ' Number of CAAB operators         ', NCAAB
        WRITE(6,*) ' Number of CC amplitudes          ', N_CC_AMP
*
        WRITE(6,*) ' Threshold for nonsingular metric eigenvalues =',
     &  THRES_SINGU
      END IF
*. Number of spin adapted operators without the unitoperator
      NSPAM1 = NSPA - 1
      N_REF = XISPSM(IREFSM,IREFSPC)
*. Size of subspace Jacobian
      MXVEC_SBSPJA = 15
      IF(I_DO_SBSPJA.EQ.1) THEN
        WRITE(6,*) 
     &  ' Subspace Jacobian will be constructed. Max. dim of subspace ',
     &  MXVEC_SBSPJA
      END IF
*
* ==============================================
* 1 : Set up zero-order Hamiltonian in WORK(KFIFA)
* ==============================================
*
*. It is assumed that one-body density over reference resides 
*  in WORK(KRHO1)
*
      CALL COPVEC(WORK(KINT1O),WORK(KFIFA),NINT1)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The original one-body hamiltonian '
        CALL APRBLM2(WORK(KINT1O),NTOOBS,NTOOBS,NSMOB,1)
      END IF
*. Calculate zero-order Hamiltonian : use either actual density
*. or Hartree-Fock densi
      I_ACT_OR_HF = 1
      IF(I_ACT_OR_HF.EQ.1) THEN
        WRITE(6,*) ' Zero-order Hamiltonian with actual density '
        CALL FIFAM(WORK(KFIFA))
      ELSE
        WRITE(6,*) ' Zero-order Hamiltonian with zero-order density '
*. IPHGAS1 should be used to divide into H,P,V, but IPHGAS is used, so swap
        CALL ISWPVE(IPHGAS(1),IPHGAS1(1),NGAS)
*
        CALL COPVEC(WORK(KINT1O),WORK(KFIFA),NINT1)
        CALL FI(WORK(KFIFA),ECC,1)
        WRITE(6,*) ' FI before zeroing : '
        CALL APRBLM2(WORK(KFIFA),NTOOBS,NTOOBS,NSMOB,1)
*. And clean up
        CALL ISWPVE(IPHGAS,IPHGAS1,NGAS)
*. zero offdiagonal elements 
        IF(I_DO_EI.EQ.0) THEN 
          CALL ZERO_OFFDIAG_BLM(WORK(KFIFA),NSMOB,NTOOBS,1)
        END IF
      END IF
*
      IF(NTEST.GE.00) THEN
        WRITE(6,*) ' One-body zero-order Hamiltonian '
        CALL APRBLM2(WORK(KFIFA),NTOOBS,NTOOBS,NSMOB,1)
      END IF
*. Scratch space for CI - has already been allocated in EI approach
      IF(I_DO_EI.EQ.0) THEN
        CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
        KVEC1P = KVEC1
        KVEC2P = KVEC2
      END IF
*
* =====================================================================
* Obtain metric matrix and nonsingular set of operators in WORK(KLXMAT)
* =====================================================================
*. Some additional scratch, dominated by two complete matrices !!
*. And a few working vectors 
      CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'VCC1  ')
      CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'VCC2  ')
      CALL MEMMAN(KLVCC3,N_CC_AMP,'ADDL  ',2,'VCC3  ')
      CALL MEMMAN(KLVCC4,N_CC_AMP,'ADDL  ',2,'VCC4  ')
      CALL MEMMAN(KLRHS ,N_CC_AMP,'ADDL  ',2,'RHS   ')
      CALL MEMMAN(KLC1  ,N_CC_AMP,'ADDL  ',2,'C1    ')
      CALL MEMMAN(KLC1O ,N_CC_AMP,'ADDL  ',2,'C1    ')
      CALL MEMMAN(KLC_REF,N_REF   ,'ADDL  ',2,'C_REF  ')
      CALL MEMMAN(KLI_REF,N_REF   ,'ADDL  ',1,'I_REF  ')
      IF(I_DO_SBSPJA.EQ.1) THEN
        LSBSPJA = 5*MXVEC_SBSPJA**2 + 2*MXVEC_SBSPJA
        CALL MEMMAN(KLSBSPJA,LSBSPJA,'ADDL  ',2,'SBSPJA')
      ELSE
        LSBSPJA = 0
        KLSBSPJA = 1
      END IF
*. Identify the unit  operator i.e. the operator with 
*. zero creation and annihilation operators
      IDOPROJ = 1
      IF(IDOPROJ.EQ.1) THEN
        CALL GET_SPOBTP_FOR_EXC_LEVEL(0,WORK(KLCOBEX_TP),NSPOBEX_TP,
     &       NUNIOP,IUNITP,WORK(KLSOX_TO_OX))
*. And the position of the unitoperator in the list of SPOBEX operators
        WRITE(6,*) ' NUNIOP, IUNITP = ', NUNIOP,IUNITP
        IF(NUNIOP.EQ.0) THEN
          WRITE(6,*) ' Unitoperator not found in exc space '
          WRITE(6,*) ' I will proceed without projection '
          IDOPROJ = 0
        ELSE
          IUNIOP = IFRMR(WORK(KLIBSOBEX),1,IUNITP)
          IF(NTEST.GE.100) WRITE(6,*) ' IUNIOP = ', IUNIOP
        END IF
      END IF
*
* We will iterate over optimization of internal and external 
* parts of the CC wavefunction, allowed number of iteration
*. 
*. Flag for iterative calculation
      IF(I_DO_EI.EQ.1) THEN
        I_IT_OR_DIR_IN_EXT   = 1
      ELSE
        I_IT_OR_DIR_IN_EXT   = 1
      END IF
*. Will we allow relaxation of coefficients defining reference
*. state
      I_RELAX_INT = 1
*. Will direct or iterative methods be used for relaxing
*. reference coefficients
      I_IT_OR_DIR_IN_RELAX = 1
*. Space for external correlation vector
      CALL MEMMAN(KLTEXT,N_CC_AMP,'ADDL  ',2,'T_EXT ')
*
*. Initial  T_EXT : zero or readin
*
      IF(IRESTRT_IC.EQ.0) THEN
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLTEXT),ZERO,NCAAB)
*. Store inital guess on unit 54 in CAAB form 
        CALL VEC_TO_DISC(WORK(KLTEXT),NCAAB,1,-1,LUSC54)
      ELSE 
        WRITE(6,*) ' T_ext restarted from  LU54'
        CALL VEC_FROM_DISC(WORK(KLTEXT),NCAAB,1,-1,LUSC54)
        WRITE(6,*) 'T_EXT read in '
      END IF
*. Allocate vectors for CI behind the curtain 
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
*
      IF(IRESTRT_IC.EQ.1) THEN
*. Copy old CI coefficients for reference space to LUC
        CALL REWINO(LUC)
        CALL COPVCD(LUSC54,LUC,WORK(KVEC1),0,-1)
        WRITE(6,*) ' Internal coefs copied from LUSC54'
      END IF

*
      MAXITG = MAXITM
      CONVER =.FALSE.
      CONVERL =.FALSE.
*. Convergence threshold for norm of vectorfunction
      VTHRES = 1.0D-11
      DO IT_IE = 1, MAXITG
        IDUM = 0
*
* ===============================================
*. Optimize T for current internal coefficients 
* ===============================================
*
C?      WRITE(6,*)  ' ITREFSPC before call to ICCC ', ITREFSPC

        IF(IT_IE.EQ.1) THEN
          INI_IT = 1
        ELSE 
          INI_IT = 0
        END IF
        IF(IT_IE.EQ.MAXITG) THEN
          IFIN_IT = 1
        ELSE 
          IFIN_IT = 0
        END IF
*. use DIIS/CROP to accelerate
        IDIIS = 2
*. Use approach  where internal and external parts are 
*. optimized simultaneously.
        ISIMULT = 1
*
*. In the calculation of the MRCC vector function 3 spaces 
*. will be used
* 1 : IREFSPC : Space of !0>
* 2 : IT2REFSPC : Space of T!0>
* 3 : ITREFSPC : Largest space needed in the calculation of e(-T) H e(T)
*. In the following it will be assumed that IT2REFSPC is the space BEFORE
*. ITREFSPC
        IT2REFSPC = ITREFSPC 
        IT2REFSPC = ITREFSPC - 1
C?      WRITE(6,*) ' After Mod: ITREFSPC, IT2REFSPC=',
C?   &                          ITREFSPC, IT2REFSPC
C?          WRITE(6,*) ' Space for T !0> : ', IT2REFSPC
*. Readin C_REF
        CALL REWINO(LUC)
        CALL FRMDSCN(WORK(KLC_REF),-1,-1,LUC)
*
        
        I_DO_COMP = 0
        IF(I_DO_COMP.EQ.1) THEN
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
            WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
          WRITE(6,*) ' Note: Complete matrix flag activated'
        END IF
*
        I_REDO_INT = 1
*
        I_CAAB_OR_ORT = 2
        IF(I_CAAB_OR_ORT.EQ.1) THEN
          CALL ICCC_OPT_SIMULT(
     &          IREFSPC,ITREFSPC,IT2REFSPC,I_SPIN_ADAPT,
     &          NROOT,WORK(KLTEXT),C_0,INI_IT,IFIN_IT,
     &          WORK(KVEC1),WORK(KVEC2),IDIIS,
     &          WORK(KLC_REF),N_REF,I_DO_COMP,CONVERL,VTHRES,
     &          I_REDO_INT,EFINAL,VNFINAL,CONVER,
     &          WORK(KLSBSPJA),MXVEC_SBSPJA,I_FIX_INTERNAL)
        ELSE
          CALL ICCC_OPT_SIMULT_ONB(
     &          IREFSPC,ITREFSPC,IT2REFSPC,I_SPIN_ADAPT,
     &          NROOT,WORK(KLTEXT),C_0,INI_IT,IFIN_IT,
     &          WORK(KVEC1),WORK(KVEC2),IDIIS,
     &          WORK(KLC_REF),N_REF,I_DO_COMP,CONVERL,VTHRES,
     &          I_REDO_INT,EFINAL,VNFINAL,CONVER,
     &          WORK(KLSBSPJA),MXVEC_SBSPJA,I_FIX_INTERNAL)
        END IF
*. transfer new C_REF to file LUC
        CALL ISTVC2(WORK(KLI_REF),0,1,N_REF)
        CALL REWINO(LUC)
        CALL WRSVCD(LUC,-1,WORK(KVEC1),WORK(KLI_REF),
     &          WORK(KLC_REF),N_REF,N_REF,LUDIA,1)
*. Save current T_ext in CAAB form and CI coefs on LUSC54
        CALL VEC_TO_DISC(WORK(KLTEXT),NCAAB,1,-1,LUSC54)
        CALL WRSVCD(LUSC54,-1,WORK(KVEC1),WORK(KLI_REF),
     &          WORK(KLC_REF),N_REF,N_REF,LUDIA,1)
        REWIND(LUSC54)
        IF(CONVER) GOTO 1001
*
        IF(ISIMULT.EQ.0.AND.I_RELAX_INT.EQ.1) THEN
* ============================================================
*. Relax coefficients of internal/reference/zero-order state
* ============================================================
*
*. Three vectors are actually allocated and kept in ICCC_COMPLETE..
*. so these could and should be reused 
           CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'ICCREL')
*
           IF(I_IT_OR_DIR_IN_RELAX.EQ.2) THEN
*
*. Construct complete matrices and diagonalize
*
*. Space for H and S in zero-order space
             N_REF = XISPSM(IREFSM,IREFSPC)
             CALL MEMMAN(KLH_REF,N_REF**2,'ADDL  ',2,'H_REF  ')
             CALL MEMMAN(KLS_REF,N_REF**2,'ADDL  ',2,'S_REF  ')
*
C     ICCC_RELAX_REFCOEFS_COM(T_EXT,H_REF,N_REF,
C    &           NCOMMU,VEC1,VEC2,IREFSPC,ITREFSPC,
C    &           ECORE,C_REF_OUT,IREFROOT)
             CALL ICCC_RELAX_REFCOEFS_COM(WORK(KLTEXT),
     &            WORK(KLH_REF),N_REF,NCOMMU_E,WORK(KVEC1),
     &            WORK(KVEC2),
     &            IREFSPC,ITREFSPC,ECORE,WORK(KLC_REF),NROOT)
*. transfer new reference vector to DISC
             CALL ISTVC2(WORK(KLI_REF),0,1,N_REF)
C  WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
             CALL REWINO(LUC)
             CALL WRSVCD(LUC,-1,WORK(KVEC1),WORK(KLI_REF),
     &            WORK(KLC_REF),N_REF,N_REF,LUDIA,1)
           ELSE 
             WRITE(6,*) ' Iterative ICCC not working yet '
           END IF
*.         ^ End of switch direct/iterative methods for reference
*.         relaxation
           CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'ICCREL')
        END IF
*.      ^ End of reference coefs should be relaxed
      END DO
*.    ^ End of loop over Internal/external correlation iterations
 1001 CONTINUE
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICCC  ')
      RETURN
      END 
      SUBROUTINE MRCC_VECFNC(CCVECFNC,T,NCOMMU,I_APPROX_HCOM,
     &           IREFSPC,ITREFSPC,IT2REFSPC,CCVECFNCI)
*
* Obtain external and internal parts of the MRCC vector function 
*
* External part : 
* ================
*
* <0!tau^{\dagger} exp(-T) H exp(T) !0>. 
*
*. Internal part 
* ================
*
* <J! exp(-T) H exp(T) !0>
*
* Input and output vectors  are in CAAB basis.
*. The commutator  exp(-T) H exp(T) is terminated after NCOMMU commutators
* (initial version using CI behind the curtains)
*
* Jeppe Olsen, August 2005
*              Latest modification : September 2005, IT2REFSPC added
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
*. Specific input
      DIMENSION T(*)
*. Output
      DIMENSION CCVECFNC(*),CCVECFNCI(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from MRCC_VECFNC'
        WRITE(6,*) ' -----------------------'
        WRITE(6,*) ' IREFSPC,ITREFSPC, IT2REFSPC =',
     &               IREFSPC,ITREFSPC, IT2REFSPC
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MRCCVF')
*
* 1 : Obtain exp(-T) H exp(T)  !0> and save on LUHC
*
C          EMNTHETO(T,LUOUT,NCOMMU,IREFSPC,ITREFSPC)
      IF(I_APPROX_HCOM.EQ.0) THEN
        CALL EMNTHETO(T,LUC,LUHC,NCOMMU,IREFSPC,ITREFSPC,IT2REFSPC)
      ELSE
*. Exact calculation of all terms with upto NCOMMU-1 commutators
        CALL EMNTHETO(T,LUC,LUHC,NCOMMU-1,IREFSPC,ITREFSPC,IT2REFSPC)
*. and add contribution from highest commutaror
*. At the moment FULL Hamiltonian is used for testing 
COLD    WRITE(6,*) ' Note : Full Hamiltonian is used in highest commu'
*. Use zero-order Hamiltonian stored in 
        I12 = 1
        CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
        CALL TCOM_H_N(T,LUC,LUHC,NCOMMU,IREFSPC,ITREFSPC,IT2REFSPC,1)
C            TCOM_H_N(T,LUINI,LUUT,NCOMMU,IREFSPC,ITREFSPC,IT2REFSPC,IAC)
        I12 = 2
        CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
      END IF
*
* 2 : Obtain  <0!tau^{\dagger} exp(-T) H exp(T) !0> = <LUC!tau^{\dagger}|LUHC>
*
      ICSPC = IREFSPC
      ISSPC = IT2REFSPC
C     WRITE(6,*) ' IREFSPC, IT2REFSPC =', IREFSPC, IT2REFSPC
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Vector on LUC '
        CALL WRTVCD(WORK(KVEC1P),LUC,1,-1)
        WRITE(6,*) ' Vector on LUHC '
        CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
      END IF
*
      ZERO = 0.0D0
      CALL SETVEC(CCVECFNC,ZERO,N_CC_AMP)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,CCVECFNC,2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 'CCVECFNC right after SIGDEN_CC'
        CALL WRTMAT(CCVECFNC,1,N_CC_AMP,1,N_CC_AMP)
      END IF
      
*
* 3 : Contract  exp(-T) H exp(T) |0> to reference space and save on LUHC
*     to obtain internal part of MRCC vector function 
*
      CALL EXPCIV(IREFSM,IT2REFSPC,LUHC,IREFSPC,LUSC34,-1,
     /            LUSC35,1,1,IDC,0)
      CALL REWINO(LUHC)
      CALL FRMDSCN(CCVECFNCI,-1,-1,LUHC)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input T-coefficients '
        CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
        WRITE(6,*) ' MRCC Vector function, external part  '
        CALL WRTMAT(CCVECFNC,1,N_CC_AMP,1,N_CC_AMP)
        WRITE(6,*) 'first element of MRCC Vector function,internal part'
        WRITE(6,*) ' (before subtracting E-term )'
        CALL WRTMAT(CCVECFNCI,1,1,1,1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'MRCCVF')
*
     
      RETURN
      END 
      SUBROUTINE EMNTHETO(T,LUINI,LUUT,NCOMMU,IREFSPC,ITREFSPC,
     &                    IT2REFSPC)
*
*. Obtain on LUOUT exp(-T) H exp(T)  !0>, truncated after NCOMMU commutators
*. Input in CAAB basis
*  Output on LUOT in SD basis 
*. LUUT should differ from scratch files used below, one possible choice is LUHC
*. Scratch files in use : LUSC1, LUSC2, LUSC3, LUSC34
*. Jeppe Olsen, August 2005
*
* The three spaces : IREFSPC : Space of !0>
*                    ITREFSPC : Largest space required for the calculation of
*                               exp(-T) H exp(T)  !0>
*                    IT2REFSPC : Space for T !0>
*. Final vector is delivered in space IT2REFSPC
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
*
*. Specific input
      DIMENSION T(*)
*. We are after Sum(N=0,Ncommu,i=0,N)(-1)^(N-I) 1/N! T^(N-I) H T^I |0>
*. So realize the calculation as a double loop
*
      NTEST = 00
      IF(NTEST.GE.10) THEN 
        WRITE(6,*) ' exp(-T) H Exp(T) |0> will be constructed '
        WRITE(6,*) ' Input T-coefficients '
        CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
        WRITE(6,*) ' EMNTHETO: IREFSPC, ITREFSPC, IT2REFSPC =',
     &  IREFSPC, ITREFSPC, IT2REFSPC
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EMNTH ')
*
* LUINI : Initial expansion |0>
* LUSC1 : T^I |0>
* LUSC2  : H T^I |0>
* LUSC3 : T^N-I H T^I |0>
*
      ONE = 1.0D0
*
      DO I = 0, NCOMMU
        ICSPC = ITREFSPC
        ISSPC = ITREFSPC
        IF(I.EQ.0) THEN
*. Expand |0> in IREFSPC on LUINI to ITREFSPC on LUSC1
           CALL EXPCIV(IREFSM,IREFSPC,LUINI,ITREFSPC,LUSC1,-1,
     /                   LUSC34,1,0,IDC,NTEST)
C       EXPCIV(ISM,ISPCIN,LUIN,
C    &                  ISPCUT,LUUT,LBLK,
C    &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
        ELSE 
*T^(I-1)|0> => T^I |0>
         CALL REWINO(LUSC1)
         CALL REWINO(LUSC2)
         CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUSC34,T,1)
         CALL COPVCD(LUSC34,LUSC1,WORK(KVEC1P),1,-1)
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' T^I |0> for I = ',I
          CALL WRTVCD(WORK(KVEC1P),LUSC1,1,-1)
        END IF
*. Calculate H T^I |0> and save on LUSC2
*. Space of H T^I |0> may be reduced to IT2REFSPC
        ICSPC = ITREFSPC
        ISSPC = IT2REFSPC
C?      WRITE(6,*) ' MV7 will be called with ISSPC=IT2REFSPC'
        CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUSC2,0,0)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Output from MV7'
          CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
        END IF
*. Compress Sigma-vector to space IT2REFSPC
C      WRITE(6,*) ' sigma vector will be contracted to IT2REFSPC'
C         CALL EXPCIV(1,ITREFSPC,LUSC2,IT2REFSPC,LUSC3,-1,
C    &                   LUSC34,1,1,IDC,NTEST)
        
*. C space may now also be restricted to IT2REFSPC
        ISSPC = IT2REFSPC
        ICSPC = IT2REFSPC
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' H T^I |0> for I = ',I
           CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
         END IF
        DO NMI = 0, NCOMMU-I
          IF(NMI.EQ.0) THEN
*. Just copy H T^I |0> to LUSC3
           CALL COPVCD(LUSC2,LUSC3,WORK(KVEC1P),1,-1)
          ELSE 
*. Calculate T^(N-I) H T^I |0> and save on LUSC3
           REWIND(LUSC3)
           REWIND(LUSC34)
           CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC34,T,1)
           CALL COPVCD(LUSC34,LUSC3,WORK(KVEC1P),1,-1)
          END IF
          IF(NTEST.GE.1000) THEN 
            WRITE(6,*) '  T^(N-I) H T^I for N-I and I ', NMI,I
            CALL WRTVCD(WORK(KVEC1P),LUSC3,1,-1)
          END IF
* We are now ready to add (-1)^(N-I) 1/N! T^(N-I) H T^I |0> to result vector
          N = NMI  + I
          IF(NMI.EQ.0) THEN 
            XNMIFAC = 1.0D0
          ELSE 
            XNMIFAC = XFAC(NMI)
          END IF
          IF(I.EQ.0) THEN
            XIFAC = 1.0D0
          ELSE 
            XIFAC = XFAC(I)
          END IF
          IF(MOD(NMI,2).EQ.0) THEN
           FACTOR = 1.0D0/(XNMIFAC*XIFAC)
          ELSE 
           FACTOR = -1.0D0/(XNMIFAC*XIFAC)
          END IF
*. First contribution : Just copy (factor is 1)
          IF(N.EQ.0) THEN
            CALL COPVCD(LUSC3,LUUT,WORK(KVEC1P),1,-1)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Initial vector copied to LUUT '
              CALL WRTVCD(WORK(KVEC1P),LUUT,1,-1)
            END IF
          ELSE 
* add : LUUT = LUUT + FACTOR*LUSC3
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
           CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),FACTOR,ONE,LUSC3,LUUT,
     &                 LUSC34,1,-1)
           CALL COPVCD(LUSC34,LUUT,WORK(KVEC1P),1,-1)
           IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' LUUT opdated for I, NMI = ', I,NMI
             CALL WRTVCD(WORK(KVEC1P),LUUT,1,-1)
           END IF
          END IF
        END DO
*       ^ End of loop over NMI
      END DO
*     ^ End of loop over I
*
*. Test Contract from ITREFSPC to IT2REFSPC, save on LUSC34
*
C?    WRITE(6,*) ' Output vector will be contracted to IT2REFSPC'
C?         CALL EXPCIV(1,ITREFSPC,LUUT,IT2REFSPC,LUSC1,-1,
C?   &                   LUSC34,1,1,IDC,NTEST)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' exp(-T) H exp(T) |0> '
        CALL WRTVCD(WORK(KVEC1P),LUUT,1,-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EMNTH ')
      RETURN
      END 
      SUBROUTINE COM_JMRCC(T,NCOMMU,I_APPROX_HCOM,XJ,VCC1,VCC2,VCC3,
     &                  VCC4,
     &                  N_CC_AMP,NSPAM1,NNONSING,IREFSPC,ITREFSPC,
     &                  XNONSING)
*
* Construct - by finite difference - the MRCC Jacobian for current 
* set of amplitudes 
*
* For the finite difference the following form is used 
*
* F' = (8*F(DELTA)-8*F(-DELTA)-E(2*DELTA)+E(2*DELTA))/(12*DELTA)
*
* The Jacobian will be returned in the Nonsingular basis as 
* defined by XNONSING.
*
* Jeppe Olsen, Aug. 2005
*
* Latest modification : Sept 2005, New form of call to MRCC_VECFNC
*
*
      INCLUDE 'wrkspc.inc'
      REAL*8 INPRDD
* 
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
*. Input 
      DIMENSION T(*), XNONSING(NSPAM1,NNONSING)
*. T is on input assumed to be in CAAB basis !
*. Output
      DIMENSION XJ(NNONSING,NNONSING)
*. Scratch
      DIMENSION VCC1(*),VCC2(*),VCC3(*),VCC4(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'COM_JM')
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' COM_JMRCC speaking '
         WRITE(6,*) ' IREFSPC, ITREFSPC = ', IREFSPC,ITREFSPC
      END IF
*. CC vector function at point of expansion in VCC2
      CALL MRCC_VECFNC(VCC2,T,NCOMMU,I_APPROX_HCOM,IREFSPC,ITREFSPC,
     &     ITREFSPC,VCC2(1+N_CC_AMP))
*. Transform to SPA basis and save in VCC1
      CALL REF_CCV_CAAB_SP(VCC2,VCC1,VCC3,1)
*. and to orthonormal basis, save in VCC1
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(XNONSING,VCC2,VCC1,NSPAM1,NNONSING,1)
*. Dimension of space in which S or H is constructed 
      DELTA = 0.0001D0
      DELTAM = -DELTA
      DELTA2= 2.0D0*DELTA
      DELTA2M = -2.0D0*DELTA
*
      ONE = 1.0D0
      ONEM = -1.0D0
      EIGHT = 8.0D0
      EIGHTM = -8.0D0
C     DO I = 1, NSPAM1
      DO I = 1, NNONSING
       IF(NTEST.GE.10)
     & WRITE(6,*) ' Jacobian will be constructed, column = ', I
*. Transform I'th direction to CAAB basis and save in VCC1
       CALL REF_CCV_CAAB_SP(VCC1,XNONSING(1,I),VCC2,2)
* ===================
* a : 8*vecfnc(Delta)
* ===================
*. ( T + delta Xnonsing(*,I)) in VCC2
       CALL VECSUM(VCC2,VCC1,T,DELTA,ONE,N_CC_AMP)
*. Vecfnc( T + delta Xnonsing(*,I)) in VCC3
       CALL MRCC_VECFNC(VCC3,VCC2,NCOMMU,I_APPROX_HCOM,IREFSPC,ITREFSPC,
     &     ITREFSPC,VCC3(1+N_CC_AMP))
*. Transform to SPA  basis and save in VCC2
       CALL REF_CCV_CAAB_SP(VCC3,VCC2,VCC4,1)
C             REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
*. and to orthonormal basis, save in VCC3
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
       CALL MATVCC(XNONSING,VCC2,VCC3,NSPAM1,NNONSING,1)
*. Save 8*Vecfnc(Delta*X(I)) in XJ(1,I)
       CALL COPVEC(VCC3,XJ(1,I),NNONSING)
       CALL SCALVE(XJ(1,I),EIGHT,NNONSING)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' XJ(1,I), first term '
         CALL WRTMAT(XJ(1,I),1,NSPAM1,1,NSPAM1)
       END IF
* ===================
* b : 8*vecfnc(-Delta)
* ===================
*. ( T - delta Xnonsing(*,I)) in VCC2
       CALL VECSUM(VCC2,VCC1,T,DELTAM,ONE,N_CC_AMP)
*. Vecfnc( T - delta Xnonsing(*,I)) in VCC3
       CALL MRCC_VECFNC(VCC3,VCC2,NCOMMU,I_APPROX_HCOM,IREFSPC,ITREFSPC,
     &     ITREFSPC,VCC3(1+N_CC_AMP))
*. Transform to SPA  basis and save in VCC2
       CALL REF_CCV_CAAB_SP(VCC3,VCC2,VCC4,1)
C             REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
*. and to orthonormal basis, save in VCC3
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
       CALL MATVCC(XNONSING,VCC2,VCC3,NSPAM1,NNONSING,1)
*. Save 8*Vecfnc(Delta*X(I)) in XJ(1,I)
       CALL VECSUM(XJ(1,I),XJ(1,I),VCC3,ONE,EIGHTM,NNONSING)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' XJ(1,I), second term '
         CALL WRTMAT(XJ(1,I),1,NSPAM1,1,NSPAM1)
       END IF
* ===================
* c : vecfnc(2*Delta)
* ===================
*. ( T +2*delta Xnonsing(*,I)) in VCC2
       CALL VECSUM(VCC2,VCC1,T,DELTA2,ONE,N_CC_AMP)
*. Vecfnc( T +2*delta Xnonsing(*,I)) in VCC3
       CALL MRCC_VECFNC(VCC3,VCC2,NCOMMU,I_APPROX_HCOM,IREFSPC,ITREFSPC,
     &     ITREFSPC,VCC3(1+N_CC_AMP))
*. Transform to SPA  basis and save in VCC2
       CALL REF_CCV_CAAB_SP(VCC3,VCC2,VCC4,1)
C             REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
*. and to orthonormal basis, save in VCC3
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
       CALL MATVCC(XNONSING,VCC2,VCC3,NSPAM1,NNONSING,1)
*. add -Vecfnc(2Delta*X(I)) in XJ(1,I)
       CALL VECSUM(XJ(1,I),XJ(1,I),VCC3,ONE,ONEM,NNONSING)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*)  ' XJ(1,I), third term '
         CALL WRTMAT(XJ(1,I),1,NSPAM1,1,NSPAM1)
       END IF
* ===================
* d : vecfnc(-2*Delta)
* ===================
*. ( T - 2*delta Xnonsing(*,I)) in VCC2
       CALL VECSUM(VCC2,VCC1,T,DELTA2M,ONE,N_CC_AMP)
*. Vecfnc( T +2*delta Xnonsing(*,I)) in VCC3
       CALL MRCC_VECFNC(VCC3,VCC2,NCOMMU,I_APPROX_HCOM,IREFSPC,ITREFSPC,
     &     ITREFSPC,VCC3(1+N_CC_AMP))
*. Transform to SPA  basis and save in VCC2
       CALL REF_CCV_CAAB_SP(VCC3,VCC2,VCC4,1)
C             REF_CCV_CAAB_SP(VEC_CAAB,VEC_SP,VEC_SCR,IWAY)
*. and to orthonormal basis, save in VCC3
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
       CALL MATVCC(XNONSING,VCC2,VCC3,NSPAM1,NNONSING,1)
*. add Vecfnc(-2Delta*X(I)) in XJ(1,I)
       CALL VECSUM(XJ(1,I),XJ(1,I),VCC3,ONE,ONE,NNONSING)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' XJ(1,I), Fourth term '
         CALL WRTMAT(XJ(1,I),1,NSPAM1,1,NSPAM1)
       END IF
*. and scale 
       FACTOR = 1.0D0/(12.0D0*DELTA)
       CALL SCALVE(XJ(1,I),FACTOR,NNONSING)
      END DO
*     ^ End of loop over nonsingular modes
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Constructed Jacobian matrix '
        WRITE(6,*) ' ==================== '
        CALL WRTMAT(XJ,NNONSING,NNONSING,NNONSING,NNONSING)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'COM_JM')
*
      RETURN
      END 
      SUBROUTINE LINSOL_FROM_LUCOMP(XL,XU,RHS,X,NDIM,SCR1)
*
* Solve linear set of equations from matrix given in LU decomposition
*
*  L U X = RHS
*
* is solved in two steps
*
* 1)  L Y = RHS
* 2)  U X = Y
*
* Jeppe Olsen, Aug. 2005
*
* LU are given in the form defined by routine LULU
*
*    L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
*    U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION XL(NDIM*(NDIM+1)/2), XU(NDIM*(NDIM+1)/2), RHS(NDIM)
*. Output
      DIMENSION X(NDIM)
*. Scratch 
      DIMENSION SCR1(NDIM)
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' LINSOL_FROM_LUCOM speaking '
      END IF
*
* 1 : L Y = RHS by forward substitution  and store in SCR1
*
      DO I = 1, NDIM
*. sum(J=1,I-1) L(I,J)Y(J)  
        SUM = 0.0D0
        DO J = 1, I-1
          SUM = SUM + XL(I*(I-1)/2+J)*SCR1(J)
        END DO
        SCR1(I) = (RHS(I)-SUM)/XL(I*(I-1)/2+I)
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Solution to L Y = RHS '
        CALL WRTMAT(SCR1,1,NDIM,1,NDIM)
      END IF
*
* 2 : U X = Y by backwards substitution 
*
      DO I = NDIM, 1, -1
*. sum(J=I+1,NDIM) U(I,J)*X(J)
        SUM = 0.0D0
        DO J = I+1, NDIM
          SUM = SUM + XU(J*(J-1)/2+I)*X(J)
        END DO
        X(I) = (SCR1(I)-SUM)/XU(I*(I-1)/2+I)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' RHS '
        CALL WRTMAT(RHS,1,NDIM,1,NDIM)
        WRITE(6,*) ' Solution to set of linear equations '
        CALL WRTMAT(X,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END 
      SUBROUTINE ICCC_RELAX_REFCOEFS_COM(T_EXT,H_REF,N_REF,
     &           NCOMMU,VEC1,VEC2,IREFSPC,ITREFSPC,
     &           ECORE,C_REF_OUT,IREFROOT)
*
*
* Relax internal coefficients for MRCC wave function 
* Initial version generating complete matrices 
*
* The wave-function is given as 
*
* |MRCC > = exp(T) |0 >
*
* and we want to solve the equations
*
* sum_J <I!exp(-T)H exp(T)!J> C(J) = E C(J)
*
* ( note that the metric disappears )
*
*. Jeppe Olsen, August 2005
* NOTE : ONLY PROGRAMMED FOR LOWEST ROOT - Easy to modify ...
      INCLUDE 'wrkspc.inc'
      REAL*8 INPRDD, INPROD
*
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
*. Input : in CAAB form 
      DIMENSION T_EXT(*)
*. Output
      DIMENSION H_REF(N_REF,N_REF)
      DIMENSION C_REF_OUT(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'CC_REL')
*
      WRITE(6,*) ' Code has should be modified to new MRCC vecfnc '
      STOP ' Code has should be modified to new MRCC vecfnc '
*
      NTEST = 10
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' Reoptimization of internal coefficients'
        WRITE(6,*) ' =======================================' 
        IF(IDIIS.EQ.1) THEN
          WRITE(6,*) ' DIIS acceleration will be used '
        END IF
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' IREFSPC, ITREFSPC ', IREFSPC,ITREFSPC
        WRITE(6,*) ' IREFROOT = ', IREFROOT
      END IF
*
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
*
      DO J = 1, N_REF
        IF(NTEST.GE.10) WRITE(6,*) ' Column J = ', J
*. Place |J> on LUSC1 
       CALL REWINO(LUSC36)
       CALL REWINO(LUDIA)
C  WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
       ONE = 1.0D0
       CALL WRSVCD(LUSC36,-1,VEC1,J,ONE,1,N_REF,LUDIA,1)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Input vector on LUSC36'
         CALL WRTVCD(VEC1,LUSC36,1,-1)
       END IF
*. 
*
*. Obtain exp(-T) H exp(T) |J>  on LUHC
C     EMNTHETO(T,LUINI,LUOUT,NCOMMU,IREFSPC,ITREFSPC)
       CALL EMNTHETO(T_EXT,LUSC36,LUHC,NCOMMU,IREFSPC,ITREFSPC,ITREFSPC)
*. Contract  exp(-T) H exp(T) |J> to reference space and save on LUHC
       CALL EXPCIV(IREFSM,ITREFSPC,LUHC,IREFSPC,LUSC34,-1,
     /             LUSC35,1,1,IDC,0)
*. and read in - the J'th column of H_REF has been constructed
       CALL REWINO(LUHC)
       CALL FRMDSCN(H_REF(1,J),-1,-1,LUHC)
C      FRMDSCN(VEC,NREC,LBLK,LU)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Effective H-matrix in reference space '
        CALL WRTMAT(H_REF,N_REF,N_REF,N_REF,N_REF)
      END IF
*
** And diagonalize 
*
C       EIGGMTN(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,Z,W,SCR)
      CALL MEMMAN(KLEIGVA_R,N_REF   ,'ADDL  ',2,'EIGVAR')
      CALL MEMMAN(KLEIGVA_I,N_REF   ,'ADDL  ',2,'EIGVAI')
      CALL MEMMAN(KLEIGVC_R,N_REF**2,'ADDL  ',2,'EIGVCR')
      CALL MEMMAN(KLEIGVC_I,N_REF**2,'ADDL  ',2,'EIGVCI')
      CALL MEMMAN(KLZ,N_REF**2,'ADDL  ',2,'Z_SCR ')
      CALL MEMMAN(KLW,N_REF   ,'ADDL  ',2,'W_SCR ')
      CALL MEMMAN(KLSCR    ,2*N_REF   ,'ADDL  ',2,'EIGSCR')
      CALL EIGGMTN(H_REF,N_REF,WORK(KLEIGVA_R),WORK(KLEIGVA_I),
     &             WORK(KLEIGVC_R),WORK(KLEIGVC_I),
     &             WORK(KLZ),WORK(KLW),WORK(KLSCR))
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Real and imaginary parts of eigenvalues '
        DO I = 1, N_REF
          WRITE(6,*) I,WORK(KLEIGVA_R-1+I),WORK(KLEIGVA_I-1+I)
        END DO
       END IF
*. Lowest eigenvalue - should really be eigenvalue IREFROOT - here 
*. are the bits of codes that should be generalized to general roots
      IMIN = 1
      EIGMIN = WORK(KLEIGVA_R-1+1)
      DO I = 2, N_REF
        IF( WORK(KLEIGVA_R-1+I).LT.EIGMIN) THEN
          EIGMIN = WORK(KLEIGVA_R-1+I)
          IMIN = I
         END IF
      END DO
      WRITE(6,*) ' Root with lowest energy ', IMIN,EIGMIN
      IF(WORK(KLEIGVA_I-1+IMIN).NE.0.0D0) THEN
        WRITE(6,*) ' Warning : Complex eigenvalue '
        WRITE(6,*) ' Real and imaginary parts ', 
     &  WORK(KLEIGVA_R-1+IMIN),WORK(KLEIGVA_I-1+IMIN)
        STOP ' Complex eigenvalue '
      END IF
*. Copy the coefficients of root IROOT to C_REF_OUT  
      CALL COPVEC(WORK(KLEIGVC_R+(IMIN-1)*N_REF),C_REF_OUT,N_REF)
*. Ensure standard normalization
      XNORM = SQRT(INPROD(C_REF_OUT,C_REF_OUT,N_REF))
      FACTOR = 1.0D0/XNORM
      CALL SCALVE(C_REF_OUT,FACTOR,N_REF)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated coefficients of reference state'
        CALL WRTMAT(C_REF_OUT,1,N_REF,1,N_REF)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'CC_REL')
*
      RETURN
      END 
      SUBROUTINE HEFF_INT_TV_ICCC(T_EXT,N_REF,
     &           NCOMMU,IAPROX_HCOM,VEC1,VEC2,IREFSPC,ITREFSPC,
     &           IT2REFSPC,ECORE,C_REF,S_REF)
*
*. Calculate Heff times vector in reference space for ICCC
*
*. S_REF = <I!exp(-T)H exp(T)|0>
*. where |0> is defined by C_REF
*
*
*. Jeppe Olsen, August 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
*. Input : in CAAB form 
      DIMENSION T_EXT(*)
      DIMENSION C_REF(*)
*. Output
      DIMENSION S_REF(N_REF)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*. Files in use pt : LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUHC
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'HEFFCC')
*
      NTEST = 0
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' Calculation of gradient for reference dets '
        WRITE(6,*) ' ===========================================' 
        WRITE(6,*) ' IREFSPC, ITREFSPC ', IREFSPC,ITREFSPC
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input C_REF '
        CALL WRTMAT(C_REF,1,N_REF,1,N_REF)
      END IF
*
      ICSPC = IREFSPC
      ISSPC = ITREFSPC
*
*. transfer new reference vector to file LUSC34 - use S_REF as integer scratch
*. and LUDIA as form 
      CALL ISTVC2(S_REF,0,1,N_REF)
C  WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
      CALL REWINO(LUSC34)
      CALL REWINO(LUDIA)
      CALL WRSVCD(LUSC34,-1,VEC1,S_REF,
     &     C_REF,N_REF,N_REF,LUDIA,1)
C  WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
*. Obtain exp(-T) H exp(T) |0>  on LUHC
      IF(IAPROX_HCOM.EQ.0) THEN
*. No approximations in highest commutator
        CALL EMNTHETO(T_EXT,LUSC34,LUHC,NCOMMU,IREFSPC,ITREFSPC,
     &                IT2REFSPC)
      ELSE 
        CALL EMNTHETO(T_EXT,LUSC34,LUHC,NCOMMU-1,IREFSPC,ITREFSPC,
     &                IT2REFSPC)
*. PT full Hamiltonian is used for testing
        CALL TCOM_H_N(T_EXT,LUSC34,LUHC,NCOMMU,IREFSPC,
     &                ITREFSPC,IT2REFSPC,1)
      END IF
*. Contract  exp(-T) H exp(T) |0> to reference space and save on LUHC
      CALL EXPCIV(IREFSM,ITREFSPC,LUHC,IREFSPC,LUSC34,-1,
     /            LUSC35,1,1,IDC,0)
      CALL REWINO(LUHC)
      CALL FRMDSCN(S_REF,-1,-1,LUHC)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Heff times vector in internal space  '
        CALL WRTMAT(S_REF,1,N_REF,1,N_REF)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'HEFFCC')
*
      RETURN
      END 
      SUBROUTINE ICCC_OPT_SIMULT(
     &        IREFSPC,ITREFSPC,IT2REFSPC,I_SPIN_ADAPT,
     &        IREFROOT,T_EXT,C_0,INI_IT,IFIN_IT,VEC1,VEC2,IDIIS,
     &        C_REF,N_REF,I_DO_COMP,CONVERL,VTHRES,I_REDO_INT,
     &        EFINAL,VNFINAL,CONVERG,SCR_SBSPJA,MXVEC_SBSPJA)

*
* Master routine for Internal Contraction Coupled Cluster 
*
* It is assumed that the excitation manifold produces 
* states that are orthogonal to the reference so 
* no projection is carried out
*
* Routine is allowed to leave without turning the lights off,
* i.e. leave routine with all allocations and marks intact.
*: Thus : Allocations are only done if INI_IT = 1
*        Deallocations are only done if IFIN_IT = 1
*
*. Preconditioners are only calculated if INI_IT = 1
*
* IF I_REDO_INT = 1, the internal states are recalculated at start
*
* IF IDIIS.EQ.1, DIIS is used
*         .EQ.2, CROP is used to accelerate convergence 
* 
*
* Jeppe Olsen, Aug. 2005, modified aug 2009 - also in Washington
*              Redo of internal states: Sept. 2009 in Sicily
*              Subspace Jacobian added: Oct. 2009
*
*
*. for DIIS units LUSC37 and LUSC36 will be used for storing vectors
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Temporary  array for debugging
      REAL*8 XNORM_EI(1000)
*
      LOGICAL CONVERL,CONVERG
*. Converl: is local iterative procedure for given internal states converged
*. converg: is global iterative procedure converged
      REAL*8
     &INPROD
*. Input and Output : Coefficients of internal and external correlation 
      DIMENSION T_EXT(*), C_REF(*)
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
      COMMON/CLOCAL2/KVEC1,KVEC2,MXCJ,
     & KLVCC1,KLVCC2,KLVCC3,KLVCC4,KLVCC5,KLSMAT,KLXMAT,KLJMAT,KLU,KLL,
     & NSING,NNONSING,KLCDIIS,KLC_INT_DIA,KLDIA,KLVCC6,KLVCC7,KLVCC8,
     & NVECP,NVEC,KLA_CROP,KLSCR_CROP
*. Scratch for CI behind the curtain 
      DIMENSION VEC1(*),VEC2(*)
*. Scratch for subspace Jacobian
      DIMENSION SCR_SBSPJA(*)
*. Threshold for convergence of norm of Vectorfuntion

C     WRITE(6,*) ' ICCC_OPT_SIMULT: I_DO_COMP =', I_DO_COMP
C     WRITE(6,*) ' ICCC_OPT_SIMULT: MAXIT,MAXITM =', MAXIT,MAXITM
      WRITE(6,*) ' ICCC_OPT_SIMULT: I_DO_SBSPJA, MXVEC_SBSPJA = ', 
     &                              I_DO_SBSPJA, MXVEC_SBSPJA
*. Number of Spin adapted functions (and NCAAB for a check)
      NSPA = N_ZERO_EI 
      NCAAB = NDIM_EI
      WRITE(6,*) ' NCAAB og NDIM_EI = ', NCAAB, NDIM_EI
*. We will not include the unit-operator so  ???
      NSPAM1 = NSPA - 1
*. Different adresses of the unit op
      IF(I_DO_EI.EQ.0) THEN
        IUNI_AD = 1
      ELSE
        IUNI_AD = NCAAB 
      END IF
*. Freeze internal expansion
CM    I_FIX_INTERNAL = 0
*. Project on nonredundant space
      I_DO_PROJ_NR = 1
*. For file access
      LBLK = -1
      NTEST = 5
      IF(NTEST.GE.2) THEN
      WRITE(6,*) 
     &  ' Simultaneous optimization of internal and external parts '
        WRITE(6,*) 
     &  ' ========================================================='
        WRITE(6,*)
        WRITE(6,*) ' Reference space is ', IREFSPC
        WRITE(6,*) ' Space for evaluating general operators  ', ITREFSPC
        WRITE(6,*) ' Space for T times reference space  ', IT2REFSPC
        WRITE(6,*) ' Number of parameters in CAAB basis ', 
     &             N_CC_AMP
        WRITE(6,*) ' Number of parameters in spincoupled/ort basis ', 
     &             NSPA
        WRITE(6,*) ' Number of coefficients  in internal space ', N_REF
        WRITE(6,*) ' INI_IT, IFIN_IT = ', INI_IT, IFIN_IT
        WRITE(6,*) ' Max. number microiterations per macro ', MAXIT
        WRITE(6,*) ' Max. number of macroiterations        ', MAXITM
        WRITE(6,*) ' Number of vectors allowed in subspace ', MXCIVG
        WRITE(6,*) ' Number of vectors allowed in initial subspace ', 
     &               MXVC_I
        IF(IDIIS.EQ.1) THEN
          WRITE(6,*)' DIIS optimization'
        ELSE IF (IDIIS.EQ.2) THEN
          WRITE(6,*)' CROP optimization'
        END IF
*
        IF(I_DO_PROJ_NR.EQ.1) THEN
          WRITE(6,*) ' Redundant directions projected out'
        ELSE
          WRITE(6,*) ' No projection of redundant directions'
        END IF
*
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial T_ext-amplitudes '
        CALL WRTMAT(T_EXT,1,N_CC_AMP,1,N_CC_AMP)
        WRITE(6,*) ' Initial C_int-amplitudes '
        CALL WRTMAT(C_REF,1,N_REF,1,N_REF)
      END IF
*. Allowed number of iterations
      NNEW_MAX = MAXIT
      MAXITL = NNEW_MAX
*
      NVAR = N_CC_AMP + N_REF
      IF(INI_IT.EQ.1) THEN
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICC_CM')
        CALL MEMMAN(KLVCC1,NVAR,'ADDL  ',2,'VCC1  ')
        CALL MEMMAN(KLVCC2,NVAR,'ADDL  ',2,'VCC2  ')
        CALL MEMMAN(KLVCC3,NVAR,'ADDL  ',2,'VCC3  ')
        CALL MEMMAN(KLVCC4,NVAR,'ADDL  ',2,'VCC4  ')
        CALL MEMMAN(KLVCC5,NVAR,'ADDL  ',2,'VCC5  ')
        CALL MEMMAN(KLVCC6,2*NVAR,'ADDL  ',2,'VCC6  ')
*. Just a few extra to be on the safe side when programming EI
*. approach
        CALL MEMMAN(KLVCC7,NVAR,'ADDL  ',2,'VCC5  ')
        CALL MEMMAN(KLVCC8,NVAR,'ADDL  ',2,'VCC5  ')
*. Complete matrices for external part, three used pt
        LEN = NSPA**2
        IF(I_DO_COMP.EQ.1) THEN
          CALL MEMMAN(KLSMAT,LEN,'ADDL  ',2,'SMAT  ')
          CALL MEMMAN(KLXMAT,LEN,'ADDL  ',2,'XMAT  ')
          CALL MEMMAN(KLJMAT,LEN,'ADDL  ',2,'JMAT  ')
*. Storage for LU decomposition of J
          LEN = NSPA*(NSPA+1)/2
          CALL MEMMAN(KLL,LEN,'ADDL  ',2,'L     ')
          CALL MEMMAN(KLU,LEN,'ADDL  ',2,'U     ')
        ELSE
*. Space for diagonal- space is allocated also for CI part.
          CALL MEMMAN(KLDIA,NVAR+1,'ADDL  ',2,'DIAORT')
        END IF
*. Space for DIIS/CROP
        IF(IDIIS.EQ.1) THEN
          CALL MEMMAN(KLCDIIS,MAXITL,'ADDL ',2,'CDIIS ') 
        ELSE IF(IDIIS.EQ.2) THEN
          CALL MEMMAN(KLA_CROP,MXCIVG*(MXCIVG+1)/2,'ADDL  ',2,'A_CROP')
          LEN_SCR_CROP = 3*MXCIVG*MXCIVG + 3*MAX(MXCIVG,NVAR)
          CALL MEMMAN(KLSCR_CROP,LEN_SCR_CROP,'ADDL  ',2,'S_CROP')
C?        WRITE(6,*) ' KLA_CROP,KLSCR_CROP, a =', KLA_CROP,KLSCR_CROP
        END IF
*. Space Diagonal for internal part
        CALL MEMMAN(KLC_INT_DIA,N_REF,'ADDL ',2,'C_DIA ')
      END IF
*.    ^ End if INI_IT.EQ.1
*
*======================================
* 0: Redo internal states if required
* =====================================
*
      IF(I_REDO_INT.EQ.1) THEN
        CALL GET_INTERNAL_STATES(N_EXTOP_TP,N_INTOP_TP,
     &     WORK(KLSOBEX),WORK(KL_N_INT_FOR_EXT),WORK(KL_IB_INT_FOR_EXT),
     &     WORK(KL_I_INT_FOR_EXT),WORK(KL_NDIM_IN_SE),
     &     WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &     WORK(KL_X1_INT_EI_FOR_SE), WORK(KL_X2_INT_EI_FOR_SE),
     &     WORK(KL_SG_INT_EI_FOR_SE),WORK(KL_S_INT_EI_FOR_SE),
     &     WORK(KL_IBX1_INT_EI_FOR_SE), WORK(KL_IBX2_INT_EI_FOR_SE),
     &     WORK(KL_IBSG_INT_EI_FOR_SE),WORK(KL_IBS_INT_EI_FOR_SE),
     &     WORK(KL_X2L_INT_EI_FOR_SE),
     &     I_IN_TP,I_INT_OFF,I_EXT_OFF)
*
C IMNNMX(IVEC,NDIM,MINMAX)
        N_INT_MAX = IMNMX(WORK(KL_N_INT_FOR_SE),N_EXTOP_TP*NSMOB,2)
*. Largest number of zero-order states of given sym and external type
        N_ORTN_MAX = IMNMX(WORK(KL_N_ORTN_FOR_SE),N_EXTOP_TP*NSMOB,2)
        WRITE(6,*) ' N_INT_MAX, N_ORTN_MAX = ', N_INT_MAX, N_ORTN_MAX
*. Largest transformation block 
        N_XEO_MAX = N_INT_MAX*N_ORTN_MAX
        IF(NTEST.GE.10) 
     &  WRITE(6,*) ' Largest (EL,ORTN) block = ', N_XEO_MAX
*. Number of zero-order states - does now include the unit-operator
        N_ZERO_EI = N_ZERO_ORDER_STATES(WORK(KL_N_ORTN_FOR_SE),
     &             WORK(KL_NDIM_EX_ST),N_EXTOP_TP,1)
        NSPA = N_ZERO_EI
       IF(NTEST.GE.10) WRITE(6,*) 
     & ' Number of zero-order states with sym 1 = ', N_ZERO_EI
      END IF
*
* ============================================================
* 1 : Prepare preconditioners for external and internal parts 
* ============================================================
*
* --------------------
*. 1a : External part 
* --------------------
*
*. Identify the unit  operator i.e. the operator with 
*. zero creation and annihilation operators
      IDOPROJ = 0
*. Construct metric (once again ..)
*. Prepare the routines used in COM_SH
*. Not used here
      C_0X = 0.0D0
      KLTOPX = -1
*. Used
      NREFX = N_REF
      IREFSPCX = IREFSPC
*. Space to be used for evaluating metric : If T = 0, then IT2REFSPC is sufficient
      ITREFSPCX = ITREFSPC
      ITREFSPCX = IT2REFSPC
*
      NCAABX = N_CC_AMP
      NSPAX = NSPA
      IPROJSPCX = IREFSPC
*. Unitoperator in SPA order ... Please check ..
      IUNIOPX = 0
      IF(I_DO_COMP.EQ.1) THEN
*. Set up or read in complete matrices
        IF(INI_IT.EQ.1.AND.IREADSJ.EQ.0) THEN
          CALL COM_SH(WORK(KLSMAT),WORK(KLSMAT),WORK(KLVCC1),
     &         WORK(KLVCC2),
     &         WORK(KLVCC3),VEC1,VEC2,
     &         N_CC_AMP,IREFSPC,IT2REFSPC,LUC,LUHC,LUSC1,LUSC2,
     &         IDOPROJ,IUNIOP,1,0,1,I_DO_EI,NSPA,0,0,0)
*. ELiminate part referring to unit operator
          CALL TRUNC_MAT(WORK(KLSMAT),NSPA,NSPA,NSPAM1,NSPAM1)
          CALL GET_ON_BASIS2(WORK(KLSMAT),NSPAM1,NSING,
     &              WORK(KLXMAT),WORK(KLVCC1),WORK(KLVCC2),THRES_SINGU)
          WRITE(6,*) ' Number of singularities in S ', NSING
          NNONSING = NSPAM1 - NSING
*. Write to LU_SJ
          CALL REWINO(LU_SJ)
          WRITE(LU_SJ) NSING,NNONSING
          WRITE(LU_SJ) (WORK(KLXMAT-1+IJ),IJ=1,NSPAM1*NNONSING)
        ELSE
*. Read in transformation  matrix from LU_SJ
          CALL REWINO(LU_SJ)
          READ(LU_SJ) NSING,NNONSING
          READ(LU_SJ) (WORK(KLXMAT-1+IJ),IJ=1,NSPAM1*NNONSING)
        END IF
*       ^ End of switch whether complete metrix should read or calc
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Transformation matrix to nonsingular basis '
          CALL WRTMAT(WORK(KLXMAT),NSPAM1,NNONSING,NSPAM1,
     &               NNONSING)
        END IF
*
        IF(INI_IT.EQ.1.AND.IREADSJ.EQ.0) THEN
*. Construct exact or approximate Jacobian
          IF(NCOMMU_J.EQ.1) THEN
*. I assume that the  space before ITREFSPC contains T*IREFSPC 
           ITREFSPC_L = ITREFSPC - 1
           WRITE(6,*) ' Space used for approximate J ', ITREFSPC_L
*. Jacobian independent of T, so use T = 0 for simplicity
           ZERO = 0.0D0
           CALL SETVEC(WORK(KLVCC6),ZERO,N_CC_AMP)
           CALL COM_JMRCC(WORK(KLVCC6),NCOMMU_J,I_APPROX_HCOM_J,
     &          WORK(KLJMAT),WORK(KLVCC1),WORK(KLVCC2), WORK(KLVCC3),
     &          WORK(KLVCC4),N_CC_AMP,NSPAM1,NNONSING,IREFSPC,
     &          ITREFSPC_L,WORK(KLXMAT) )
          ELSE 
*. More than one commutator, so J depends on T
           CALL COM_JMRCC(T_EXT,NCOMMU_J,I_APPROX_HCOM_J,
     &          WORK(KLJMAT),WORK(KLVCC1),WORK(KLVCC2), WORK(KLVCC3),
     &          WORK(KLVCC4),N_CC_AMP,NSPAM1,NNONSING,IREFSPC,
     &          ITREFSPC,WORK(KLXMAT) )
          END IF
*         ^ End if more than one commutator
          WRITE(LU_SJ) (WORK(KLJMAT-1+IJ),IJ=1,NNONSING*NNONSING)
*. Rewind to flush buffer
          CALL REWINO(LU_SJ)
        ELSE
*. Read Approximate Jacobian in from LU_SJ
          READ(LU_SJ) (WORK(KLJMAT-1+IJ),IJ=1,NNONSING*NNONSING)
        END IF
*       ^ End if matrix should be constructed or read in
        I_ADD_SHIFT = 0
        IF(I_ADD_SHIFT.EQ.1) THEN
*. Add a shift to the diagonal of J
          SHIFT = 10.0D0
          WRITE(6,*) ' A shift will be added to initial Jacobian'
          WRITE(6,'(A,E14.7)') ' Value of shift = ', SHIFT
          CALL ADDDIA(WORK(KLJMAT),SHIFT,NNONSING,0)
        END IF
*       ^ End if shift should be added
*
        I_DIAG_J = 0
        IF(I_DIAG_J.EQ.1) THEN
*. Obtain eigenvalues of approximate Jacobian
*. S-matrix is not used anymore to use this space for 
*. diagonalization 
         WRITE(6,*) ' Approximate Jacobian will be diagonalized '
         CALL COPVEC(WORK(KLJMAT),WORK(KLSMAT),NNONSING*NNONSING)
         CALL EIGGMT3(WORK(KLSMAT),NNONSING,WORK(KLVCC1),WORK(KLVCC2),
     &                XDUM,XDUM,XDUM,WORK(KLVCC3),WORK(KLVCC6),1,0)
         WRITE(6,*) ' Real and imaginary part of eigenvalues of J '
         WRITE(6,*) ' ========================================== '
         CALL WRT_2VEC(WORK(KLVCC1),WORK(KLVCC2),NNONSING)
        END IF
*. Obtain LU-Decomposition of Jacobian 
        CALL LULU(WORK(KLJMAT),WORK(KLL),WORK(KLU),NNONSING)
      ELSE
        IF(INI_IT.EQ.1) THEN
*. Complete matrix is not constructed, rather just a diagonal
*. Obtain diagonal of H 
C         GET_DIAG_H0_EI(DIAG,I_IN_TP)
          CALL GET_DIAG_H0_EI(WORK(KLDIA))
*. The last element in KLDIA is the zero-order energy(without core)
          E0 = WORK(KLDIA-1+N_ZERO_EI)
          IF(NTEST.GE.0)
     &    WRITE(6,*) ' Zero-order energy without core term ', E0
*. To get diagonal approximation to J, subtract E0
          DO I = 1, N_ZERO_EI
           WORK(KLDIA-1+I) = WORK(KLDIA-1+I) - E0
          END DO
*. The last term in KLDIA corresponds to the zero-order state.
*. This will not contribute, but to eliminate errors occuring 
*. from dividing by zero do
*. Checl for diagonal values close to zero, and shift these
C         MODDIAG(H0DIAG,NDIM,XMIN)
          WORK(KLDIA-1+N_ZERO_EI) = 300656.0
          XMIN = 0.2D0
          CALL MODDIAG(WORK(KLDIA),N_ZERO_EI,XMIN)
*. And save on LU_SJ
          CALL VEC_TO_DISC(WORK(KLDIA),N_ZERO_EI-1,1,LBLK,LU_SJ)
*. test norm of the E-blocks of diagonal
          IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Norm of various E-blocks of diagonal'
          CALL NORM_T_EI(WORK(KLDIA),2,1,XNORM_EI,1)
          END IF
C NORM_T_EI(T,IEO,ITSYM,XNORM_EI,IPRT)
          IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Diagonal J-approx in ort. zero-order basis'
           CALL WRTMAT(WORK(KLDIA),1,N_ZERO_EI,1,N_ZERO_EI)
          END IF
        END IF
*.      ^ End if it was first iteration
      END IF
*     ^ End of complete or diagonal matrix should be set up
*
* ---------------------
*. 1b : internal part  - constructed in all its.. no problem
* ---------------------
*
      CALL REWINO(LUDIA)
      CALL FRMDSCN(WORK(KLC_INT_DIA),-1,-1,LUDIA)
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Diagonal preconditioner for internal correlation'
         CALL WRTMAT(WORK(KLC_INT_DIA),1,N_REF,1,N_REF)
      END IF
*
      IF(IDIIS.EQ.1.OR.(IDIIS.EQ.2.AND.INI_IT.EQ.1)) THEN
        CALL REWINO(LUSC37)
        CALL REWINO(LUSC36)
      END IF
*. Ensure proper defs
      I12 = 2
      ICSM = IREFSM
      ISSM = IREFSM
      IF(NTEST.GE.100)
     &  WRITE(6,*) ' After const of precond: ITREFSPC, IT2REFSPC =',
     &  ITREFSPC, IT2REFSPC
*
C?    WRITE(6,*) ' KINT before entering optimization'
C?    CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,1)
*. Loop over iterations 
      WRITE(6,*)
      WRITE(6,*) ' -------------------------- '
      WRITE(6,*) ' Entering optimization part ' 
      WRITE(6,*) ' -------------------------- '
      WRITE(6,*)
*. Number of vectors in initial space for DIIS/CROP optimization
      IF(INI_IT.EQ.1) THEN
        NVECP = 0
        NVEC  = 0
      END IF
*. (If INI_IT .ne. 0, MXVC_I vectors from previous macro are used)
      IF(I_DO_SBSPJA.EQ.1) THEN
*. Initialize files that will be used for subspace Jacobian)
        WRITE(6,*) ' LU_CCVECT,LU_CCVECF, LU_CCVECFL = ',
     &               LU_CCVECT,LU_CCVECF, LU_CCVECFL
        CALL REWINO(LU_CCVECT)
        CALL REWINO(LU_CCVECF)
        CALL REWINO(LU_CCVECFL)
      END IF
      DO IT = 1, NNEW_MAX
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 
          WRITE(6,*) ' Information for iteration ', IT
          WRITE(6,*) 
        END IF
        IF(IT.EQ.1) THEN
          MXVC_SUB = MXVC_I
        ELSE
          MXVC_SUB = MXCIVG
        END IF
*
*
* ==================================================================
*. Construct vectorfunction/gradient for external and internal parts
* ==================================================================
*
*. CC vector function for external part  in VCC5 
C?      WRITE(6,*) ' NCAAB before MRCC.. ', NCAAB
        CALL MRCC_VECFNCN(WORK(KLVCC5),T_EXT,
     &       IREFSPC,ITREFSPC,IT2REFSPC,WORK(KLVCC5+N_CC_AMP),
     &       C_REF, N_REF,I_DO_PROJ_NR, 
     &       E_INT,E_EXT,ECORE,1,1)
*
C?      WRITE(6,*) ' Jeppe has asked med to analyze gradient '
C?      CALL ANA_GENCC(WORK(KLVCC5),1)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) 
     &    ' The CC vector function  including internal part'
          CALL WRTMAT(WORK(KLVCC5),1,N_CC_AMP+N_REF,1,N_CC_AMP+N_REF)
        END IF
        IF(NTEST.GE.10) WRITE(6,'(A,I4,1E22.15)')
     &  ' It, Energy from external and internal ', IT, E_EXT + ECORE,
     &        E_INT+ECORE
        VCFNORM_EXT =SQRT(INPROD(WORK(KLVCC5),WORK(KLVCC5),NCAAB))
        VCFNORM_INT = SQRT(
     &  INPROD(WORK(KLVCC5+N_CC_AMP),WORK(KLVCC5+N_CC_AMP),
     &                N_REF)) 
*. Update energy and residual norms
        VNFINAL = VCFNORM_EXT+VCFNORM_INT
        E = E_INT 
        EFINAL = E_INT + ECORE
*. Converged?
        IF(VCFNORM_EXT+VCFNORM_INT.LE.VTHRES) THEN
*. Local iterative procedure converged
          CONVERL = .TRUE.
*. Is global procedure also converged?
          IF((I_REDO_INT.NE.1            ) .OR.
     &       (I_REDO_INT.EQ.1.AND.IT.EQ.1)) THEN
             CONVERG = .TRUE.
          END IF
          WRITE(6,*) ' Iterative procedure converged'
          WRITE(6,'(A,I4,E22.15,2E12.5)')
     &  ' It, energy ,  vecfnc_ext, vecfnc_int ', 
     &    IT, E + ECORE, VCFNORM_EXT, VCFNORM_INT
          GOTO 1001
        END IF
*       ^ End if local procedure is converged
*
* ======================================================================
*. Save vectorfunction in form that will be used in later subspace opt.
* ======================================================================
*
*
        IF(I_DO_SBSPJA.EQ.1) THEN
*. Save Vectorfunction and change in vectorfunction 
*. in EO form if subspace Jacobian is in use
*. Vecfunc in CAAB in VCC5 to Vecfunc in EI in VCC2
*. zero-order state is not to be included
          N_ZERO_EIM = N_ZERO_EI - 1
          CALL TRANS_CAAB_ORTN(WORK(KLVCC5),WORK(KLVCC2),1,1,2,
     &         WORK(KLVCC7),1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Vector function in EI basis '
            CALL WRTMAT(WORK(KLVCC2),1,N_ZERO_EIM,1,N_ZERO_EIM)
          END IF
          IF(IT.GE.2)  THEN
*. Read previous vectorfunction in VCC7 from CCVECFL
            CALL VEC_FROM_DISC(WORK(KLVCC7),N_ZERO_EIM,1,LBLK,
     &           LU_CCVECFL)
            ONE = 1.0D0
            ONEM =-1.0D0
*. Store in VCC7: Delta V  = Vecfnc(ITER) - Vecfnc(ITER-1)
            CALL VECSUM(WORK(KLVCC7),WORK(KLVCC7),WORK(KLVCC2),
     &                  ONEM,ONE,N_ZERO_EIM)
*. Add CCVF(X_{i+1})-CCVF(X_{i}) as vector IT-1 in FILE LU_CCVECF
            CALL SKPVCD(LU_CCVECF,IT-2,WORK(KLVCC6),1,LBLK)
            CALL VEC_TO_DISC(WORK(KLVCC7),N_ZERO_EIM,0,LBLK,LU_CCVECF)
          END IF
*. Save current vector-function in EO form in LU_CCVECFL
          CALL VEC_TO_DISC(WORK(KLVCC2),N_ZERO_EIM,1,LBLK,LU_CCVECFL)
        END IF
*       ^ End if subspace method in use
*
* ========================================================
* Diis/CROP/SBSPJA based on current and previous vectors 
* ========================================================
*
*. Vectors are stored in CAAB basis - not the smartest- Oh yes it was-
*. helps a lot that a common simple basis is used and not  a
*. specific nonsingular basis!
*
        IF(IDIIS.EQ.1.OR.IDIIS.EQ.2) THEN
*. It is assumed that DIIS left the file at end of file 
*. T_ext,C_int on LUSC37, VECFNC on LUSC36
          CALL COPVEC(T_EXT,WORK(KLVCC1),NCAAB)
          CALL COPVEC(C_REF,WORK(KLVCC1+NCAAB),N_REF)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Combined T_ext, C_int coefficients '
            CALL WRTMAT(WORK(KLVCC1),1,NVAR,1,NVAR)
          END IF
          CALL VEC_TO_DISC(WORK(KLVCC1),NVAR,0,-1,LUSC37)
          CALL VEC_TO_DISC(WORK(KLVCC5),NVAR,0,-1,LUSC36)
        END IF
*. We have now a number of vectors in LUSC36, find combination with lowest 
*. norm 
*. DIIS:
        IF(IDIIS.EQ.1) THEN
*. Simple DIIS with no restart
          CALL DIIS_SIMPLE(LUSC36,IT,NVAR,WORK(KLCDIIS))
*. Obtain combination of parameters given in CDIIS
          CALL MVCSMD(LUSC37,WORK(KLCDIIS),LUSC39,LUSC38,
     &                WORK(KLVCC1),WORK(KLVCC2),IT,1,-1)
          CALL VEC_FROM_DISC(WORK(KLVCC1),NVAR,1,-1,LUSC39)
          CALL COPVEC(WORK(KLVCC1),T_EXT,NCAAB)
          CALL COPVEC(WORK(KLVCC1+NCAAB),C_REF,N_REF)
*. Calculate new vectorfunction in VCC5 for T_EXT  and C_INT using sums 
          CALL MVCSMD(LUSC36,WORK(KLCDIIS),LUSC39,LUSC38,
     &                WORK(KLVCC1),WORK(KLVCC2),IT,1,-1)
          CALL VEC_FROM_DISC(WORK(KLVCC5),NVAR,1,-1,LUSC39)
        ELSE IF(IDIIS.EQ.2) THEN
*. CROP:
*. The CROP version of DIIS
*. Matrices are reconstructed in each IT
          IDIRDEL = 1
          NVEC = NVEC + 1
C     CROP(NVEC,NVECP,MXNVEC,NDIM,LUE,LUP,A,
C    &                EOUT,POUT,SCR,LUSCR,IDIRDEL)
*. Note: NVECP is number of vectors for which subspace matrix 
*. has been constructed and saved- CROP updates this
          CALL CROP(NVEC,NVECP,MXVC_SUB,NVAR,LUSC36,LUSC37,
     &         WORK(KLA_CROP),
     &         WORK(KLVCC5),WORK(KLVCC1),WORK(KLSCR_CROP),LUSC39,
     &         IDIRDEL)
*Change of T-coefs
          ONE = 1.0D0
          ONEM = -1.0D0
          CALL VECSUM(WORK(KLVCC1),WORK(KLVCC1),T_EXT,ONE,ONEM,NCAAB)
*. Check if change is to large..
          XNORM = SQRT(INPROD(WORK(KLVCC1),WORK(KLVCC1),NCAAB))
          WRITE(6,*) ' Norm of CROP-correction ', XNORM
          XNORM_MAX = 0.5D0
          I_DO_SCALE = 1
          IF(XNORM.GT.XNORM_MAX.AND.I_DO_SCALE.EQ.1) THEN
            WRITE(6,*) 
     &      ' CROPStep is scaled: from and to to ', XNORM,XNORM_MAX
            FACTOR = XNORM_MAX/XNORM
            CALL SCALVE(WORK(KLVCC1),FACTOR,NCAAB)
            CALL VECSUM(T_EXT,T_EXT,WORK(KLVCC1),ONE,ONE,NCAAB)
          END IF
C         CALL COPVEC(WORK(KLVCC1+NCAAB),C_REF,N_REF)
*.        NOTE: If CI-coefs are changed, they should be renormalized!!
        END IF
*.      ^ End of DIIS/CROP should be used 
        VCFNORM = SQRT(INPROD(WORK(KLVCC5),WORK(KLVCC5),NVAR))
        IF(NTEST.GE.10) WRITE(6,'(A,I4,1E12.5)')
     &  ' From DIIS/CROP : It, norm of approx vecfnc  ',
     &  IT,  VCFNORM 
*
* ===================================================================
* Obtain new direction by applying preconditioners to approx vecfunc
* ===================================================================
*
* --------------
* External part
* --------------
*
*. EI- Approach: Transform Vecfunc to Orthonormal basis, 
*  multiply with diagonal transform result back to CAAB basis
*. Vectorfunction
*. Vecfunc in CAAB in VCC5 to Vecfunc in EI in VCC2
        CALL COPVEC(WORK(KLVCC5),WORK(KLVCC6),NDIM_EI)
        CALL TRANS_CAAB_ORTN(WORK(KLVCC6),WORK(KLVCC2),1,1,2,
     &       WORK(KLVCC7),1)
C            TRANS_CAAB_ORTN(T_CAAB,T_ORTN,ITSYM,ICO,ILR,SCR,
C    &       ICOCON)
          WRITE(6,*) ' Norm of various E-blocks of Vecfnc'
          CALL NORM_T_EI(WORK(KLVCC2),2,1,XNORM_EI,1)
C NORM_T_EI(T,IEO,ITSYM,XNORM_EI,IPRT)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Vectorfunction i ort zero-order basis'
          CALL WRTMAT(WORK(KLVCC2),1,N_ZERO_EI,1,N_ZERO_EI)
        END IF
*
        IF(I_DO_SBSPJA.EQ.0) THEN
*·  New direction = -diag-1 * Vecfunc
          DO I = 1, N_ZERO_EI
            WORK(KLVCC2-1+I) = - WORK(KLVCC2-1+I)/WORK(KLDIA-1+I)
          END DO
*. And no correction for the zero-order state
          WORK(KLVCC2-1+IUNI_AD) = 0.0D0
          WRITE(6,*) ' Norm of various E-blocks of step'
          CALL NORM_T_EI(WORK(KLVCC2),2,1,XNORM_EI,1)
        ELSE
*. Use subspace Jacobian to solve equations
*. Multiply current CC vector function with approximate Jacobian
*. to obtain new step
          NSBSPC_VEC = IT-1
          MAXVEC = MXVEC_SBSPJA
          CALL APRJAC_TV(NSBSPC_VEC,LU_CCVECFL,LUSC41,LU_CCVECT,
     &                   LU_CCVECF,LU_SJ,WORK(KLVCC6),WORK(KLVCC7),
     &                   SCR_SBSPJA,N_ZERO_EIM,LUSC43,LUSC44,
     &                   MAXVEC)
C              APRJAC_TV(NVEC,LUIN,LUOUT,LUVEC,LUJVEC,
C    &                   LUJDIA,VEC1,VEC2,SCR,N_CC_AMP,LUSCR,LUSCR2,
C    &                   MAXVEC)
*. The new correction vector is now residing in LUSC41,
*. Fetch and multiply with -1
          CALL VEC_FROM_DISC(WORK(KLVCC2),N_ZERO_EIM,1,LBLK,LUSC41)
          ONEM = -1.D0
          CALL SCALVE(WORK(KLVCC2),ONEM,N_ZERO_EIM)
*. And no correction for the zero-order state
          WORK(KLVCC2-1+IUNI_AD) = 0.0D0
*. Add step to LU_CCVECT for future use
          CALL SKPVCD(LU_CCVECT,IT-1,WORK(KLVCC6),1,LBLK)
          CALL VEC_TO_DISC(WORK(KLVCC2),N_ZERO_EIM,0,LBLK,LU_CCVECT)
        END IF
*.      ^ End if subspace Jacobian used for generating new step
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' direction in ort zero-order basis'
          CALL WRTMAT(WORK(KLVCC2),1,N_ZERO_EI,1,N_ZERO_EI)
        END IF
*. Dir in EI in VCC2 to Dir in CAAB in VCC1
        CALL TRANS_CAAB_ORTN(WORK(KLVCC1),WORK(KLVCC2),1,2,2,
     &         WORK(KLVCC6),2)
        IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Direction in EI approach, CAAB basis'
          CALL WRTMAT(WORK(KLVCC1),1,NDIM_EI,1,NDIM_EI)
        END IF
*. Norm of change
        XNORM_CAAB = SQRT(INPROD(WORK(KLVCC1),WORK(KLVCC1),N_CC_AMP))
        IF(NTEST.GE.10) WRITE(6,*) ' Norm of correction ', XNORM_CAAB
        XNORM_MAX = 0.5D0
        I_DO_SCALE = 1
        IF(XNORM_CAAB.GT.XNORM_MAX.AND.I_DO_SCALE.EQ.1) THEN
          WRITE(6,*) 
     &    ' Step is scaled: from and to to ', XNORM_CAAB,XNORM_MAX
          FACTOR = XNORM_MAX/XNORM_CAAB
          CALL SCALVE(WORK(KLVCC1),FACTOR,N_CC_AMP)
          XNORM_CAAB = XNORM_MAX
          IF(I_DO_SBSPJA.EQ.1) THEN
*. Well, step was scaled, read in EI form of step and scale this
            CALL SKPVCD(LU_CCVECT,IT-2,WORK(KLVCC2),1,LBLK)
            CALL VEC_FROM_DISC(WORK(KLVCC2),N_ZERO_EIM,0,LBLK,LU_CCVECT)
            CALL SCALVE(WORK(KLVCC2),FACTOR,N_ZERO_EIM)
            CALL SKPVCD(LU_CCVECT,IT-2,WORK(KLVCC2),1,LBLK)
            CALL VEC_TO_DISC(WORK(KLVCC2),N_ZERO_EIM,0,LBLK,LU_CCVECT)
          END IF
        END IF
*. And update the T-coefficients
        ONE = 1.0D0
        CALL VECSUM(T_EXT,T_EXT,WORK(KLVCC1),ONE,ONE,N_CC_AMP)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated T-coefficients in CAAB basis '
          CALL WRTMAT(T_EXT,1,N_CC_AMP,1,N_CC_AMP)
        END IF
*
* --------------
* Internal part
* --------------
*
        IF(N_REF.EQ.1) THEN
          C_REF(1) = 1
          XNORM_CI = 0.0D0
        ELSE
          DO I = 1, N_REF
           XNORM_CI = 0.0D0
           IF(ABS(WORK(KLC_INT_DIA-1+I)-E).GE.1.0D-10) THEN
             DELTA = - WORK(KLVCC5+NCAAB-1+I)/(WORK(KLC_INT_DIA-1+I)-E)
             XNORM_CI = XNORM_CI + DELTA**2
             C_REF(I) = C_REF(I)  + DELTA
           END IF
          END DO
        END IF
        XNORM_CI = SQRT(XNORM_CI)
        WRITE(6,'(A)')
     &  ' It, Energy,  vecfn_ext, vecfn_int, step_ext, step_int: ' 
        WRITE(6,'(I4,1X,E22.15,2x,4(2X,E12.5))')
     &    IT, E + ECORE, VCFNORM_EXT, VCFNORM_INT, XNORM_CAAB, XNORM_CI
*. And normalize the internal part
        CNORM2 = INPROD(C_REF,C_REF,N_REF)
        FACTOR = 1.0D0/SQRT(CNORM2)
        CALL SCALVE(C_REF,FACTOR,N_REF)
*. Write new C_ref to file LUC - used by vector function 
        CALL ISTVC2(WORK(KLVCC2),0,1,N_REF)
        CALL REWINO(LUC)
        CALL WRSVCD(LUC,-1,VEC1,WORK(KLVCC2),
     &              C_REF,N_REF,N_REF,LUDIA,1)
*
      END DO
*     ^ End of loop over iterations
 1001 CONTINUE
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Info from T optimization ', IREFROOT
        WRITE(6,*) ' Updated amplitudes '
        CALL WRTMAT(T_EXT,1,NCAAB,1,NCAAB)
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Analysis of external amplitudes'
        CALL ANA_GENCC(T_EXT,1)
      END IF
*
      IF(IFIN_IT.EQ.1.OR.CONVERG) 
     &CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICC_CMP')
      RETURN
      END 
      SUBROUTINE TCOM_H_N(T,LUINI,LUUT,NCOMMU,IREFSPC,ITREFSPC,
     &           IT2REFSPC,IAC)
*
* Obtain 1/NCOMMU! * NCOMMU-fold commutator of T with H
*
*. Input in CAAB basis
*  Output on LUOT in SD basis 
*. LUUT should differ from scratch files used below, one possible choice is LUHC
*. Scratch files in use : LUSC1, LUSC2, LUSC3, LUSC34
*. Jeppe Olsen, August 2005, Drinking coffee in the early morning at Red Roof Inn in Washington with Jette
*
* IAC = 1 : Add results to LUUT
* IAC = 2 : copy result to LUUT
*
   
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
*
*. Specific input
      DIMENSION T(*)
*. Calculated as sum_I (-1)^(NCOMMU-I) 1/(I!(NCOMMU-1)!) T^(N-I) H T^I |0>
*. So realize the calculation as a loop over I
*
      NTEST = 000
      IF(NTEST.GE.10) THEN 
        WRITE(6,*) ' Task : 1/NCOMMU! times [H,T],T], ... ]]] |0> '
        WRITE(6,*) ' Ncommu = ', NCOMMU
        WRITE(6,*) ' Input T-coefficients '
        CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
        WRITE(6,'(A,3I3)') ' TCOM.., IREFSPC, IT2REFSPC, IAC = ',  
     &                               IREFSPC, IT2REFSPC, IAC
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TCOMHN')
*
* LUINI : Initial expansion |0>
* LUSC1 : T^I |0>
* LUSC2  : H T^I |0>
* LUSC3 : T^N-I H T^I |0>
*
      ONE = 1.0D0
*
      DO I = 0, NCOMMU
        ICSPC = ITREFSPC
        ISSPC = ITREFSPC
C?      WRITE(6,*) ' I = ', I
        IF(I.EQ.0) THEN
*. Expand |0> in IREFSPC on LUINI to ITREFSPC on LUSC1
           CALL EXPCIV(IREFSM,IREFSPC,LUINI,ITREFSPC,LUSC1,-1,
     /                   LUSC34,1,0,IDC,NTEST)
C?         WRITE(6,*) ' After EXPCIV'
        ELSE 
*T^(I-1)|0> => T^I |0> on LUSC1
         CALL REWINO(LUSC1)
         CALL REWINO(LUSC2)
         CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUSC34,T,1)
C?       WRITE(6,*) ' After SIGDEN_CC'
         CALL COPVCD(LUSC34,LUSC1,WORK(KVEC1P),1,-1)
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' T^I |0> for I = ',I
          CALL WRTVCD(WORK(KVEC1P),LUSC1,1,-1)
        END IF
*. Calculate H T^I |0> and save on LUSC2
*. Space of H T^I |0> may be reduced to IT2REFSPC
        ICSPC = ITREFSPC
        ISSPC = IT2REFSPC
C?      WRITE(6,*) ' Before MV7 '
        CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUSC2,0,0)
C?      WRITE(6,*) ' After MV7 '
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' H T^I |0> for I = ',I
           CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
        END IF
*. C space may now also be restricted to IT2REFSPC
        ISSPC = IT2REFSPC
        ICSPC = IT2REFSPC
*. Calculate  T^(NOMMU-I)H T^I on LUSC3
        CALL COPVCD(LUSC2,LUSC3,WORK(KVEC1P),1,-1)
        DO J = 1, NCOMMU-I
C?        WRITE(6,*) ' J = ', J
*. Calculate T * T^(J-1) H T^I |0> and save on LUSC3
          REWIND(LUSC3)
          REWIND(LUSC34)
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC34,T,1)
C?        WRITE(6,*) 'After SIGDEN_CC, 2 '
          CALL COPVCD(LUSC34,LUSC3,WORK(KVEC1P),1,-1)
C?        WRITE(6,*) ' After second COPVCD '
          IF(NTEST.GE.1000) THEN 
            WRITE(6,*) '  T^(J) H T^I for J and I ', J,I
            CALL WRTVCD(WORK(KVEC1P),LUSC3,1,-1)
          END IF
        END DO
C?      WRITE(6,*) ' After J loop '
*. Add (-1)**(NCOMMU-I)1/(NCOMMU-I)!/I! T^(NCOMMU-I) H T^I |0>
        IF(NCOMMU-I.EQ.0) THEN 
          XNMIFAC = 1.0D0
        ELSE 
          XNMIFAC = XFAC(NCOMMU-I)
        END IF
        IF(I.EQ.0) THEN
          XIFAC = 1.0D0
        ELSE 
          XIFAC = XFAC(I)
        END IF
        IF(MOD(NCOMMU-I,2).EQ.0) THEN
         FACTOR = 1.0D0/(XNMIFAC*XIFAC)
        ELSE 
         FACTOR = -1.0D0/(XNMIFAC*XIFAC)
        END IF
*. First contribution : Add or copy
        IF(I.EQ.0) THEN
          IF(IAC.EQ.2) THEN
C                SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
              CALL SCLVCD(LUSC3,LUUT,FACTOR,WORK(KVEC1P),1,-1)
          ELSE 
C?          WRITE(6,*) ' Before VECSMD'
            CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),FACTOR,ONE,LUSC3,
     &           LUUT,LUSC34,1,-1)
C?          WRITE(6,*) ' After VECSMD'
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
              CALL COPVCD(LUSC34,LUUT,WORK(KVEC1P),1,-1)
          END IF
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Initial vector scaled to LUUT '
            CALL WRTVCD(WORK(KVEC1P),LUUT,1,-1)
          END IF
        ELSE 
* add : LUUT = LUUT + FACTOR*LUSC3
          CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),FACTOR,ONE,LUSC3,LUUT,
     &                LUSC34,1,-1)
          CALL COPVCD(LUSC34,LUUT,WORK(KVEC1P),1,-1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' LUUT opdated for I, NCOMMU-I = ', I,NCOMMU-I
            CALL WRTVCD(WORK(KVEC1P),LUUT,1,-1)
          END IF
        END IF
      END DO
*     ^ End of loop over I
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' 1/NCOMMU! [[[H,T,],T..]] |0> (n-fold commutator)'
        CALL WRTVCD(WORK(KVEC1P),LUUT,1,-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TCOMHN')
      RETURN
      END 
      SUBROUTINE GET_GENOP_INFO(NOBEX_TP,IOBEX_TP,NOCCLS,
     &           IOBEX_TP_TO_OCCLS,
     &           KLCOBEX_TP,KLAOBEX_TP,NSPOBEX_TP,
     &           MXSPOXL,KLSOBEX,KLSOX_TO_OX,KIBSOX_FOR_OX,KNSOX_FOR_OX,
     &           KISOX_FOR_OX,KLLSOBEX,KLIBSOBEX,KLSPOBEX_AC,
     &           KIBSOX_FOR_OCCLS,KNSOX_FOR_OCCLS,KISOX_FOR_OCCLS,
     &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &           LEN_T_VEC,MSCOMB_CC,MX_TBLK_AS,
     &           NAOBEX_TP,NBOBEX_TP,KLAOBEX,KLBOBEX,
     &           MAXLENA,MAXLENB,MAXLEN_I1)
*
*. Generate information for general operators as defined by the 
*  NOBEX_TP excitationtypes in IOBEX_TP
*
* Jeppe Olsen, September 05
*
* For working with more than one set of general operators
* 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
C     INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
*. Input
       INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
       INTEGER IOBEX_TP_TO_OCCLS(NOBEX_TP)
*
      NTEST = 10
      IF(NTEST.GE.5) THEN
         WRITE(6,*)
         WRITE(6,*) ' Generation of general operator information '
         WRITE(6,*) ' ========================================== '
         WRITE(6,*)
         WRITE(6,*) ' Orbital excitations : '
C             WRT_ORBEX_LIST(IOBOX,NOBEX,NGAS)
         CALL WRT_ORBEX_LIST(IOBEX_TP,NOBEX_TP,NGAS)
      END IF
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
*. Number of creation and annihilation operators per op
      CALL MEMMAN(KLCOBEX_TP,NOBEX_TP,'ADDL ',1,'COBEX ')
      CALL MEMMAN(KLAOBEX_TP,NOBEX_TP,'ADDL ',1,'AOBEX ')
      CALL GET_NCA_FOR_ORBOP(NOBEX_TP,IOBEX_TP,
     &     WORK(KLCOBEX_TP),WORK(KLAOBEX_TP),NGAS)
*. Number of spinorbital excitations
      IZERO = 0
      MXSPOXL = 0
      IACT_SPC = 0
      IAAEXC_TYP = 3
      IREFSPCX = 0
      MSCOMB_CC = 0
      CALL OBEX_TO_SPOBEX(1,IOBEX_TP,WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TP,NGAS,
     &     NOBPT,0,IZERO,IAAEXC_TYP,IACT_SPC,IPRCC,IDUMMY,
     &     MXSPOXL,IDUMMY,IDUMMY,IDUMMY,NAEL,NBEL,IREFSPCX)
*. And the actual spinorbital excitations
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
*. First SOX of given OX ( including zero operator )
      CALL MEMMAN(KIBSOX_FOR_OX,NOBEX_TP,'ADDL  ',1,'IBSOXF')
*. Number of SOX's for given OX
      CALL MEMMAN(KNSOX_FOR_OX,NOBEX_TP,'ADDL  ',1,'IBSOXF')
*. SOX for given OX
      CALL MEMMAN(KISOX_FOR_OX,NSPOBEX_TP,'ADDL  ',1,'IBSOXF')
*
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,IPRCC,
     &     WORK(KLSOX_TO_OX),MXSPOXL,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPCX)
*
*. Mapping spinorbital excitations => occupation classes
      CALL MEMMAN(KIBSOX_FOR_OCCLS,NOCCLS,'ADDL  ',1,'IBSXOC')
      CALL MEMMAN(KNSOX_FOR_OCCLS,NOCCLS,'ADDL  ',1,' NSXOC')
      CALL MEMMAN(KISOX_FOR_OCCLS,NSPOBEX_TPE,'ADDL  ',1,' ISXOC')
C       SPOBEX_FOR_OCCLS(
C    &           IEXTP_TO_OCCLS,NOCCLS,ISOX_TO_OX,NSOX,
C    &           NSOX_FOR_OCCLS,ISOX_FOR_OCCLS,IBSOX_FOR_OCCLS)
      CALL SPOBEX_FOR_OCCLS(WORK(KEX_TO_OC),NOCCLS,WORK(KLSOX_TO_OX),
     &     NSPOBEX_TPE,WORK(KNSOX_FOR_OCCLS),WORK(KISOX_FOR_OCCLS),
     &     WORK(KIBSOX_FOR_OCCLS))
*
* Dimension and offsets of IC operators
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TP,'ADDL  ',1,'SPOBAC')
*. ALl spinorbital excitations are initially active
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TPE)
*
      ITOP_SM = 1
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ITOP_SM,
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MSCOMB_CC,MX_TBLK_AS,
     &     WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &     NTCONF,IPRCC)
      N_CC_AMP = LEN_T_VEC
      WRITE(6,*) ' Number of IC parameters ', N_CC_AMP
      WRITE(6,*) ' Dimension of the various types '
      CALL IWRTMA(WORK(KLLSOBEX),1,NSPOBEX_TP,1,NSPOBEX_TP)
*
      MX_ST_TSOSO_MX = MX_ST_TSOSO
      MX_ST_TSOSO_BLK_MX = MX_ST_TSOSO_BLK
      MX_TBLK_MX = MX_TBLK
      MX_TBLK_AS_MX = MX_TBLK_AS
      LEN_T_VEC_MX =  LEN_T_VEC
*. Some more scratch etc
*. Alpha- and beta-excitations constituting the spinorbital excitations
*. Number 
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      CALL MEMMAN(KLAOBEX,LENA,'ADDL  ',2,'IAOBEX')
      CALL MEMMAN(KLBOBEX,LENB,'ADDL  ',2,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     0,NAOBEX_TP,NBOBEX_TP,WORK(KLAOBEX),WORK(KLBOBEX))
*. Max dimensions of CCOP !KSTR> = !ISTR> maps
*. For alpha excitations
      IATP = 1
      IOCTPA = IBSPGPFTP(IATP)
      NOCTPA = NSPGPFTP(IATP)
      CALL LEN_GENOP_STR_MAP(
     &     NAOBEX_TP,WORK(KLAOBEX),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IBTP = 2
      IOCTPB = IBSPGPFTP(IBTP)
      NOCTPB = NSPGPFTP(IBTP)
      CALL LEN_GENOP_STR_MAP(
     &     NBOBEX_TP,WORK(KLBOBEX),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(NTEST.GE.5) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
      RETURN
      END
*
      SUBROUTINE WRT_ORBEX_LIST(IOBOX,NOBEX,NGAS)
*
* Print NOBEX orbital excitations given in IOBEX
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOBEX(2*NGAS,NOBEX)
*
      DO JOBEX = 1, NOBEX
        WRITE(6,*) ' Orbital excitation ', JOBEX
        CALL WRT_ORBEX(IOBEX(1,JOBEX),NGAS)
      END DO
*
      RETURN
      END 
      SUBROUTINE WRT_ORBEX(IOBEX,NGAS)
*
* Print orbital excitation 
*
      INCLUDE 'implicit.inc'
      INTEGER IOBEX(NGAS,2)
*
      WRITE(6,'(A,16I3)') ' Crea for each GASpace : ', 
     &                     (IOBEX(I,1),I=1,NGAS)
      WRITE(6,'(A,16I3)') ' Anni for each GASpace : ', 
     &                     (IOBEX(I,2),I=1,NGAS)
*
      RETURN
      END 
      SUBROUTINE GET_ON_BASIS2(S,NVEC,NSING,X,SCRVEC1,SCRVEC2,
     &           THRES_SINGU)
*
* NVEC vectors with overlap matrix S are given.
* Obtain transformation matrix to orthonormal basis
*
* NSING is the number of singularities obtained 
* If there are singularities, the nonsingular transformation 
* os obtained as a NVEC x (NVEC-NSING) matrix in X 
* First vectors. The eigenvectors corresponding to the 
* singular eigenvectors are lost. 
*
*
* Jeppe Olsen, Palermo, oct 2002
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION S(NVEC*NVEC)
*. Output
      DIMENSION X(NVEC*NVEC)
*. Local scratch
      DIMENSION SCRVEC1(*), SCRVEC2(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GET_ON_BASIS speaking '
        WRITE(6,*) ' Input overlap matrix '
        CALL WRTMAT(S,NVEC,NVEC,NVEC,NVEC)
      END IF
*1 : Diagonalize S and save eigenvalues in SCRVEC1
      CALL COPVEC(S,X,NVEC*NVEC)
C          DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
      CALL DIAG_SYMMAT_EISPACK(X,SCRVEC1,SCRVEC2,NVEC,IRETURN)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Eigenvalues of metric '
        CALL WRTMAT(SCRVEC1,1,NVEC,1,NVEC)
      END IF
*2 : Count number of nonsingularities
      NNONSING = 0
      THRES = 1.0D-14
      DO I = 1, NVEC
        IF(ABS(SCRVEC1(I)).GT.THRES) THEN
          NNONSING = NNONSING + 1
          IF(I.NE.NNONSING) THEN
            SCRVEC1(NNONSING) = SCRVEC1(I)
            CALL COPVEC(X((I-1)*NVEC+1), X((NNONSING-1)*NVEC+1),NVEC)
          END IF
        END IF
      END DO
      NSING = NVEC - NNONSING
*2 : Rearrange so the nonsingular
*    eigenvectors and eigenvalues are  the first parts of X and 
*    SCRVEC1
CE    ISING = 0
CE    INONSING = 0
CE    DO I = 1, NVEC
CE      IF(ABS(SCRVEC1(I)) .GT. THRES) THEN
*. A nonsingular eigenpair
CE        INONSING = INONSING + 1
CE        ITO = INONSING
CE      ELSE 
*. A singular eigenpair
CE        ISING = ISING + 1
CE        ITO = ISING + NNONSING
CE      END IF
CE      IF(ITO.NE.I) THEN
CE        SCRVEC1(ITO) = SCRVEC1(I)
CE        CALL COPVEC(X((I-1)*NVEC+1), X((ITO-1)*NVEC+1),NVEC)
CE      END IF
CE    END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Nonsingular eigenvalues of metric '
        CALL WRTMAT(SCRVEC1,1,NNONSING,1,NNONSING)
      END IF
*3 : Construct orthonormal basis using 
*  X = U sigma^{-1/2}, 
*  where U are the nonsingular 
*. eigenvectors of S and sigma are the corresponding eigenvalues
      DO I = 1, NNONSING
        SCALE = 1/SQRT(SCRVEC1(I))
        IBX = (I-1)*NVEC+1
        CALL SCALVE(X(IBX),SCALE,NVEC)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformation matrix to nonsingular basis '
        CALL WRTMAT(X,NVEC,NNONSING,NVEC,NNONSING)
      END IF
*
      RETURN
      END 
C              PRECTV(VEC1,VEC2,E,LUDIAM,LUDIAS)
      SUBROUTINE H0_EI_TV(VECIN,VECOUT,E,LUDIA,LUDIAS,VECSCR)
*
* A vector, VECIN, is given in the zero-order basis. 
* Multiply with inverse diagonal of LUDIA
*
*. Jeppe Olsen, Sicily sept. 2009
*
      INCLUDE 'implicit.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'cshift.inc'
      REAL*8 INPROD
*
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT(*)
*. Scratch
      DIMENSION VECSCR(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'H0EITV')
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Information from H0_EI_TV '
      END IF
*
      VECIN_ORT= INPROD(VECIN,VECIN,N_ZERO_EI-1)
*. read in approximate (and unshifted) Jacobian in VECSCR
      CALL VEC_FROM_DISC(VECSCR,N_ZERO_EI,1,-1,LUDIA)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Diagonal read in '
        CALL WRTMAT(VECSCR,1,N_ZERO_EI,1,N_ZERO_EI)
      END IF
      E0 = VECSCR(N_ZERO_EI)
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' EREFX, E, E0 = ', EREFX,E,E0
      END IF
*·  New direction = - Vecfunc/(diag - e)
      DO I = 1, N_ZERO_EI - 1
       VECOUT(I) = -VECIN(I)/(VECSCR(I) - E0)
      END DO
*. And the final element- corresponding to the zero-order state
      IF(ABS(EREFX-E).GT.1.0D-10) THEN
        VECOUT(N_ZERO_EI) = -VECIN(N_ZERO_EI)/(EREFX-E)
      ELSE
        VECOUT(N_ZERO_EI) = 0.0D0
      END IF
*
      VECOUT_ORT= INPROD(VECOUT,VECOUT,N_ZERO_EI-1)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' VECIN_0, VECIN_ORT = ', 
     &               VECIN(N_ZERO_EI),VECIN_ORT
        WRITE(6,*) ' VECOUT_0, VECOUT_ORT = ', 
     &               VECOUT(N_ZERO_EI),VECOUT_ORT
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' direction in ort zero-order basis'
        CALL WRTMAT(VECOUT,1,N_ZERO_EI,1,N_ZERO_EI)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'H0EITV')
      RETURN
      END
      SUBROUTINE LUCIA_ICCI(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                      EFINAL,CONVER,VNFINAL)
*
* Master routine for Internal Contraction CI
*
* LUCIA_IC is assumed to have been called to do the 
* preperatory work for working with internal contraction
*
* Jeppe Olsen, October 2009 (as separate routine)
*
C     INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
      LOGICAL CONVER,CONVER_INT,CONVER_EXT
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cei.inc'
*. Transfer common block for communicating with H_EFF * vector routines
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
*. Transfer block for communicating zero order energy to 
*. routien for performing H0-E0 * vector
      INCLUDE 'cshift.inc'
*
      CHARACTER*6 ICTYP
      EXTERNAL MTV_FUSK, STV_FUSK
      EXTERNAL H_S_EFF_ICCI_TV,H_S_EXT_ICCI_TV
      EXTERNAL HOME_SD_INV_T_ICCI
      EXTERNAL H0_EI_TV
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICCI  ')
      NTEST = 10
      WRITE(6,*)
      WRITE(6,*) ' ===================='
      WRITE(6,*) ' ICCI section entered '
      WRITE(6,*) ' ===================='
      WRITE(6,*)
*
      IF(IEI_VERSION.EQ.0) THEN
        I_DO_EI = 0
      ELSE
        I_DO_EI = 1
      END IF
*
      IF(I_DO_EI.EQ.1) THEN
       WRITE(6,*) ' EI approach in use'
      ELSE
       WRITE(6,*) ' Partial spin-adaptation in use'
      END IF
*

      WRITE(6,*) ' Energy of reference state ', EREF
*. Number of parameters with and without spinadaptation
      IF(I_DO_EI.EQ.0) THEN
        CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
      ELSE
*. zero-particle operator is included in N_ZERO_EI
        NSPA = N_ZERO_EI
*. Note: NCAAB includes unitop
        NCAAB = NDIM_EI
      END IF
      IF(I_DO_EI.EQ.0) THEN
          WRITE(6,*) ' Number of spin-adapted operators ', NSPA
      ELSE
          WRITE(6,*) ' Number of orthonormal zero-order states',
     &    N_ZERO_EI
      END IF
      WRITE(6,*) ' Number of CAAB operators         ', NCAAB
*. Number of spin adapted operators without the unitoperator
      I_IT_OR_DIR = 1
      IF(I_IT_OR_DIR.EQ.2) THEN
        WRITE(6,*) ' Explicit construction of all matrices'
      ELSE
        WRITE(6,*) ' Iterative solution of equations'
      END IF
*
      I_RELAX_INT = 1
* 
*
      N_REF = XISPSM(IREFSM,IREFSPC)
*. Space for external correlation vector
      CALL MEMMAN(KLTEXT,NCAAB,'ADDL  ',2,'T_EXT ')
*. Initial  guess to T_EXT: just a 1 for the zero order state
      IF(IRESTRT_IC.EQ.0) THEN 
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLTEXT),ZERO,NSPA)
        WORK(KLTEXT-1+NSPA) = 1.0D0
*. Store inital guess on unit 54
        CALL VEC_TO_DISC(WORK(KLTEXT),NSPA,1,-1,LUSC54)
      END IF
*
      CONVER =.FALSE.
      CONVER_INT = .FALSE.
      CONVER_EXT = .FALSE.
      I12 = 2
*
      MAXIT_MACRO = MAXITM
      MAXITL  = MAXIT
      MAXVECL = MXCIV
      WRITE(6,'(A,2I4)') 
     &' Allowed number of outer and inner iterations', 
     &  MAXIT_MACRO, MAXITL
*. Convergence will be defined as energy change
      I_ER_CONV = 1
*. There is no external converence threshold for linear equations,
*. just use sqrt of energythreshold
      THRES_R = SQRT(THRES_E)
      DO IT_IE = 1, MAXIT_MACRO
*
        IF(NTEST.GE.1) THEN
          WRITE(6,*)
          WRITE(6,*) ' ------------------------------------------'
          WRITE(6,*) ' Information from outer iteration ', IT_IE
          WRITE(6,*) ' ------------------------------------------'
          WRITE(6,*)
        END IF
        IDUM = 0
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'COMP_M')
*. Start by obtaining set of internal states
        I_REDO_ZERO = 1
        IF(I_DO_EI.EQ.1.AND.I_REDO_ZERO.EQ.1) THEN
          WRITE(6,*) ' Zero-order states recalculated'
          CALL GET_INTERNAL_STATES_OUTER
          N_INT_MAX = IMNMX(WORK(KL_N_INT_FOR_SE),N_EXTOP_TP*NSMOB,2)
*. Largest number of zero-order states of given sym and external type
          N_ORTN_MAX = IMNMX(WORK(KL_N_ORTN_FOR_SE),N_EXTOP_TP*NSMOB,2)
          WRITE(6,*) ' N_INT_MAX, N_ORTN_MAX = ', N_INT_MAX, N_ORTN_MAX
*. Largest transformation block 
          N_XEO_MAX = N_INT_MAX*N_ORTN_MAX
          IF(NTEST.GE.5) WRITE(6,*) ' Largest (EL,ORTN) block = ', 
     &    N_XEO_MAX
*. Number of zero-order states - does now include the unit-operator
          N_ZERO_EI = N_ZERO_ORDER_STATES(WORK(KL_N_ORTN_FOR_SE),
     &                WORK(KL_NDIM_EX_ST),N_EXTOP_TP,1)
          NSPA = N_ZERO_EI
        END IF
*
* ======================================================
*. Coefficients for external correlation for root NROOT
* ======================================================
        IF(NTEST.GE.0) THEN
           WRITE(6,*)
           WRITE(6,*) ' Optimization of external correlation part'
           WRITE(6,*) ' .........................................'
           WRITE(6,*)
        END IF
*
*. Prepare transfer common block used for H(ICCI) * v, S(ICCI) * v ( also used for constructing H,S)
*. Not used here 
        C_0X = 0.0D0
        KLTOPX = -1
*. Used 
        NREFX = N_REF
        IREFSPCX = IREFSPC
        ITREFSPCX = ITREFSPC
        NCAABX = N_CC_AMP
        NSPAX = NSPA
        IPROJSPCX = IREFSPC
*. Unitoperator in SPA order ... Please check ..
        IUNIOPX = NSPA
        IF (I_IT_OR_DIR.EQ.2 ) THEN
*. Construct matrices explicit and diagonalize
*. Not used here 
          C_0X = 0.0D0
          KLTOPX = -1
*. Used 
          NREFX = N_REF
          IREFSPCX = IREFSPC
          ITREFSPCX = ITREFSPC
          NCAABX = N_CC_AMP
          NSPAX = NSPA
          IPROJSPCX = IREFSPC
          CALL ICCI_COMPLETE_MAT2(IREFSPC,ITREFSPC,I_SPIN_ADAPT,
     &         NROOT,WORK(KLTEXT),C_0,E_EXTOP)

          EFINAL = E_EXTOP
          CONVER_EXT = .TRUE.
          VNFINAL_EXT = 0.0D0
        ELSE 
*. Iterative approach to solving ICCI equations ....
*. Currently : no preconditioning and no elimination of singularities 
*              ( Yes, I am still an optimist ( or desperate ))
          NTESTL = 10
*. Space for CI behind the curtain 
CMOVED    CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
CMOVED    KVEC1P = KVEC1
CMOVED    KVEC2P = KVEC2
*. Allocate space for iterative solver 
          CALL MEMMAN(KL_EXTVEC1,NCAAB,'ADDL  ',2,'EXTVC1')
          CALL MEMMAN(KL_EXTVEC2,NCAAB,'ADDL  ',2,'EXTVC2')
          CALL MEMMAN(KL_EXTVEC3,NCAAB,'ADDL  ',2,'EXTVC3')
*         ^ KLEXTVEC3 is also used as scratch in reformat 
          CALL MEMMAN(KL_EXTVEC4,NCAAB,'ADDL  ',2,'EXTVC3')
*
          CALL MEMMAN(KL_RNRM,MAXITL*NROOT,'ADDL  ',2,'RNRM  ')
          CALL MEMMAN(KL_EIG ,MAXITL*NROOT,'ADDL  ',2,'EIG   ')
          CALL MEMMAN(KL_FINEIG,NROOT,'ADDL  ',2,'FINEIG')
*
          CALL MEMMAN(KL_APROJ,MAXVECL**2,'ADDL  ',2,'APROJ ')
          CALL MEMMAN(KL_SPROJ,MAXVECL**2,'ADDL  ',2,'SPROJ ')
          CALL MEMMAN(KL_AVEC ,MAXVECL**2,'ADDL  ',2,'AVEC  ')
          LLWORK = 5*MAXVECL**2 + 2*MAXVECL
          CALL MEMMAN(KL_WORK ,LLWORK   ,'ADDL  ',2,'WORK  ')
          CALL MEMMAN(KL_AVEC ,MAXVECL**2,'ADDL  ',2,'AVECP ')
          CALL MEMMAN(KL_AVECP,MAXVECL**2,'ADDL  ',2,'AVECP ')
*. Obtain diagonal of H and S
          I_DO_PRE_IN_EXT = 0
          IF(I_DO_PRE_IN_EXT.EQ.1) THEN
           IF(I_DO_EI.EQ.0) THEN
             CALL GET_HS_DIA(WORK(KL_EXTVEC3),WORK(KL_EXTVEC4),
     &            1,1,1,WORK(KL_EXTVEC1),WORK(KL_EXTVEC2),
     &              WORK(KVEC1),WORK(KVEC2),IREFSPC,ITREFSPC,
     &            IUNIOPX,NSPA,0,IDUM,IDUM)
           ELSE
*. EI approach
             CALL GET_DIAG_H0_EI(WORK(KL_EXTVEC3))
*. clean up
             I12 = 2
*. States are normalized, so
             ONE = 1.0D0
             CALL SETVEC(WORK(KL_EXTVEC4),ONE,NSPA)
           END IF
          ELSE
           ONE = 1.0D0
           CALL SETVEC(WORK(KL_EXTVEC3),ONE,NSPA)
           CALL SETVEC(WORK(KL_EXTVEC4),ONE,NSPA)
          END IF
*. And write diagonal to disc as single record files
          CALL VEC_TO_DISC(WORK(KL_EXTVEC3),NSPA,1,-1,LUSC53)
          CALL VEC_TO_DISC(WORK(KL_EXTVEC4),NSPA,1,-1,LUSC51)
*. (LUSC51 is not used)
          IF(IRESTRT_IC.EQ.1) THEN
*. Copy old CI coefficients for reference space to LUC
            CALL COPVCD(LUEXC,LUC,WORK(KVEC1),1,-1)
          END IF
          DO IMAC = 1, 1
* LUSC53 is LU_DIAH, LUSC51 is LU_DIAS
*. 2 implies that advanced preconditioner is called 
*- Save reference energy for use with diagonal preconditioner
            EREFX = EREF
*
            IF(IT_IE.GT.1) THEN
              I_ENFORCE_COLD_START = 0
              IF(I_ENFORCE_COLD_START.EQ.1) THEN
                WRITE(6,*) ' Enforced start with Text = 0'
                ZERO = 0.0D0
                CALL SETVEC(WORK(KLTEXT),ZERO,NSPA)
                WORK(KLTEXT-1+NSPA) = 1.0D0
                CALL VEC_TO_DISC(WORK(KLTEXT),NSPA,1,-1,LUSC54)
              ELSE
*. Use the previous coefficients to start. 
                T_CAAB_NORM =
     &          SQRT(INPROD(WORK(KLTEXT),WORK(KLTEXT),NCAAB))
                WRITE(6,*) ' Norm of T in CAAB basis before MINGENEIG',
     &          T_CAAB_NORM
                WRITE(6,*) ' T(zero-op) in CAAB basis ', 
     &          WORK(KLTEXT-1+NCAAB)
*. Transform to zero-order basis- used in MINGENEIG
                CALL TRANS_CAAB_ORTN(WORK(KLTEXT),WORK(KL_EXTVEC1),
     &                               1,1,2,WORK(KL_EXTVEC3),2) 
*. Test back-transformation to CAAB basis
                CALL TRANS_CAAB_ORTN(WORK(KL_EXTVEC4),WORK(KL_EXTVEC1),
     &                               1,2,2,WORK(KL_EXTVEC3),2) 
                T_CAAB_NORM2 =
     &          SQRT(INPROD(WORK(KL_EXTVEC4),WORK(KL_EXTVEC4),NCAAB))
                WRITE(6,*) ' Norm of T in CAAB basis backtransformed',
     &          T_CAAB_NORM2
*. End of test
                T_ORT_NORM =
     &          SQRT(INPROD(WORK(KL_EXTVEC1),WORK(KL_EXTVEC1),NSPA))
                WRITE(6,*) ' Norm of T in Ort basis before MINGENEIG',
     &          T_ORT_NORM
                WRITE(6,*) ' T(zero-op) in ort basis ', 
     &          WORK(KL_EXTVEC1-1+NSPA)
                CALL VEC_TO_DISC(WORK(KL_EXTVEC1),NSPA,1,-1,LUSC54)
              END IF
            END IF
*           ^ End if not first IE-iteration 
*
            I12 = 2
            IF(I_DO_EI.EQ.0) THEN
              IPREC_FORM = 1
              SHIFT =  0.0D0
              CALL MINGENEIG(H_S_EXT_ICCI_TV,HOME_SD_INV_T_ICCI,
     &             IPREC_FORM,THRES_E,THRES_R,I_ER_CONV,
     &             WORK(KL_EXTVEC1),WORK(KL_EXTVEC2),WORK(KL_EXTVEC3),
     &             LUSC54, LUSC37,
     &             WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &             NSPA,LUSC38,LUSC39,LUSC40,LUSC53,LUSC51,LUSC52,
     &             NROOT,MAXVECL,NROOT,WORK(KL_APROJ),
     &             WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &             NTESTL,SHIFT,WORK(KL_AVECP),I_DO_PRE_IN_EXT,
     &             CONVER_EXT,E_EXTOP,VNFINAL_EXT)
            ELSE
              IPREC_FORM = 2
              SHIFT = 0.0D0
              CALL MINGENEIG(H_S_EXT_ICCI_TV,H0_EI_TV,
     &             IPREC_FORM,THRES_E,THRES_R,I_ER_CONV,
     &             WORK(KL_EXTVEC1),WORK(KL_EXTVEC2),WORK(KL_EXTVEC3),
     &             LUSC54, LUSC37,
     &             WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &             NSPA,LUSC38,LUSC39,LUSC40,LUSC53,LUSC51,LUSC52,
     &             NROOT,MAXVECL,NROOT,WORK(KL_APROJ),
     &             WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &             NTESTL,SHIFT,WORK(KL_AVECP),I_DO_PRE_IN_EXT,
     &             CONVER_EXT,E_EXTOP,VNFINAL_EXT)
            END IF
           EFINAL = E_EXTOP
          END DO
*         ^ End of loop over reset eigenvalue problem
          CALL VEC_FROM_DISC(WORK(KL_EXTVEC1),NSPA,1,-1,LUSC54)
*
          T_ORT_NORM =
     &    SQRT(INPROD(WORK(KL_EXTVEC1),WORK(KL_EXTVEC1),NSPA))
          WRITE(6,*) ' Norm of T in Ort basis after MINGENEIG',
     &    T_ORT_NORM
          C_0 = WORK(KL_EXTVEC1-1+NSPA)
*. And reform to CAAB basis and store in KLTEXT
          IF(I_DO_EI.EQ.0) THEN
            CALL REF_CCV_CAAB_SP(WORK(KLTEXT),WORK(KL_EXTVEC1),
     &                       WORK(KL_EXTVEC3),2) 
          ELSE
            CALL TRANS_CAAB_ORTN(WORK(KLTEXT),WORK(KL_EXTVEC1),1,2,2,
     &                            WORK(KL_EXTVEC3),2) 
          END IF
          T_CAAB_NORM =
     &    SQRT(INPROD(WORK(KLTEXT),WORK(KLTEXT),NCAAB))
          WRITE(6,*) ' Norm of T in CAAB basis after MINGENEIG',
     &    T_CAAB_NORM
*
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' coefficient of zero-order state ', C_0
            WRITE(6,*) ' Analysis of external amplitudes in CAAB basis'
            CALL ANA_GENCC(WORK(KLTEXT),1)
          END IF
    
        END IF
*       ^ End of switch direct/iterative approach for T_EXT
        IF(I_RELAX_INT.EQ.1) THEN
* ============================================================
*. Relax coefficients of internal/reference/zero-order state 
* ============================================================
*
        IF(NTEST.GE.0) THEN
           WRITE(6,*)
           WRITE(6,*) ' Optimization of internal correlation part'
           WRITE(6,*) ' .........................................'
           WRITE(6,*)
        END IF
COLD       CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
COLD       KVEC1P = KVEC1
COLD       KVEC2P = KVEC2
*
           IF(I_IT_OR_DIR.EQ.2) THEN
*
*. Construct complete matrices and diagonalize
*
*. Space for H and S in zero-order space 
             CALL MEMMAN(KLH_REF,N_REF**2,'ADDL  ',2,'H_REF  ')
             CALL MEMMAN(KLS_REF,N_REF**2,'ADDL  ',2,'S_REF  ')
             CALL MEMMAN(KLC_REF,N_REF   ,'ADDL  ',2,'C_REF  ')
             CALL MEMMAN(KLI_REF,N_REF   ,'ADDL  ',1,'I_REF  ')
*
             CALL ICCI_RELAX_REFCOEFS_COM(WORK(KLTEXT),NSPA,
     &            WORK(KLH_REF),
     &            WORK(KLS_REF),N_REF,WORK(KVEC1),WORK(KVEC2),1,
     &            IREFSPC,ITREFSPC,C_0,ECORE,WORK(KLC_REF),NROOT,
     &            NCAAB,E_INTOP)
             CONVER_INT =.TRUE.
             VNFINAL_INT = 0.0D0
             EFINAL = E_INTOP
*. transfer new reference vector to DISC
             CALL ISTVC2(WORK(KLI_REF),0,1,N_REF)
C  WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
             CALL REWINO(LUC)
             CALL WRSVCD(LUC,-1,WORK(KVEC1),WORK(KLI_REF),
     &            WORK(KLC_REF),N_REF,N_REF,LUDIA,1)
           ELSE 
*. Use iterative methods to reoptimize reference coefficients
             MAXITL = MAXIT
             MAXVEC = MXCIV
*
             CALL MEMMAN(KL_REFVEC1,N_REF,'ADDL  ',2,'REFVC1')
             CALL MEMMAN(KL_REFVEC2,N_REF,'ADDL  ',2,'REFVC2')
             CALL MEMMAN(KL_REFVEC3,N_REF,'ADDL  ',2,'REFVC3')
*
             CALL MEMMAN(KL_RNRM,MAXIT*NROOT,'ADDL  ',2,'RNRM  ')
             CALL MEMMAN(KL_EIG ,MAXIT*NROOT,'ADDL  ',2,'EIG   ')
             CALL MEMMAN(KL_FINEIG,NROOT,'ADDL  ',2,'FINEIG')
*
             CALL MEMMAN(KL_APROJ,MAXVEC**2,'ADDL  ',2,'APROJ ')
             CALL MEMMAN(KL_SPROJ,MAXVEC**2,'ADDL  ',2,'SPROJ ')
             CALL MEMMAN(KL_AVEC ,MAXVEC**2,'ADDL  ',2,'AVEC  ')
             LLWORK = 5*MAXVEC**2 + 2*MAXVEC
             CALL MEMMAN(KL_WORK ,LLWORK   ,'ADDL  ',2,'WORK  ')
             CALL MEMMAN(KL_AVEC ,MAXVEC**2,'ADDL  ',2,'AVECP ')
             CALL MEMMAN(KL_AVECP,MAXVEC**2,'ADDL  ',2,'AVECP ')
*
* Well, there is pt a conflict between the form of files 
* in mingeneig and in the general CI programs
*. In MINGENEIG all vectors are single record files, whereas
*  the vectors are multirecord files in the general LUCIA 
* world. Reformatting is therefore required..
*. LUC is LUC
*. LUSC36 is LUDIA
*. LUSC51 is LUDIAS
*
*. Reform LUC to single record file
             CALL REWINO(LUC)
             CALL FRMDSCN(WORK(KL_REFVEC1),-1,-1,LUC)
             CALL REWINO(LUC)
             CALL VEC_TO_DISC(WORK(KL_REFVEC1),N_REF,1,-1,LUC)
*. Reform LUDIA to single record file on LUSC36
             CALL REWINO(LUDIA)
             CALL FRMDSCN(WORK(KL_REFVEC1),-1,-1,LUDIA)
             CALL VEC_TO_DISC(WORK(KL_REFVEC1),N_REF,1,-1,LUSC36)
*. Write diagonal of S as unit mat as single vector file
             ONE = 1.0D0
             CALL SETVEC(WORK(KL_REFVEC1),ONE,N_REF)
             CALL VEC_TO_DISC(WORK(KL_REFVEC1),N_REF,1,-1,LUSC51)
*. (LUSC51 is not used)
*
* As preconditioners, the standard CI diagonal and the 
* unit diagonal will be used for H and S, respectively.
* This is fine if the T operator is not too large...
*
*. Prepare transfer common block for communicating with
*. matrix-vector routines
C            C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX
             C_0X = C_0
             KLTOPX = KLTEXT
             NREFX = N_REF
             IREFSPCX = IREFSPC
             ITREFSPCX = ITREFSPC
             NCAABX = N_CC_AMP
             NSPAX = NSPA
*. Unitoperator in SPA order ... Please check ..
             IUNIOPX = NSPA
*
             SHIFT = 0.0D0
             CALL MINGENEIG( H_S_EFF_ICCI_TV,HOME_SD_INV_T_ICCI,1,
     &            THRES_E,THRES_R,I_ER_CONV,
     &            WORK(KL_REFVEC1),WORK(KL_REFVEC2),WORK(KL_REFVEC3),
     &            LUC, LUSC37,
     &            WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &            N_REF,LUSC38,LUSC39,LUSC40,LUSC36,LUSC51,LUSC52,
     &            NROOT,MXCIV,NROOT,WORK(KL_APROJ),
     &            WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &            NTESTL,SHIFT,WORK(KL_AVECP),1,
     &            CONVER_INT,E_INTOP,VNFINAL_INT)
                  EFINAL = E_INTOP
C                 MINGENEIG(MTV,STV,
C    &                VEC1,VEC2,VEC3,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
C    &                NVAR,
C    &                LU3,LU4,LU5,LUDIAM,LUDIAS,LUS,NROOT,MAXVEC,
C    &                NINVEC,
C    &                APROJ,AVEC,SPROJ,WORK,IPRT,EIGSHF,AVECP,I_DO_PRECOND)
*
*. Read new eigenvector from LUC
             CALL REWINO(LUC)
             CALL FRMDSCN(WORK(KL_REFVEC1),-1,-1,LUC)
* The eigenvector is normalized with respect to the <i!T+P P T|j>
*. metric, normalize with standard unit metrix
             XNORM = INPROD(WORK(KL_REFVEC1),WORK(KL_REFVEC1),N_REF)
             FACTOR = 1.0D0/SQRT(XNORM)
             CALL SCALVE(WORK(KL_REFVEC1),FACTOR,N_REF)
*. And write to disc in a form suitable for the other parts of LUCIA
             CALL ISTVC2(WORK(KL_REFVEC2),0,1,N_REF)
             CALL REWINO(LUC)
             CALL REWINO(LUDIA)
             CALL WRSVCD(LUC,-1,WORK(KVEC1P),WORK(KL_REFVEC2),
     &                   WORK(KL_REFVEC1),N_REF,N_REF,LUDIA,1)
             IF(NTEST.GE.100) THEN
               WRITE(6,*) ' New reference coefficients '
               CALL WRTVCD(WORK(KVEC1P),LUC,1,-1)
             END IF
           END IF 
*.         ^ End of switch direct/iterative methods for reference relaxation 
        END IF
*.      ^ End of reference coefs should be relaxed
        CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'COMP_M')
        IF(CONVER_INT.AND.CONVER_EXT.AND.
     &     ABS(E_INTOP-E_EXTOP).LE.THRES_E) CONVER = .TRUE.
        IF(CONVER) GOTO 1001
      END DO
 1001 CONTINUE
*
      IF(MAXIT_MACRO.GT.0) THEN
       IF(NTEST.GE.10) THEN
        WRITE(6,*) ' coefficient of zero-order state ', C_0
        WRITE(6,*) 
     &  ' Analysis of final external amplitudes in CAAB basis'
        CALL ANA_GENCC(WORK(KLTEXT),1)
       END IF
*
       VNFINAL = VNFINAL_INT + VNFINAL_EXT
       WRITE(6,*) ' VNFINAL_INT, VNFINAL_EXT =', 
     &              VNFINAL_INT,VNFINAL_EXT
*. Print the final coefs ..
C?     CALL VEC_FROM_DISC(WORK(KL_EXTVEC1),NSPA,1,-1,LUSC54)
C?     WRITE(6,*) ' Final list of IC-coefficients '
C?     CALL WRTMAT(WORK(KL_EXTVEC1),NSPA,1,NSPA,1)
      END IF ! There were iterations to analyze
      
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICCI  ')
      RETURN
      END 
      SUBROUTINE GET_INTERNAL_STATES_OUTER
*
* Outer routine for obtaining set of orthonormal internal states
*
*. For hiding an ugly parameter list..
*
*. Jeppe Olsen, Oct. 2009
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
*
      WRITE(6,*) ' GET_INTERNAL..., I_INT_HAM = ', I_INT_HAM
      CALL GET_INTERNAL_STATES(N_EXTOP_TP,N_INTOP_TP,
     &     WORK(KLSOBEX),WORK(KL_N_INT_FOR_EXT),WORK(KL_IB_INT_FOR_EXT),
     &     WORK(KL_I_INT_FOR_EXT),WORK(KL_NDIM_IN_SE),
     &     WORK(KL_N_ORTN_FOR_SE),WORK(KL_N_INT_FOR_SE),
     &     WORK(KL_X1_INT_EI_FOR_SE), WORK(KL_X2_INT_EI_FOR_SE),
     &     WORK(KL_SG_INT_EI_FOR_SE),WORK(KL_S_INT_EI_FOR_SE),
     &     WORK(KL_IBX1_INT_EI_FOR_SE), WORK(KL_IBX2_INT_EI_FOR_SE),
     &     WORK(KL_IBSG_INT_EI_FOR_SE),WORK(KL_IBS_INT_EI_FOR_SE),
     &     WORK(KL_X2L_INT_EI_FOR_SE),
     &     I_IN_TP,I_INT_OFF,I_EXT_OFF) 
*
      RETURN 
      END
      SUBROUTINE MRCC_VECFNCN(CCVECFNC,T,
     &           IREFSPC,ITREFSPC,IT2REFSPC,CCVECFNCI,C_REF,N_REF,
     &           I_DO_PROJ_NR,E_INT,E_EXT,ECORE,I_INI_CO,I_FIN_CO)
*
* Obtain external and internal parts of the MRCC vector function 
*
*. Version allowing various forms of input and output and 
*. includes calculation of internal part for NCOMMU_E .ne N_COMMU_V
*
* I_INI_CO = 1 => Initial guess is in CAAB basis, 
*          = 2 => Initial guess is in Orthornormal basis
* I_FIN_CO = 1 => Final guess is in CAAB basis, 
*          = 2 => Final guess is in Orthornormal basis
*
* Jeppe Olsen, Feb. 20, 2010 from MRCC_VECFNC
*
* Unclean: Internal CI-coefficients are handled
* borh through LUC  and C_REF...
*
* External part: 
* ================
*
* <0!tau^{\dagger} exp(-T) H exp(T) !0>. 
*. The commutator  exp(-T) H exp(T) is terminated after NCOMMU_V commutators
*
*. Internal part:
* ================
*
* <J! exp(-T) H exp(T) - E !0>
*. The commutator  exp(-T) H exp(T) is terminated after NCOMMU_E commutators
*
* (initial version using CI behind the curtains)
*
*
      INCLUDE 'wrkspc.inc'
      REAL*8
     &INPROD
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cei.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      DIMENSION C_REF(N_REF)
*. Specific input
      DIMENSION T(*)
*. Output
      DIMENSION CCVECFNC(*),CCVECFNCI(*)
*
      NTEST = 05
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from MRCC_VECFNCN'
        WRITE(6,*) ' -----------------------'
        WRITE(6,*) ' IREFSPC,ITREFSPC, IT2REFSPC =',
     &               IREFSPC,ITREFSPC, IT2REFSPC
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MRCCVF')
*
      CALL MEMMAN(KLVCC1,N_CC_AMP,'ADDL  ',2,'LCCVC1')
      CALL MEMMAN(KLVCC2,N_CC_AMP,'ADDL  ',2,'LCCVC2')
*
      IF(I_INI_CO.EQ.2) THEN
*. Initial guess is in orthonormal basis, change to CAAB basis
*. Dir in EI in T to Dir in CAAB in VCC1
        CALL TRANS_CAAB_ORTN(WORK(KLVCC1),T,1,2,2,
     &         WORK(KLVCC2),2)
      ELSE
        CALL COPVEC(T,WORK(KLVCC1),N_CC_AMP)
      END IF
*
* 1 : Obtain exp(-T) H exp(T)  !0> and save on LUHC
*
C          EMNTHETO(T,LUOUT,NCOMMU,IREFSPC,ITREFSPC)
      IF(I_APPROX_HCOM_V.EQ.0) THEN
        CALL EMNTHETO(WORK(KLVCC1),LUC,LUHC,NCOMMU_V,IREFSPC,ITREFSPC,
     &                IT2REFSPC)
      ELSE
*. Exact calculation of all terms with upto NCOMMU_V-1 commutators
        CALL EMNTHETO(WORK(KLVCC1),LUC,LUHC,NCOMMU_V-1,IREFSPC,ITREFSPC,
     &                IT2REFSPC)
*. and add contribution from highest commutator
*. Use zero-order Hamiltonian stored in 
        I12 = 1
        CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
        CALL TCOM_H_N(WORK(KLVCC1),LUC,LUHC,NCOMMU_V,IREFSPC,ITREFSPC,
     &               IT2REFSPC,1)
C            TCOM_H_N(T,LUINI,LUUT,NCOMMU,IREFSPC,ITREFSPC,IT2REFSPC,IAC)
        I12 = 2
        CALL SWAPVE(WORK(KINT1),WORK(KFIFA),NINT1)
      END IF
*
* 2 : Obtain  <0!tau^{\dagger} exp(-T) H exp(T) !0> = <LUC!tau^{\dagger}|LUHC>
*
      ICSPC = IREFSPC
      ISSPC = IT2REFSPC
C     WRITE(6,*) ' IREFSPC, IT2REFSPC =', IREFSPC, IT2REFSPC
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Vector on LUC '
        CALL WRTVCD(WORK(KVEC1P),LUC,1,-1)
        WRITE(6,*) ' Vector on LUHC '
        CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
      END IF
*
      ZERO = 0.0D0
      CALL SETVEC(CCVECFNC,ZERO,N_CC_AMP)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUHC,CCVECFNC,2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 'CCVECFNC right after SIGDEN_CC'
        CALL WRTMAT(CCVECFNC,1,N_CC_AMP,1,N_CC_AMP)
      END IF
      XN_CAAB = INPROD(CCVECFNC,CCVECFNC,N_CC_AMP-1)
      WRITE(6,*) ' Norm of CCVEC in CAAB basis = ', XN_CAAB
*
*. 2.5. Project redundant directions out if requested
      IF(I_DO_PROJ_NR.EQ.1) THEN
        IF(NTEST.GE.5)
     &  WRITE(6,*) ' Redundant directions projected out in MRCC...'
C              PROJ_TO_NONRED(VECIN,VECOUT,ITSYM,VECSCR)
        CALL PROJ_TO_NONRED(CCVECFNC,WORK(KLVCC1),1,WORK(KLVCC2))
        CALL COPVEC(WORK(KLVCC1),CCVECFNC,N_CC_AMP)
      END IF
*. The energy obtained from the external vectorfunction
      E_EXT = CCVECFNC(N_CC_AMP)
      IF(NTEST.GE.5)
     &WRITE(6,*) ' Energy from external part of vecfnc ', E_EXT
C    &WRITE(6,*) ' Energy from external part of vecfnc ', E_EXT+ECORE
*. And clear element corresponding to N_CC_AMP- not really part of 
*. vectorfunction
      CCVECFNC(N_CC_AMP) = 0.0D0
*
*. 2.6: Transform if required vector function to orthonormal basis
      IF(I_FIN_CO.EQ.2) THEN
*. Vecfunc in CAAB in VCC5 to Vecfunc in EI in VCC2
*. zero-order state is not to be included
        N_ZERO_EIM = N_ZERO_EI - 1
        CALL TRANS_CAAB_ORTN(CCVECFNC,WORK(KLVCC1),1,1,2,
     &                       WORK(KLVCC2),1)
        CALL COPVEC(WORK(KLVCC1),CCVECFNC,N_ZERO_EIM)
*. To be sure..
        CCVECFNC(N_ZERO_EI) = 0.0D0
      END IF
*
* 3 : Contract  exp(-T) H exp(T) |0> to reference space and save on LUHC
*     to obtain part of internal part of MRCC vector function 
*
      IF((NCOMMU_E.NE.NCOMMU_V.AND.
     &  .NOT.(NCOMMU_E.EQ.4.AND.NCOMMU_V.GT.4)) .OR.
     &        I_APPROX_HCOM_V.NE.I_APPROX_HCOM_E) THEN 
*. Recalculate Internal part of MRCC vector function 
        IF(NTEST.GE.10)
     &  WRITE(6,*) ' Internal part of vector-function recalculated'
        CALL HEFF_INT_TV_ICCC(T,N_REF,NCOMMU_E,I_APPROX_HCOM_E,
     &  dbl_mb(VEC1P),dbl_mb(KVEC2P),IREFSPC,ITREFSPC,IT2REFSPC,
     &  0.0D0,C_REF,CCVECFNCI)
      ELSE 
        CALL EXPCIV(IREFSM,IT2REFSPC,LUHC,IREFSPC,LUSC34,-1,
     /              LUSC35,1,1,IDC,0)
        CALL REWINO(LUHC)
        CALL FRMDSCN(CCVECFNCI,-1,-1,LUHC)
      END IF
*. Energy from internal part
      E_INT = INPROD(C_REF,CCVECFNCI,N_REF)
      IF(NTEST.GE.5)
     &WRITE(6,*) ' Energy from internal part of vecfnc ', E_INT
C    &WRITE(6,*) ' Energy from internal part of vecfnc ', E_INT+ECORE
*. And the internal vector function
      ONE = 1.0D0
      FACTOR = -E_INT
      CALL VECSUM(CCVECFNCI,CCVECFNCI,C_REF,ONE,FACTOR,N_REF)
*. Zero internal if requested
*  - after all the work... - could be done in a more elegant way...
      IF(I_FIX_INTERNAL.EQ.1) THEN
*. set internal gradient to zero
        ZERO = 0.0D0
        CALL SETVEC(CCVECFNCI,ZERO,N_REF)
        WRITE(6,*) ' Internal gradient set to zero '
      END IF
*
      IF(NTEST.GE.100) THEN
*
        IF(I_INI_CO.EQ.1) THEN
          WRITE(6,*) ' Input T-coefficients in CAAB basis'
          CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
        ELSE
          WRITE(6,*) ' Input T-coefficients in ortn. basis'
          CALL WRTMAT(T,1,N_ZERO_EI,1,N_ZERO_EI)
        END IF
*
        IF(I_FIN_CO.EQ.1) THEN
          WRITE(6,*) ' MRCC Vector function, external part (CAAB) '
          CALL WRTMAT(CCVECFNC,1,N_CC_AMP,1,N_CC_AMP)
        ELSE
          WRITE(6,*) ' MRCC Vector function, external part (ortn) '
          CALL WRTMAT(CCVECFNC,1,N_ZERO_EI,1,N_ZERO_EI)
        END IF
*
        WRITE(6,*) 'MRCC Vector function,internal part'
        CALL WRTMAT(CCVECFNCI,1,N_REF,1,N_REF)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'MRCCVF')
*
     
      RETURN
      END 
      SUBROUTINE ICCC_OPT_SIMULT_ONB(
     &        IREFSPC,ITREFSPC,IT2REFSPC,I_SPIN_ADAPT,
     &        IREFROOT,T_EXT,C_0,INI_IT,IFIN_IT,VEC1,VEC2,IDIIS,
     &        C_REF,N_REF,I_DO_COMP,CONVERL,VTHRES,I_REDO_INT,
     &        EFINAL,VNFINAL,CONVERG,SCR_SBSPJA,MXVEC_SBSPJA)

*
* Master routine for Internal Contraction Coupled Cluster 
*
* It is assumed that the excitation manifold produces 
* states that are orthogonal to the reference so 
* no projection is carried out
*
* Routine is allowed to leave without turning the lights off,
* i.e. leave routine with all allocations and marks intact.
*: Thus : Allocations are only done if INI_IT = 1
*        Deallocations are only done if IFIN_IT = 1
*
*. Preconditioners are only calculated if INI_IT = 1
*
* IF I_REDO_INT = 1, the internal states are recalculated at start
*
* IF IDIIS.EQ.1, DIIS is used
*         .EQ.2, CROP is used to accelerate convergence 
* 
*
* Jeppe Olsen, Aug. 2005, modified aug 2009 - also in Washington
*              Redo of internal states: Sept. 2009 in Sicily
*              Subspace Jacobian added: Oct. 2009
*              ONB version: March 2010
*
* ONB: Orthonormal basis version: all calc in zero-order basis
*
*. for DIIS units LUSC37 and LUSC36 will be used for storing vectors
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cei.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
*. Temporary  array for debugging
      REAL*8 XNORM_EI(1000), XJ1(1000),XJ2(1000)
*
      LOGICAL CONVERL,CONVERG
*. Converl: is local iterative procedure for given internal states converged
*. converg: is global iterative procedure converged
      REAL*8
     &INPROD,INPRDD
*. Input and Output : Coefficients of internal and external correlation 
      DIMENSION T_EXT(*), C_REF(*)
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
      COMMON/CLOCAL2/KVEC1,KVEC2,MXCJ,
     & KLVCC1,KLVCC2,KLVCC3,KLVCC4,KLVCC5,KLSMAT,KLXMAT,KLJMAT,KLU,KLL,
     & NSING,NNONSING,KLCDIIS,KLC_INT_DIA,KLDIA,KLVCC6,KLVCC7,KLVCC8,
     & NVECP,NVEC,KLA_CROP,KLSCR_CROP
*. Scratch for CI behind the curtain 
      DIMENSION VEC1(*),VEC2(*)
*. Scratch for subspace Jacobian
      DIMENSION SCR_SBSPJA(*)
*. Threshold for convergence of norm of Vectorfuntion

C     WRITE(6,*) ' ICCC_OPT_SIMULT: I_DO_COMP =', I_DO_COMP
C     WRITE(6,*) ' ICCC_OPT_SIMULT: MAXIT,MAXITM =', MAXIT,MAXITM
      WRITE(6,*) ' ICCC_OPT_SIMULT: I_DO_SBSPJA, MXVEC_SBSPJA = ', 
     &                              I_DO_SBSPJA, MXVEC_SBSPJA
      NCAAB = NDIM_EI
      WRITE(6,*) ' NCAAB og NDIM_EI = ', NCAAB, NDIM_EI
*. We will not include the unit-operator so  ???
*. Project on nonredundant space
      I_DO_PROJ_NR = 0
*. For file access
      LBLK = -1
      NTEST = 5
      IF(NTEST.GE.2) THEN
      WRITE(6,*) 
     &  ' Simultaneous optimization of internal and external parts '
        WRITE(6,*) 
     &  ' ========================================================='
        WRITE(6,*)
        WRITE(6,*) ' CROP/DIIS performed in ortn. zero-order basis'
        WRITE(6,*) ' Reference space is ', IREFSPC
        WRITE(6,*) ' Space for evaluating general operators  ', ITREFSPC
        WRITE(6,*) ' Space for T times reference space  ', IT2REFSPC
        WRITE(6,*) ' Number of parameters in CAAB basis ', 
     &             N_CC_AMP
        WRITE(6,*) ' Number of parameters in spincoupled/ort basis ', 
     &             NSPA
        WRITE(6,*) ' Number of coefficients  in internal space ', N_REF
        WRITE(6,*) ' INI_IT, IFIN_IT = ', INI_IT, IFIN_IT
        WRITE(6,*) ' Max. number microiterations per macro ', MAXIT
        WRITE(6,*) ' Max. number of macroiterations        ', MAXITM
        WRITE(6,*) ' Number of vectors allowed in subspace ', MXCIVG
        WRITE(6,*) ' Number of vectors allowed in initial subspace ', 
     &               MXVC_I
        IF(IDIIS.EQ.1) THEN
          WRITE(6,*)' DIIS optimization'
        ELSE IF (IDIIS.EQ.2) THEN
          WRITE(6,*)' CROP optimization'
        END IF
*
        IF(I_DO_PROJ_NR.EQ.1) THEN
          WRITE(6,*) ' Redundant directions projected out'
        ELSE
          WRITE(6,*) ' No projection of redundant directions'
        END IF
*
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial T_ext-amplitudes '
        CALL WRTMAT(T_EXT,1,N_CC_AMP,1,N_CC_AMP)
        WRITE(6,*) ' Initial C_int-amplitudes '
        CALL WRTMAT(C_REF,1,N_REF,1,N_REF)
      END IF
*. Allowed number of iterations
      NNEW_MAX = MAXIT
      MAXITL = NNEW_MAX
*
      NVAR_CAAB = N_CC_AMP + N_REF
      IF(INI_IT.EQ.1) THEN
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICC_CM')
        CALL MEMMAN(KLVCC1,NVAR_CAAB,'ADDL  ',2,'VCC1  ')
        CALL MEMMAN(KLVCC2,NVAR_CAAB,'ADDL  ',2,'VCC2  ')
        CALL MEMMAN(KLVCC3,NVAR_CAAB,'ADDL  ',2,'VCC3  ')
        CALL MEMMAN(KLVCC4,NVAR_CAAB,'ADDL  ',2,'VCC4  ')
        CALL MEMMAN(KLVCC5,NVAR_CAAB,'ADDL  ',2,'VCC5  ')
        CALL MEMMAN(KLVCC6,2*NVAR_CAAB,'ADDL  ',2,'VCC6  ')
*. Just a few extra to be on the safe side when programming EI
*. approach
        CALL MEMMAN(KLVCC7,NVAR_CAAB,'ADDL  ',2,'VCC5  ')
        CALL MEMMAN(KLVCC8,NVAR_CAAB,'ADDL  ',2,'VCC5  ')
        CALL MEMMAN(KLDIA,NVAR_CAAB,'ADDL  ',2,'DIAORT')
*. Space for DIIS/CROP
        IF(IDIIS.EQ.1) THEN
          CALL MEMMAN(KLCDIIS,MAXITL,'ADDL ',2,'CDIIS ') 
        ELSE IF(IDIIS.EQ.2) THEN
          CALL MEMMAN(KLA_CROP,MXCIVG*(MXCIVG+1)/2,'ADDL  ',2,'A_CROP')
          LEN_SCR_CROP = 3*MXCIVG*MXCIVG + 3*MAX(MXCIVG,NVAR_CAAB)
          CALL MEMMAN(KLSCR_CROP,LEN_SCR_CROP,'ADDL  ',2,'S_CROP')
C?        WRITE(6,*) ' KLA_CROP,KLSCR_CROP, a =', KLA_CROP,KLSCR_CROP
        END IF
*. Space Diagonal for internal part
        CALL MEMMAN(KLC_INT_DIA,N_REF,'ADDL ',2,'C_DIA ')
      END IF
*.    ^ End if INI_IT.EQ.1
*
*======================================
* 0: Redo internal states if required
* =====================================
*
      IF(I_REDO_INT.EQ.1) THEN
        CALL GET_INTERNAL_STATES_OUTER
        N_INT_MAX = IMNMX(WORK(KL_N_INT_FOR_SE),N_EXTOP_TP*NSMOB,2)
*. Largest number of zero-order states of given sym and external type
        N_ORTN_MAX = IMNMX(WORK(KL_N_ORTN_FOR_SE),N_EXTOP_TP*NSMOB,2)
        WRITE(6,*) ' N_INT_MAX, N_ORTN_MAX = ', N_INT_MAX, N_ORTN_MAX
*. Largest transformation block 
        N_XEO_MAX = N_INT_MAX*N_ORTN_MAX
        IF(NTEST.GE.10) 
     &  WRITE(6,*) ' Largest (EL,ORTN) block = ', N_XEO_MAX
*. Number of zero-order states - does now include the unit-operator
        N_ZERO_EI = N_ZERO_ORDER_STATES(WORK(KL_N_ORTN_FOR_SE),
     &             WORK(KL_NDIM_EX_ST),N_EXTOP_TP,1)
        NVAR = N_ZERO_EI + N_REF
        NSPA = N_ZERO_EI
        NSPAM1 = NSPA - 1
*. Adresses of the unit op
        IUNI_AD = N_ZERO_EI
       IF(NTEST.GE.10) WRITE(6,*) 
     & ' Number of zero-order states with sym 1 = ', N_ZERO_EI
      END IF
*
*. Memory for complete matrices can now be defined
*. Complete matrices for external part, three used pt
      IF(INI_IT.EQ.1.AND.I_DO_COMP.EQ.1) THEN
        LEN = N_ZERO_EI**2
        CALL MEMMAN(KLSMAT,LEN,'ADDL  ',2,'SMAT  ')
        CALL MEMMAN(KLXMAT,LEN,'ADDL  ',2,'XMAT  ')
        CALL MEMMAN(KLJMAT,LEN,'ADDL  ',2,'JMAT  ')
*. Storage for LU decomposition of J
        LEN = N_ZERO_EI*(N_ZERO_EI+1)/2
          CALL MEMMAN(KLL,LEN,'ADDL  ',2,'L     ')
          CALL MEMMAN(KLU,LEN,'ADDL  ',2,'U     ')
        ELSE
*. Space for diagonal- space is allocated also for CI part.
        END IF
*
* ============================================================
* 1 : Prepare preconditioners for external and internal parts 
* ============================================================
*
* --------------------
*. 1a : External part 
* --------------------
*
*. Identify the unit  operator i.e. the operator with 
*. zero creation and annihilation operators
      IDOPROJ = 0
*. Construct metric (once again ..)
*. Prepare the routines used in COM_SH
*. Not used here
      C_0X = 0.0D0
      KLTOPX = -1
*. Used
      NREFX = N_REF
      IREFSPCX = IREFSPC
*. Space to be used for evaluating metric : If T = 0, then IT2REFSPC is sufficient
      ITREFSPCX = ITREFSPC
      ITREFSPCX = IT2REFSPC
*
      NCAABX = N_CC_AMP
      NSPAX = N_ZERO_EI
      IPROJSPCX = IREFSPC
*. Unitoperator in SPA order ... Please check ..
      IUNIOPX = 0
*
      NVAR_EXT = N_ZERO_EI - 1
      IF(I_DO_COMP.EQ.1) THEN
*
*. Set up or read in Jacobian in orthonormal basis
*
        IF(INI_IT.EQ.1.AND.IREADSJ.EQ.0) THEN
*. Construct exact or approximate Jacobian
          IF(NCOMMU_J.EQ.1) THEN
*. I assume that the  space before ITREFSPC contains T*IREFSPC 
           ITREFSPC_L = ITREFSPC - 1
           WRITE(6,*) ' Space used for approximate J ', ITREFSPC_L
*. Do not include zero-order state
           INCLUDE0 = 0
           CALL COM_JAC_1COM(IREFSPC,IT2REFSPC,WORK(KLJMAT),INCLUDE0)
          ELSE 
*. More than one commutator, so J depends on T
           CALL COM_JMRCC(T_EXT,NCOMMU_J,I_APPROX_HCOM_J,
     &          WORK(KLJMAT),WORK(KLVCC1),WORK(KLVCC2), WORK(KLVCC3),
     &          WORK(KLVCC4),N_CC_AMP,NSPAM1,N_ZERO_EI,IREFSPC,
     &          ITREFSPC,WORK(KLXMAT) )
          END IF
*         ^ End if more than one commutator
          WRITE(LU_SJ) (WORK(KLJMAT-1+IJ),IJ=1,NVAR_EXT*NVAR_EXT)
*. Rewind to flush buffer
          CALL REWINO(LU_SJ)
        ELSE
*. Read Approximate Jacobian in from LU_SJ
          CALL REWINO(LU_SJ)
          READ(LU_SJ) (WORK(KLJMAT-1+IJ),IJ=1,NVAR_EXT*NVAR_EXT)
        END IF
*       ^ End if matrix should be constructed or read in
        I_ADD_SHIFT = 0
        IF(I_ADD_SHIFT.EQ.1) THEN
*. Add a shift to the diagonal of J
          SHIFT = 10.0D0
          WRITE(6,*) ' A shift will be added to initial Jacobian'
          WRITE(6,'(A,E14.7)') ' Value of shift = ', SHIFT
          CALL ADDDIA(WORK(KLJMAT),SHIFT,NVAR_EXT,0)
        END IF
*       ^ End if shift should be added
*
        I_DIAG_J = 0
        IF(I_DIAG_J.EQ.1) THEN
*. Obtain eigenvalues of approximate Jacobian
*. S-matrix is not used anymore to use this space for 
*. diagonalization 
         WRITE(6,*) ' Approximate Jacobian will be diagonalized '
         CALL COPVEC(WORK(KLJMAT),WORK(KLSMAT),NVAR_EXT*NVAR_EXT)
         CALL EIGGMT3(WORK(KLSMAT),NVAR_EXT,WORK(KLVCC1),WORK(KLVCC2),
     &                XDUM,XDUM,XDUM,WORK(KLVCC3),WORK(KLVCC6),1,0)
         WRITE(6,*) ' Real and imaginary part of eigenvalues of J '
         WRITE(6,*) ' ========================================== '
         CALL WRT_2VEC(WORK(KLVCC1),WORK(KLVCC2),NVAR_EXT)
        END IF
*. Obtain LU-Decomposition of Jacobian 
        CALL LULU(WORK(KLJMAT),WORK(KLL),WORK(KLU),NVAR_EXT)
      ELSE
        IF(INI_IT.EQ.1) THEN
*. Complete matrix is not constructed, rather just a diagonal
*. Obtain diagonal of H 
          CALL GET_DIAG_H0_EI(WORK(KLDIA))
*. The last element in KLDIA is the zero-order energy
          E0 = WORK(KLDIA-1+N_ZERO_EI)
          IF(NTEST.GE.0)
     &    WRITE(6,*) ' Zero-order energy  ', E0
*. To get diagonal approximation to J, subtract E0
          DO I = 1, N_ZERO_EI
           WORK(KLDIA-1+I) = WORK(KLDIA-1+I) - E0
          END DO
*. The last term in KLDIA corresponds to the zero-order state.
*. This will not contribute, but to eliminate errors occuring 
*. from dividing by zero 
          WORK(KLDIA-1+N_ZERO_EI) = 300656.0
*. Check for diagonal values close to zero, and shift these
          XMIN = 0.2D0
          CALL MODDIAG(WORK(KLDIA),N_ZERO_EI,XMIN)
C              MODDIAG(H0DIAG,NDIM,XMIN)
*. And save on LU_SJ
          CALL VEC_TO_DISC(WORK(KLDIA),N_ZERO_EI-1,1,LBLK,LU_SJ)
*. test norm of the E-blocks of diagonal
          WRITE(6,*) ' Norm of various E-blocks of diagonal'
          CALL NORM_T_EI(WORK(KLDIA),2,1,XNORM_EI,1)
C              NORM_T_EI(T,IEO,ITSYM,XNORM_EI,IPRT)
          IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Diagonal J-approx in ort. zero-order basis'
           CALL WRTMAT(WORK(KLDIA),1,N_ZERO_EI,1,N_ZERO_EI)
          END IF
        END IF
*.      ^ End if it was first iteration
      END IF
*     ^ End of complete or diagonal matrix should be set up
*
* ---------------------
*. 1b : internal part  - Fetch in all macroiterations
* ---------------------
*
      CALL REWINO(LUDIA)
      CALL FRMDSCN(WORK(KLC_INT_DIA),-1,-1,LUDIA)
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Diagonal preconditioner for internal correlation'
         CALL WRTMAT(WORK(KLC_INT_DIA),1,N_REF,1,N_REF)
      END IF
*
      IF(IDIIS.EQ.1.OR.(IDIIS.EQ.2.AND.INI_IT.EQ.1)) THEN
        CALL REWINO(LUSC37)
        CALL REWINO(LUSC36)
      END IF
*. Ensure proper defs
      I12 = 2
      ICSM = IREFSM
      ISSM = IREFSM
      IF(IUSE_PH.EQ.1) THEN
        CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
      END IF
*
      IF(NTEST.GE.100)
     &  WRITE(6,*) ' After const of precond: ITREFSPC, IT2REFSPC =',
     &  ITREFSPC, IT2REFSPC
*
*. Transformation of T from CAAB to orthonormal basis should 
*. initialize procedure
      CALL TRANS_CAAB_ORTN(T_EXT,WORK(KLVCC1),1,1,2,
     &         WORK(KLVCC2),2)
      CALL COPVEC(WORK(KLVCC1),T_EXT,N_ZERO_EI)
      XTNORM_INI = SQRT(INPROD(T_EXT,T_EXT,N_ZERO_EI))
      WRITE(6,*) ' Norm of initial T-vector', XTNORM_INI
*
*. Loop over iterations 
      WRITE(6,*)
      WRITE(6,*) ' -------------------------- '
      WRITE(6,*) ' Entering optimization part ' 
      WRITE(6,*) ' -------------------------- '
      WRITE(6,*)
*. Number of vectors in initial space for DIIS/CROP optimization
      IF(INI_IT.EQ.1) THEN
        NVECP = 0
        NVEC  = 0
      END IF
*. (If INI_IT .ne. 0, MXVC_I vectors from previous macro are used)
      IF(I_DO_SBSPJA.EQ.1) THEN
*. Initialize files that will be used for subspace Jacobian)
        WRITE(6,*) ' LU_CCVECT,LU_CCVECF, LU_CCVECFL = ',
     &               LU_CCVECT,LU_CCVECF, LU_CCVECFL
        CALL REWINO(LU_CCVECT)
        CALL REWINO(LU_CCVECF)
        CALL REWINO(LU_CCVECFL)
      END IF
*
      DO IT = 1, NNEW_MAX
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 
          WRITE(6,*) ' Information for iteration ', IT
          WRITE(6,*) 
        END IF
        IF(IT.EQ.1) THEN
          MXVC_SUB = MXVC_I
        ELSE
          MXVC_SUB = MXCIVG
        END IF
*
*
* ==================================================================
*. Construct vectorfunction/gradient for external and internal parts
* ==================================================================
*
*. CC vector function for external part  in VCC5 
C?      WRITE(6,*) ' NCAAB before MRCC.. ', NCAAB
        CALL MRCC_VECFNCN(WORK(KLVCC5),T_EXT,
     &       IREFSPC,ITREFSPC,IT2REFSPC,WORK(KLVCC5+N_CC_AMP),
     &       C_REF, N_REF,I_DO_PROJ_NR, 
     &       E_INT,E_EXT,ECORE,2,2)
        CALL COPVEC(WORK(KLVCC5+N_CC_AMP),WORK(KLVCC5+N_ZERO_EI),
     &              N_REF)
*
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' Norm of various E-blocks of Vecfnc'
            CALL NORM_T_EI(WORK(KLVCC5),2,1,XNORM_EI,1)
          END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) 
     &    ' The CC vector function  including internal part'
          CALL WRTMAT(WORK(KLVCC5),1,NVAR,1,NVAR,1)
        END IF
        IF(NTEST.GE.10) WRITE(6,'(A,I4,2E22.15)')
     &  ' It, Energy from external and internal ', IT, E_EXT ,
     &        E_INT
C    &  ' It, Energy from external and internal ', IT, E_EXT + ECORE,
C    &        E_INT+ECORE
        VCFNORM_EXT =SQRT(INPROD(WORK(KLVCC5),WORK(KLVCC5),N_ZERO_EI))
        VCFNORM_INT = SQRT(
     &  INPROD(WORK(KLVCC5+N_ZERO_EI),WORK(KLVCC5+N_ZERO_EI),N_REF)) 
*. Update energy and residual norms
        VNFINAL = VCFNORM_EXT+VCFNORM_INT
        E = E_INT 
        EFINAL = E_INT 
*. Converged?
        IF(VCFNORM_EXT+VCFNORM_INT.LE.VTHRES) THEN
*. Local iterative procedure converged
          CONVERL = .TRUE.
*. Is global procedure also converged?
          IF((I_REDO_INT.NE.1            ) .OR.
     &       (I_REDO_INT.EQ.1.AND.IT.EQ.1)) THEN
             CONVERG = .TRUE.
          END IF
          WRITE(6,*) ' Iterative procedure converged'
          WRITE(6,'(A,I4,E22.15,2E12.5)')
     &  ' It, energy ,  vecfnc_ext, vecfnc_int ', 
     &    IT, E, VCFNORM_EXT, VCFNORM_INT
          GOTO 1001
        END IF
*       ^ End if local procedure is converged
*
* ======================================================================
*. Save vectorfunction in form that will be used in later subspace opt.
* ======================================================================
*
*
        IF(I_DO_SBSPJA.EQ.1) THEN
*
* Has not been bebugged for Zero-order states
*. Save Vectorfunction and change in vectorfunction 
*. if subspace Jacobian is in use
          N_ZERO_EIM = N_ZERO_EI - 1
          IF(IT.GE.2)  THEN
*. Read previous vectorfunction in VCC7 from CCVECFL
            CALL VEC_FROM_DISC(WORK(KLVCC7),N_ZERO_EIM,1,LBLK,
     &           LU_CCVECFL)
            ONE = 1.0D0
            ONEM =-1.0D0
*. Store in VCC7: Delta V  = Vecfnc(ITER) - Vecfnc(ITER-1)
            CALL VECSUM(WORK(KLVCC7),WORK(KLVCC5),WORK(KLVCC2),
     &                  ONEM,ONE,N_ZERO_EIM)
*. Add CCVF(X_{i+1})-CCVF(X_{i}) as vector IT-1 in FILE LU_CCVECF
            CALL SKPVCD(LU_CCVECF,IT-2,WORK(KLVCC6),1,LBLK)
            CALL VEC_TO_DISC(WORK(KLVCC7),N_ZERO_EIM,0,LBLK,LU_CCVECF)
          END IF
*. Save current vector-function in EO form in LU_CCVECFL
          CALL VEC_TO_DISC(WORK(KLVCC5),N_ZERO_EIM,1,LBLK,LU_CCVECFL)
        END IF
*       ^ End if subspace method in use
*
* ========================================================
* Diis/CROP/SBSPJA based on current and previous vectors 
* ========================================================
*
* Subspace is in this version saved in orthonormal basis
*
        IF(IDIIS.EQ.1.OR.IDIIS.EQ.2) THEN
*. It is assumed that DIIS left the file at end of file 
*. T_ext,C_int on LUSC37, VECFNC on LUSC36
          CALL COPVEC(T_EXT,WORK(KLVCC1),N_ZERO_EI)
          CALL COPVEC(C_REF,WORK(KLVCC1+N_ZERO_EI),N_REF)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Combined T_ext, C_int coefficients '
            CALL WRTMAT(WORK(KLVCC1),1,NVAR,1,NVAR)
          END IF
          CALL VEC_TO_DISC(WORK(KLVCC1),NVAR,0,-1,LUSC37)
          CALL VEC_TO_DISC(WORK(KLVCC5),NVAR,0,-1,LUSC36)
        END IF
*. We have now a number of vectors in LUSC36, find combination with lowest 
*. norm 
*. DIIS:
        IF(IDIIS.EQ.1) THEN
*. Simple DIIS with no restart
          CALL DIIS_SIMPLE(LUSC36,IT,NVAR,WORK(KLCDIIS))
*. Obtain combination of parameters given in CDIIS
          CALL MVCSMD(LUSC37,WORK(KLCDIIS),LUSC39,LUSC38,
     &                WORK(KLVCC1),WORK(KLVCC2),IT,1,-1)
          CALL VEC_FROM_DISC(WORK(KLVCC1),NVAR,1,-1,LUSC39)
          CALL COPVEC(WORK(KLVCC1),T_EXT,N_ZERO_EI)
          CALL COPVEC(WORK(KLVCC1+N_ZERO_EI),C_REF,N_REF)
*. Calculate new vectorfunction in VCC5 for T_EXT  and C_INT using sums 
          CALL MVCSMD(LUSC36,WORK(KLCDIIS),LUSC39,LUSC38,
     &                WORK(KLVCC1),WORK(KLVCC2),IT,1,-1)
          CALL VEC_FROM_DISC(WORK(KLVCC5),NVAR,1,-1,LUSC39)
        ELSE IF(IDIIS.EQ.2) THEN
*. CROP:
*. The CROP version of DIIS
*. Matrices are reconstructed in each IT
          IDIRDEL = 1
          NVEC = NVEC + 1
*. Note: NVECP is number of vectors for which subspace matrix 
*. has been constructed and saved- CROP updates this
*. Obtain improved amplitudes in VCC1, improved vectorfunction in VCC4
          CALL CROP(NVEC,NVECP,MXVC_SUB,NVAR,LUSC36,LUSC37,
     &         WORK(KLA_CROP),
     &         WORK(KLVCC4),WORK(KLVCC1),WORK(KLSCR_CROP),LUSC39,
     &         IDIRDEL)
C     CROP(NVEC,NVECP,MXNVEC,NDIM,LUE,LUP,A,
C    &                EOUT,POUT,SCR,LUSCR,IDIRDEL)
*Change of T-coefs 
          ONE = 1.0D0
          ONEM = -1.0D0
          CALL VECSUM(WORK(KLVCC1),WORK(KLVCC1),T_EXT,ONE,ONEM,
     &                N_ZERO_EI)
*. Update of external coefficients
*. Check if change is to large..
          XNORM = SQRT(INPROD(WORK(KLVCC1),WORK(KLVCC1),N_ZERO_EI))
          WRITE(6,*) ' Norm of CROP external correction ', XNORM
          XNORM_MAX = 0.5D0
          I_DO_SCALE = 1
          IF(XNORM.GT.XNORM_MAX.AND.I_DO_SCALE.EQ.1) THEN
            WRITE(6,*) 
     &      ' CROPStep is scaled: from and to to ', XNORM,XNORM_MAX
            FACTOR = XNORM_MAX/XNORM
            CALL VECSUM(T_EXT,T_EXT,WORK(KLVCC1),ONE,FACTOR,N_ZERO_EI)
*. Well, if change in parameters was reduced, then change in 
*. vector function should also be reduced
* VEC5 = VEC5 + Factor*(vec4-vec5) = (1-factor)vec5 + factor*vec4
            FACTOR5 = 1.0D0-FACTOR
            FACTOR4 = FACTOR
            CALL VECSUM(WORK(KLVCC5),WORK(KLVCC5),WORK(KLVCC4),
     %                  FACTOR5,FACTOR4, N_ZERO_EI)
          ELSE
            CALL VECSUM(T_EXT,T_EXT,WORK(KLVCC1),ONE,ONE,N_ZERO_EI)
            CALL COPVEC(WORK(KLVCC4),WORK(KLVCC5),N_ZERO_EI)
          END IF
*. And update internal (CI-)coefficients 
          CALL COPVEC(WORK(KLVCC1+N_ZERO_EI),C_REF,N_REF)
          XNORM  = INPROD(C_REF,C_REF,N_REF)
          FACTOR = 1.0D0/SQRT(XNORM)
          FACTOR = 1.0D0
          WRITE(6,*) ' No normalization of C_REF in CROP'
          CALL SCALVE(C_REF,FACTOR,N_REF)
*. And scale CI-vector function
          CALL COPVEC(WORK(KLVCC4+N_ZERO_EI),WORK(KLVCC5+N_ZERO_EI),
     &                N_REF)
          CALL SCALVE(WORK(KLVCC5+N_ZERO_EI),FACTOR,N_REF)
        END IF
*.      ^ End of DIIS/CROP should be used 
        VCFNORM = SQRT(INPROD(WORK(KLVCC5),WORK(KLVCC5),NVAR))
        IF(NTEST.GE.5) WRITE(6,'(A,I4,1E12.5)')
     &  ' From DIIS/CROP : It, norm of approx vecfnc  ',
     &  IT,  VCFNORM 
*
* ===================================================================
* Obtain new direction by applying preconditioners to approx vecfunc
* ===================================================================
*
* --------------
* External part
* --------------
*
*  multiply with diagonal transform 
*. Vectorfunction
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' Norm of various E-blocks of apr Vecfnc'
            CALL NORM_T_EI(WORK(KLVCC5),2,1,XNORM_EI,1)
          END IF
*
        IF(I_DO_COMP.EQ.1) THEN
*
*. Complete matrix approximation to J in use
*
*. Solve Linear equations J Delta = - Vecfnc, store solution in VCC1
          ONEM = -1.0D0
          CALL SCALVE(WORK(KLVCC5),ONEM,NVAR_EXT)
          CALL LINSOL_FROM_LUCOMP(WORK(KLL),WORK(KLU),WORK(KLVCC5),
     &         WORK(KLVCC1),NVAR_EXT,WORK(KLVCC2))

*. And no correction for the zero-order state
            WORK(KLVCC1-1+IUNI_AD) = 0.0D0
        ELSE 
*
*. Complete matrices not in use..
*
          IF(I_DO_SBSPJA.EQ.0) THEN
*·  New direction = -diag-1 * Vecfunc
            DO I = 1, N_ZERO_EI
              WORK(KLVCC1-1+I) = - WORK(KLVCC5-1+I)/WORK(KLDIA-1+I)
            END DO
*. And no correction for the zero-order state
            WORK(KLVCC1-1+IUNI_AD) = 0.0D0
            IF(NTEST.GE.10) THEN
              WRITE(6,*) ' Norm of various E-blocks of step'
              CALL NORM_T_EI(WORK(KLVCC1),2,1,XNORM_EI,1)
            END IF
          ELSE
*. Use subspace Jacobian to solve equations
*. Multiply current CC vector function with approximate Jacobian
*. to obtain new step
            NSBSPC_VEC = IT-1
            MAXVEC = MXVEC_SBSPJA
            CALL APRJAC_TV(NSBSPC_VEC,LU_CCVECFL,LUSC41,LU_CCVECT,
     &                     LU_CCVECF,LU_SJ,WORK(KLVCC6),WORK(KLVCC7),
     &                     SCR_SBSPJA,N_ZERO_EIM,LUSC43,LUSC44,
     &                     MAXVEC)
C                APRJAC_TV(NVEC,LUIN,LUOUT,LUVEC,LUJVEC,
C    &                     LUJDIA,VEC1,VEC2,SCR,N_CC_AMP,LUSCR,LUSCR2,
C    &                     MAXVEC)
*. The new correction vector is now residing in LUSC41,
*. Fetch and multiply with -1
            CALL VEC_FROM_DISC(WORK(KLVCC1),N_ZERO_EIM,1,LBLK,LUSC41)
            ONEM = -1.D0
            CALL SCALVE(WORK(KLVCC1),ONEM,N_ZERO_EIM)
*. And no correction for the zero-order state
            WORK(KLVCC1-1+IUNI_AD) = 0.0D0
*. Add step to LU_CCVECT for future use
            CALL SKPVCD(LU_CCVECT,IT-1,WORK(KLVCC6),1,LBLK)
            CALL VEC_TO_DISC(WORK(KLVCC1),N_ZERO_EIM,0,LBLK,LU_CCVECT)
          END IF
*.        ^ End if subspace Jacobian used for generating new step
        END IF
*       ^ End of switch between complete matrices and not complete
*       matrices
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' direction in ort zero-order basis'
          CALL WRTMAT(WORK(KLVCC1),1,N_ZERO_EI,1,N_ZERO_EI)
        END IF
*. Norm of change
        XNORM = SQRT(INPROD(WORK(KLVCC1),WORK(KLVCC1),N_ZERO_EI))
        IF(NTEST.GE.10) WRITE(6,*) ' Norm of correction ', XNORM
        XNORM_MAX = 0.5D0
        I_DO_SCALE = 1
        IF(XNORM.GT.XNORM_MAX.AND.I_DO_SCALE.EQ.1) THEN
          WRITE(6,*) 
     &    ' Step is scaled: from and to to ', XNORM,XNORM_MAX
          FACTOR = XNORM_MAX/XNORM
          CALL SCALVE(WORK(KLVCC1),FACTOR,N_ZERO_EI)
          XNORM = XNORM_MAX
          IF(I_DO_SBSPJA.EQ.1) THEN
*. Well, step was scaled, read in EI form of step and scale this
            CALL SKPVCD(LU_CCVECT,IT-2,WORK(KLVCC2),1,LBLK)
            CALL VEC_FROM_DISC(WORK(KLVCC2),N_ZERO_EIM,0,LBLK,LU_CCVECT)
            CALL SCALVE(WORK(KLVCC2),FACTOR,N_ZERO_EIM)
            CALL SKPVCD(LU_CCVECT,IT-2,WORK(KLVCC2),1,LBLK)
            CALL VEC_TO_DISC(WORK(KLVCC2),N_ZERO_EIM,0,LBLK,LU_CCVECT)
          END IF
        END IF
*. And update the T-coefficients
        ONE = 1.0D0
        CALL VECSUM(T_EXT,T_EXT,WORK(KLVCC1),ONE,ONE,N_ZERO_EI)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated T-coefficients in ortn. basis '
          CALL WRTMAT(T_EXT,1,N_ZERO_EI,1,N_ZERO_EIP)
        END IF
*
* --------------
* Internal part
* --------------
*
        IF(N_REF.EQ.1) THEN
          C_REF(1) = 1
          XNORM_CI = 0.0D0
        ELSE
          DO I = 1, N_REF
           XNORM_CI = 0.0D0
           IF(ABS(WORK(KLC_INT_DIA-1+I)-E).GE.1.0D-10) THEN
             DELTA = 
     &       - WORK(KLVCC5+N_ZERO_EI-1+I)/(WORK(KLC_INT_DIA-1+I)-E)
             XNORM_CI = XNORM_CI + DELTA**2
             C_REF(I) = C_REF(I)  + DELTA
           END IF
          END DO
        END IF
        XNORM_CI = SQRT(XNORM_CI)
        WRITE(6,'(A)')
     &  ' It, Energy,  vecfn_ext, vecfn_int, step_ext, step_int: ' 
        WRITE(6,'(I4,1X,E22.15,2x,4(2X,E12.5))')
     &    IT, E, VCFNORM_EXT, VCFNORM_INT, XNORM, XNORM_CI
*. And normalize the internal part
        CNORM2 = INPROD(C_REF,C_REF,N_REF)
        FACTOR = 1.0D0/SQRT(CNORM2)
        CALL SCALVE(C_REF,FACTOR,N_REF)
*. Write new C_ref to file LUC - used by vector function 
        CALL ISTVC2(WORK(KLVCC2),0,1,N_REF)
        CALL REWINO(LUC)
        CALL WRSVCD(LUC,-1,VEC1,WORK(KLVCC2),
     &              C_REF,N_REF,N_REF,LUDIA,1)
*
      END DO
*     ^ End of loop over iterations
 1001 CONTINUE
*
*. Transformation of T to CAAB from orthonormal basis 
*. finalize  procedure
      CALL TRANS_CAAB_ORTN(WORK(KLVCC1),T_EXT,1,2,2,
     &         WORK(KLVCC2),2)
      CALL COPVEC(WORK(KLVCC1),T_EXT,NCAAB)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Info from T optimization ', IREFROOT
        WRITE(6,*) ' Updated amplitudes '
        CALL WRTMAT(T_EXT,1,NCAAB,1,NCAAB)
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Analysis of external amplitudes'
        CALL ANA_GENCC(T_EXT,1)
      END IF
*
      IF(IFIN_IT.EQ.1.OR.CONVERG) 
     &CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICC_CMP')
      RETURN
      END 
      SUBROUTINE COM_JAC_1COM(IREFSPC,IT2REFSPC,XJ,INCLUDE0)
*
*. Obtain in the orthonormal EI basis, 
*  the complete one-commutator approximation to Jacobian:
*  XJ(I,J) = <0!O+(I)[H,O(J)]|0>
*
* If INCLUDE0 = 1, then the zero-order state is included in Jacobian
*
* The spaces: IREFSPC : Space of !0>
*             IT2REFSPC : Space for T !0>
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cei.inc'
*. Output
      DIMENSION XJ(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'COMJC1')
*
      CALL MEMMAN(KLVCC1,NDIM_EI,'ADDL  ',2,'LVEC1 ')
      CALL MEMMAN(KLVCC2,NDIM_EI,'ADDL  ',2,'LVEC2 ')
      CALL MEMMAN(KLVCC3,NDIM_EI,'ADDL  ',2,'LVEC3 ')
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' --------------------------------- '
        WRITE(6,*) ' COM_JAC_1COM reporting to service '
        WRITE(6,*) ' --------------------------------- '
        WRITE(6,*)
      END IF
*
      IF(INCLUDE0.EQ.1) THEN
        NVAR = N_ZERO_EI
      ELSE
        NVAR = N_ZERO_EI - 1
      END IF
*
*. Part 1: <0| O(+)i H O j|0>
*

      ZERO = 0.0D0
      ONE = 1.0D0
      ONEM = -1.0D0
      WRITE(6,*) 'N_ZERO_EI = ', N_ZERO_EI

      DO J = 1, NVAR
       IF(NTEST.GE.10) WRITE(6,*) ' Part I, J =', J
       CALL SETVEC(WORK(KLVCC1),ZERO,N_ZERO_EI)
       WORK(KLVCC1-1+J) = 1.0D0
*. transform to CAAB basis
*. Dir in EI in T to Dir in CAAB in VCC1
        CALL TRANS_CAAB_ORTN(WORK(KLVCC2),WORK(KLVCC1),1,2,2,
     &         WORK(KLVCC3),2)
        CALL COPVEC(WORK(KLVCC2),WORK(KLVCC1),NDIM_EI)
* O(j) |0> on LUSC34
        ICSPC = IREFSPC
        ISSPC = IT2REFSPC
        CALL REWINO(LUSC34)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUSC34,
     &  WORK(KLVCC1),1)
*. Space of H T^I |0> may be reduced to IT2REFSPC
*. H O(j) |0>
        ICSPC = IT2REFSPC
        ISSPC = IT2REFSPC
        CALL REWINO(LUSC34)
        CALL REWINO(LUSC2)
        CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUSC34,LUSC2,0,0)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Output from MV7'
          CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
        END IF
*. The density <0|o+(CAAB) H O(j)|0>
        ZERO = 0.0D0
        ICSPC = IREFSPC
        ISSPC = IT2REFSPC
        CALL SETVEC(WORK(KLVCC1),ZERO,N_CC_AMP)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUSC2,
     &  WORK(KLVCC1),2)
*. And transform to obtain  <0|o+(i) H O(j)|0>
*. Vecfunc in CAAB in VCC1 to Vecfunc in ortn in VCC2
        CALL TRANS_CAAB_ORTN(WORK(KLVCC1),WORK(KLVCC2),1,1,2,
     &                       WORK(KLVCC3),1)
        CALL COPVEC(WORK(KLVCC2),XJ((J-1)*NVAR+1),NVAR)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The matrix <0|O+(i) H O(j)|0> '
        CALL WRTMAT(XJ,NVAR,NVAR,NVAR,NVAR)
      END IF
*
*. Part 2: -<0| O(+)iO j H|0>
*
*. H |0> on LUSC2
      ICSPC = IREFSPC
      ISSPC = IT2REFSPC
      CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUC,LUSC2,0,0)
      DO J = 1, NVAR
       IF(NTEST.GE.10) WRITE(6,*) ' Part II, J =', J
* O j H|0>
        CALL SETVEC(WORK(KLVCC1),ZERO,N_ZERO_EI)
        WORK(KLVCC1-1+J) = 1.0D0
*. transform to CAAB basis
*. Dir in ortn in VCC1 to Dir in CAAB in VCC2
        CALL TRANS_CAAB_ORTN(WORK(KLVCC2),WORK(KLVCC1),1,2,2,
     &         WORK(KLVCC3),2)
        CALL COPVEC(WORK(KLVCC2),WORK(KLVCC1),NDIM_EI)
*. O(j) H |0>
        ISSPC = IT2REFSPC
*. ISSPC kan reduceres til IREFSPC
        ICSPC = IT2REFSPC
        CALL REWINO(LUSC34)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC2,LUSC34,
     &  WORK(KLVCC1),1)
*. The density <0|o+(CAAB) O(j) H|0>
        ZERO = 0.0D0
        ICSPC = IREFSPC
*. ISSPC kan reduceres til IREFSPC
        ISSPC = IT2REFSPC
        CALL SETVEC(WORK(KLVCC1),ZERO,N_CC_AMP)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUC,LUSC34,
     &  WORK(KLVCC1),2)
*. And transform to obtain  <0|o+(i) O(j) H|0>
*. Vecfunc in CAAB in VCC1 to Vecfunc in ortn in VCC2
        CALL TRANS_CAAB_ORTN(WORK(KLVCC1),WORK(KLVCC2),1,1,2,
     &                       WORK(KLVCC3),1)
        CALL VECSUM(XJ((J-1)*NVAR+1),XJ((J-1)*NVAR+1),
     &              WORK(KLVCC2),ONE,ONEM,NVAR)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The matrix <0|O+(i) [H, O(j)]|0> '
        CALL WRTMAT(XJ,NVAR,NVAR,NVAR,NVAR)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'COMJC1')
      RETURN
      END 
      SUBROUTINE LUCIA_GICCI(
     &           ICTYP,EREF,EFINAL,CONVER,VNFINAL)
*
* Master routine for General Internal Contraction CI
* (alowing more than one external operators)
*
* LUCIA_IC is assumed to have been called to do the 
* preperatory work for working with internal contraction
*
* Jeppe Olsen, March 2010 for the Zurich tensor meeting
*
* Last modifications; Oct. 27, 2012; Jeppe Olsen; aligning..
*
C     INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 
     &INPROD, INPRDD
      LOGICAL CONVER,CONVER_INT,CONVER_EXT
C     INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cei.inc'
*. Transfer common block for communicating with H_EFF * vector routines
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
      INCLUDE 'gicci.inc'
*.Pointers for the external correlation operators
*.Number of parameters in the various spaces
*. Transfer block for communicating zero order energy to 
*. routien for performing H0-E0 * vector
      INCLUDE 'cshift.inc'
*
      CHARACTER*6 ICTYP
      EXTERNAL MTV_FUSK, STV_FUSK
      EXTERNAL H_S_EFF_ICCI_TV,H_S_EXT_ICCI_TV
      EXTERNAL H_S_EFF_GICCI_TV,H_S_EXT_GICCI_TV
      EXTERNAL HOME_SD_INV_T_ICCI
      EXTERNAL H0_EI_TV
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GICCI ')
      NTEST = 10
      WRITE(6,*)
      WRITE(6,*) ' ======================='
      WRITE(6,*) '  GICCI section entered '
      WRITE(6,*) ' ======================='
      WRITE(6,*)
*
      IF(IEI_VERSION.EQ.0) THEN
        I_DO_EI = 0
      ELSE
        I_DO_EI = 1
      END IF
*
      IF(I_DO_EI.EQ.1) THEN
       WRITE(6,*) ' EI approach in use'
      ELSE
       WRITE(6,*) ' Partial spin-adaptation in use'
      END IF
*. Notes
*
* In the initial version of this approach, a CI calculation typically
* preceeded the internal contraction calculations. In the GICCI approach
* T-operators are used for all correlation.
*
* The wavefunction is therefore: |0> = t_s(T(n)T(n-1)...T(1)|ref>
*                                          +T(n-1)...T(1)|ref>
*                                  .....
*                                          +|ref>)
*                                   
*. So space I is the initial HF or CAS space (|ref>)
*
*  
*. Transfer information on spaces
      NTEXC_GX  = NTEXC_G
      DO IEX = 1, NTEXC_G
       IPTCSPC_GX(IEX) = IPTCSPC_G(IEX)
       ITCSPC_GX(IEX) = ITCSPC_G(IEX)
      END DO
      
      IREFSPC = 1
      WRITE(6,*) ' Energy of reference state ', EREF
*
*. Information about the various CI spaces
*
      NCAAB_MX = 0
      NCAAB_TOT = 0
      NSPA_TOT = 0
      DO IEX = 1,  NTEXC_G
*. Prepare 
       CALL PREPARE_FOR_IEX(IEX)
*. Number of parameters with and without spinadaptation
       IF(I_DO_EI.EQ.0) THEN
         CALL NSPA_FOR_EXP_FUSK(NSPA,NCAAB)
         NCAAB_FOR_IEX(IEX) = NCAAB
         NSPA_FOR_IEX(IEX) = NSPA
         NCAAB_MX = MAX(NCAAB_MX,NCAAB)
         NSPA_MX = MAX(NSPA_MX,NSPA)
         NSPA_TOT = NSPA_TOT + NSPA
         NCAAB_TOT = NCAAB_TOT + NCAAB
       ELSE
*. Not updated pt
*. zero-particle operator is included in N_ZERO_EI
         NSPA = N_ZERO_EI
*. Note: NCAAB includes unitop
         NCAAB = NDIM_EI
       END IF
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Information about External operators '
        WRITE(6,*) ' ------------------------------------ '
        WRITE(6,*)
        WRITE(6,*) ' Operator   NCAAB    NSPA  '
        WRITE(6,*) '---------------------------'
        DO IEX = 1, NTEXC_G
         WRITE(6,'(3X,I3,3X,I8,3X,I8)')
     &   IEX, NCAAB_FOR_IEX(IEX),NSPA_FOR_IEX(IEX)
       END DO
      END IF
      I_IT_OR_DIR = 1
      IF(I_IT_OR_DIR.EQ.2) THEN
        WRITE(6,*) ' Explicit construction of all matrices'
      ELSE
        WRITE(6,*) ' Iterative solution of equations'
      END IF
      I_RELAX_INT = 0
      IF(I_RELAX_INT.EQ.1) THEN
        WRITE(6,*) ' Expansion of |ref> will be reoptimized '
      ELSE
        WRITE(6,*) ' Expansion of |ref> will be not be reoptimized '
      END IF
*. Space for CI behind the curtain 
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
      KVEC1P = KVEC1
      KVEC2P = KVEC2
* Allocate space and define pointers for two complete 
* external operators in the CAAB basis 
      CALL MEMMAN(KTEX_FOR_IEX(1),NCAAB_TOT+1,
     &  'ADDL  ',2,'T_EXT ')
      CALL MEMMAN(KTEXP_FOR_IEX(1),NCAAB_TOT+1,
     &  'ADDL  ',2,'T_EXT ')

      DO IEX = 2, NTEXC_G + 1
        KTEX_FOR_IEX(IEX) = KTEX_FOR_IEX(IEX-1)+NSPA_FOR_IEX(IEX-1)
        KTEXP_FOR_IEX(IEX) = KTEXP_FOR_IEX(IEX-1)+NSPA_FOR_IEX(IEX-1)
      END DO
*. And a vector that can hold the expansion for any given IEX_G
      CALL MEMMAN(KLTACT,NCAAB_MX,'ADDL  ',2,'TACT  ')
*
      N_REF = XISPSM(IREFSM,IREFSPC)
*. Initial  guess to T_EXT: Just the reference state: 
*  Zeroes in all T and coefficient one for the reference
      IF(IRESTRT_IC.EQ.0) THEN 
        ZERO = 0.0D0
        DO IEX = 1, NTEXC_G
          NSPA = NSPA_FOR_IEX(IEX)
          KLTEXT = KTEX_FOR_IEX(IEX) 
          CALL SETVEC(WORK(KLTEXT),ZERO,NSPA)
        END DO
*. And the coefficient for the reference state
        WORK(KTEX_FOR_IEX(NTEXC_G+1)) = 1.0D0
C            WRT_GICCI_VEC(KTEX)
C?      WRITE(6,*) ' TEX as set '
C?      CALL WRT_GICCI_VEC(KTEX_FOR_IEX)
C?      WRITE(6,*) ' KTEX_FOR_IEX(1), KTEX_FOR_IEX(NTEXC_G+1) =',
C?   &               KTEX_FOR_IEX(1), KTEX_FOR_IEX(NTEXC_G+1)
*. Store inital guess on unit 54
C     GIC_VEC_TO_DISC(KTEX,LEN_TEX,NTEX_G,IREW,LU)
        CALL GIC_VEC_TO_DISC(KTEX_FOR_IEX,NSPA_FOR_IEX,NTEXC_G,
     &                       1,LUSC54)
      END IF
*
      CONVER =.FALSE.
      CONVER_INT = .FALSE.
      CONVER_EXT = .FALSE.
      I12 = 2
      MAXIT_MACRO = MAXITM
*. Convergence will be defined as energy change
      I_ER_CONV = 1
*. There is no external converence threshold for residual
*. just use sqrt of energythreshold
      THRES_R = SQRT(THRES_E)
      DO IT_IE = 1, MAXIT_MACRO
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'COMP_M')
*
        IF(NTEST.GE.1) THEN
          WRITE(6,*)
          WRITE(6,*) ' ------------------------------------------'
          WRITE(6,*) ' Information from outer iteration ', IT_IE
          WRITE(6,*) ' ------------------------------------------'
          WRITE(6,*)
        END IF
        IDUM = 0
*. In iteration IT_IE, the internal operators upto and including T(IT_IE)
* are reoptimed
*
        ITMAX = MIN(IT_IE,NTEXC_G)
        WRITE(6,*) ' Number of operators to be optimized ', ITMAX
* and loop over the various T-operators to be optimized
        DO ITACT = 1,  ITMAX
          WRITE(6,*)
          WRITE(6,*) 
     &    ' Information about optimization of operator ', ITACT
          WRITE(6,*) 
     &    ' .........................................'
          WRITE(6,*)
*. Prepare for calculation in this space
          CALL PREPARE_FOR_IEX(ITACT)
*. Number of parameters with and without spinadaptation
          NCAAB = NCAAB_FOR_IEX(ITACT) 
          NSPA = NSPA_FOR_IEX(ITACT) 
*
          IF (I_IT_OR_DIR.EQ.2 ) THEN
*
* --------------------------------------------
*. Construct matrices explicit and diagonalize
* --------------------------------------------
*
            CALL ICCI_COMPLETE_MAT2(IREFSPC,ITREFSPC,I_SPIN_ADAPT,
     &           NROOT,WORK(KLTEXT),C_0,E_EXTOP)
            EFINAL = E_EXTOP
            CONVER_EXT = .TRUE.
            VNFINAL_EXT = 0.0D0
          ELSE 
*
*.------------------------------------------------
* Iterative methods used to solve GICCI equations
*.------------------------------------------------
*
*. Currently : no preconditioning and no elimination of singularities 
*              ( Yes, I am still an optimist ( or desperate ))
            NTESTL = 10
            MAXITL  = MAXIT
            MAXVECL = MXCIV
*. Jeppe Playing around
CD          IF(ITACT.EQ.1) THEN
CD           MAXITL = 2
CD           DO I = 1, 100
CD             WRITE(6,*) ' MAXITL = 2 for ITACT = 1 set by Jeppe !!'
CD           END DO
CD          END IF
*- End of Jeppe playing around
*. Allocate space for iterative solver 
            CALL MEMMAN(KL_EXTVEC1,NCAAB,'ADDL  ',2,'EXTVC1')
            CALL MEMMAN(KL_EXTVEC2,NCAAB,'ADDL  ',2,'EXTVC2')
            CALL MEMMAN(KL_EXTVEC3,NCAAB,'ADDL  ',2,'EXTVC3')
            CALL MEMMAN(KL_EXTVEC4,NCAAB,'ADDL  ',2,'EXTVC3')
*
            CALL MEMMAN(KL_RNRM,MAXITL*NROOT,'ADDL  ',2,'RNRM  ')
            CALL MEMMAN(KL_EIG ,MAXITL*NROOT,'ADDL  ',2,'EIG   ')
            CALL MEMMAN(KL_FINEIG,NROOT,'ADDL  ',2,'FINEIG')
*
            CALL MEMMAN(KL_APROJ,MAXVECL**2,'ADDL  ',2,'APROJ ')
            CALL MEMMAN(KL_SPROJ,MAXVECL**2,'ADDL  ',2,'SPROJ ')
            CALL MEMMAN(KL_AVEC ,MAXVECL**2,'ADDL  ',2,'AVEC  ')
            LLWORK = 5*MAXVECL**2 + 2*MAXVECL
            CALL MEMMAN(KL_WORK ,LLWORK   ,'ADDL  ',2,'WORK  ')
            CALL MEMMAN(KL_AVEC ,MAXVECL**2,'ADDL  ',2,'AVECP ')
            CALL MEMMAN(KL_AVECP,MAXVECL**2,'ADDL  ',2,'AVECP ')
*. Obtain diagonal of H and S
            I_DO_PRE_IN_EXT = 0
            IF(I_DO_PRE_IN_EXT.EQ.1) THEN
*. Generate non-trivial preconditioner
             IF(I_DO_EI.EQ.0) THEN
               CALL GET_HS_DIA(WORK(KL_EXTVEC3),WORK(KL_EXTVEC4),
     &              1,1,1,WORK(KL_EXTVEC1),WORK(KL_EXTVEC2),
     &                WORK(KVEC1),WORK(KVEC2),IREFSPC,ITREFSPC,
     &              IUNIOPX,NSPA,0,IDUM,IDUM)
             ELSE
*. EI approach
               CALL GET_DIAG_H0_EI(WORK(KL_EXTVEC3))
*. clean up
               I12 = 2
*. States are normalized, so
               ONE = 1.0D0
               CALL SETVEC(WORK(KL_EXTVEC4),ONE,NSPA)
             END IF
            ELSE
*. Generate trivial preconditioner
             ONE = 1.0D0
             CALL SETVEC(WORK(KL_EXTVEC3),ONE,NSPA)
             CALL SETVEC(WORK(KL_EXTVEC4),ONE,NSPA)
            END IF
*. And write diagonal to disc as single record files
            CALL VEC_TO_DISC(WORK(KL_EXTVEC3),NSPA,1,-1,LUSC53)
            CALL VEC_TO_DISC(WORK(KL_EXTVEC4),NSPA,1,-1,LUSC51)
*. (LUSC51 is not used)
            IF(IRESTRT_IC.EQ.1) THEN
*. Copy old CI coefficients for reference space to LUC
              CALL COPVCD(LUEXC,LUC,WORK(KVEC1),1,-1)
            END IF
*. Obtain current amplitudes for TACT and save in LUSC34
C                GIC_VEC_FROM_DISC(KTEX,LEN_TEX,NTEX_G,IREW,LU)
C?          WRITE(6,*) ' Before GIC_VEC_FROM... '
            CALL GIC_VEC_FROM_DISC(KTEX_FOR_IEX,NSPA_FOR_IEX,NTEXC_G,
     &                             1,LUSC54)
C?          WRITE(6,*) ' After GIC_VEC_FROM... T read in '
C?          CALL WRT_GICCI_VEC(KTEX_FOR_IEX)
*
            C0 = WORK(KTEX_FOR_IEX(1)-1+NSPA_TOT+1) 
C?          WRITE(6,*) ' coefficient of ref before MINGENEIG', 
C?   &                 C0
            CALL COPVEC(WORK(KTEX_FOR_IEX(ITACT)),WORK(KLTACT),NSPA-1)
*. Coefficient for constant part of expansion (independent of T(IACT))
            WORK(KLTACT-1+NSPA) = 1.0D0
C?          WRITE(6,*) ' KLTACT, KLTACT-1+NSPA+1=',
C?   &                   KLTACT, KLTACT-1+NSPA+1
C?          WRITE(6,*) ' WORK(KLTACT) as defined'
C?          CALL WRTMAT(WORK(KLTACT),1,NSPA,1,NSPA)
*. and save amplitudes
            CALL VEC_TO_DISC(WORK(KLTACT),NSPA,1,-1,LUSC34)
            DO IMAC = 1, 1
* LUSC53 is LU_DIAH, LUSC51 is LU_DIAS, LUSC36 is LUC where 
* eigenvector is stored
*. 2 implies that advanced preconditioner is called 
*- Save reference energy for use with diagonal preconditioner
              EREFX = EREF
*
C?            WRITE(6,*) ' I_DO_EI = ', I_DO_EI
              I12 = 2
              IF(I_DO_EI.EQ.0) THEN
                IPREC_FORM = 1
                SHIFT = 0.0D0
                CALL MINGENEIG(H_S_EXT_GICCI_TV,HOME_SD_INV_T_ICCI,
     &               IPREC_FORM,THRES_E,THRES_R,I_ER_CONV,
     &               WORK(KL_EXTVEC1),WORK(KL_EXTVEC2),WORK(KL_EXTVEC3),
     &               LUSC34, LUSC37,
     &               WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &               NSPA,LUSC38,LUSC39,LUSC40,LUSC53,LUSC51,LUSC52,
     &               NROOT,MAXVECL,NROOT,WORK(KL_APROJ),
     &               WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &               NTESTL,SHIFT,WORK(KL_AVECP),I_DO_PRE_IN_EXT,
     &               CONVER_EXT,E_EXTOP,VNFINAL_EXT)
              ELSE
                IPREC_FORM = 2
                CALL MINGENEIG(H_S_EXT_GICCI_TV,H0_EI_TV,
     &               IPREC_FORM,THRES_E,THRES_R,I_ER_CONV,
     &               WORK(KL_EXTVEC1),WORK(KL_EXTVEC2),WORK(KL_EXTVEC3),
     &               LUSC34, LUSC37,
     &               WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &               NSPA,LUSC38,LUSC39,LUSC40,LUSC53,LUSC51,LUSC52,
     &               NROOT,MAXVECL,NROOT,WORK(KL_APROJ),
     &               WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &               NTESTL,SHIFT,WORK(KL_AVECP),I_DO_PRE_IN_EXT,
     &               CONVER_EXT,E_EXTOP,VNFINAL_EXT)
              END IF
             EFINAL = E_EXTOP
            END DO
*           ^ End of loop over reset eigenvalue problem
*. Update T-coefficients on LU54
            CALL GIC_VEC_FROM_DISC(KTEX_FOR_IEX,NSPA_FOR_IEX,
     &           NTEXC_G,1,LUSC54)
            CALL VEC_FROM_DISC(WORK(KLTACT),NSPA,1,-1,LUSC34)
C                UPDATE_GICCI_VEC(KTEX,I_EX_ACT,TACTVEC,ISCALE)
            CALL UPDATE_GICCI_VEC(KTEX_FOR_IEX,ITACT,WORK(KLTACT),1)
*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Updated T-coefficients to be written '
              CALL WRT_GICCI_VEC(KTEX_FOR_IEX)
C                  WRT_GICCI_VEC(KTEX)
            END IF
*
            CALL GIC_VEC_TO_DISC(KTEX_FOR_IEX,NSPA_FOR_IEX,NTEXC_G,
     &           1,LUSC54)
*. Test: construct wave function
            CALL GET_GICCI_0(KTEX_FOR_IEX,LUSC38,LUC,LUSC39,LUSC40)
            XNORM0 = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUSC38,LUSC38,
     &                      1,-1)
C?          WRITE(6,*) ' Square norm of |0> after MINGENEIG', XNORM0
C   GET_GICCI_0(KTEXG,LUOUT,LUC,LUSC2,LUSC3)
            C_0 = WORK(KTEX_FOR_IEX(NTEXC_G+1))
*. And the current T(ACT)
            CALL COPVEC(WORK(KTEX_FOR_IEX(ITACT)),WORK(KLTACT),NSPA)
            IF(I_DO_EI.EQ.0) THEN
              CALL PREPARE_FOR_IEX(ITACT)
              CALL REF_CCV_CAAB_SP(WORK(KL_EXTVEC1),WORK(KLTACT),
     &             WORK(KL_EXTVEC3),2) 
            ELSE
              CALL TRANS_CAAB_ORTN(WORK(KL_EXTVEC1),WORK(KLTACT),1,2,2,
     &             WORK(KL_EXTVEC3),2) 
            END IF
            T_CAAB_NORM =
     &      SQRT(INPROD(WORK(KL_EXTVEC1),WORK(KL_EXTVEC1),NCAAB))
            WRITE(6,*) ' Norm of T in CAAB basis after MINGENEIG',
     &      T_CAAB_NORM
*
            IF(NTEST.GE.10) THEN
              WRITE(6,*) ' coefficient of zero-order state ', C_0
              WRITE(6,*) 
     &        ' Analysis of external amplitudes in CAAB basis'
              CALL ANA_GENCC(WORK(KL_EXTVEC1),1)
            END IF
          END IF
*         ^ End of switch direct/iterative approach for T_EXT
         END DO
*        ^ End of loop over Operators to be optimized in this outer
*        iteration
    
        VNFINAL_INT = 0.0D0
        IF(I_RELAX_INT.EQ.1) THEN
* ============================================================
*. Relax coefficients of internal/reference/zero-order state 
* ============================================================
*
        IF(NTEST.GE.0) THEN
           WRITE(6,*)
           WRITE(6,*) ' Optimization of internal correlation part'
           WRITE(6,*) ' .........................................'
           WRITE(6,*)
        END IF
           CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
           KVEC1P = KVEC1
           KVEC2P = KVEC2
*
           IF(I_IT_OR_DIR.EQ.2) THEN
*
*. Construct complete matrices and diagonalize
*
*. Space for H and S in zero-order space 
             CALL MEMMAN(KLH_REF,N_REF**2,'ADDL  ',2,'H_REF  ')
             CALL MEMMAN(KLS_REF,N_REF**2,'ADDL  ',2,'S_REF  ')
             CALL MEMMAN(KLC_REF,N_REF   ,'ADDL  ',2,'C_REF  ')
             CALL MEMMAN(KLI_REF,N_REF   ,'ADDL  ',1,'I_REF  ')
*
             CALL ICCI_RELAX_REFCOEFS_COM(WORK(KLTEXT),NSPA,
     &            WORK(KLH_REF),
     &            WORK(KLS_REF),N_REF,WORK(KVEC1),WORK(KVEC2),1,
     &            IREFSPC,ITREFSPC,C_0,ECORE,WORK(KLC_REF),NROOT,
     &            NCAAB,E_INTOP)
             CONVER_INT =.TRUE.
             VNFINAL_INT = 0.0D0
             EFINAL = E_INTOP
*. transfer new reference vector to DISC
             CALL ISTVC2(WORK(KLI_REF),0,1,N_REF)
C  WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
             CALL REWINO(LUC)
             CALL WRSVCD(LUC,-1,WORK(KVEC1),WORK(KLI_REF),
     &            WORK(KLC_REF),N_REF,N_REF,LUDIA,1)
           ELSE 
*. Use iterative methods to reoptimize reference coefficients
             MAXITL = MAXIT
             MAXVEC = MXCIV
*
             CALL MEMMAN(KL_REFVEC1,N_REF,'ADDL  ',2,'REFVC1')
             CALL MEMMAN(KL_REFVEC2,N_REF,'ADDL  ',2,'REFVC2')
             CALL MEMMAN(KL_REFVEC3,N_REF,'ADDL  ',2,'REFVC3')
*
             CALL MEMMAN(KL_RNRM,MAXIT*NROOT,'ADDL  ',2,'RNRM  ')
             CALL MEMMAN(KL_EIG ,MAXIT*NROOT,'ADDL  ',2,'EIG   ')
             CALL MEMMAN(KL_FINEIG,NROOT,'ADDL  ',2,'FINEIG')
*
             CALL MEMMAN(KL_APROJ,MAXVEC**2,'ADDL  ',2,'APROJ ')
             CALL MEMMAN(KL_SPROJ,MAXVEC**2,'ADDL  ',2,'SPROJ ')
             CALL MEMMAN(KL_AVEC ,MAXVEC**2,'ADDL  ',2,'AVEC  ')
             LLWORK = 5*MAXVEC**2 + 2*MAXVEC
             CALL MEMMAN(KL_WORK ,LLWORK   ,'ADDL  ',2,'WORK  ')
             CALL MEMMAN(KL_AVEC ,MAXVEC**2,'ADDL  ',2,'AVECP ')
             CALL MEMMAN(KL_AVECP,MAXVEC**2,'ADDL  ',2,'AVECP ')
*
* Well, there is pt a conflict between the form of files 
* in mingeneig and in the general CI programs
*. In MINGENEIG all vectors are single record files, whereas
*  the vectors are multirecord files in the general LUCIA 
* world. Reformatting is therefore required..
*. LUC is LUC
*. LUSC36 is LUDIA
*. LUSC51 is LUDIAS
*
*. Reform LUC to single record file
             CALL REWINO(LUC)
             CALL FRMDSCN(WORK(KL_REFVEC1),-1,-1,LUC)
             CALL REWINO(LUC)
             CALL VEC_TO_DISC(WORK(KL_REFVEC1),N_REF,1,-1,LUC)
*. Reform LUDIA to single record file on LUSC36
             CALL REWINO(LUDIA)
             CALL FRMDSCN(WORK(KL_REFVEC1),-1,-1,LUDIA)
             CALL VEC_TO_DISC(WORK(KL_REFVEC1),N_REF,1,-1,LUSC36)
*. Write diagonal of S as unit mat as single vector file
             ONE = 1.0D0
             CALL SETVEC(WORK(KL_REFVEC1),ONE,N_REF)
             CALL VEC_TO_DISC(WORK(KL_REFVEC1),N_REF,1,-1,LUSC51)
*. (LUSC51 is not used)
*
* As preconditioners, the standard CI diagonal and the 
* unit diagonal will be used for H and S, respectively.
* This is fine if the T operator is not too large...
*
*. Prepare transfer common block for communicating with
*. matrix-vector routines
C            C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX
             C_0X = C_0
             KLTOPX = KLTEXT
             NREFX = N_REF
             IREFSPCX = IREFSPC
             ITREFSPCX = ITREFSPC
             NCAABX = N_CC_AMP
             NSPAX = NSPA
*. Unitoperator in SPA order ... Please check ..
             IUNIOPX = NSPA
*
             NTESTL = 10
             CALL MINGENEIG( H_S_EFF_ICCI_TV,HOME_SD_INV_T_ICCI,1,
     &            THRES_E,THRES_R,I_ER_CONV,
     &            WORK(KL_REFVEC1),WORK(KL_REFVEC2),WORK(KL_REFVEC3),
     &            LUC, LUSC37,
     &            WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &            N_REF,LUSC38,LUSC39,LUSC40,LUSC36,LUSC51,LUSC52,
     &            NROOT,MXCIV,NROOT,WORK(KL_APROJ),
     &            WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &            NTESTL,ECORE,WORK(KL_AVECP),1,
     &            CONVER_INT,E_INTOP,VNFINAL_INT)
                  E_FINAL = E_INTOP
C                 MINGENEIG(MTV,STV,
C    &                VEC1,VEC2,VEC3,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
C    &                NVAR,
C    &                LU3,LU4,LU5,LUDIAM,LUDIAS,LUS,NROOT,MAXVEC,
C    &                NINVEC,
C    &                APROJ,AVEC,SPROJ,WORK,IPRT,EIGSHF,AVECP,I_DO_PRECOND)
*
*. Read new eigenvector from LUC
             CALL REWINO(LUC)
             CALL FRMDSCN(WORK(KL_REFVEC1),-1,-1,LUC)
* The eigenvector is normalized with respect to the <i!T+P P T|j>
*. metric, normalize with standard unit metrix
             XNORM = INPROD(WORK(KL_REFVEC1),WORK(KL_REFVEC1),N_REF)
             FACTOR = 1.0D0/SQRT(XNORM)
             CALL SCALVE(WORK(KL_REFVEC1),FACTOR,N_REF)
*. And write to disc in a form suitable for the other parts of LUCIA
             CALL ISTVC2(WORK(KL_REFVEC2),0,1,N_REF)
             CALL REWINO(LUC)
             CALL REWINO(LUDIA)
             CALL WRSVCD(LUC,-1,WORK(KVEC1P),WORK(KL_REFVEC2),
     &                   WORK(KL_REFVEC1),N_REF,N_REF,LUDIA,1)
             IF(NTEST.GE.100) THEN
               WRITE(6,*) ' New reference coefficients '
               CALL WRTVCD(WORK(KVEC1P),LUC,1,-1)
             END IF
           END IF 
*.         ^ End of switch direct/iterative methods for reference relaxation 
        END IF
*.      ^ End of reference coefs should be relaxed
        CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'COMP_M')
        IF(CONVER_INT.AND.CONVER_EXT.AND.
     &     ABS(E_INTOP-E_EXTOP).LE.THRES_E) CONVER = .TRUE.
        IF(CONVER) GOTO 1001
      END DO
 1001 CONTINUE
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' coefficient of zero-order state ', C_0
        WRITE(6,*) 
     &  ' Analysis of final external amplitudes in CAAB basis'
        CALL ANA_GENCC(WORK(KLTEXT),1)
      END IF
*
      VNFINAL = VNFINAL_INT + VNFINAL_EXT
      WRITE(6,*) ' VNFINAL_INT, VNFINAL_EXT =', 
     &             VNFINAL_INT,VNFINAL_EXT
*.    ^ End of loop over Internal/external correlation iterations
*. Print the final coefs ..
C?    CALL VEC_FROM_DISC(WORK(KL_EXTVEC1),NSPA,1,-1,LUSC54)
C?    WRITE(6,*) ' Final list of IC-coefficients '
C?    CALL WRTMAT(WORK(KL_EXTVEC1),NSPA,1,NSPA,1)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GICCI ')
      RETURN
      END 
      SUBROUTINE LUCIA_IC(IREFSPC,ITREFSPC,ICTYP,EREF,I_DO_CUMULANTS,
     &                    EFINAL,CONVER,VNFINAL)
*
*
* Master routine for internally contracted CI calculations,
* Fall 02 version 
*
* Allowing CAS as well as RAS and MRSDCI references -I hope
*
* Jeppe Olsen, September 02
*
* Last modification; Oct. 21, 2012; Jeppe Olsen; error in defining NSPOBEX_TPE corrected
*
* Also used for generating cumulant matrices 
* 
      INCLUDE 'wrkspc.inc'
      REAL*8 
     &INPROD
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
*. Transfer common block for communicating with H_EFF * vector routines
      COMMON/COM_H_S_EFF_ICCI_TV/
     &       C_0X,KLTOPX,NREFX,IREFSPCX,ITREFSPCX,NCAABX,
     &       IUNIOPX,NSPAX,IPROJSPCX
*. A bit of local scratch
      DIMENSION ICASCR(MXPNGAS)
      CHARACTER*6 ICTYP
      LOGICAL CONVER
*
      EXTERNAL MTV_FUSK, STV_FUSK
      EXTERNAL H_S_EFF_ICCI_TV,H_S_EXT_ICCI_TV
      EXTERNAL HOME_SD_INV_T_ICCI
*. Test of new transformer
C?    CALL tranma_lm_test
C?    STOP ' Enforced stop after  tranma_lm_test '
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ICCI  ')
*. I will play with spinadaptation in this routine so 
*. It is probably not working of I_SPIN_ADAPT = 0 is used !!!
      IF(I_DO_CUMULANTS.EQ.0) THEN 
         I_SPIN_ADAPT = 1
      ELSE
         I_SPIN_ADAPT = 0
      END IF
*
      NTEST = 10
      IF(NTEST.GE.5) THEN
        IF(I_DO_CUMULANTS.EQ.0) THEN 
         WRITE(6,*)
         WRITE(6,*) ' Internal contracted section entered '
         WRITE(6,*) ' ==================================== '
         WRITE(6,*)
         WRITE(6,*) '  Symmetri of reference vector ' , IREFSM 
         WRITE(6,*) '  Space of Reference vector ', IREFSPC
         WRITE(6,*) '  Space of Internal contracted vector ', ITREFSPC
         WRITE(6,*)
         WRITE(6,*) ' Parameters defining internal contraction '
         WRITE(6,*) '       Min excitation rank  ', ICEXC_RANK_MIN
         WRITE(6,*) '       Max excitation rank  ', ICEXC_RANK_MAX
         WRITE(6,'(A,A)') ' Form of calculation  ', ICTYP
         IF(ICEXC_INT.EQ.1) THEN
           WRITE(6,*) 
     &   ' Internal (ina->ina, sec->sec) excitations allowed'
         ELSE
           WRITE(6,*) 
     &   ' Internal (ina->ina, sec->sec) excitations not allowed'
         END IF
         WRITE(6,*) 
     &   '  Largest number of vectors in iterative supspace ', MXCIV
         WRITE(6,*) 
     &   '  Largest initial number of vectors in iterative supspace ',
     &     MXVC_I
         IF(IRESTRT_IC.EQ.1) THEN
           WRITE(6,*) ' Restarted calculation : '
           WRITE(6,*) '      IC coefficients  read from LUSC54'
           WRITE(6,*) '      CI for reference read from LUSC54 '
         END IF
        ELSE
         WRITE(6,*) ' Cumulants will be calculated upto order ',
     &              ICUMULA
        END IF
*
      END IF
*
      IDUM = 0
*. Divide orbital spaces into inactive, active, secondary using 
*. space 1
      CALL CC_AC_SPACES(1,IREFTYP)
C     CC_AC_SPACES(ISPC,IREFTYP)
*
*. Orbital excitations to work in reference state
*
*. Number of orbital excitations
C     IC_ORBOP(IWAY,NIC_ORBOP,IC_ORBOP,MX_OP_RANK,MN_OP_RANK,
C    &                   IONLY_EXCOP)
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      IF(ICEXC_INT.EQ.1) THEN
         IONLY_EXCOP = 0
      ELSE
         IONLY_EXCOP = 1
      END IF
      IF(I_DO_CUMULANTS.EQ.0) THEN
*. Normal internal contracted run - unit operator included
       IADD_UNI = 1
       CALL GEN_IC_ORBOP(1,NOBEX_TP,IDUMMY,
     &               ICEXC_RANK_MAX,ICEXC_RANK_MIN,
     &               IONLY_EXCOP,IREFSPC,ITREFSPC,IADD_UNI,
     &               IPRSTR)
*. and the orbital excitations
       CALL MEMMAN(KOBEX_TP,2*NGAS*NOBEX_TP,'ADDL ',2,'IC_OBX')
       KLOBEX = KOBEX_TP
       CALL GEN_IC_ORBOP(2,NOBEX_TP,WORK(KOBEX_TP),
     &               ICEXC_RANK_MAX,ICEXC_RANK_MIN,
     &               IONLY_EXCOP,IREFSPC,ITREFSPC,IADD_UNI,
     &               IPRSTR)
       NOBEX_TPE = NOBEX_TP+1
      ELSE
*. Cumulant calculation 
C     GEN_IC_IN_ORBSPC(IWAY,NIC_ORBOP,IC_ORBOP,MX_OP_NUM,
C    &                               IORBSPC)
*. Identify the active space ( determined in CC_AC_SPACES)
       NACT_SPC = 0
       DO IGAS = 1, NGAS
         IF(IHPVGAS(IGAS).EQ.3) THEN
          IACTSPC = IGAS
          NACT_SPC = NACT_SPC + 1
         END IF
       END DO
       IF(NACT_SPC.GT.1) THEN
         WRITE(6,*) ' More than one active space in cumulant expansion'
         WRITE(6,*) ' Cumulant code currently assumes one active space '
         STOP ' More than one active space for cumulant calculation '
       END IF
       IF(NACT_SPC.EQ.0) THEN
         WRITE(6,*) ' No active space '
         WRITE(6,*) ' Cumulant matrices only calculated in active space'
         WRITE(6,*) ' I am therefore finished and stop '
         STOP ' Zero active space for cumulant calculation '
       END IF
       CALL GEN_IC_IN_ORBSPC(1,NOBEX_TP,IDUMMY,ICUMULA,IACTSPC)
*. and the orbital excitations
       CALL MEMMAN(KOBEX_TP,2*NGAS*NOBEX_TP,'ADDL ',2,'IC_OBX')
       KLOBEX = KOBEX_TP
       CALL GEN_IC_IN_ORBSPC(2,NOBEX_TP,WORK(KLOBEX),ICUMULA,IACTSPC)
       NOBEX_TPE = NOBEX_TP+1
      END IF
*
      IF(I_SPIN_ADAPT.EQ.1) THEN
*
*. Excitation operators will be spin adapted
*
        DO JOBEX_TP = 1, NOBEX_TP
C?        WRITE(6,*) ' Constructing CA confs for JOBEX_TP = ', JOBEX_TP
*. Integer arrays for creation and annihilation part 
          CALL ICOPVE2(WORK(KOBEX_TP),1+(JOBEX_TP-1)*2*NGAS,2*NGAS,
     &                  ICASCR)
          NOP_C = IELSUM(ICASCR,NGAS)
          NOP_A = IELSUM(ICASCR(1+NGAS),NGAS)
          NOP_CA = NOP_C + NOP_A
          CALL GET_CA_CONF_FOR_ORBEX(ICASCR,ICASCR(1+NGAS),
     &         NCOC_FSM(1,JOBEX_TP),NAOC_FSM(1,JOBEX_TP),
     &         IBCOC_FSM(1,JOBEX_TP),IBAOC_FSM(1,JOBEX_TP),
     &         KCOC(JOBEX_TP),KAOC(JOBEX_TP),
     &         KZC(JOBEX_TP),KZA(JOBEX_TP),
     &         KCREO(JOBEX_TP),KAREO(JOBEX_TP))
C?        WRITE(6,*) ' NCOC_FSM and NAOC_FSM after GET_CA ... '
C?        CALL IWRTMA(NCOC_FSM,1,NSMST,1,NSMST)
C?        CALL IWRTMA(NAOC_FSM,1,NSMST,1,NSMST)
         
*. Offsets in CA block for given symmetry of creation occ
C IOFF_SYMBLK_MAT(NSMST,NA,NB,ITOTSM,IOFF,IRESTRICT
          CALL IOFF_SYMBLK_MAT(NSMST,NCOC_FSM(1,JOBEX_TP),
     &         NAOC_FSM(1,JOBEX_TP),1,IBCAOC_FSM(1,JOBEX_TP),0)
C                           NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
          NCAOC(JOBEX_TP) = NDIM_1EL_MAT(1,NCOC_FSM(1,JOBEX_TP),
     &                      NAOC_FSM(1,JOBEX_TP),NSMST,0)
*. And the actual configurations 
          CALL MEMMAN(KCAOC(JOBEX_TP),NOP_CA*NCAOC(JOBEX_TP),'ADDL  ',
     &                2,'CA_OC ')
C     GET_CONF_FOR_ORBEX(NCOC_FSM,NAOC_FSM,ICOC,IAOC,
C    &           NOP_C,NOP_A, IBCOC_FSM,IBAOC_FSM,NSMST,IOPSM,
C    &           ICAOC)
          CALL GET_CONF_FOR_ORBEX(
     &         NCOC_FSM(1,JOBEX_TP),NAOC_FSM(1,JOBEX_TP),
     &         WORK(KCOC(JOBEX_TP)),WORK(KAOC(JOBEX_TP)),
     &         NOP_C, NOP_A,
     &         IBCOC_FSM(1,JOBEX_TP),IBAOC_FSM(1,JOBEX_TP),
     &         NSMST,1,WORK(KCAOC(JOBEX_TP)) )
        END DO
      END IF
*. Number of creation and annihilation operators per op
      CALL MEMMAN(KLCOBEX_TP,NOBEX_TPE,'ADDL ',1,'COBEX ')
      CALL MEMMAN(KLAOBEX_TP,NOBEX_TPE,'ADDL ',1,'AOBEX ')
      CALL GET_NCA_FOR_ORBOP(NOBEX_TP,WORK(KOBEX_TP),
     &     WORK(KLCOBEX_TP),WORK(KLAOBEX_TP),NGAS)
*. Number of spinorbital excitations
      IZERO = 0
      MXSPOX = 0
      IACT_SPC = 0
      IAAEXC_TYP = 3
      IREFSPCX = 0
      MSCOMB_CC = 0
      CALL OBEX_TO_SPOBEX(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TPE,NGAS,
     &     NOBPT,0,IZERO ,IAAEXC_TYP,IACT_SPC,IPRCC,IDUMMY,
     &     MXSPOX,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPCX)
*CJO, Oct, 21, 2012, start
C     NSPOBEX_TPE = NSPOBEX_TP + 1
      NSPOBEX_TP = NSPOBEX_TPE - 1
*CJO, Oct, 21, 2012, end
*. And the actual spinorbital excitations
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. First SOX of given OX ( including zero operator )
      CALL MEMMAN(KIBSOX_FOR_OX,NOBEX_TPE,'ADDL  ',1,'IBSOXF')
*. Number of SOX's for given OX
      CALL MEMMAN(KNSOX_FOR_OX,NOBEX_TPE,'ADDL  ',1,'IBSOXF')
*. SOX for given OX
      CALL MEMMAN(KISOX_FOR_OX,NSPOBEX_TPE,'ADDL  ',1,'IBSOXF')
*
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TPE,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,IPRCC,
     &     WORK(KLSOX_TO_OX),MXSPOX,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPCX)
      IF(I_DO_CUMULANTS.EQ.0) THEN
*
* A bit of info on prototype-excitations
*
*. Number of prototype-excitations
C      NPROTO_CA(NOBEX_TP,IOBEX_TP,NGAS)
       NPROTO_CA_EX = NPROTO_CA(NOBEX_TP,WORK(KOBEX_TP),NGAS)
*. And  info on the prototypes 
       CALL MEMMAN(K_MX_DLB_C,NOBEX_TP,'ADDL  ',2,'MXDB_C')
       CALL MEMMAN(K_MX_DLB_A,NOBEX_TP,'ADDL  ',2,'MXDB_A')
       CALL MEMMAN(K_IB_PROTO,NOBEX_TP,'ADDL  ',2,'IB_PRO')
       CALL MEMMAN(K_NCOMP_FOR_PROTO,NPROTO_CA_EX,'ADDL  ',2,
     &             'NCO_PR')
       CALL INFO2_FOR_PROTO_CA(
     &       NOBEX_TP,WORK(KOBEX_TP),WORK(KISOX_FOR_OX),
     &       WORK(KNSOX_FOR_OX),WORK(KIBSOX_FOR_OX),
     &       WORK(KLSOBEX),NGAS,
     &       WORK(K_IB_PROTO),WORK(K_MX_DLB_C),WORK(K_MX_DLB_A),
     &       WORK(K_NCOMP_FOR_PROTO),NPROTO_CA_EX)
C      INFO2_FOR_PROTO_CA(
C    &            NOBEX_TP,IOBEX_TP,ISOX_FOR_OX,NSOX_FOR_OX,IBSOX_FOR_OX,
C    &            ISPOBEX_TP,NGAS,
C    &            IB_PROTO_CA, MX_DBL_C_CA, MX_DBL_A_CA,
C    &            NCOMP_FOR_PROTO_CA,NPROTO_CA)
      END IF
*
* Dimension and offsets of IC operators
*
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TPE,'ADDL  ',1,'SPOBAC')
      CALL MEMMAN(KLSPOBEX_FRZ,NSPOBEX_TPE,'ADDL  ',1,'SPOBAC')
*. ALl spinorbital excitations are initially active
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TPE)
*. And none are frozen
      IZERO = 0
      CALL ISETVC(WORK(KLSPOBEX_FRZ),IZERO,NSPOBEX_TPE)
*
      ITOP_SM = 1
C?    WRITE(6,*) ' IREFSPC before IDIM.. ', IREFSPC
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TPE,ITOP_SM,
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MSCOMB_CC,MX_TBLK_AS,
     &     WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &     NTCONF,IPRCC)
      N_CC_AMP = LEN_T_VEC
      WRITE(6,*) ' Number of IC parameters ', N_CC_AMP
      WRITE(6,*) ' Dimension of the various types '
      CALL IWRTMA(WORK(KLLSOBEX),1,NSPOBEX_TP,1,NSPOBEX_TP)
*
      MX_ST_TSOSO_MX = MX_ST_TSOSO
      MX_ST_TSOSO_BLK_MX = MX_ST_TSOSO_BLK
      MX_TBLK_MX = MX_TBLK
      MX_TBLK_AS_MX = MX_TBLK_AS
      LEN_T_VEC_MX =  LEN_T_VEC
*. Some more scratch etc
*. Alpha- and beta-excitations constituting the spinorbital excitations
*. Number 
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      CALL MEMMAN(KLAOBEX,LENA,'ADDL  ',2,'IAOBEX')
      CALL MEMMAN(KLBOBEX,LENB,'ADDL  ',2,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     0,NAOBEX_TP,NBOBEX_TP,WORK(KLAOBEX),WORK(KLBOBEX))
*. Max dimensions of CCOP !KSTR> = !ISTR> maps
*. For alpha excitations
      IATP = 1
      IOCTPA = IBSPGPFTP(IATP)
      NOCTPA = NSPGPFTP(IATP)
      CALL LEN_GENOP_STR_MAP(
     &     NAOBEX_TP,WORK(KLAOBEX),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IBTP = 2
      IOCTPB = IBSPGPFTP(IBTP)
      NOCTPB = NSPGPFTP(IBTP)
      CALL LEN_GENOP_STR_MAP(
     &     NBOBEX_TP,WORK(KLBOBEX),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(NTEST.GE.5) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
*. Space for old fashioned CI behind the curtain
*. For calculations without EI VEC1, VEC2, VEC3 have not been defined, do this.
* There must be inserted a check to see if EI calculation is called or move 
* allocation
        CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
        KVEC1P = KVEC1
        KVEC2P = KVEC2
      IF(I_DO_CUMULANTS.EQ.1) THEN 
*. 1 : construct standard density matrices 
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CUMULA')
*. Space for old fashioned CI behind the curtain
*. For calculations without EI VEC1, VEC2, VEC3 have not been defined, do this.
* There must be inserted a check
        CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
        KVEC1P = KVEC1
        KVEC2P = KVEC2
*. and space for the reduced density matrices/cumulants
        WRITE(6,*) ' IREFSPC = ', IREFSPC
        ICSPC = IREFSPC
        ISSPC = IREFSPC
        CALL MEMMAN(KLCUMULANTS,N_CC_AMP,'ADDL  ',2,'CUMULA')
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLCUMULANTS),ZERO,N_CC_AMP)
*. And an independent copy of the reference vector
        CALL COPVCD(LUC,LUSC1,WORK(KVEC1),1,-1)
*. Calculate reduced density matrices 
        CALL SIGDEN_CC(WORK(KVEC1),WORK(KVEC2),LUC,LUSC1,
     &                 WORK(KLCUMULANTS),2)
*. And reform to cumulant expansion
        WRITE(6,*) ' RDM => Cumulant reformer will be called '
        CALL REFORM_RDM_TO_CUMULANTS(WORK(KLCUMULANTS),WORK(KLSOBEX),
     &       WORK(KLLSOBEX))
C     REFORM_RDM_TO_CUMULANTS(CUMULANTS,ISPOBEX_TP,LSOBEX_TP)
        CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CUMULA')
      END IF
*
      IF(I_SPIN_ADAPT.EQ.1) THEN
*. Generate maps CAAB excitations to CA .ie. the spinorbital
*. excitations belonging to the various orbital excitations
         DO JOBEX = 1, NOBEX_TP
*. Number of spinorbital excitations belonging to this orbital 
*, excitation type
           NSOX = IFRMR(WORK(KNSOX_FOR_OX),1,JOBEX)
           IBSOX = IFRMR(WORK(KIBSOX_FOR_OX),1,JOBEX)
           NCAAB = IGATSUM(WORK(KLLSOBEX),WORK(KISOX_FOR_OX),
     &                     IBSOX,NSOX)
           WRITE(6,*) ' JOBEX, NSOX, IBSOX, NCAAB = ',
     &                  JOBEX, NSOX, IBSOX, NCAAB 
           NCA   = NCAOC(JOBEX)
C                          IGATSUM(IVEC,IGAT,IOFF,NELMNT)
           NOP_C = IFRMR(WORK(KLCOBEX_TP),1,JOBEX)
           NOP_CA = 2*NOP_C
           WRITE(6,*) ' NOP_CA = ', NOP_CA
*
*. Allocate space 
* KICAAB_FOR_CA_OP : The CA CB AA AB operators for each CAAB
           LEN = NOP_CA*NCAAB
           CALL MEMMAN(KICAAB_FOR_CA_OP(JOBEX),LEN,'ADDL  ',2,'ICAABO') 
* KICAAB_FOR_CA_NUM : A number for each CAAB
           LEN = NCAAB
           CALL MEMMAN(KICAAB_FOR_CA_NUM(JOBEX),LEN,'ADDL  ',2,'ICAABN')
*.KLCAAB_FOR_CA : Length of CA CB AA AB for each CAAB
           LEN = 4*NCAAB
           CALL MEMMAN(KLCAAB_FOR_CA(JOBEX),LEN,'ADDL  ',2,'LCAAB ')
*.KNCAAB_FOR_CA : A length for each CA
           LEN = NCA
           CALL MEMMAN(KNCAAB_FOR_CA(JOBEX),LEN,'ADDL  ',2,'NCAAB ')
*.KIBCAAB_FOR_CA : First CAAB for given CA
           LEN = NCA
           CALL MEMMAN(KIBCAAB_FOR_CA(JOBEX),LEN,'ADDL  ',2,'IBCAAB')
*
           CALL CAAB_TO_CA_OC(1,WORK(KLSOBEX),WORK(KLOBEX),JOBEX,
     &          WORK(KISOX_FOR_OX),WORK(KIBSOX_FOR_OX),
     &          WORK(KNSOX_FOR_OX),WORK(KLIBSOBEX),
     &          MX_ST_TSOSO_BLK_MX,NOP_CA,
     &          WORK(KZC(JOBEX)),WORK(KZA(JOBEX)),WORK(KCREO(JOBEX)),
     &          WORK(KAREO(JOBEX)),WORK(KCAOC(JOBEX)),
     &          IBCAOC_FSM(1,JOBEX),NCOC_FSM(1,JOBEX),
     &          WORK(KIBCAAB_FOR_CA(JOBEX)),
     &          WORK(KICAAB_FOR_CA_OP(JOBEX)),
     &          WORK(KICAAB_FOR_CA_NUM(JOBEX)),
     &          WORK(KLCAAB_FOR_CA(JOBEX)),
     &          WORK(KNCAAB_FOR_CA(JOBEX)),NCA,NCAAB,
     &          WORK(K_NCOMP_FOR_PROTO) )
     
        END DO
      IF(NTEST.GE.100) CALL WRITE_CAAB_CONFM
*. Construct reorder array, CONF => CAAB order
      CALL MEMMAN(KLREORDER_CAAB,N_CC_AMP,'ADDL  ',1,'RECAAB')
      CALL GEN_REORDER_CAABM(WORK(KLREORDER_CAAB))
C     GEN_REORDER_CAABM(ICAAB_REO)
*
* Construct matrices for Spinadaptation 
*
      CALL PROTO_SPIN_MAT
*. Number of SPA and CAAB excitations per orbital excitation type
      CALL DIM_FOR_OBEXTP
C          DIM_FOR_OBEXTP
      END IF
*     ^ End if spinadaptation 
*
* Call routines for explicit construction of matrices 
* and complete diagonalizations
*
      I_ANALYZE_SING = 0
      IF(I_ANALYZE_SING.EQ.1) THEN 
*. Check single excitation like operators for singularities
        CALL SXLIKE_SING(IREFSPC,ITREFSPC,NSXLIKE,I_SPIN_ADAPT)
C?      WRITE(6,*) ' Enforced stop after SXLIKE_SING '
C?      STOP       ' Enforced stop after SXLIKE_SING '
*. Still checking singularities : Find singularities in SX and a+p a+h a ah ah
* space
*
        WRITE(6,*) 
     &  ' singularities in space spanned by sx,a+pa+ha h a h,a+pa+papah'
        WRITE(6,*) ' =================================================='
        ICASCR(1) = 1
        ICASCR(2) = 2
        ICASCR(3) = 4
        CALL SING_IN_OCCLS(IREFSPC,ITREFSPC,ICASCR,3)
*
        WRITE(6,*) 
     &  ' singularities in space spanned by a+pa+ha h a h,a+pa+papah'
        WRITE(6,*) ' =================================================='
        ICASCR(1) = 2
        ICASCR(2) = 4
        CALL SING_IN_OCCLS(IREFSPC,ITREFSPC,ICASCR,2)
*
        WRITE(6,*) ' singularities in space spanned by sx, a+pa+hahah'
        WRITE(6,*) ' ================================================'
        ICASCR(1) = 1
        ICASCR(2) = 2
        CALL SING_IN_OCCLS(IREFSPC,ITREFSPC,ICASCR,2)
*
        WRITE(6,*) ' singularities in space spanned by SX, a+pa+papah '
        WRITE(6,*) ' ================================================='
        ICASCR(1) = 1
        ICASCR(2) = 4
        CALL SING_IN_OCCLS(IREFSPC,ITREFSPC,ICASCR,2)
*
        WRITE(6,*) ' singularities in space spanned by a+pa+ha h a h '
        WRITE(6,*) ' ================================================'
        ICASCR(1) = 2
        CALL SING_IN_OCCLS(IREFSPC,ITREFSPC,ICASCR,1)
*
        WRITE(6,*) ' singularities in space spanned by  a+pa+pa p a h '
        WRITE(6,*) ' ================================================='
        ICASCR(1) = 4
        CALL SING_IN_OCCLS(IREFSPC,ITREFSPC,ICASCR,1)
*
        WRITE(6,*)  ' Enforced stop After checking singularities '
        STOP ' Enforced stop After checking singularities '
      END IF
*
*. Analyze singularities in SX-space by diagonaling
*. the various 2-e spin-densities
C     CALL GET_SING_IN_SX_SPACE(IREFSPC)
C     GET_SING_IN_SX_SPACE
      IF(ICTYP(1:4).EQ.'ICCI') THEN
*
*                    ==============================
*                    Internal contracted CI section 
*                    ==============================
*
* Solve Internal contracted CI problem 
        CALL LUCIA_ICCI(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                 EFINAL,CONVER,VNFINAL)
*
      ELSE IF(ICTYP(1:4).EQ.'ICPT') THEN
*
*                    ==========================================
*                    Internal contracted Perturbation expansion 
*                    ==========================================
*
        CALL LUCIA_ICPT(IREFSPC,ITREFSPC,ICTYP,EREF,
     &                 EFINAL,CONVER,VNFINAL)
*
      ELSE IF(ICTYP(1:4).EQ.'ICCC') THEN
* Internal contracted coupled cluster 
*
*                    ======================================
*                    Internal contracted Coupled Cluster 
*                    =======================================
*
        CALL LUCIA_ICCC(IREFSPC,ITREFSPC,ICTYP,EREF,EFINAL,
     &                  CONVER,VNFINAL)
      END IF
*
*.
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ICCI  ')
*
      RETURN
      END
      SUBROUTINE GET_TEX_INFO(
     &           IREFSPC,ITREFSPC,
     &           MX_ST_TSOSO, MX_ST_TSOSO_BLK, MX_TBLK,  MX_TBLK_AS)
*
* Generate all information about orbital and spin-orbital excitations
* Information is stored in scalars in CTCC
* 
*
*. Jeppe Olsen, collecting and restructuring for GICCI etc.
*. March 27, 2010
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc' 
      INCLUDE 'cprnt.inc'
*. Controlling print flag: IPRSTR
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      IF(ICEXC_INT.EQ.1) THEN
         IONLY_EXCOP = 0
      ELSE
         IONLY_EXCOP = 1
      END IF
*
      IADD_UNI = 1
      IDUM = 0
      CALL GEN_IC_ORBOP(1,NOBEX_TP,IDUM,
     &              ICEXC_RANK_MAX,ICEXC_RANK_MIN,
     &              IONLY_EXCOP,IREFSPC,ITREFSPC,IADD_UNI,
     &               IPRSTR)
*. and the orbital excitations
      CALL MEMMAN(KOBEX_TP,2*NGAS*NOBEX_TP,'ADDL ',2,'IC_OBX')
      KLOBEX = KOBEX_TP
      CALL GEN_IC_ORBOP(2,NOBEX_TP,WORK(KOBEX_TP),
     &              ICEXC_RANK_MAX,ICEXC_RANK_MIN,
     &              IONLY_EXCOP,IREFSPC,ITREFSPC,IADD_UNI,
     &               IPRSTR)
      NOBEX_TPE = NOBEX_TP+1
*. Number of creation and annihilation operators per op
      CALL MEMMAN(KLCOBEX_TP,NOBEX_TPE,'ADDL ',1,'COBEX ')
      CALL MEMMAN(KLAOBEX_TP,NOBEX_TPE,'ADDL ',1,'AOBEX ')
      CALL GET_NCA_FOR_ORBOP(NOBEX_TP,WORK(KOBEX_TP),
     &     WORK(KLCOBEX_TP),WORK(KLAOBEX_TP),NGAS)
*. Number of spinorbital excitations
      IZERO = 0
      MXSPOX = 0
      IACT_SPC = 0
      IAAEXC_TYP = 3
      IREFSPCX = 0
      MSCOMB_CC = 0
      CALL OBEX_TO_SPOBEX(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUM,NSPOBEX_TPE,NGAS,
     &     NOBPT,0,IZERO ,IAAEXC_TYP,IACT_SPC,IPRSTR,IDUM,
     &     MXSPOX,IDUM,
     &     IDUM,IDUM,NAEL,NBEL,IREFSPCX)
      NSPOBEX_TP = NSPOBEX_TPE 
*. And the actual spinorbital excitations
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. First SOX of given OX ( including zero operator )
      CALL MEMMAN(KIBSOX_FOR_OX,NOBEX_TPE,'ADDL  ',1,'IBSOXF')
*. Number of SOX's for given OX
      CALL MEMMAN(KNSOX_FOR_OX,NOBEX_TPE,'ADDL  ',1,'IBSOXF')
*. SOX for given OX
      CALL MEMMAN(KISOX_FOR_OX,NSPOBEX_TPE,'ADDL  ',1,'IBSOXF')
*. KLSOBEX,KIBSOX_FOR_OX,KNSOX_FOR_OX,KISOX_FOR_OX,
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TPE,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,IPRSTR,
     &     WORK(KLSOX_TO_OX),MXSPOX,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPCX)
C?    WRITE(6,*) 'ISOX_FOR_OX after OBEX_TO.....'
C?    CALL IWRTMA(WORK(KISOX_FOR_OX),1,NSPOBEX_TP,1,NSPOBEX_TP)
* Dimension and offsets of IC operators
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TPE,'ADDL  ',1,'SPOBAC')
      CALL MEMMAN(KLSPOBEX_FRZ,NSPOBEX_TPE,'ADDL  ',1,'SPOBAC')
*. KLLSOBEX, KLIBSOBEX, KLSPOBEX_AC, KLSPOBEX_FRZ
*. ALl spinorbital excitations are initially active
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TPE)
*. And none are frozen
      IZERO = 0
      CALL ISETVC(WORK(KLSPOBEX_FRZ),IZERO,NSPOBEX_TPE)
*
      ITOP_SM = 1
*. Dimension of blocks of CC and of total expansion
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ITOP_SM,
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MSCOMB_CC,MX_TBLK_AS,
     &     WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &     NTCONF,IPRSTR)
      N_CC_AMP = LEN_T_VEC
      WRITE(6,*) ' Number of IC parameters ', N_CC_AMP
      IF(IPRSTR.GE.5) THEN
        WRITE(6,*) ' Dimension of the various types '
        CALL IWRTMA(WORK(KLLSOBEX),1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*  MX_ST_TSOSO, MX_ST_TSOSO_BLK, MX_TBLK,  MX_TBLK_AS,
      MX_ST_TSOSO_MX = MX_ST_TSOSO
      MX_ST_TSOSO_BLK_MX = MX_ST_TSOSO_BLK
      MX_TBLK_MX = MX_TBLK
      MX_TBLK_AS_MX = MX_TBLK_AS
      LEN_T_VEC_MX =  LEN_T_VEC
*. Alpha- and beta-excitations constituting the spinorbital excitations
*. Number 
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      CALL MEMMAN(KLAOBEX,LENA,'ADDL  ',2,'IAOBEX')
      CALL MEMMAN(KLBOBEX,LENB,'ADDL  ',2,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     0,NAOBEX_TP,NBOBEX_TP,WORK(KLAOBEX),WORK(KLBOBEX))
*. Max dimensions of CCOP !KSTR> = !ISTR> maps
*. For alpha excitations
      IATP = 1
      IOCTPA = IBSPGPFTP(IATP)
      NOCTPA = NSPGPFTP(IATP)
      CALL LEN_GENOP_STR_MAP(
     &     NAOBEX_TP,WORK(KLAOBEX),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IBTP = 2
      IOCTPB = IBSPGPFTP(IBTP)
      NOCTPB = NSPGPFTP(IBTP)
      CALL LEN_GENOP_STR_MAP(
     &     NBOBEX_TP,WORK(KLBOBEX),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(NTEST.GE.5) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
      RETURN
      END
      SUBROUTINE TRANSFER_T_OFFSETS(I_FT_GLOBAL,IEX_G)
*
*. Transfer offsets for  T-operators between specific and 
*. general arrays for offsets and lengths
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
*
*. Jeppe Olsen, March 2009
*
*. Last modification; Oct. 27, 2012; Jeppe Olsen; NSPOBEX_TPE added
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ---------------------------'
        WRITE(6,*) ' Entering TRANSFER_T_OFFSETS'
        WRITE(6,*) ' ---------------------------'
        WRITE(6,*)
        WRITE(6,*) ' I_FT_GLOBAL, IEX_G =', I_FT_GLOBAL,IEX_G
      END IF
*
      IF(I_FT_GLOBAL.EQ.2) THEN
*. Write information to permanent arrays
        NOBEX_TP_G(IEX_G) = NOBEX_TP
C?      WRITE(6,*) ' NOBEX_TP_G, NOBEX_TP = ',
C?   &               NOBEX_TP_G(IEX_G), NOBEX_TP
        KOBEX_TP_G(IEX_G) = KOBEX_TP
C?      WRITE(6,*) ' KOBEX_TP_G(IEX_G), KOBEX_TP (a) ',
C?   &               KOBEX_TP_G(IEX_G), KOBEX_TP
        KLCOBEX_TP_G(IEX_G) = KLCOBEX_TP
        KLAOBEX_TP_G(IEX_G) = KLAOBEX_TP
        NSPOBEX_TP_G(IEX_G) = NSPOBEX_TP
        KLSOBEX_G(IEX_G) =    KLSOBEX
        KIBSOX_FOR_OX_G(IEX_G) = KIBSOX_FOR_OX
        KNSOX_FOR_OX_G(IEX_G) = KNSOX_FOR_OX
        KISOX_FOR_OX_G(IEX_G) = KISOX_FOR_OX
        KLSOX_TO_OX_G(IEX_G) = KLSOX_TO_OX
C?      WRITE(6,*) 'KISOX_FOR_OX_G, KISOX_FOR_OX(a)'
C?      WRITE(6,*) KISOX_FOR_OX_G(IEX_G),KISOX_FOR_OX
        KLLSOBEX_G(IEX_G) = KLLSOBEX
        KLIBSOBEX_G(IEX_G) = KLIBSOBEX
        KLSPOBEX_AC_G(IEX_G) = KLSPOBEX_AC
        KLSPOBEX_FRZ_G(IEX_G) = KLSPOBEX_FRZ
        N_CC_AMP_G(IEX_G) = N_CC_AMP
        NAOBEX_TP_G(IEX_G) = NAOBEX_TP
        NBOBEX_TP_G(IEX_G) = NBOBEX_TP
        KLAOBEX_G(IEX_G) = KLAOBEX
        KLBOBEX_G(IEX_G) = KLBOBEX 
      ELSE
        NOBEX_TP = NOBEX_TP_G(IEX_G) 
C?      WRITE(6,*) ' NOBEX_TP_G, NOBEX_TP = ',
C?   &               NOBEX_TP_G(IEX_G), NOBEX_TP
        KOBEX_TP = KOBEX_TP_G(IEX_G) 
C?      WRITE(6,*) ' KOBEX_TP_G(IEX_G), KOBEX_TP (b) ',
C?   &               KOBEX_TP_G(IEX_G), KOBEX_TP
        KLCOBEX_TP = KLCOBEX_TP_G(IEX_G) 
        KLAOBEX_TP = KLAOBEX_TP_G(IEX_G) 
        NSPOBEX_TP = NSPOBEX_TP_G(IEX_G) 
        NSPOBEX_TPE = NSPOBEX_TP
        KLSOBEX = KLSOBEX_G(IEX_G) 
        KIBSOX_FOR_OX = KIBSOX_FOR_OX_G(IEX_G) 
        KNSOX_FOR_OX = KNSOX_FOR_OX_G(IEX_G) 
        KISOX_FOR_OX = KISOX_FOR_OX_G(IEX_G)
        KLSOX_TO_OX = KLSOX_TO_OX_G(IEX_G)
C?      WRITE(6,*) 'KISOX_FOR_OX_G, KISOX_FOR_OX(b)'
C?      WRITE(6,*) KISOX_FOR_OX_G(IEX_G),KISOX_FOR_OX
        KLLSOBEX = KLLSOBEX_G(IEX_G)
        KLIBSOBEX = KLIBSOBEX_G(IEX_G)
        KLSPOBEX_AC = KLSPOBEX_AC_G(IEX_G)
        KLSPOBEX_FRZ = KLSPOBEX_FRZ_G(IEX_G)
        N_CC_AMP = N_CC_AMP_G(IEX_G)
        NAOBEX_TP = NAOBEX_TP_G(IEX_G)
        NBOBEX_TP = NBOBEX_TP_G(IEX_G)
        KLAOBEX = KLAOBEX_G(IEX_G)
        KLBOBEX = KLBOBEX_G(IEX_G)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_SP_INFO
*
*. Information in partial spin-adaptation of excitation operators
*  Information is stored in specific arrays in corbex, ctcc. glbbas
*
*. Jeppe Olsen, march 27, 2010
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
*
      DIMENSION ICASCR(MXPNGAS)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ----------------------------'
        WRITE(6,*) ' Information from GET_SP_INFO'
        WRITE(6,*) ' ----------------------------'
        WRITE(6,*)
      END IF
*
      DO JOBEX_TP = 1, NOBEX_TP
        IF(NTEST.GE.100) 
     &  WRITE(6,*) ' Constructing CA confs for JOBEX_TP = ', JOBEX_TP
*. Integer arrays for creation and annihilation part 
        CALL ICOPVE2(WORK(KOBEX_TP),1+(JOBEX_TP-1)*2*NGAS,2*NGAS,
     &                  ICASCR)
        NOP_C = IELSUM(ICASCR,NGAS)
        NOP_A = IELSUM(ICASCR(1+NGAS),NGAS)
        NOP_CA = NOP_C + NOP_A
        CALL GET_CA_CONF_FOR_ORBEX(ICASCR,ICASCR(1+NGAS),
     &       NCOC_FSM(1,JOBEX_TP),NAOC_FSM(1,JOBEX_TP),
     &       IBCOC_FSM(1,JOBEX_TP),IBAOC_FSM(1,JOBEX_TP),
     &       KCOC(JOBEX_TP),KAOC(JOBEX_TP),
     &       KZC(JOBEX_TP),KZA(JOBEX_TP),
     &       KCREO(JOBEX_TP),KAREO(JOBEX_TP))
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' NCOC_FSM and NAOC_FSM after GET_CA ... '
          CALL IWRTMA(NCOC_FSM,1,NSMST,1,NSMST)
          CALL IWRTMA(NAOC_FSM,1,NSMST,1,NSMST)
        END IF
*. Offsets in CA block for given symmetry of creation occ
C IOFF_SYMBLK_MAT(NSMST,NA,NB,ITOTSM,IOFF,IRESTRICT
        CALL IOFF_SYMBLK_MAT(NSMST,NCOC_FSM(1,JOBEX_TP),
     &       NAOC_FSM(1,JOBEX_TP),1,IBCAOC_FSM(1,JOBEX_TP),0)
C                           NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
        NCAOC(JOBEX_TP) = NDIM_1EL_MAT(1,NCOC_FSM(1,JOBEX_TP),
     &                    NAOC_FSM(1,JOBEX_TP),NSMST,0)
*. And the actual configurations 
        CALL MEMMAN(KCAOC(JOBEX_TP),NOP_CA*NCAOC(JOBEX_TP),'ADDL  ',
     &              2,'CA_OC ')
C     GET_CONF_FOR_ORBEX(NCOC_FSM,NAOC_FSM,ICOC,IAOC,
C    &           NOP_C,NOP_A, IBCOC_FSM,IBAOC_FSM,NSMST,IOPSM,
C    &           ICAOC)
        CALL GET_CONF_FOR_ORBEX(
     &       NCOC_FSM(1,JOBEX_TP),NAOC_FSM(1,JOBEX_TP),
     &       WORK(KCOC(JOBEX_TP)),WORK(KAOC(JOBEX_TP)),
     &       NOP_C, NOP_A,
     &       IBCOC_FSM(1,JOBEX_TP),IBAOC_FSM(1,JOBEX_TP),
     &       NSMST,1,WORK(KCAOC(JOBEX_TP)) )
      END DO
*
* A bit of info on prototype-excitations
*
*. Number of prototype-excitations
C     NPROTO_CA(NOBEX_TP,IOBEX_TP,NGAS)
      NPROTO_CA_EX = NPROTO_CA(NOBEX_TP,WORK(KOBEX_TP),NGAS)
*. And  info on the prototypes 
      CALL MEMMAN(K_MX_DLB_C,NOBEX_TP,'ADDL  ',2,'MXDB_C')
      CALL MEMMAN(K_MX_DLB_A,NOBEX_TP,'ADDL  ',2,'MXDB_A')
      CALL MEMMAN(K_IB_PROTO,NOBEX_TP,'ADDL  ',2,'IB_PRO')
      CALL MEMMAN(K_NCOMP_FOR_PROTO,NPROTO_CA_EX,'ADDL  ',2,
     &            'NCO_PR')
      CALL INFO2_FOR_PROTO_CA(
     &      NOBEX_TP,WORK(KOBEX_TP),WORK(KISOX_FOR_OX),
     &      WORK(KNSOX_FOR_OX),WORK(KIBSOX_FOR_OX),
     &      WORK(KLSOBEX),NGAS,
     &      WORK(K_IB_PROTO),WORK(K_MX_DLB_C),WORK(K_MX_DLB_A),
     &      WORK(K_NCOMP_FOR_PROTO),NPROTO_CA_EX)
C?    WRITE(6,*) ' After INFO2'
C     INFO2_FOR_PROTO_CA(
C    &           NOBEX_TP,IOBEX_TP,ISOX_FOR_OX,NSOX_FOR_OX,IBSOX_FOR_OX,
C    &           ISPOBEX_TP,NGAS,
C    &           IB_PROTO_CA, MX_DBL_C_CA, MX_DBL_A_CA,
C    &           NCOMP_FOR_PROTO_CA,NPROTO_CA)
*
*
*. Generate maps CAAB excitations to CA .ie. the spinorbital
*. excitations belonging to the various orbital excitations
*
      DO JOBEX = 1, NOBEX_TP
*. Number of spinorbital excitations belonging to this orbital 
*, excitation type
        NSOX = IFRMR(WORK(KNSOX_FOR_OX),1,JOBEX)
        IBSOX = IFRMR(WORK(KIBSOX_FOR_OX),1,JOBEX)
        NCAAB = IGATSUM(WORK(KLLSOBEX),WORK(KISOX_FOR_OX),
     &                  IBSOX,NSOX)
        NCA   = NCAOC(JOBEX)
C                       IGATSUM(IVEC,IGAT,IOFF,NELMNT)
        NOP_C = IFRMR(WORK(KLCOBEX_TP),1,JOBEX)
        NOP_CA = 2*NOP_C
*
*. Allocate space 
* KICAAB_FOR_CA_OP : The CA CB AA AB operators for each CAAB
        LEN = NOP_CA*NCAAB
        CALL MEMMAN(KICAAB_FOR_CA_OP(JOBEX),LEN,'ADDL  ',2,'ICAABO') 
* KICAAB_FOR_CA_NUM : A number for each CAAB
        LEN = NCAAB
        CALL MEMMAN(KICAAB_FOR_CA_NUM(JOBEX),LEN,'ADDL  ',2,'ICAABN')
*.KLCAAB_FOR_CA : Length of CA CB AA AB for each CAAB
        LEN = 4*NCAAB
        CALL MEMMAN(KLCAAB_FOR_CA(JOBEX),LEN,'ADDL  ',2,'LCAAB ')
*.KNCAAB_FOR_CA : A length for each CA
        LEN = NCA
        CALL MEMMAN(KNCAAB_FOR_CA(JOBEX),LEN,'ADDL  ',2,'NCAAB ')
*.KIBCAAB_FOR_CA : First CAAB for given CA
        LEN = NCA
        CALL MEMMAN(KIBCAAB_FOR_CA(JOBEX),LEN,'ADDL  ',2,'IBCAAB')
*
        CALL CAAB_TO_CA_OC(1,WORK(KLSOBEX),WORK(KOBEX_TP),JOBEX,
     &       WORK(KISOX_FOR_OX),WORK(KIBSOX_FOR_OX),
     &       WORK(KNSOX_FOR_OX),WORK(KLIBSOBEX),
     &       MX_ST_TSOSO_BLK_MX,NOP_CA,
     &       WORK(KZC(JOBEX)),WORK(KZA(JOBEX)),WORK(KCREO(JOBEX)),
     &       WORK(KAREO(JOBEX)),WORK(KCAOC(JOBEX)),
     &       IBCAOC_FSM(1,JOBEX),NCOC_FSM(1,JOBEX),
     &       WORK(KIBCAAB_FOR_CA(JOBEX)),
     &       WORK(KICAAB_FOR_CA_OP(JOBEX)),
     &       WORK(KICAAB_FOR_CA_NUM(JOBEX)),
     &       WORK(KLCAAB_FOR_CA(JOBEX)),
     &       WORK(KNCAAB_FOR_CA(JOBEX)),NCA,NCAAB,
     &       WORK(K_NCOMP_FOR_PROTO) )
C?    WRITE(6,*) ' After CAAB_TO'
     
      END DO
*
      IF(NTEST.GE.100) CALL WRITE_CAAB_CONFM
*. Construct reorder array, CONF => CAAB order
C?    WRITE(6,*) ' N_CC_AMP before GEN_REORDER... ', N_CC_AMP
      CALL MEMMAN(KLREORDER_CAAB,N_CC_AMP,'ADDL  ',1,'RECAAB')
      CALL GEN_REORDER_CAABM(WORK(KLREORDER_CAAB))
C     GEN_REORDER_CAABM(ICAAB_REO)
*. Number of SPA and CAAB excitations per orbital excitation type
      CALL DIM_FOR_OBEXTP
C          DIM_FOR_OBEXTP
*
      RETURN
      END
      SUBROUTINE TRANSFER_SPIN_OFFSETS(I_FT_GLOBAL,IEX_G)
*
* Transfer from (I_FT_GLOBAL=1) or to (I_FT_GLOBAL=2) 
* global arrays from specific/actual arrays
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'corbex.inc'
      INCLUDE 'glbbas.inc'
      COMMON/PROTO_SP_MAT/NSPA_FOP(6),NCAAB_FOP(6),IB_FOP(6),XTRA(100),
     &                    NSPA_FOP_G(6,MXPCYC),NCAAB_FOP_G(6,MXPCYC),
     &                    IB_FOP_G(6,MXPCYC)
*
      IF(I_FT_GLOBAL.EQ.2) THEN
        DO IOBEX_TP = 1, NOBEX_TP
         CALL ICOPVE(NCOC_FSM(1,IOBEX_TP),NCOC_FSM_G(1,IOBEX_TP,IEX_G),
     &        NSMST)
         CALL ICOPVE(NAOC_FSM(1,IOBEX_TP),NAOC_FSM_G(1,IOBEX_TP,IEX_G),
     &        NSMST)
         CALL ICOPVE(IBCOC_FSM(1,IOBEX_TP),
     &               IBCOC_FSM_G(1,IOBEX_TP,IEX_G),NSMST)
         CALL ICOPVE(IBAOC_FSM(1,IOBEX_TP),
     &               IBAOC_FSM_G(1,IOBEX_TP,IEX_G),NSMST)
         KCOC_G(IOBEX_TP,IEX_G) = KCOC(IOBEX_TP)
         KAOC_G(IOBEX_TP,IEX_G) = KAOC(IOBEX_TP)
         KZC_G(IOBEX_TP,IEX_G) = KZC(IOBEX_TP)
         KZA_G(IOBEX_TP,IEX_G) = KZA(IOBEX_TP)
         KCREO_G(IOBEX_TP,IEX_G) = KCREO(IOBEX_TP)
         KAREO_G(IOBEX_TP,IEX_G) = KAREO(IOBEX_TP)
         CALL ICOPVE(IBCAOC_FSM(1,IOBEX_TP),
     &               IBCAOC_FSM_G(1,IOBEX_TP,IEX_G),NSMST)
         NCAOC_G(IOBEX_TP,IEX_G) = NCAOC(IOBEX_TP)
         KCAOC_G(IOBEX_TP,IEX_G) = KCAOC(IOBEX_TP)
*
         KICAAB_FOR_CA_NUM_G(IOBEX_TP,IEX_G) = 
     &   KICAAB_FOR_CA_NUM(IOBEX_TP)
         KICAAB_FOR_CA_OP_G(IOBEX_TP,IEX_G) =
     &   KICAAB_FOR_CA_OP(IOBEX_TP)
         KLCAAB_FOR_CA_G(IOBEX_TP,IEX_G) = KLCAAB_FOR_CA(IOBEX_TP)
         KNCAAB_FOR_CA_G(IOBEX_TP,IEX_G) = KNCAAB_FOR_CA(IOBEX_TP)
         KIBCAAB_FOR_CA_G(IOBEX_TP,IEX_G) = KIBCAAB_FOR_CA(IOBEX_TP)
         NSPA_FOR_OCCLS_G(IOBEX_TP,IEX_G) = NSPA_FOR_OCCLS(IOBEX_TP)
         NCAAB_FOR_OCCLS_G(IOBEX_TP,IEX_G) = NCAAB_FOR_OCCLS(IOBEX_TP)
         IBSPA_FOR_OCCLS_G(IOBEX_TP,IEX_G) = IBSPA_FOR_OCCLS(IEX_G)
        END DO
*
        K_NCOMP_FOR_PROTO_G(IEX_G) = K_NCOMP_FOR_PROTO
        K_MX_DLB_C_G(IEX_G) = K_MX_DLB_C
        K_MX_DLB_A_G(IEX_G) = K_MX_DLB_A
        K_IB_PROTO_G(IEX_G) = K_IB_PROTO
        KLREORDER_CAAB_G(IEX_G) = KLREORDER_CAAB
*
        MAXNDET = 6
        CALL ICOPVE(NSPA_FOP,NSPA_FOP_G(1,IEX_G),MAXNDET)
        CALL ICOPVE(NCAAB_FOP,NCAAB_FOP_G(1,IEX_G),MAXNDET)
        CALL ICOPVE(IB_FOP,IB_FOP_G(1,IEX_G),MAXNDET)
       ELSE
*.  From general to specific/actual 
        DO IOBEX_TP = 1, NOBEX_TP
         CALL ICOPVE(NCOC_FSM_G(1,IOBEX_TP,IEX_G),NCOC_FSM(1,IOBEX_TP),
     &        NSMST)
         CALL ICOPVE(NAOC_FSM_G(1,IOBEX_TP,IEX_G),NAOC_FSM(1,IOBEX_TP),
     &        NSMST)
         CALL ICOPVE(IBCOC_FSM_G(1,IOBEX_TP,IEX_G),
     &        IBCOC_FSM(1,IOBEX_TP),NSMST)
         CALL ICOPVE(IBAOC_FSM_G(1,IOBEX_TP,IEX_G),
     &        IBAOC_FSM(1,IOBEX_TP),NSMST)
         KCOC(IOBEX_TP) = KCOC_G(IOBEX_TP,IEX_G) 
         KAOC(IOBEX_TP) = KAOC_G(IOBEX_TP,IEX_G) 
         KZC(IOBEX_TP) = KZC_G(IOBEX_TP,IEX_G) 
         KZA(IOBEX_TP) = KZA_G(IOBEX_TP,IEX_G) 
         KCREO(IOBEX_TP) = KCREO_G(IOBEX_TP,IEX_G) 
         KAREO(IOBEX_TP) = KAREO_G(IOBEX_TP,IEX_G) 
         CALL ICOPVE(IBCAOC_FSM_G(1,IOBEX_TP,IEX_G),
     &        IBCAOC_FSM(1,IOBEX_TP),NSMST)
         NCAOC(IOBEX_TP) = NCAOC_G(IOBEX_TP,IEX_G) 
         KCAOC(IOBEX_TP) = KCAOC_G(IOBEX_TP,IEX_G) 
*
         KICAAB_FOR_CA_NUM(IOBEX_TP) =
     &   KICAAB_FOR_CA_NUM_G(IOBEX_TP,IEX_G)   
         KICAAB_FOR_CA_OP =
     &   KICAAB_FOR_CA_OP_G(IOBEX_TP,IEX_G) 
         KLCAAB_FOR_CA(IOBEX_TP) = KLCAAB_FOR_CA_G(IOBEX_TP,IEX_G) 
         KNCAAB_FOR_CA(IOBEX_TP) = KNCAAB_FOR_CA_G(IOBEX_TP,IEX_G) 
         KIBCAAB_FOR_CA(IOBEX_TP) = KIBCAAB_FOR_CA_G(IOBEX_TP,IEX_G) 
         NSPA_FOR_OCCLS(IOBEX_TP) = NSPA_FOR_OCCLS_G(IOBEX_TP,IEX_G) 
         NCAAB_FOR_OCCLS(IOBEX_TP) = NCAAB_FOR_OCCLS_G(IOBEX_TP,IEX_G) 
         IBSPA_FOR_OCCLS(IOBEX_TP) = IBSPA_FOR_OCCLS_G(IOBEX_TP,IEX_G) 
        END DO
*
        K_NCOMP_FOR_PROTO = K_NCOMP_FOR_PROTO_G(IEX_G) 
        K_MX_DLB_C = K_MX_DLB_C_G(IEX_G) 
        K_MX_DLB_A = K_MX_DLB_A_G(IEX_G) 
        K_IB_PROTO =  K_IB_PROTO_G(IEX_G) 
        KLREORDER_CAAB = KLREORDER_CAAB_G(IEX_G) 
*
        MAXNDET = 6
        CALL ICOPVE(NSPA_FOP_G(1,IEX_G),NSPA_FOP,MAXNDET)
        CALL ICOPVE(NCAAB_FOP_G(1,IEX_G),NCAAB_FOP,MAXNDET)
        CALL ICOPVE(IB_FOP_G(1,IEX_G),IB_FOP,MAXNDET)
      END IF
*
      RETURN
      END
      SUBROUTINE PREPARE_FOR_IEX(IEX)
*
*. Prepare setup for calculation with general excitation operator IEX
*
*. Jeppe Olsen, on the way to Zurick, march 2010
*
      INCLUDE 'implicit.inc'
*
      I_FT_GLOBAL = 1
      CALL TRANSFER_T_OFFSETS(I_FT_GLOBAL,IEX)
      CALL TRANSFER_SPIN_OFFSETS(I_FT_GLOBAL,IEX)
*
      RETURN
      END
      SUBROUTINE GIC_VEC_TO_DISC(KTEX,LEN_TEX,NTEX_G,IREW,LU)
*
* put a GIC vector to DISC
*
*. Jeppe Olsen, Billund on the way to Zurich, march 2010
*
      INCLUDE 'wrkspc.inc'
*. Input: pointers to start and length of each TEX
      INTEGER KTEX(NTEX_G), LEN_TEX(NTEX_G)
*
      NTEST = 00
      IF(NTEST.GE.100) WRITE(6,*) ' Entering GIC_VEC_TO_DISC'
      IF(IREW.EQ.1) CALL REWINO(LU)
*
      DO IEX = 1, NTEX_G
C?      WRITE(6,*) ' Record to be written ', IEX
        KP = KTEX(IEX)
        LEN = LEN_TEX(IEX)
        CALL VEC_TO_DISC(WORK(KP),LEN,-1,-1,LU)
      END DO
      KP = KTEX(NTEX_G+1)
      LEN = 1
C?    WRITE(6,*) ' Reference coefficient written', WORK(KP)
      CALL VEC_TO_DISC(WORK(KP),LEN,-1,-1,LU)
C?    IF(NTEST.GE.100) WRITE(6,*) ' Leaving GIC_VEC_TO_DISC'
*
      RETURN
      END
      SUBROUTINE GIC_VEC_FROM_DISC(KTEX,LEN_TEX,NTEX_G,IREW,LU)
*
* Read a GIC vector to DISC
*
*. Jeppe Olsen, Billund on the way to Zurich, march 2010
*
      INCLUDE 'wrkspc.inc'
*. Input: pointers to start and length of each TEX
      INTEGER KTEX(NTEX_G), LEN_TEX(NTEX_G)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Entering GIC_VEC_FROM_DISC'
        WRITE(6,*) ' IREW, LU = ', IREW, LU
      END IF
*
      IF(IREW.EQ.1) CALL REWINO(LU)
*
      DO IEX = 1, NTEX_G
C?      WRITE(6,*) ' Record to be read ', IEX
        KP = KTEX(IEX)
        LEN = LEN_TEX(IEX)
        CALL VEC_FROM_DISC(WORK(KP),LEN,-1,-1,LU)
C?      WRITE(6,*) ' Record read '
      END DO
*. And the coefficient of the reference state
      KP = KTEX(NTEX_G+1)
      LEN = 1
      CALL VEC_FROM_DISC(WORK(KP),LEN,-1,-1,LU)
C?    WRITE(6,*) ' coefficient read in', WORK(KP)
*
      IF(NTEST.GE.100) WRITE(6,*) ' Leaving GIC_VEC_FROM_DISC'
      RETURN
      END
      SUBROUTINE H_S_EXT_GICCI_TV(VECIN,VECOUT_H,VECOUT_S,
     &                           I_DO_H,I_DO_S)
*
*. Obtain gradient of general GICCI vector function for
*  active operator ITACT (given in gicci)
*
* The current set of T-parameters are stored at KTEX_FOR_IEX
*
* The input is the T-coefficients for the active operators
* The remaining operators are accessed through KTEX. 
* KTEX is also updated with the coefficients in VECIN
*
*
* If(I_DO_H.eq.1) vecout_h(i): (I = ITACT)
*     <L|O(i,I)|R>
*     <F(I)!H!0'>
* where
*     |R> = T(I-1)...T(1)|ref>
*     |L> = P(I)(H|0'> + O+(I+1)H|0'> + .... + O+(N)...O(I+1)H|0'>)
*     |F(I)> = (1 + O(1) + O(2)O(1) + ... + O(I-1)...O(1)|ref>
*
* if(I_DO_S.eq.1) vecout_s(i) : 
*     <L'|O(i,I)|R>
*     <F(I)!0'>
* where
*     |R> = O(I-1)...O(1)|ref>
*     |L'> = P(I) (|0'> + O+(I+1)|0'> + .... + O+(N)...O(I+1)|0'>)
*
*  where O(J) as usual is a combination of a projection operator
* and a two-electron operator
*
*  O(J) = P(J) T(J)
*  P(J) projects on a space (ITCSPC(J)) and projects a space out
*  (IPTCSPC(J))
*
* <0!0> is assumed normalized
*
* Vecin is supposed to be delivered in SPA basis (if I_DO_EI = 0)
* or in the Zeroorder basis (if I_DO_EI = 1)
*
* Jeppe Olsen, March 2010 for the Zurich conference
*
      INCLUDE 'wrkspc.inc'
      REAL*8
     &INPRDD
      INCLUDE 'clunit.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc'
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT_H(*), VECOUT_S(*)
*. For transfer of data
      INCLUDE 'gicci.inc'
      NTEST = 00
*
      NSPA = NSPA_FOR_IEX(ITACT)
      NCAAB = NCAAB_FOR_IEX(ITACT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '---------------------------------'
        WRITE(6,*) ' Reporting from  H_S_EXT_GICCI_TV '
        WRITE(6,*) '---------------------------------'
        WRITE(6,*)
        WRITE(6,*) ' ITACT = ', ITACT
        WRITE(6,*) ' I_DO_H, I_DO_S =', I_DO_H, I_DO_S
        WRITE(6,*) ' NSPA, NCAAB = ', NSPA, NCAAB
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input vector for active operator'
        CALL WRTMAT(VECIN,1,NSPA,1,NSPA)
        WRITE(6,*) ' The current set of T-parameters'
        CALL WRT_GICCI_VEC(KTEX_FOR_IEX)
C            WRT_GICCI_VEC(KTEX)
      END IF
      
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HSG_TV')
      CALL MEMMAN(KL_VIC1,NCAAB_MX+1,'ADDL  ',2,'VIC1  ')
      CALL MEMMAN(KL_VIC2,NCAAB_MX+1,'ADDL  ',2,'VIC2  ')
*
*.  Obtain GICCI vector |0'> corresponding to set of coefficients
*   for active operator 
*. Obtain T-coefficients for |0'> in KTEXP_FOR_IEX
COLD  CALL COPVEC(WORK(KTEX_FOR_IEX(1)),WORK(KTEXP_FOR_IEX(1)),NSPA_TOT)
COLD  WORK(KTEXP_FOR_IEX(NTEXC_GX+1)) = WORK(KTEX_FOR_IEX(NTEXC_GX+1))
C     UPDATE_GICCI_VEC(KTEX,I_EX_ACT,TACTVEC,ISCALE)
COLD  CALL UPDATE_GICCI_VEC(KTEXP_FOR_IEX,ITACT,VECIN,1)
C     GET_GICCI_DELTA(KTEXG,IACT,TACT,LUC,LUOUT,LUSC2,
C    &                         LUSC3)
*- Obtain |0'> on LUSC1 using LUSC2 and LUSC3 as scratch
COLD  CALL GET_GICCI_0(KTEXP_FOR_IEX,LUSC1,LUC,LUSC35,LUSC2,LUSC3)
C     CALL GET_GICCI_DELTA(KTEX_FOR_IEX,IACT,TACT,LUC,LUSC1,
C    &                     LUSC2,LUSC3)
C     GET_GICCI_DELTA(KTEXG,IACT,TACT,LUC,LUOUT,LUSC2,
C    &                         LUSC3)
*
      IF(I_DO_H.EQ.1) THEN
*
* ================
*. Hamiltonian terms 
* ================
*
* If(I_DO_H.eq.1) vecout_h(i) :
*     <L|O(i,I)|R>
*     <F(I)!H!0>
* where
*     |R> = O(I-1)...O(1)|ref>
*     |L> = P(I)(H|0'> + O+(I+1)H|0'> + .... + O+(N)...O(I+1)H|0'>)
*
*. 1.05: Obtain |0'> on LUSC1 using LUSC2 and LUSC3 as scratch
      CALL GET_GICCI_DELTA(KTEX_FOR_IEX,ITACT,VECIN,LUC,LUSC1,
     &                     LUSC2,LUSC3)
      XNORM0P = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUSC1,1,-1)
      IF(NTEST.GE.5) WRITE(6,*) ' Square norm of |0''> ', XNORM0P
*. 1: Obtain |L> on LUSC2
*. For simplicity evrything is calculated in the largest space
      ICSPC = ITCSPC_GX(NTEXC_GX)
      ISSPC = ITCSPC_GX(NTEXC_GX)
*
*. 1.1: H|0'> on LUHC
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input to MV7 '
        CALL WRTVCD(WORK(KVEC1P),LUSC1,1,-1)
      END IF
      CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUHC,0,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Result of MV7'
        CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
      END IF
      DHD = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUHC,1,-1)
      IF(NTEST.GE.5) WRITE(6,*) ' <Delta 0|H|Delta 0> = ', DHD
*
*. 1.2: Obtain |L> on LUSC2, start with H|0'>
*     |L> = P(I)(H|0'> + O+(I+1)H|0'> + .... + O+(N)...O(I+1)H|0'>)
*
      CALL COPVCD(LUHC,LUSC2,WORK(KVEC1P),1,-1)
      ICSPC = ITCSPC_GX(NTEXC_GX)
      ISSPC = ITCSPC_GX(NTEXC_GX)
      DO IEX = ITACT+1, NTEXC_GX
C?       WRITE(6,*) ' IEX = ', IEX
*. obtain O+(ITACT+1) ... O+(IEX)H|0'> on LUSC3
        CALL COPVCD(LUHC,LUSC3,WORK(KVEC1P),1,-1)
        DO ISUB = 0, IEX-ITACT-1
          JEX = IEX-ISUB
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IEX, ISUB, JEX =', IEX, ISUB, JEX
          CALL PREPARE_FOR_IEX(JEX)
*. Obtain T(JEX) amplitudes in CAAB basis in KL_VIC2
          KP = KTEX_FOR_IEX(JEX)
          CALL REF_CCV_CAAB_SP(WORK(KL_VIC2),WORK(KP),
     &          WORK(KL_VIC1),2) 
*. Conjugate amplitudes 
          CALL CONJ_CCAMP(WORK(KL_VIC2),1,WORK(KL_VIC1))
*. and conjugate spinorbital classes
          CALL CONJ_T
          CALL REWINO(LUSC3)
          CALL REWINO(LUSC35)
*. Start by projection- conjugated operator, copy result back to LUSC3
          IPROJSPC = IPTCSPC_GX(JEX)
          IF(IPROJSPC.NE.0) THEN
            LUSCX = -1
            CALL REWINO(LUSC3)
            CALL REWINO(LUSC35)
            CALL EXTR_CIV(IREFSM,ISSPC,LUSC3,IPROJSPC,2,
     &                    LUSC35,-1,LUSCX,1,1,IDC,NTEST)
          END IF
          CALL REWINO(LUSC3)
          CALL REWINO(LUSC35)
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC35,
     &               WORK(KL_VIC1),1)
          CALL COPVCD(LUSC35,LUSC3,WORK(KVEC1P),1,-1)
*. Clean up by conjugating classes back to original
          CALL CONJ_T
        END DO
*. and add to LUSC2
        ONE = 1.0D0
*  VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),ONE,ONE,LUSC2,LUSC3,
     &              LUSC35,1,-1)
        CALL COPVCD(LUSC35,LUSC2,WORK(KVEC1P),-1,-1)
      END DO
*. And project for active op
      IPROJSPC = IPTCSPC_GX(ITACT)
      IF(IPROJSPC.NE.0) THEN
        LUSCX = -1
        CALL REWINO(LUSC2)
        CALL REWINO(LUSC35)
        CALL EXTR_CIV(IREFSM,ISSPC,LUSC2,IPROJSPC,2,
     & LUSC35,-1,LUSCX,1,1,IDC,NTEST)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The L-vector '
        CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
      END IF
*
*.2     |R> = O(I-1)...O(1)|ref> on LUSC3
*
*. Expand [ref>
      ICSPC = ITCSPC_GX(NTEXC_GX)
      ISSPC = ITCSPC_GX(NTEXC_GX)
      CALL EXPCIV(IREFSM,1,LUC,ISSPC,LUSC3,-1,
     &            LUSC35,1,0,IDC,NTEST)
*
      DO IEX = 1, ITACT-1
* T(IEX) LUSC3 on LUSC35
        CALL PREPARE_FOR_IEX(IEX)
        KP = KTEX_FOR_IEX(IEX)
        CALL REF_CCV_CAAB_SP(WORK(KL_VIC2),WORK(KP),
     &          WORK(KL_VIC1),2) 
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC35,
     &               WORK(KL_VIC2),1)
*. P(IEX)T(IEX) LUSC3 on LUSC3
        IPROJSPC = IPTCSPC_GX(IEX)
        IF(IPROJSPC.EQ.0) THEN
*. Just copy
          CALL COPVCD(LUSC35, LUSC3,WORK(KVEC1P),1,-1)
        ELSE
          CALL EXTR_CIV(IREFSM,ISSPC,LUSC35,IPROJSPC,2,
     &                    LUSC3,-1,LUSCX,1,0,IDC,NTEST)
        END IF
*
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The R-vector '
        CALL WRTVCD(WORK(KVEC1P),LUSC3,1,-1)
      END IF
*. We are now ready to calculate obtain the density <L!O(mu,ITACT)|R>
      CALL PREPARE_FOR_IEX(ITACT)
      ZERO = 0.0D0
      NCAAB = NCAAB_FOR_IEX(ITACT)
      NSPA = NSPA_FOR_IEX(ITACT)
      CALL SETVEC(WORK(KL_VIC1),ZERO,NCAAB)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC2,
     &               WORK(KL_VIC1),2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The Sigma vector in the CAAB basis '
        CALL WRTMAT(WORK(KL_VIC1),1,NCAAB,1,NCAAB)
      END IF
*. And reform to SPA basis
      CALL REF_CCV_CAAB_SP(WORK(KL_VIC1),VECOUT_H,WORK(KL_VIC2),1) 
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The Sigma vector in the SPA basis '
        CALL WRTMAT(VECOUT_H,1,NSPA,1,NSPA)
      END IF
*. 2. Obtain on LUSC1 |F(I)> = (C_0 + T(1) + T(2)T(1) + ... + T(I-1)...T(1)|ref>
C     GET_GICCI_EXP(KTEXG,IEX_MAX,LUC,LUOUT,LUSC2,LUSC3)
      CALL GET_GICCI_EXP(KTEX_FOR_IEX,ITACT-1,LUC,LUSC1,LUSC2,LUSC3)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The F(I) vector '
        CALL WRTVCD(WORK(KVEC1P),LUSC1,1,-1)
      END IF
*.2.1 and <F(I)|H|0'>
      FIH0P = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUHC,1,-1)
      IF(NTEST.GE.1000) WRITE(6,*) ' FIH0P = ', FIH0P
      VECOUT_H(NSPA) = FIH0P
      END IF
*     ^ End of Hamiltonian terms were to be calculated
      IF(I_DO_S.EQ.1) THEN
*
* ================
*. Overlap terms 
* ================
*
* vecout_S(i) :
*     <L'|O(i,I)|R>
*     <F(I)!0'>
* where
*     |R> = O(I-1)...O(1)|ref>
*     |L'> = P(I)(|0'> + O+(I+1)|0'> + .... + O+(N)...O(I+1)|0'>)
*
*. 3. Obtain |L'> 
*
*. 3.05: Obtain |0'> on LUSC1 using LUSC2 and LUSC3 as scratch
C     CALL GET_GICCI_0(KTEXP_FOR_IEX,LUSC1,LUC,LUSC35,LUSC2,LUSC3)
      CALL GET_GICCI_DELTA(KTEX_FOR_IEX,ITACT,VECIN,LUC,LUSC1,
     &                     LUSC2,LUSC3)
C?    WRITE(6,*) ' After GET_GICCI_DELTA'
*
*. 3.1: Obtain |L'> on LUSC2, start with |0'>
*
      CALL COPVCD(LUSC1,LUSC2,WORK(KVEC1P),1,-1)
      ICSPC = ITCSPC_GX(NTEXC_GX)
      ISSPC = ITCSPC_GX(NTEXC_GX)
      DO IEX = ITACT+1, NTEXC_GX
*. obtain O+(ITACT+1) ... O+(IEX)|0'> on LUSC3
        CALL COPVCD(LUSC1,LUSC3,WORK(KVEC1P),1,-1)
        DO ISUB = 0, IEX-ITACT-1
          JEX = IEX-ISUB
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IEX, ISUB, JEX =', IEX, ISUB, JEX
          CALL PREPARE_FOR_IEX(JEX)
*. Obtain T(JEX) amplitudes in CAAB basis in KL_VIC2
          KP = KTEX_FOR_IEX(JEX)
          CALL REF_CCV_CAAB_SP(WORK(KL_VIC2),WORK(KP),
     &          WORK(KL_VIC1),2) 
*. Conjugate amplitudes 
          CALL CONJ_CCAMP(WORK(KL_VIC2),1,WORK(KL_VIC1))
*. and conjugate spinorbital classes
          CALL CONJ_T
          CALL REWINO(LUSC3)
          CALL REWINO(LUSC35)
*. Start by projection- conjugated operator, copy result back to LUSC3
          IPROJSPC = IPTCSPC_GX(JEX)
          IF(IPROJSPC.NE.0) THEN
            LUSCX = -1
            CALL REWINO(LUSC3)
            CALL REWINO(LUSC35)
            CALL EXTR_CIV(IREFSM,ISSPC,LUSC3,IPROJSPC,2,
     &                    LUSC35,-1,LUSCX,1,1,IDC,NTEST)
          END IF
          CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC35,
     &               WORK(KL_VIC1),1)
          CALL COPVCD(LUSC35,LUSC3,WORK(KVEC1P),1,-1)
*. Clean up by conjugating classes back to original
          CALL CONJ_T
        END DO
*. and add to LUSC2
        ONE = 1.0D0
*  VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),ONE,ONE,LUSC2,LUSC3,
     &              LUSC35,1,-1)
        CALL COPVCD(LUSC35,LUSC2,WORK(KVEC1P),-1,-1)
      END DO
*. And project for active op
      IPROJSPC = IPTCSPC_GX(ITACT)
      IF(IPROJSPC.NE.0) THEN
        LUSCX = -1
        CALL REWINO(LUSC2)
        CALL REWINO(LUSC35)
        CALL EXTR_CIV(IREFSM,ISSPC,LUSC2,IPROJSPC,2,
     & LUSC35,-1,LUSCX,1,1,IDC,NTEST)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The L(prime)-vector '
        CALL WRTVCD(WORK(KVEC1P),LUSC2,1,-1)
      END IF
C?    WRITE(6,*) ' After L(prime)'
*     |R> = O(I-1)...O(1)|ref> on LUSC3
*. Expand [ref>
      ICSPC = ITCSPC_GX(NTEXC_GX)
      ISSPC = ITCSPC_GX(NTEXC_GX)
      CALL EXPCIV(IREFSM,1,LUC,ISSPC,LUSC3,-1,
     &            LUSC35,1,0,IDC,NTEST)
      DO IEX = 1, ITACT-1
* T(IEX) LUSC3 on LUSC35
        CALL PREPARE_FOR_IEX(IEX)
        KP = KTEX_FOR_IEX(IEX)
        CALL REF_CCV_CAAB_SP(WORK(KL_VIC2),WORK(KP),
     &          WORK(KL_VIC1),2) 
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC35,
     &               WORK(KL_VIC2),1)
*. P(IEX)T(IEX) LUSC3 on LUSC3
        IPROJSPC = IPTCSPC_GX(IEX)
        IF(IPROJSPC.EQ.0) THEN
*. Just copy
          CALL COPVCD(LUSC35, LUSC3,WORK(KVEC1P),1,-1)
        ELSE
          CALL EXTR_CIV(IREFSM,ISSPC,LUSC35,IPROJSPC,2,
     &                    LUSC3,-1,LUSCX,1,0,IDC,NTEST)
        END IF
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The R-vector( for S) '
        CALL WRTVCD(WORK(KVEC1P),LUSC3,1,-1)
      END IF
*. We are now ready to calculate obtain the density <L'!O(mu,ITACT)|R>
      CALL PREPARE_FOR_IEX(ITACT)
      ZERO = 0.0D0
      NCAAB = NCAAB_FOR_IEX(ITACT)
      CALL SETVEC(WORK(KL_VIC1),ZERO,NCAAB)
      CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC3,LUSC2,
     &               WORK(KL_VIC1),2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The S-vector before REF_CCV '
        CALL WRTMAT(WORK(KL_VIC1),1,NCAAB,1,NCAAB)
      END IF
*. And reform to SPA basis
      CALL REF_CCV_CAAB_SP(WORK(KL_VIC1),VECOUT_S,WORK(KL_VIC2),1) 
      NSPA = NSPA_FOR_IEX(ITACT)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The S-vector after REF_CCV '
        CALL WRTMAT(VECOUT_S,1,NSPA,1,NSPA)
      END IF
*. 4. Obtain on LUHC |F(I)> = (C_0 + O(1) + O(2)O(1) + ... + O(I-1)...O(1)|ref>
C     GET_GICCI_EXP(KTEXG,IEX_MAX,LUC,LUOUT,LUSC2,LUSC3)
      CALL GET_GICCI_EXP(KTEX_FOR_IEX,ITACT-1,LUC,LUHC,LUSC2,LUSC3)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The F(I) vector '
        CALL WRTVCD(WORK(KVEC1P),LUHC,1,-1)
      END IF
*.4.1 and <F(I)|0>
      FI0P = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUSC1,LUHC,1,-1)
      IF(NTEST.GE.1000) WRITE(6,*) ' FI0P = ', FI0P
      VECOUT_S(NSPA) = FI0P
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Direct ICCI, external part '
        WRITE(6,*) ' Input vector '
        CALL WRTMAT(VECIN,1,NSPA,1,NSPA)
        IF(I_DO_H.EQ.1) THEN
          WRITE(6,*) ' H(ICCI) times input vector '
          CALL WRTMAT(VECOUT_H,1,NSPA,1,NSPA)
        END IF
        IF(I_DO_S.EQ.1) THEN
          WRITE(6,*) ' S(ICCI) times input vector '
          CALL WRTMAT(VECOUT_S,1,NSPA,1,NSPA)
        END IF
      END IF
* 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HSG_TV')
      RETURN
      END
      SUBROUTINE GET_GICCI_0(KTEXG,LUOUT,LUC,LUSC2,LUSC3)
*
* Obtain GICCI wavefunction as defined by amplitudes in WORK(KTEXG)
* and save in LUOUT
*
      INCLUDE 'wrkspc.inc'
      DIMENSION KTEXG(MXPCYC)
      INCLUDE 'gicci.inc'
*
C?    SCALE = WORK(KTEXG(NTEXC_GX+1))
C?    WRITE(6,*) ' scale from GET_GICCI =', SCALE
C?    WRITE(6,*) ' LUOUT, LUC, LUSC, LUSC2, LUSC3 =',
C?   &             LUOUT, LUC, LUSC, LUSC2, LUSC3
      CALL GET_GICCI_EXP(KTEXG,NTEXC_GX,LUC,LUOUT,LUSC2,LUSC3)
*
      RETURN
      END
      SUBROUTINE GET_GICCI_EXP(KTEXG,IEX_MAX,LUC,LUOUT,LUSC2,LUSC3)
*
* Obtain on LUOUT GICCI expansion of wavefunction i
* to excitation operator IEX_MAX:
*
* |GICCI> = C_0|ref> + O_1|ref> + O_2 O_1|ref> + .... 
*         + O_IEX_MAX ...O_1|ref>
*
*. For the set of GICCI coefficients in WORK(KTEXG)
*
*. Jeppe Olsen, Zurich, march 2010
*
      INCLUDE 'wrkspc.inc'
      REAL*8
     &INPRDD
C     INCLUDE 'clunit.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc'
*. Offsets to the individual excitation vectors
      INTEGER KTEXG(MXPCYC)
      INCLUDE 'gicci.inc'
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' -----------------------------'
        WRITE(6,*) ' Reporting from GET_GICCI_EXP '
        WRITE(6,*) ' -----------------------------'
        WRITE(6,*)
        WRITE(6,*) ' Excitations are included upto ', IEX_MAX
        WRITE(6,*) ' LUC, LUSC2, LUSC3, LUOUT =', 
     &               LUC, LUSC2, LUSC3, LUOUT
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTGIC0')
      CALL MEMMAN(KLVEC1,NCAAB_MX,'ADDL  ',2,'LVEC1 ')
      CALL MEMMAN(KLVEC2,NCAAB_MX,'ADDL  ',2,'LVEC2 ')
* reference vector is on LUC
*
*. Initialize Ref on LUOUT (|0>)
*             Ref on LUSC2 (|S_0>)
*.
*
      ICSPC = ITCSPC_GX(NTEXC_GX)
      ISSPC = ITCSPC_GX(NTEXC_GX)
*
      CALL REWINO(LUC)
      CALL REWINO(LUOUT)
*. expand reference to complete space
      CALL EXPCIV(IREFSM,1,LUC,ITCSPC_GX(NTEXC_GX),LUOUT,-1,
     &            LUSC2,1,0,IDC,NTEST)
      CALL COPVCD(LUOUT,LUSC2,WORK(KVEC1P),1,-1)
*
   
*. Iterate
      DO IEX = 1, IEX_MAX
        IF(NTEST.GE.1000) WRITE(6,*) ' IEX, ICSPC =', IEX,ICSPC
C            PREPARE_FOR_IEX(IEX)
        CALL PREPARE_FOR_IEX(IEX)
*. Obtain in KLVEC1 T(IEX) in CAAB basis
        CALL REF_CCV_CAAB_SP(WORK(KLVEC1),WORK(KTEXG(IEX)),
     &  WORK(KLVEC2),2) 
        NSPA_L = NSPA_FOR_IEX(IEX)
        NCAAB_L = NCAAB_FOR_IEX(IEX)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' CAAB and SPA expansion of T(IEX)-vector'
          CALL WRTMAT(WORK(KLVEC1),1,NCAAB_L,1,NCAAB_L)
          CALL WRTMAT(WORK(KTEXG(IEX)),1,NSPA_L,1,NSPA_L)
        END IF
        
*. |S_I> = O_I|S_I-1> on LUSC3
        CALL REWINO(LUSC2)
        CALL REWINO(LUSC3)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC2,LUSC3,
     &               WORK(KLVEC1),1)
*. Project space IPTCSCP(IEX) out
        IF(IPTCSPC_GX(IEX).EQ.0) THEN
*. No projections, transfer |S_I> to LUSC2
          CALL COPVCD(LUSC3,LUSC2,WORK(KVEC1P),1,-1)
        ELSE
*. Project space IPTCSCP(IEX) out
          IPROJSPC = IPTCSPC_GX(IEX)
*. T |vecin> on LUSC3 => P T |vecin> on LUSC2
*. No scratch file is needed for 1 root
          LUSCX = -1
          CALL REWINO(LUSC2)
          CALL REWINO(LUSC3)
          CALL EXTR_CIV(IREFSM,ISSPC,LUSC3,IPROJSPC,2,
     &                    LUSC2,-1,LUSCX,1,0,IDC,NTEST)
C              EXTR_CIV(ISM,ISPCIN,LUIN,
C    &                  ISPCX,IEX_OR_DE,LUUT,LBLK,
C    &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
        END IF
*. Add |S_I> to |0>
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        ONE = 1.0D0
        CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),ONE,ONE,LUOUT,LUSC2,LUSC3,
     &              1,-1)
        CALL COPVCD(LUSC3,LUOUT,WORK(KVEC1P),1,-1)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Result after operator ', IEX
          CALL WRTVCD(WORK(KVEC1P),LUOUT,1,-1)
        END IF
      END DO
*. We are now only missing to change the  coefficient of the 
*  reference state to C_0
      C_0 = WORK(KTEXG(NTEXC_GX+1))
C?    WRITE(6,*) ' C_0 in GET_GICCI', C_0
      ONE = 1.0D0
      FACTOR = C_0 - 1.0D0 
      CALL EXPCIV(IREFSM,1,LUC,ITCSPC_GX(NTEXC_GX),LUSC3,-1,
     &            LUSC2,1,0,IDC,NTEST)
      CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),ONE, FACTOR,LUOUT,LUSC3,
     &            LUSC2,1,-1)
      CALL COPVCD(LUSC2,LUOUT,WORK(KVEC1P),1,-1)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Final GICCI vector '
        CALL WRTVCD(WORK(KVEC1P),LUOUT,1,-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTGIC0')
*
      RETURN
      END
      SUBROUTINE UPDATE_GICCI_VEC(KTEX,I_EX_ACT,TACTVEC,ISCALE)
*
*  Modify the collected GICCI vector by coefficients in TACTVEC
*  which  are the coefficient for excitation I_EX_ACT
* and a coefficient for the operators preceeding I_EX_ACT
*
*. The coefficients in TACTVEC is in the SPA basis
*
* ISCALE is inactive
*
*. Jeppe Olsen, March 2010
*
      INCLUDE 'wrkspc.inc'
      INTEGER KTEXG(MXPCYC)
      INCLUDE 'gicci.inc'
*. Input
      DIMENSION TACTVEC(*),KTEX(MXPCYC)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from UPDATE_GICI_VEC'
        WRITE(6,*) ' ---------------------------'
        WRITE(6,*) ' Active excitation operator: ', I_EX_ACT
      END IF
*
*. The update: 
*. ===========
*
* I = I_EX_ACT:
* I = 1:
* -----
* C_0(new) = delta_0 C_0
* T_1(new) = delta
* T_J(new) = T_J for J> 1
*
* I > 1:
* ------
* C_0(new) = delta_0 C_0
* T_1(new) = T_1*delta_0
* T_I(new) = delta/delta_0
* T_J(new) = T_J for J neq 1,I
*
      NSPA = NSPA_FOR_IEX(I_EX_ACT)
      NSPA1 = NSPA_FOR_IEX(1)
      KP = KTEX(I_EX_ACT)
      K1 = KTEX(1)
      KREF = KTEX(1)-1+NSPA_TOT+1
      DELTA_0 = TACTVEC(NSPA)
      IF(NTEST.GE.100) WRITE(6,*) 
     & ' NSPA, KP, KREF DELTA_0 =', NSPA,KP, KREF, DELTA_0
*. Updated coefficient for reference state
      WORK(KREF) = DELTA_0*WORK(KREF)
*. Active excitations
      CALL COPVEC(TACTVEC,WORK(KP),NSPA-1)
      IF(I_EX_ACT.NE.1) THEN
        FACTOR = 1.0D0/DELTA_0
        CALL SCALVE(WORK(KP),FACTOR,NSPA-1)
      END IF
*. First excitationvector
      IF(I_EX_ACT.NE.1) THEN
        CALL SCALVE(WORK(K1),DELTA_0,NSPA1-1)
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Updated T_GICCI vector'
        CALL WRT_GICCI_VEC(KTEX)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_GICCI_VEC(KTEX)
* Write GICCI vector with coefficent KTEX and specifications 
* defined in  COM_H_S_EFF_GICCI_TV
*
      INCLUDE 'wrkspc.inc'
      INTEGER KTEX(MXPCYC)
      INCLUDE 'gicci.inc'
*
      DO IEX = 1, NTEXC_GX
        WRITE(6,*) ' Excitation operator number', IEX
        KP = KTEX(IEX)
        NSPA = NSPA_FOR_IEX(IEX)
        CALL WRTMAT(WORK(KP),1,NSPA,1,NSPA)
      END DO
      WRITE(6,*) ' Coefficient of reference =', 
     &            WORK(KTEX(1)-1+NSPA_TOT+1)
*
      RETURN  
      END
      SUBROUTINE GET_GICCI_DELTA(KTEXG,IACT,TACT,LUC,LUOUT,LUSC2,
     &                         LUSC3)
*
* Obtain on LUOUT the correction to the GICCI vector defined by 
* TACT and KTEXG
*
*
* |GICCI> = Delta*(C_0|ref> + O_1|ref> + ... O_(IACT-1)... O_2 O_1|ref> )
*         + O_IACT O_(IACT-1)....O_1|ref>
*         + O_(IACT+1) O_IACT .... O_1|ref>
*         + .....  
*         + O_IEX_MAX ...O_1|ref>
*
*. For O(I, I.NE. IACT) the coefficients in WORK(KTEXG) are used
*  whereas Delta and O(IACT) are defined by TACT
*
*. Jeppe Olsen, Aarhus, april 2010
*
      INCLUDE 'wrkspc.inc'
      REAL*8
     &INPRDD
C     INCLUDE 'clunit.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc'
*. Offsets to the individual excitation vectors
      INTEGER KTEXG(MXPCYC)
*. And active vector
      DIMENSION TACT(*)
      INCLUDE 'gicci.inc'
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' -----------------------------'
        WRITE(6,*) ' Reporting from GET_GICCI_DELTA '
        WRITE(6,*) ' -----------------------------'
        WRITE(6,*)
        WRITE(6,*) ' Active excitation ', IACT
        WRITE(6,*) ' LUC, LUSC2, LUSC3, LUOUT =', 
     &               LUC, LUSC2, LUSC3, LUOUT
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTGIDE')
      CALL MEMMAN(KLVEC1,NCAAB_MX,'ADDL  ',2,'LVEC1 ')
      CALL MEMMAN(KLVEC2,NCAAB_MX,'ADDL  ',2,'LVEC2 ')
* reference vector is on LUC
*
*. Initialize C_0 Ref on LUOUT (|0>)
*                 Ref on LUSC2 (|S_0>)
*.
*
      ICSPC = ITCSPC_GX(NTEXC_GX)
      ISSPC = ITCSPC_GX(NTEXC_GX)
*
      CALL REWINO(LUC)
      CALL REWINO(LUOUT)
      CALL REWINO(LUSC2)
*. expand reference to complete space
      CALL EXPCIV(IREFSM,1,LUC,ITCSPC_GX(NTEXC_GX),LUSC2,-1,
     &            LUOUT,1,0,IDC,NTEST)
      C_0 = WORK(KTEXG(NTEXC_GX+1))
C?    WRITE(6,*) ' C_0 in GET_GICCI', C_0
      CALL SCLVCD(LUSC2,LUOUT,C_0,WORK(KVEC1P),1,-1)
   
*. Iterate
      DO IEX = 1, NTEXC_GX
        IF(NTEST.GE.1000) WRITE(6,*) ' IEX, ICSPC =', IEX,ICSPC
        CALL PREPARE_FOR_IEX(IEX)
C            PREPARE_FOR_IEX(IEX)
        NSPA_L = NSPA_FOR_IEX(IEX)
        NCAAB_L = NCAAB_FOR_IEX(IEX)
*
        IF(IEX.EQ.IACT) THEN
*. Scale (C_0 + O_1 + O_2O_1 + ... O_(IACT-1)... O(1))|ref> with delta
          DELTA = TACT(NSPA_FOR_IEX(IACT))
          IF(NTEST.GE.1000) WRITE(6,*) ' DELTA = ', DELTA
          CALL SCLVCD(LUOUT,LUSC3,DELTA,WORK(KVEC1P),1,-1)
          CALL COPVCD(LUSC3,LUOUT,WORK(KVEC1P),1,-1)
        END IF
*. Obtain in KLVEC1 T(IEX) in CAAB basis
        IF(IEX.NE.IACT) THEN
          CALL REF_CCV_CAAB_SP(WORK(KLVEC1),WORK(KTEXG(IEX)),
     &    WORK(KLVEC2),2) 
        ELSE
          CALL REF_CCV_CAAB_SP(WORK(KLVEC1),TACT,
     &    WORK(KLVEC2),2) 
        END IF
*. Zero coef for unit op
        WORK(KLVEC1) = 0.0D0
*
        IF(NTEST.GE.10000) THEN
          WRITE(6,*) ' CAAB and SPA expansion of T(IEX)-vector'
          CALL WRTMAT(WORK(KLVEC1),1,NCAAB_L,1,NCAAB_L)
          WRITE(6,*)
          IF(IEX.NE.IACT) THEN
            CALL WRTMAT(WORK(KTEXG(IEX)),1,NSPA_L,1,NSPA_L)
          ELSE
            CALL WRTMAT(TACT,1,NSPA_L,1,NSPA_L)
          END IF
        END IF
        
*. |S_I> = O_I|S_I-1> on LUSC3
        CALL REWINO(LUSC2)
        CALL REWINO(LUSC3)
        CALL SIGDEN_CC(WORK(KVEC1P),WORK(KVEC2P),LUSC2,LUSC3,
     &               WORK(KLVEC1),1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' The  unprojected |S_I> '
          CALL WRTVCD(WORK(KVEC1P),LUSC3,1,-1)
        END IF
*. Project space IPTCSCP(IEX) out
        IF(IPTCSPC_GX(IEX).EQ.0) THEN
*. No projections, transfer |S_I> to LUSC2
          CALL COPVCD(LUSC3,LUSC2,WORK(KVEC1P),1,-1)
        ELSE
*. Project space IPTCSCP(IEX) out
          IPROJSPC = IPTCSPC_GX(IEX)
*. T |vecin> on LUSC3 => P T |vecin> on LUSC2
*. No scratch file is needed for 1 root
          LUSCX = -1
          CALL REWINO(LUSC2)
          CALL REWINO(LUSC3)
          CALL EXTR_CIV(IREFSM,ISSPC,LUSC3,IPROJSPC,2,
     &                    LUSC2,-1,LUSCX,1,0,IDC,NTEST)
C              EXTR_CIV(ISM,ISPCIN,LUIN,
C    &                  ISPCX,IEX_OR_DE,LUUT,LBLK,
C    &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
        END IF
*. Add |S_I> to |0>
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        ONE = 1.0D0
        CALL VECSMD(WORK(KVEC1P),WORK(KVEC2P),ONE,ONE,LUOUT,LUSC2,LUSC3,
     &              1,-1)
        CALL COPVCD(LUSC3,LUOUT,WORK(KVEC1P),1,-1)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Result after operator ', IEX
          CALL WRTVCD(WORK(KVEC1P),LUOUT,1,-1)
        END IF
*
      END DO
*     ^ End of loop over excitation operators
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Final GICCI_DELTA vector '
        CALL WRTVCD(WORK(KVEC1P),LUOUT,1,-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTGIDE')
*
      RETURN
      END
    

c $Id$
