      SUBROUTINE MAT_P_MATT(A,B,NR,NC,COEF)
*
* A(I,J) = A(I,J) + Coef*B(J,I)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output
      DIMENSION A(NR,NC)
*. Input
      DIMENSION B(NC,NR)
*
      DO J = 1, NC
      DO I = 1, NR
        A(I,J) = A(I,J) + COEF*B(J,I)
      END DO
      END DO
*
      RETURN
      END
      SUBROUTINE CON_BLOCKS(IATP,IBTP,JATP,JBTP,IASM,IBSM,JASM,JBSM,
     &                      IOCTPA,IOCTPB,ICONSPA,ICONSPB,
     &                      NOCTPA,NOCTPB,MXEXC,IH_OCC_CON,
     &                      INTERACT)
*
* Does CI blocks IATP IBTP interact with blocks JATP JBTP
*
*. Input
* ======
* IATP IBTP JATP JBTP : Supergroups, relative numbers
* IOCTPA, IOBTPB : Offset for type
* ICONSPA, ICONSPB : Connection matrices giving exciation 
*                    level between two string types
* MXEXC : Largest excitation level
* IH_OCC_CON : = 1 => Use only occupation conserving part of
*                     Hamiltonian
*
*. Output
*. INTERACT : =1 => The two blocks does interact 
* Jeppe Olsen, April 99
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
#include "global.fh"
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*
      INTERACT = 0
      IA_EXC = ICONSPA(IATP,JATP)
      IB_EXC = ICONSPB(IBTP,JBTP)
      IF(IH_OCC_CON.EQ.0) THEN
*. Usual one- or two- electron operator
        IF(MXEXC.EQ.1) THEN
          IF((IA_EXC.LE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM).OR.
     &       (IB_EXC.LE.1.AND.IATP.EQ.JATP.AND.IASM.EQ.JASM)    )
     &        INTERACT = 1
        ELSE IF(MXEXC.EQ.2) THEN
          IF((IA_EXC.LE.1.AND.IB_EXC.LE.1)                  .OR.
     &       (IA_EXC.EQ.2.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM).OR.
     &       (IB_EXC.EQ.2.AND.IATP.EQ.JATP.AND.IASM.EQ.JASM)    )
     &        INTERACT = 1
        END IF
      ELSE
*. Orbital conserving part of  Hamiltonian
        IF(IA_EXC.EQ.IB_EXC .AND. IB_EXC.LE.1 ) THEN
          IATP_ABS = IATP + IOCTPA-1
          IBTP_ABS = IBTP + IOCTPB-1
          JATP_ABS = JATP + IOCTPA-1
          JBTP_ABS = JBTP + IOCTPB-1
*. Find Orb space where alpha strings differ
          IPGAS = 0
          IMGAS = 0
          DO IGAS = 1, NGAS
            IAEL = NELFSPGP(IGAS,IATP_ABS)
            JAEL = NELFSPGP(IGAS,JATP_ABS)
            IF(IAEL-JAEL.EQ.1) IPGAS = IGAS
            IF(IAEL-JAEL.EQ.-1)IMGAS = IGAS 
          END DO
          IF(IPGAS.NE.0) THEN
            IPDIF = NELFSPGP(IPGAS,IBTP_ABS)-NELFSPGP(IPGAS,JBTP_ABS) 
          ELSE
            IPDIF = 0 
          END IF
*. corresponding differences in beta
          IF(IMGAS.NE.0) THEN
            IMDIF = NELFSPGP(IMGAS,IBTP_ABS)-NELFSPGP(IMGAS,JBTP_ABS) 
          ELSE
            IMDIF = 0
          END IF
          IF(IPGAS.EQ.0.AND.IMGAS.EQ.0) INTERACT = 1
          IF(IPGAS.NE.0.AND.IMGAS.NE.0) THEN
            IF(IPDIF.EQ.-1.AND.IMDIF.EQ.1) INTERACT = 1
          END IF
        END IF
      END IF
*
C     IF((IATP.EQ.18.AND.JATP.EQ.4.OR.IATP.EQ.4.AND.JATP.EQ.18).AND.
C    &   (IBTP.EQ.18.AND.JBTP.EQ.4.OR.IBTP.EQ.4.AND.JBTP.EQ.18)) THEN
C         WRITE(6,*) ' In CON .. '
C         WRITE(6,*) ' IA_EXC, IB_EXC, INTERACT = ',
C    &                 IA_EXC, IB_EXC, INTERACT
C     END IF
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
       WRITE(6,'(a,10i3)') ' PROC ', ga_nodeid(), IATP,IBTP,JATP,JBTP,
     &       IASM, IBSM, JASM, JBSM, INTERACT
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CONBLOCKS '
        WRITE(6,*) ' PROCESSOR ', ga_nodeid(), 'IS SPEAKING'
        WRITE(6,*) ' IATP IBTP JATP JBTP ',IATP,IBTP,JATP,JBTP
        WRITE(6,*) ' IASM IBSM, JASM, JBSM = ', IASM, IBSM, JASM, JBSM
        WRITE(6,*) ' IH_OCC_CON, INTERACT = ', IH_OCC_CON,INTERACT
      END IF
*
      RETURN
      END 
      REAL*8 FUNCTION INPRODB(VEC1,VEC2,NBLK,LBLK,I0BLK)
*
* Inner products between blocked vectors with check of 
* zero blocks. Zero blocks are flagged by a unit entry 
* in I0BLK or a negative block length 
*
*Jeppe Olsen, May 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC1(*),VEC2(*), LBLK(NBLK),I0BLK(NBLK)
*
      REAL*8 INPROD
*
      X = 0.0D0
      IOFF = 1
      DO IBLK = 1, NBLK
        LLBLK = LBLK(IBLK)
C?      WRITE(6,*) ' INPRODB IBLK LLBLK IOFF ',IBLK,LLBLK,IOFF
        IF(I0BLK(IBLK).EQ.0.AND.LLBLK.GT.0) THEN
          X = X + INPROD(VEC1(IOFF),VEC2(IOFF),LLBLK)
C?        WRITE(6,*) ' Vec1 and Vec2 blocks '
C?        CALL WRTMAT(VEC1(IOFF),1,LLBLK,1,LLBLK)
C?        CALL WRTMAT(VEC2(IOFF),1,LLBLK,1,LLBLK)
C?        write(6,*) ' Updated x', X
        END IF
        IF(LLBLK.GT.0) IOFF = IOFF + LLBLK
      END DO
*
      INPRODB = X
*
      RETURN
      END
      SUBROUTINE PART_CIV2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                  NSMST,MXLNG,IOCOC,ISMOST,
     &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                  ISIMSYM)
*
* Jeppe Olsen
*
* Last update : May 1999 : ISIMSYM added
*
* Partition a CI vector into batches of blocks. 
* The length of a batch must be atmost MXLNG 
* If ISIMSYM.EQ.1, TTS blocks that differs only in symmetry are not
*                  split.
*
* IF ICOMP. eq. 1 the complete civector is constructed 
*
*
* Compared to PART_CIV, the NOOS arrays have been eliminated.
* They are becoming the size defining arrays - atleast at
* the laptop
*
*. Output 
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start 
*          of batch
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta sym
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta 
*   IBATCH(5,*) : Offset of block with respect to start of batch in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of batch in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form                   
*   IBATCH(8,*) : Length of block, packed form 
*
*
*
* Jeppe Olsen, August 1995   
*
cVOG  IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'wrkspc.inc'
      include 'strinp.inc'
#include "global.fh"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 0
      if (newciblks.eq.1) then
       IF(NTEST.GE.100.and.ga_nodeid().eq.0) write(*,*) 'You have work!'
       stop
      else
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     PART_CIV2     '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' ISMOST array '
        CALL IWRTMA(ISMOST,1,NSMST,1,NSMST) 
        WRITE(6,*) ' IBLTP array '
        CALL IWRTMA(IBLTP,1,NSMST,1,NSMST) 
        WRITE(6,*) ' NSSOA, NSSOB '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
        WRITE(6,*) 'ISIMSYM, ICOMP = ', ISIMSYM,ICOMP
        WRITE(6,*) ' MXLNG = ', MXLNG
      END IF
*. Tired of compiler warnings
      INCLUDE = 0
      LBLOCKP = 0
*
*. block 1 
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
*. New order : ISM,IB,IA (leftmost inner loop )
        IF(ISM.LT.NSMST) THEN
          ISM = ISM + 1
        ELSE
          ISM = 1
          IF(IB.LT.NOCTPB) THEN
            IB = IB + 1
          ELSE
            IB = 1
            IF(IA.LT.NOCTPA) THEN
              IA = IA + 1
            ELSE
              IFINI = 1
            END IF
          END IF
        END IF
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
*. Size of TT block ( all symmetries)
      LBLOCK_AS = 0
      IF(ISIMSYM.EQ.1 .AND. ISM. EQ. 1 ) THEN
        DO IASM = 1, NSMST
          IBSM = ISMOST(IASM)
          NSTA = NSSOA(IASM,IA)
          NSTB = NSSOB(IBSM,IB)
          IF(IBLTP(IASM).EQ.0) GOTO 99
          IF(IBLTP(IASM).EQ.2.AND.IA.LT.IB) GOTO 99
          LBLOCK_AS = LBLOCK_AS + NSTA*NSTB
   99   CONTINUE
        END DO
        INCLUDE = 0
C?      WRITE(6,*) ' IA IB LBLOCK_AS', IA,IB, LBLOCK_AS
        IF(LENGTH+LBLOCK_AS.LE.MXLNG.OR.ICOMP.EQ.1) INCLUDE = 1
      END IF
*. Should this block be included
      IASM = ISM
      IBSM = ISMOST(IASM)
      IF(IDC.EQ.2) THEN
        IF(IA.LT.IB) GOTO 1000
        IF(IA.EQ.IB.AND.IASM.LT.IBSM) GOTO 1000
      END IF
*. can this block be included
        NSTA = NSSOA(ISM,IA)
        NSTB = NSSOB(IBSM,IB)
        LBLOCK= NSTA*NSTB
        if (lblock.eq.0) goto 1000
        IF(IDC.EQ.1.OR.IA.GT.IB.OR.(IA.EQ.IB.AND.IASM.GT.IBSM)) THEN
          LBLOCKP = NSTA*NSTB
        ELSE IF(IDC.EQ.2.AND.IA.EQ.IB.AND.IASM.EQ.IBSM) THEN
          LBLOCKP = NSTA*(NSTA+1)/2
        END IF
        IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
cVOG      WRITE(6,*) ' IASM IBSM IA IB LBLOCKP,LBLOCK' ,    
          WRITE(6,'(a,6I5)') ' IASM IBSM IA IB LBLOCKP,LBLOCK' ,
     &                 IASM,IBSM,IA,IB,LBLOCKP,LBLOCK
        END IF
*
        IF(ISIMSYM.EQ.0) THEN
          INCLUDE = 0
          IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) INCLUDE = 1
        END IF
*
        IF(INCLUDE.EQ.1) THEN
          NBLOCK = NBLOCK + 1
          IBLOCK = IBLOCK + 1
          LBATCH(NBATCH) = LBATCH(NBATCH)+1
          IBATCH(1,IBLOCK) = IA
          IBATCH(2,IBLOCK) = IB
          IBATCH(3,IBLOCK) = ISM
          IBATCH(4,IBLOCK) = IBSM
          IBATCH(5,IBLOCK) = LENGTH+1
          IBATCH(6,IBLOCK) = LENGTHP+1
          IBATCH(7,IBLOCK) = LBLOCK     
          IBATCH(8,IBLOCK) = LBLOCKP     
          LENGTH = LENGTH + LBLOCK
          LENGTHP= LENGTHP+ LBLOCKP
          LEBATCH(NBATCH) = LENGTHP           
          GOTO 1000
        ELSE IF(ICOMP.EQ.0.AND.INCLUDE.EQ.0.AND.NBLOCK.EQ.0) THEN
          WRITE(6,*) 
     &    ' Not enough space to include a single Block'
          WRITE(6,*) ' Since I cannot procede I will stop '
          WRITE(6,*) ' Insufficient space detected in PART_CIV'
          WRITE(6,*) ' Alter GAS space or raise Buffer from ', MXLNG
          WRITE(6,*) ' Required dimension is ', LBLOCK
          WRITE(6,*) 
          WRITE(6,*) ' IASM, IA, NSTA, IBSM, IB, NSTB = ',
     &                 IASM, IA, NSTA, IBSM, IB, NSTB 
          CALL MEMCHK
          STOP 'Error in PART_CIV2'
        ELSE
*. This batch is finished, goto next batch
          GOTO 2000
        END IF
 2002 CONTINUE
*
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
C?      WRITE(6,*) 'Output from PART_CIV'
C?      WRITE(6,*) '====================='
        WRITE(6,*)
        WRITE(6,*) ' Number of batches ', NBATCH    
        IBLOCKT = 0
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lengths of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            IBLOCKT = IBLOCKT + 1
            WRITE(6,'(10X,I5,2X,4I3,4I8)') 
     &      IBLOCKT,(IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
      end if !End if newciblks
*
      RETURN
      END
      SUBROUTINE FRMDSCN2(VEC,NREC,LBLK,LU,NO_ZEROING,I_AM_ZERO)
*
* Read  VEC as multiple record file, NREC records read
* If NO_ZEROING.EQ.0 then zero blocks are not
* set to zero;  a 1 is instead flagged in the relevant block
* of I_AM_ZERO
*
      IMPLICIT REAL*8(A-H,O-Z)
*. OUtput
      DIMENSION VEC(*),I_AM_ZERO(*)
*
      IOFF = 1
      DO IREC = 1, NREC
        CALL IFRMDS(LREC,1,LBLK,LU)
        CALL FRMDSC2(VEC(IOFF),LREC,LBLK,LU,IMZERO,IAMPACK,
     &               NO_ZEROING)
        I_AM_ZERO(IREC) = IMZERO
        IOFF = IOFF + LREC
      END DO
*
      RETURN
      END
      SUBROUTINE CLS_TO_BLK(NBLOCK,IBLK_TO_CLS,ICLS_A,IBLK_A)
*
* an array ICLS_A is given for each class as well as 
* a block to class array IBLK_TO_CLS.
* Obtain ICLS_A in block form (Well, I am bit tired and my 
* pedagogical explanations can be pretty lousy even when I
* am awake !A)
*
* Jeppe Olsen, Jan. 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*.input
      DIMENSION ICLS_A(*),IBLK_TO_CLS(*)
*.output
      INTEGER IBLK_A(*)
*
      DO IBLK = 1, NBLOCK
        IBLK_A(IBLK) = ICLS_A(IBLK_TO_CLS(IBLK))
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CLS_TO BLK : IBLKS_a array '
        CALL IWRTMA(IBLK_A,1,NBLOCK,1,NBLOCK)
      END IF
*
      RETURN
      END 
      SUBROUTINE ZAP_BLOCK_VEC(LUIN,LBLK,IBLKS_A,SEGMNT,LUSCR)
*
* Zap blocks in vector in file LUIN for which IBLKS_A is zero
*
* Vector is initially constructed on LUSCR, and is copied back to 
* LUIN after use
*
* Note : Files are always rewinded 
* Packed version
*
* LBLK defines file type
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
      INTEGER IBLKS_A(*)
*
      CALL REWINE(LUIN  ,LBLK)
      CALL REWINE(LUSCR ,LBLK)
*
      IBLK = 0
*. Loop over blocks
C?      write(6,*) ' ZAP_BLOCK_VEC :  LBLK = ', LBLK
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUSCR) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUSCR)
        END IF
*
        IF( LBL .GE. 0 ) THEN
*
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
*
          IF(IBLKS_A(IBLK).EQ.1) THEN
            NO_ZEROING = 1
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &      NO_ZEROING)
          ELSE
            CALL SKPRCD2(LBL,-1,LUIN)
            IAMPACK = 1
            IMZERO = 1
          END IF
*
* Not packed
*
          IF(IAMPACK.EQ.0) THEN
            IF(IBLKS_A(IBLK).EQ.0) THEN
              ZERO = 0.0D0
              CALL SETVEC(SEGMNT,ZERO,LBL)
            END IF
            CALL TODSC(SEGMNT,LBL,KBLK,LUSCR)
          END IF
*
* Packed
*
          IF(IAMPACK.EQ.1) THEN 
            IF(IMZERO.EQ.1.OR.IBLKS_A(IBLK).EQ.0) THEN
              CALL ZERORC(LBL,LUSCR,IAMPACK)
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUSCR)
            END IF
          END IF
*
        END IF
*       ^ End if LBL .GE. 0


      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
*. And then copy back to LUIN
C          COPVCDP(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      IREW = 1
      CALL LUCIAQENTER('COPVC')
      CALL COPVCD(LUSCR,LUIN,SEGMNT,IREW,LBLK)
      CALL LUCIAQEXIT('COPVC')
*
      RETURN
      END
      SUBROUTINE CLASS_TRUNC(NCLS,ICLS_L,RCLS_L,CLS_CT,CLS_ET,CLS_C,
     &                       CLS_E,
     &                       E_CONV,ICLS_A,N_TRN_CLS,E_TRUNC,W_TRUNC)
*
* Decide which classes of parameters that can be eliminated
*
* Jeppe Olsen, Jan 97
*              March '97 updated
*              March '99 Modified to contribution per det
*                        (Previous version favored blocks 
*                         with many small dets)
*
*. Note in current version all energy contributions are 
*. positive. 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS_L(NCLS),RCLS_L(NCLS)     
      DIMENSION CLS_CT(NCLS),CLS_ET(NCLS),CLS_C(NCLS),CLS_E(NCLS)
*. Output
      INTEGER ICLS_A(NCLS)
*. Giving  the truncated classes
*  Additional output is
*              N_TRN_CLS : Number of truncated classes
*              E_TRUNC   : Estimated error of eliminating these classes
*                          (compared to expansion with only largest 
*                           coefficients included).
*                          E_TRUNC is thus the error that should be 
*                          added to the estmated error arising from
*                          eliminating small terms
*
      NTEST = 00
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' Welcome to CLASS_TRUNC '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Required threshold for convergence of energy',E_CONV
      END IF
*. Tired of compiler warnings
      E_THRES = 0.0D0
      E_THRES2 = 0.0D0
*
*. First find total energy correction and 
*. largest class contribution 
*
      E_TOT = 0.0D0
      E_CLS_MAX= 0.0D0
      E_CLS_MAX_PDET = 0.0D0
      DO JCLS = 1, NCLS
       E_TOT = E_TOT + CLS_ET(JCLS)
       E_CLS_MAX = MAX(E_CLS_MAX,CLS_ET(JCLS))
       IF(RCLS_L(JCLS).GT.0) THEN
         E_CLS_PDET = ABS(CLS_ET(JCLS))/RCLS_L(JCLS)
         E_CLS_MAX_PDET =  MAX(E_CLS_MAX_PDET,E_CLS_PDET)    
       END IF
      END DO
      WRITE(6,*) 'E_CLS_MAX_PDET = ', E_CLS_MAX_PDET
C?    WRITE(6,*) ' E_CLS_MAX = ', E_CLS_MAX
*
*. The truncation is done in two steps, first an overall 
*  threshold for deleting classes is constructed,
*  and all classes with energy contributions less than this 
*  threshold are excluded.
*. Then one takes a pass through remaining classes and
*  eliminates untill the quota is filled
*
*. Threshold for including determinants, start with E_CLS_MAX 
*. and decrease until quota is filled
*
      WRITE(6,*)
      WRITE(6,*) 
     &' Part I : Eliminating all classes below common threshold :'
      WRITE(6,*)
*
      ILOOP = 0
      STEP = 1.2
 1000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
          E_THRES =  E_CLS_MAX_PDET
        ELSE
          E_THRES = E_THRES/STEP
        END IF
        E_ELI = 0.0D0
        N_ELI = 0
        DO JCLS = 1, NCLS
          IF(ABS(CLS_ET(JCLS))/RCLS_L(JCLS).LE.E_THRES) THEN
            E_ELI = E_ELI + ABS(CLS_ET(JCLS))
            N_ELI = N_ELI + 1
          END IF
        END DO
        IF(ILOOP.EQ.100000000) THEN
          WRITE(6,*) ' Loop count exceeded 100000000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          STOP'CLASS_TRUNC: Iloop.eq.100000000'
        END If
      IF(E_ELI.GT. E_CONV) GOTO 1000
      IF(NTEST.GE.10) WRITE(6,*)
     & ' Common treshold for eliminating classes',E_THRES
*
      N_PAS_CLS = 0
      XL_PAS_CLS = 0.0D0
      N_TRN_CLS = 0
      XL_TRN_CLS = 0.D0
      N_ACT_CLS = 0
      XL_ACT_CLS = 0.0D0
*
      E_TRUNC = 0.0D0
      E_TRUNCT = 0.0D0
      W_TRUNCT= 0.0D0
      W_TRUNC = 0.0D0
*
*. Eliminate classes with energy contribution less than E_THRES
*
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).NE.0.0D0.AND
     &     .ABS(CLS_ET(JCLS))/RCLS_L(JCLS).LE. E_THRES) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          XL_TRN_CLS = XL_TRN_CLS + RCLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS) 
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS) 
          ICLS_A(JCLS) = 0
        END IF
      END DO
      IF(NTEST.GE.10) THEN
      WRITE(6,*) 
     &  ' Number of classes eliminated by common threshold',
     &   N_TRN_CLS
      WRITE(6,*) ' Energy contributions from these classes ',
     &   E_TRUNCT
      END IF
* Eliminate remaining classes until thres hold is obtained
      WRITE(6,*)
      WRITE(6,*) 
     &' Part II : Eliminate classes until total threshold is obtained'
      WRITE(6,*)
      NTRN_CLS3 = 0
      XLTRN_CLS3 = 0.0D0
      ETRN_CLS3 = 0.0D0
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).EQ.0.0D0) THEN
*. Passive class, no contribution before truncation
          N_PAS_CLS = N_PAS_CLS + 1
          XL_PAS_CLS = XL_PAS_CLS + RCLS_L(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1
     &  .AND.ABS(E_TRUNCT+CLS_ET(JCLS)).LT.E_CONV) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          XL_TRN_CLS = XL_TRN_CLS + RCLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS) 
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS) 
          NTRN_CLS3 = NTRN_CLS3 + 1
          XLTRN_CLS3 = XLTRN_CLS3 + RCLS_L(JCLS)
          ETRN_CLS3 = ETRN_CLS3 + CLS_E(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1) THEN 
*. Class is active
          N_ACT_CLS = N_ACT_CLS + 1
          XL_ACT_CLS = XL_ACT_CLS + RCLS_L(JCLS)
           ICLS_A(JCLS) = 1
        END IF
      END DO
      IF(NTEST.GE.10) THEN
      WRITE(6,*) 
     &  ' Number of classes eliminated in second part',
     &   NTRN_CLS3
      WRITE(6,*) ' Energy contributions from these classes ',
     &   ETRN_CLS3
      END IF
*. Correct for missing minus in first order correction (not here, not now )
      E_TRUNC  =  E_TRUNC 
      E_TRUNCT =  E_TRUNCT 
*
      IF(NTEST.GE.10) THEN
      WRITE(6,*) 
      WRITE(6,*) ' Result of truncation, part I + II : '
      WRITE(6,*) ' ===================================='
      WRITE(6,*) '   Estimated complete energy contribution of ',
     &            '  eliminated classes ', -E_TRUNCT
      WRITE(6,*) '   Estimated truncated energy contribution of ',
     &            '  eliminated classes ', -E_TRUNC
      WRITE(6,*) '   Truncation error in eliminated classes was ',
     &              -E_TRUNC+E_TRUNCT
      WRITE(6,*) '   Estimated truncated weight contribution of ',
     &            '  eliminated classes ', W_TRUNC
      WRITE(6,*) '   Estimated energy contribution without trunc',
     &             -E_TOT
      WRITE(6,*) '   Energy contribution of of active classes ',
     &             -E_TOT+E_TRUNCT
      END IF
       
C?    WRITE(6,*)
C?    WRITE(6,*) '  Class      Number     Dimension  '
C?    WRITE(6,*) ' ================================= '
C?    WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Passive  ', N_PAS_CLS,L_PAS_CLS
C?    WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Truncated', N_TRN_CLS,L_TRN_CLS
C?    WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Active   ', N_ACT_CLS,L_ACT_CLS
*
*. Eliminate classes with energy contributions less than E_TOT* FAC2   
*  (temporary elimination, these classes may be invoked later in 
*  the iterative sequence). Only the active classes are examined.
*
      WRITE(6,*)
      WRITE(6,*) 
     &' Part III : Eliminate classes using dynamic threshold'
      WRITE(6,*)
      FAC2 = 0.1D0
*C Changed for comparisons
      FAC2 = 0.0D0
      WRITE(6,*) 
     &' Dynamic threshold is ', FAC2, ' times expected energychange'
*
      ILOOP = 0
      STEP = 1.2
      E_TEMP_TRUNC =  ABS(E_TOT)*FAC2
 2000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
          E_THRES2 = E_CLS_MAX_PDET
        ELSE
          E_THRES2 = E_THRES2/STEP
        END IF
        E_ELI2 = 0.0D0
        N_ELI2 = 0
        DO JCLS = 1, NCLS
          IF(ICLS_A(JCLS).EQ.1.AND.
     &       ABS(CLS_ET(JCLS))/RCLS_L(JCLS) .LE. E_THRES2) THEN
            E_ELI2 = E_ELI2 + ABS(CLS_ET(JCLS))
            N_ELI2 = N_ELI2 + 1
          END IF
        END DO
        IF(ILOOP.EQ.100000000) THEN
          WRITE(6,*) ' Loop count exceeded 100000000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          STOP'CLASS_TRUNC: Iloop.eq.10000000'
        END If
      IF(E_ELI2.GT. E_TEMP_TRUNC) GOTO 2000
      IF(NTEST.GE.10) THEN
      WRITE(6,*) 
     & ' Temporary elimination of classes with total contribution'
      WRITE(6,*) ' less than ',E_TEMP_TRUNC 
      WRITE(6,*) ' gives threshold for temporary elimination',E_THRES2
      END IF
*
*. Eliminate classes with energy contribution less than E_THRES2
*
      N_TRN_CLS2 = 0
      XL_TRN_CLS2 = 0
      E_TRUNC2 = 0
      W_TRUNC2 = 0
      DO JCLS = 1, NCLS
        IF(ICLS_A(JCLS).EQ.1.AND.CLS_CT(JCLS).NE.0.0D0.AND.
     &     ABS(CLS_ET(JCLS))/RCLS_L(JCLS).LE. E_THRES2) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          N_TRN_CLS2= N_TRN_CLS2+ 1
          XL_TRN_CLS = XL_TRN_CLS + RCLS_L(JCLS)
          XL_TRN_CLS2= XL_TRN_CLS2+ RCLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNC2 = E_TRUNC2 + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS) 
          W_TRUNC2 = W_TRUNC2+ CLS_C (JCLS) 
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS) 
          ICLS_A(JCLS) = 0 
        END IF
      END DO
      IF(NTEST.GE.1) THEN
      WRITE(6,*)
      WRITE(6,'(A,F25.12)')
     &  ' Energy contributions from ALL eliminated classes ', -E_TRUNC
      WRITE(6,'(A,F25.12)')
     &  ' Norm of ALL eliminated classes                   ',
     &   SQRT(W_TRUNC)
      END IF
*
      N_ACT_CLS = N_ACT_CLS - N_TRN_CLS2
      XL_ACT_CLS = XL_ACT_CLS - XL_TRN_CLS2
*
      E_TRUNC = - E_TRUNC
      E_TRUNCT = -  E_TRUNCT 
      
      WRITE(6,*)
      WRITE(6,*) '             Number            Dimension  '
      WRITE(6,*) ' ========================================= '
      WRITE(6,'(1H ,A,5X,I5,5X,E22.15)') 
     &'Passive  ', N_PAS_CLS,XL_PAS_CLS
      WRITE(6,'(1H ,A,5X,I5,5X,E22.15)') 
     &'Truncated', N_TRN_CLS,XL_TRN_CLS
      WRITE(6,'(1H ,A,5X,I5,5X,E22.15)') 
     &'(Temp)   ', 
     &N_TRN_CLS2,XL_TRN_CLS2
      WRITE(6,'(1H ,A,5X,I5,5X,E22.15)') 
     &'Active   ', N_ACT_CLS,XL_ACT_CLS
*
      WRITE(6,*)
      WRITE(6,*) ' Information about classes '
      WRITE(6,*) ' =========================='
      WRITE(6,*)
      WRITE(6,*)
     & ' Class    Dimension      E          E(Trunc)      C',
     & '         C(Trunc)   Active'
      WRITE(6,*)
     & ' =============================================================',
     & '================'
       DO JCLS = 1, NCLS
         IF(ABS(CLS_ET(JCLS)).GT.0.0D0.OR.ABS(CLS_CT(JCLS)).GT.0.0D0)
     &   THEN
           WRITE(6,'(1H ,2X,I4,1X,E20.12,3X,4(E12.6,1X),1X,I2)') 
     &     JCLS,RCLS_L(JCLS),CLS_ET(JCLS),CLS_E(JCLS),
     &     CLS_CT(JCLS),CLS_C(JCLS),ICLS_A(JCLS)
         END IF
       END DO
*
       WRITE(6,*) ' The ICLS_A array (active classes)'
       CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
*
      RETURN
      END
      SUBROUTINE CLASS_PROD3(VEC1,VEC2,IBLOCK_OFF,NBLOCK,IBLOCK,
     &                       IBLTOCLS,NCLS,CLSVEC)
*
* Two vectors in blocked form are given.
* Find contributions of product to each occupation class
*
* Current version (PROD3 !!) uses IBLTOCLS to give relation
* between blocks and classes
      IMPLICIT REAL*8(A-H,O-Z)
*
      REAL*8 INPROD
*. Input
      DIMENSION VEC1(*),VEC2(*)
      INTEGER IBLOCK(8,*), IBLTOCLS(*)            
*. Input/output
      DIMENSION CLSVEC(*)
      
*
      IOFF = 1
      DO JBLOCK = IBLOCK_OFF,IBLOCK_OFF-1+NBLOCK
      IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        JATP = IBLOCK(1,JBLOCK)
        JBTP = IBLOCK(2,JBLOCK)
        NELMNT = IBLOCK(8,JBLOCK)      
        JCLS = IBLTOCLS(JBLOCK)
        XTERM = INPROD(VEC1(IOFF),VEC2(IOFF),NELMNT)
C?      WRITE(6,*)
C?   &  ' CLASS_PROD : CLASS and XTERM = ', JCLS,XTERM
          CLSVEC(JCLS) = CLSVEC(JCLS) + XTERM
COLD    END IF
        IOFF = IOFF + NELMNT
      END IF
      END DO
*
      NTEST = 0
      IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Updated CLSVEC '
         CALL WRTMAT(CLSVEC,1,NCLS,1,NCLS)
      END IF
*
      RETURN
      END  
      SUBROUTINE HINTV(LURHS,LUX,SHIFT,SHIFT_DIA,VECIN,VECOUT,
     &                LBLK,LUPROJ,LUPROJ2)
*
* Solve  (H+Shift)X = RHS
*
* Where H is matrix rep of operator defined by /OPER/ in space defined by
* /CANDS/
*
* If ICISTR.EQ.1 VECIN contains RHS, else RHS is assumed  on LURHS
* Output : solution is on LUX
*
* Jeppe Olsen, Winter of 1996
* 
* Jan. 98 : SHIFT_DIA added
*
c      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD , INPRDD
      LOGICAL CONVER
      INCLUDE 'cands.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'clunit.inc'
* SCRATCH files used : LUSC3,LUSC34,LUSC35,LUSC37 
*. These are not used in SIMPRT ! 
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Max number of iterations is picked from MAXCIT in crun
      EXTERNAL H0TVM
      DIMENSION ERROR(100)

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HINTV ')
*
      NTEST = 0
*
*
* 1 : Construct diagonal  on LUSC3
* ================================
*
*. Use type of H0 as type of zero order operator  in all spaces
       IF(IAPR.NE.0) THEN
*. Not exact Hamiltonian
         IF(IH0INSPC(1).EQ.1.OR.IH0INSPC(1).EQ.3) THEN
*. Mp operator
           I12 = 1
CBERT: Just swapping of pointers would do
cnw        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
cnw        CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
           ktmp = kint1
           kint1 = kfi
           kfi = ktmp
           ktmp = kint1o
           kint1o = kfio
           kfio = kint1o
         ELSE IF (IH0INSPC(1).EQ.2.OR.IH0INSPC(1).EQ.4) THEN
*. EN diagonal
           I12 = 2
         END IF
       END IF
       ECOREX = SHIFT_DIA
*
* Partitioning and blockstructure of CI vector
*
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      NTTS = MXNTTS
      CALL MEMMAN(KLCLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')    !done
      CALL MEMMAN(KLCLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')   !done
      CALL MEMMAN(KLCI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')    !done
      CALL MEMMAN(KLCIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')    !done
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'CIOIO ')   !done
      CALL IAIBCM(ISSPC,int_mb(KLCIOIO))
*
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')   !done
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,dbl_mb(KLCBLTP),int_mb(KSVST))
*. Allocate memory for diagonalization
      IF(ISIMSYM.EQ.0) THEN
        LBLOCK = MXSOOB
      ELSE
        LBLOCK = MXSOOB_AS
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
*. Batches  of C vector
      CALL PART_CIV2(IDC,dbl_mb(KLCBLTP),int_mb(KNSTSO(IATP)),
     &              int_mb(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,int_mb(KLCIOIO),
     &              ISMOST(1,ISSM),NBATCH,int_mb(KLCLBT),
     &              int_mb(KLCLEBT),int_mb(KLCI1BT),int_mb(KLCIBT),0,
     &              ISIMSYM)
*. Number of BLOCKS
        NBLOCK = IFRMR(int_mb(KLCI1BT),1,NBATCH)
     &         + IFRMR(int_mb(KLCLBT),1,NBATCH) - 1
C?      WRITE(6,*) ' HINTV : Number of blocks ', NBLOCK
        CALL GASDIAT(VECIN,LUSC3,ECOREX,ICISTR,I12,
     &               dbl_mb(KLCBLTP),NBLOCK,int_mb(KLCIBT))
C      CALL GASDIAT(VECIN,ISSM,ISSPC,LUSC3,ECOREX,ICISTR,I12)
*. Clean up time
       IF(IH0INSPC(1).EQ.1.OR.IH0INSPC(1).EQ.3) THEN
*. MP operator
CBERT: Just swapping of pointers would do
           ktmp = kint1
           kint1 = kfi
           kfi = ktmp
           ktmp = kint1o
           kint1o = kfio
           kfio = kint1o
cnw      CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
cnw      CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
       END IF
*
* 2 : Solve linear set of equations
* ==================================
*
      ZERO = 0.0D0
      IF(LBLK.GT.0 ) THEN
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        STOP       ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
*. Two vectors can be in core
*. Initial Guess on LUX   
        NDIM = LBLK
        CALL SETVEC(VECOUT,ZERO,NDIM)
        CALL REWINO(LUX)
        CALL TODSC(VECOUT,NDIM,-1,LUX)
*. Right hand side on LUSC34
        CALL REWINE(LUSC34,-1)
        CALL TODSC(VECIN,NDIM,-1,LUSC34)
*. Max number of its and convergence thresholds are picked up from
* corresponding eigenvalue info
        CONVER = .FALSE.
        TEST = SQRT(THRES_E) * SQRT(INPROD(VECIN,VECIN,NDIM))
        ILNPRT = MAX(NTEST-10,0)
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXITLE = ',MXITLE
        CALL MINGCG(MV8,LUX,LUSC34,LUSC35,LUSC3,VECIN,VECOUT,        
     &              MXIT_LOC,CONVER,TEST,SHIFT,ERROR,NDIM,
     &              LUPROJ,ILNPRT)
        CALL REWINE(LUX,-1)
        CALL FRMDSC(VECOUT,NDIM,-1,LUX,IMZERO,IAMPACK)     
*
         IF(NTEST.GE.5) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(VECOUT,1,NDIM,1,NDIM)
         END IF
*
      ELSE IF(LBLK.LE.0)   THEN
*
*. Use path allowing us to work with segments of vectors
*
*. Initial guess on LUX
        CALL SETVCD(LUSC3,LUX,VECOUT,ZERO,1,LBLK)
*. (Right hand side vector is assumed in place)
*. convergence threshold is picked up from
* corresponding eigenvalue info
        TEST = 
     &  SQRT(THRES_E) * SQRT(INPRDD(VECIN,VECOUT,LURHS,LURHS,1,-1))
        ILNPRT = NTEST
        SHIFT2 = 0.0D0
        CONVER = .FALSE.
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXIT_LOC ',MXIT_LOC
        CALL MICGCG(H0TVM,LUX,LURHS,LUSC34,LUSC35,LUSC37,LUSC3,
     &              VECIN,VECOUT,MXIT_LOC,
     &              CONVER,TEST,SHIFT2,ERROR,NDIM,LUPROJ,LUPROJ2,
     &              VFINAL,ILNPRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Solution to linear set of Equations '
          CALL WRTVCD(VECIN,LUX,1,LBLK)
        END IF
*
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HINTV ')
      RETURN
      END 
      SUBROUTINE SETVCD(LUIN,LUOUT,SEGMNT,VALUE,IREW,LBLK)
*
* Construct a vector with the same structure as LUIN 
* with values VALUE
*
* LBLK DEFINES STRUCTURE OF FILE
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
*
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF
 
*
* Loop over blocks
*
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUOUT)
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          CALL SETVEC(SEGMNT,VALUE,LBL)
          CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
      RETURN
      END



      SUBROUTINE H0TVM(VEC1,VEC2,LLUC,LLUHC)
*
* Outer routine for zero order operator + shift times vector
* 
*. Input  vector : on LLUC
*. Output fector : on LLUHC
*
* Jeppe Olsen, February 1996
*
c      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      DIMENSION VEC1(*),VEC2(*)
*
      INCLUDE 'clunit.inc'
C
*. Transfer of zero order energy
      COMMON/CENOT/E0
*. Transfer of shift 
      INCLUDE 'cshift.inc'
*. Default block parameter
      LBLK = -1 
*.  Zero order vector is assumed on LUC
      IF(IPROJ.EQ.0) THEN
       LU0 = 0
      ELSE IF (IPROJ.EQ.1) THEN
       LU0 = LUC
      ELSE
       WRITE(6,*)  ' H0TVM, Unknown IPROJ = ', IPROJ
       STOP ' H0TVM, Unknown IPROJ  '
      END IF
      LUSCR1 = LUSC40
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*)
        WRITE(6,*) '============== '
        WRITE(6,*) ' H0TVM entered '
        WRITE(6,*) '============== '
        WRITE(6,*)
        WRITE(6,*) ' LLUC LLUHC LU0 and LUSCR1 ',
     &               LLUC,LLUHC,LU0,LUSCR1
        WRITE(6,*) ' E0 , Shift : ', E0 , SHIFT 
      END IF
*. A scratch file not used by linear solver in SIMPRT : LUSCR1
      IF(SHIFT.EQ.0.0D0) THEN
        CALL H0TV2(VEC1,VEC2,LLUC,LLUHC,LU0,LUSCR1,E0,LBLK) 
      ELSE
*. H0TV on LUSCR1
        CALL H0TV2(VEC1,VEC2,LLUC,LUSCR1,LU0,LLUHC,E0,LBLK) 
*. Add shift and save on LLUHC
        ONE = 1.0D0
        CALL VECSMD(VEC1,VEC2,ONE,SHIFT,LUSCR1,LLUC,LLUHC,1,LBLK)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input and output vectors from H0TVM '
        CALL WRTVCD(VEC1,LLUC,1,LBLK)
        WRITE(6,*)
        CALL WRTVCD(VEC1,LLUHC,1,LBLK)
      END IF
*
      RETURN
      END 
      SUBROUTINE H0TV2(VEC1,VEC2,LUC,LUHC,LU0,LUSCR1,E0,LBLK)
*
* Multiply vector in LUC with H0 where H0 is defined as 
*
* H0 = (1-|0><0|) H apr  (1-|0><0>) + E0 |0><0>
*
* Where H apr is defined by call to MV7 ( with IPERTOP = 1)
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL*8 INPRDD
*. For communicating with sigma routine
      include 'oper.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
*. Tired of compiler warnings
      SC0 = 0.0D0
*
C     WRITE(6,*) ' H0TV2 : LUC, LUHC, LU0, LUSCR1 ',
C    &                     LUC, LUHC, LU0, LUSCR1
*. Overlap <C|0>
      IF(LU0.GT.0) THEN
        SC0 = INPRDD(VEC1,VEC2,LUC,LU0,1,LBLK)
*. C -  <C|0> |0> on LUSCR1
        FAC1 = 1.0D0
        FAC2 = -SC0
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUC,LU0,LUSCR1,1,LBLK)
      ELSE
        CALL COPVCD(LUC,LUSCR1,VEC1,1,LBLK)
      END IF
*. Multiply with H apr, result on LUHC
*. Place for trouble in the future.
      IF(IH0INSPC(1).NE.4) THEN
         IPERTOP = 1
      END IF
      WRITE(6,*) ' H0TV2: MV7 will be called in a few NANOSECONDS'
      CALL MV7(VEC1,VEC2,LUSCR1,LUHC,0,0)
*. Orthogonalize LUHC to LU0
      IF(LU0.NE.0) THEN
        SSIGMA0 = INPRDD(VEC1,VEC2,LUHC,LU0,1,LBLK)
        FAC1 = 1.0D0
        FAC2 = -SSIGMA0
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUHC,LU0,LUSCR1,1,LBLK)
*. and add E0 <C|0> |0>
        FAC1 = 1.0D0
        FAC2 = E0 * SC0
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUSCR1,LU0,LUHC,1,LBLK)
      ELSE
CSEPT29 CALL COPVCD(LUSCR1,LUHC,VEC1,1,LBLK)
      END IF
*.
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' results from H0TV2 '
        WRITE(6,*) ' ==================='
        write(6,*) ' SC0, SSIGMA0 ', SC0,SSIGMA0
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input vector ' 
        CALL WRTVCD(VEC1,LUC,1,LBLK)
        WRITE(6,*)
        WRITE(6,*) ' Output vector ' 
        CALL WRTVCD(VEC1,LUHC,1,LBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,
     &                    IOCA,IOCB,NGAS,MXPNGAS,INTH0SPC,NELFTP)
*
* Set up INTH0SPC : Division of CI space, so only determinants
*                   belonging to the same space  have nonvanishing 
*                   matrix elements of H0
*
* =====
* Input
* =====
*
* IH0SPC : ne. 0 : Interacting subspaces have been defined 
*          .eq.0 : Interacting subspaces not defined, let 
*                  evrything interact
* NPTSPC : Number of subspaces defined
* IOCPTSPC : Allowed occumulated occupation of each subspace
* NOCTPA :  Number of alpha occupation types
* NOCTPB : Number of beta occupation types
* IOCA : Occupation  of alpha string
* IOCB : Occupation  of beta string
*
* Jeppe Olsen, January 1996
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      DIMENSION IOCPTSPC(2,MXPNGAS,*)
      DIMENSION IOCA(MXPNGAS,*),IOCB(MXPNGAS,*)
      DIMENSION NELFTP(*)
*. Output 
      DIMENSION INTH0SPC(NOCTPA,NOCTPB)
*
      IF(IH0SPC.EQ.0) THEN
*. All interactions allowed 
        IONE = 1
        CALL ISETVC(INTH0SPC,IONE,NOCTPA*NOCTPB)
      ELSE
*. Explicit construction of matrix giving partitionning of 
*  subspaces
        IZERO = 0
        CALL ISETVC(INTH0SPC,IZERO,NOCTPA*NOCTPB)
*
        DO ISPC = 1, NPTSPC
          DO IATP = 1, NOCTPA
            DO IBTP = 1, NOCTPB  
              IAMOKAY = 1
              IEL = 0
C?            WRITE(6,*) ' ISPC IATP IBTP ', ISPC,IATP,IBTP
              DO IGAS = 1, NGAS
               IEL = IEL
     &             + NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
C?             WRITE(6,*) ' IGAS IEL ', IGAS,IEL
C?             WRITE(6,*)
C?   &          ' Limits :',IOCPTSPC(1,IGAS,ISPC),IOCPTSPC(2,IGAS,ISPC)
               IF(IEL.LT.IOCPTSPC(1,IGAS,ISPC).OR.
     &            IEL.GT.IOCPTSPC(2,IGAS,ISPC)    ) IAMOKAY = 0
              END DO
C?            WRITE(6,*) ' IAMOKAY = ', IAMOKAY
*. Allowed 
              IF(IAMOKAY.EQ.1.AND.INTH0SPC(IATP,IBTP).EQ.0) THEN
                INTH0SPC(IATP,IBTP) = ISPC
              END IF
            END DO
          END DO
        END DO
      END IF
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Output from  H0INTSPC '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        CALL IWRTMA(INTH0SPC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END  
      SUBROUTINE GET_BATCH_OF_SECOND
     &           (N2BLOCK,I2OFF,N2ELMNT,I2BLOCK,I2SM,I2SPC,
     &            N1BATCH,N1BLOCK,L1BLOCK,I1BLOCK,I1SM,I1SPC,
     &            N0BATCH,N0BLOCK,I0BLOCK,I0SPC,
     &            E0,E1,E2,SIN,LU0,LUSCR,IDC,
     &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &            IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,HAPR01T,LU1,
     &            N1BLOCKT)
*
* Construct a part (a batch of blocks ) of the second order correction 
*
*  |2> = (H0-E0) -1 {(E1-V)|1> + E2 * |0>}  - S02 |0>
*
* If LU1 .gt. 0, then it is assumed that |1> has been constructed and is 
* stored on LU1
*
* Jeppe Olsen, October 1995
*              Feb 96 : LU1 added
*
*. Note : The terms S02 |0> is  not included here, since S02 is in genral not
*. Known 
*. Usually our H0 has the form
*
* H0 = Q H apr Q + E0 P
* so the second order correction can be written
*
* |2> = Q (H apr -E0 ) -1 {(E1-H + Hapr)|1> + ( E2 - <0|H(apr)|1>) * |0>} 
*
*     
*
* =======
*  Input  
* =======
*         N2BLOCK  : Number of blocks to be calculated
*         I2OFF    : Absolute number of first block to be calculated 
*         N2ELMNT  : Number of ELEMENTS to be calculated
*         I2BLOCK  : List of all blocks in |2>
*         I2SM     : Symmetry of |2>
*         I2SPC    : Space of |2>
*         N1BATCH  : Number of batches in |1>
*         N1BLOCK  : Number of blocks in each batch of |1>
*         L1BLOCK  : Number of elements in each batch of |1>
*         I1BLOCK  : Blocks in |1>
*         E0       : reference energy
*         SIN      : assumed overlap between |1> and |0>
*         LU0      : file containing 0 
*         IDC      : packing in use
*         IOIO0    : allowed combination of alpha and beta supergroups for |0>
*         NOCTPA,NOCTPB : Number of alpha and beta types 
*         NSMST    : Number of symmetries of strings 
*         NSASO    : Number of alpha strings per sym and occupation type
*         NSBSO    : Number of beta  strings per sym and occupation type
*         N0BLOCK  : Total number of blocks in |0>
*         I0BLOCK  : Blocks in |0>
*         E1,E2    : First and second order energy corrections
*         LUSCR    : Scratch file
*         LU1      : Nonvanishing indicates that |1> is stored on LU1
*         N1BLOCKT : Total number of blocks in |1> 
*
* ======
* Output 
* ======
*
*
*     XB : Contains blocks of first order correction 
*     SOUT : Overlap between |2> as above  and |0>
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IOIO0(NOCTPA,NOCTPB)
*. BLock structure of |2>
      DIMENSION I2BLOCK(8,*)
*. Block structure of |1>
      DIMENSION I1BLOCK(8,*)
      DIMENSION N1BLOCK(*), L1BLOCK(*)
*. Block structure of |0>
      DIMENSION I0BLOCK(8,*)
*. Scratch
      DIMENSION IXBLOCK(8,*)
*
      DIMENSION SB(*),CB(*),XB(*)
* For communcating with MV7
      INCLUDE 'cands.inc'
      INCLUDE 'oper.inc'
C     COMMON/OPER/I12,IPERTOP,IAPR,MNRS1E,MXRS3E,IPART

*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_SE')
*
      WRITE(6,*) ' KVEC3 properly defined 3 block is used ' 
      stop'new implementation of kvec3'
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Information from  GET_SECOND_ORDER_CORRECTION '
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks to be constructed ', N2BLOCK
        WRITE(6,*) ' TTSS for Blocks to be constructed '
        DO IBLOCK = I2OFF, I2OFF-1+N2BLOCK
          WRITE(6,'(10X,4I3,2I8)') (I2BLOCK(II,IBLOCK),II=1,4)
        END DO
        CALL GFLUSH(6)
      END IF 
      IF(NTEST.GE.1)
     & WRITE(6,*) ' Entering GET_BATCH_OF_SECOND'
*
* ====================
* 1 : (V - E1 ) | 1 > 
* ====================
*
*. Loop over partionings of |1>
      DO NN1BATCH = 1, N1BATCH
        IF(NN1BATCH.EQ.1) THEN
          I1OFF = 1
        ELSE
          I1OFF = I1OFF + NN1BLOCK
        END IF
        NN1BLOCK = N1BLOCK(NN1BATCH)
        LL1BLOCK = L1BLOCK(NN1BATCH)
*. Obtain this batch of first order corrections
C       INSPC = I0SPC
C       IUTSPC = I1SPC
*
        ICSPC = I0SPC 
        ISSPC = I1SPC 
*
        IF(NTEST.GE.200) write(6,*) ' NN1BATCH I1OFF NN1BLOCK',
     &               NN1BATCH,I1OFF,NN1BLOCK
        CALL GFLUSH(6)
        CALL GET_BATCH_OF_FIRST
     &       (NN1BLOCK,I1OFF,LL1BLOCK,I1BLOCK,I1SM,I1SPC,
     &        E0,E1,SIN,LU0,IDC,
     &        IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &        N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2, 
     &        SB,XB,CB,SOUT2,
     &        HAPR01,LU1,N1BLOCKT)
*. First order correction returned in CB
        IF(NTEST.GE.5) WRITE(6,*)
     &   ' Returning from GET_BATCH_OF_FIRST'
* (E1-V)!|1> = (E1-H+H0+E1)|1>
*
*.1.1 : H |1(batch)>
*
*
         ICSPC = I1SPC
         ISSPC = I2SPC
*
         IF(NN1BATCH.EQ.1)  THEN
           ZERO = 0.0D0
           CALL SETVEC(SB,ZERO,N2ELMNT)
         ELSE
*. retrieve previous contributions to (V-E1) |1>
           CALL REWINO(LUSCR)
           CALL FRMDSC(SB,N2ELMNT,-1,LUSCR,IMZERO,IAMPACK)
         END IF
         CALL COPVEC(CB,XB,LL1BLOCK)
C?       WRITE(6,*) ' GET_SEC : Just before SBATCH2 '
         CALL GFLUSH(6)                                    
         I12 = 2
         IRESTR = 0
         CALL SBATCH2_FROM_CBATCH(N2BLOCK,I2OFF,NN1BLOCK,I1OFF,
     &                            I2BLOCK,I1BLOCK,XB,SB,IRESTR)               
C?       WRITE(6,*) ' GET_SEC : Home form SBATCH2 '
         CALL GFLUSH(6)                                    
*. SB contains now previous contributions to (V-E1) |1>
*. and H |1(batch)>, save it !
         CALL REWINO(LUSCR)
         CALL TODSC(SB,N2ELMNT,-1,LUSCR)
*
*
*. 1.2 : H0+E1 |1(batch)>
*
* obtain blocks in current batch of |1> this is in the space of current
* batch of |2>.
C                GET_TTS_BLK_IN_VECTOR
C    &           (NBLOCKI,IBLOCKI,IOFFI,VECI,
C    &            NBLOCKO,IBLOCKO,IOFFO,
C    &            NBLOCKIO,IBLCOKIO,VECIO )
        CALL GET_TTS_BLK_IN_VECTOR
     &          (NN1BLOCK,I1BLOCK,I1OFF,CB,
     &          N2BLOCK,I2BLOCK,I2OFF,
     &          N12BLOCK,L12BLOCK,IXBLOCK,XB)
         IF(NTEST.GE.5)
     &   WRITE(6,*) ' GET___SEC : Home form GET_TTS '
* 
        FACTOR =  E1
        ITASK = 2
        ECOREX = 0.0D0
        JPERT = 1
        CALL DIATERM_GAS(FACTOR,ITASK,XB,N12BLOCK,
     &       IXBLOCK,1,JPERT,0,0)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H0+E1 |1 > ' 
          CALL WRTTTS(XB,I2BLOCK(1,I2OFF),N2BLOCK,
     &                NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
        END IF
*. XB contains now relevant block of (H0+E1), subtract from prev + H|1(batch)>
        ONEM = -1.0D0
        CALL ADD_TTS_BLK_TO_VECTOR
     &       (N12BLOCK,IXBLOCK,1,XB,
     &        N2BLOCK,I2BLOCK,I2OFF,SB,ONEM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' (V - E1 ) | 1 > '
          CALL WRTTTS(SB,I2BLOCK(1,I2OFF),N2BLOCK,
     &                NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
        END IF
*. SB contains now sum(i1batch = 1, n1batch) (V-E1)|1(ibatch)>, save it
         CALL REWINO(LUSCR)
         CALL TODSC(SB,N2ELMNT,-1,LUSCR)
      END DO
*. It was actually (E1-V)|1> we were interested in so 
      CALL SCALVE(SB,ONEM,N2ELMNT)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Final (E1 - V ) | 1 > '
        CALL WRTTTS(SB,I2BLOCK(1,I2OFF),N2BLOCK,
     &                NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
*
*==================================
* 2. Add (E2 - <0|H(apr)|1>) * |0>
*==================================
*
*. Blocks of |0> in current batch of |2>
      CALL REDBLK(N2BLOCK,I2BLOCK,I2OFF,IOIO0,
     &            NOCTPA,NOCTPB,N20BLOCK,IXBLOCK,IX1)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from REDBLK '
*.Fetch the corresponding C blocks. 
      IREW = 1   
      ISCAL=0
      CALL GET_TTS_BATCH(CB,N20BLOCK,IXBLOCK,N0BLOCK,I0BLOCK,
     &                   NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                   IDC,LU0,IX2,IREW,ISCAL)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from GET_TTS'
      FACTOR = E2 - HAPR01T
C?    WRITE(6,*) ' GET_SEC : FACTOR ', FACTOR
      CALL ADDBLKV(SB,CB,FACTOR,
     &            N20BLOCK,IX1,I2BLOCK,I2OFF)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' (E1 - V) | 1 > + E2 | 0 >'
        CALL WRTTTS(SB,I2BLOCK(1,I2OFF),N2BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from ADDBLKV'
*
* ==============================================
* 3 : Multiply with (H0-E0) to give 
*     |2> = (H0-E0)-1 ((E1-V)|1> + E2 * |0>)
* ==============================================
      CALL COPVEC(SB,XB,N2ELMNT)
      FACTOR = - E0 
      ECOREX = 0.0D0
      JPERT = 1
      CALL DIATERM_GAS(FACTOR,1,XB,N2BLOCK,
     &            I2BLOCK,I2OFF,JPERT,0,0)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from DIATERM'
*. overlap between |2> and |0>
      I2EQ1 = 0
      SOUT = SCALAR_PRODUCT_OF_B
     &       (XB,CB,N2BLOCK,I2BLOCK,I2EQ1,N20BLOCK,IX1,I2OFF)
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' Batch of second order correction blocks '
        CALL WRTTTS(XB,I2BLOCK(1,I2OFF),N2BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
        
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_SE')
      IF(NTEST.GE.1)
     & WRITE(6,*) ' Leaving GET_BATCH_OF_SECOND'
      RETURN
      END 
      SUBROUTINE SBATCH2_FROM_CBATCH(NLBLOCK,ILOFF,NRBLOCK,IROFF,
     &                               ILBLOCK,IRBLOCK,
     &                               CB,HCB,IRESTRICT)
*
* Generate the parts of H times right batch that is in a given 
* left batch 
*
* The complete set of blocks of the vector is specified by IBLOCK
*
* The right batch is specified by IROFF,NRBLOCK
* The left  batch is specified by ILOFF,NLBLOCK
*
* the input right batch is assumed to reside in CB 
*
* The output blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order 
*
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      INCLUDE 'cands.inc'
*. Sigma blocks require
      INTEGER ILBLOCK(8,*)
*. C blocks included
      INTEGER IRBLOCK(8,*)
*
*./ORBINP/ : NACOB used
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
*. Used : NSMOB
      INCLUDE 'lucinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
*
      INCLUDE 'csmprd.inc'
c      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
c      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
c     &              ADSXA(MXPOBS,2*MXPOBS),
c     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
      NTEST = 0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBFCB ')
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1              
      IBTP = 2             
*. alpha and beta strings with an electron removed
      IATPM1 = 3 
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5 
      IBTPM2 = 6
*
*. Number of supergroups 
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups 
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')   !done
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')   !done
      CALL STSTSM(dbl_mb(KSTSTS),dbl_mb(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR           
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB,MXSTBL0)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
C?    WRITE(6,*) ' MXTSOB = ', MXTSOB
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSOOB_AS
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(IPRCIX.GE.1)
     &WRITE(6,*) ' SBATCH  : LCSBLK = ', LCSBLK
      LSCR1 = MAX(LCSBLK,LSCR1)
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1

*.vectors able to hold strings of given sym and type
C     MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
C     LSCR3 = MAX(MAXIK*MXTSOB*MXTSOB,MXSTBL0)
C     CALL MEMMAN(KI1,LSCR3       ,'ADDL  ',1,'I1    ')
C     CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
*
C     CALL MEMMAN(KI2,LSCR3       ,'ADDL  ',1,'I2    ')
C     CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
*
C     CALL MEMMAN(KI3,MAXIK*MXTSOB,'ADDL  ',1,'I3    ')
C     CALL MEMMAN(KXI3S,MAXIK*MXTSOB,'ADDL  ',2,'XI3S  ')
*
C     CALL MEMMAN(KI4,MAXIK*MXTSOB,'ADDL  ',1,'I4    ')
C     CALL MEMMAN(KXI4S,MAXIK*MXTSOB,'ADDL  ',2,'XI4S  ')
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')   !done
*. Arrays giving allowed type combinations 
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'CIOIO ')
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'SIOIO ')
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
C          IAIBCM(ICISPC,IAIB)
      CALL IAIBCM(ISSPC,int_mb(KSIOIO))
      CALL IAIBCM(ICSPC,int_mb(KCIOIO))
*. Arrays giving block type
COLD  CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',1,'SVST  ')  !done
        CALL SIGVST(int_mb(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESC(int_mb(KSIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
      IF(IPRCIX.GE.2) THEN
        WRITE(6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
        WRITE(6,*) 'SBLOCK : MXADKBLK ', MXADKBLK
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.2)  WRITE(6,*) ' LSCR2 = ', LSCR2
*
      LSCR12 = MAX(LSCR1,2*LSCR2)  
      CALL MEMMAN(KC2,LSCR12,'ADDL  ',2,'KC2   ')   !done
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')   !done
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')   !done
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')   !done
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')   !done
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')   !done
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')   !done
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')   !done
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')   !done
*
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,int_mb(KCBLTP),int_mb(KSVST))
*.Some TTS arrays 
      NOOS = NSMCI*NOCTPA*NOCTPB   
      NTTS = MXNTTS
*
      CALL MEMMAN(KOOS1,NOOS,'ADDL  ',2,'OOS1  ')   !done
      CALL MEMMAN(KOOS2,NOOS,'ADDL  ',2,'OOS2  ')   !done
      CALL MEMMAN(KOOS3,NOOS,'ADDL  ',2,'OOS3  ')   !done
      CALL MEMMAN(KOOS4,NOOS,'ADDL  ',2,'OOS4  ')   !done
      CALL MEMMAN(KOOS5,NOOS,'ADDL  ',2,'OOS5  ')   !done
*. Space for four blocks of string occupations and arrays of 
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      DO I1234 = 1, 4
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')         
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*. Place perturbation integrals over one body integrals
      IF(IPERTOP.NE.0) then !CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        ktmp=kint1
        kint1 = kfi
        kfi = ktmp
      end if
*. Scale and reformat
      IF(IDC.EQ.2) THEN
*. scale
        CALL SCDTTS(CB,IRBLOCK(1,IROFF),NRBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,2,NTEST)
        CALL SCDTTS(HCB,ILBLOCK(1,ILOFF),NLBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,2,NTEST)
*. reform
        CALL RFTTS(CB,dbl_mb(KC2),IRBLOCK(1,IROFF),NRBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,2,NTEST)
        CALL RFTTS(HCB,dbl_mb(KC2),ILBLOCK(1,ILOFF),NLBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,2,NTEST)
      END IF
*
      CALL SBLOCKS2(NLBLOCK,ILBLOCK(1,ILOFF),NRBLOCK,IRBLOCK(1,IROFF),
     &            CB,HCB,dbl_mb(KC2),
     &            int_mb(KCIOIO),ISMOST(1,ICSM),int_mb(KCBLTP),
     &            NACOB,int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &            NAEL,IATP,NBEL,IBTP,
     &            IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &            ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
     &            LSCR1,dbl_mb(KINSCR),dbl_mb(KCSCR),dbl_mb(KSSCR),
     &            SXSTSM,dbl_mb(KSTSTS),dbl_mb(KSTSTD),SXDXSX,
     &            ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &            dbl_mb(KOOS1),dbl_mb(KOOS2),dbl_mb(KOOS3),
     $            dbl_mb(KOOS4),
     &            dbl_mb(KOOS5),int_mb(KI1),dbl_mb(KXI1S),
     &            int_mb(KI2),dbl_mb(KXI2S),IDOH2,MXPOBS,int_mb(KSVST),
     &            PSSIGN,IPRDIA,LUC,ICJKAIB,dbl_mb(KCJRES),
     &            dbl_mb(KSIRES),int_mb(KI3),dbl_mb(KXI3S),
     &            int_mb(KI4),dbl_mb(KXI4S),MXSXST,MXSXBL,
     &            MOCAA,MOCAB,IAPR,IRESTRICT)
*
      IF(IDC.EQ.2) THEN    
*. reform 
        CALL RFTTS(HCB,CB,ILBLOCK(1,ILOFF),NLBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,ILBLOCK(1,ILOFF),NLBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF

*. restore order
      IF(IPERTOP.NE.0) then !CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        ktmp=kint1
        kint1 = kfi
        kfi = ktmp
      end if
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SBFCB ')
      RETURN
      END
      SUBROUTINE SBATCH_FROM_CBATCH(NLBLOCK,ILOFF,NRBLOCK,IROFF,IBLOCK,
     &                               CB,HCB,IRESTRICT)
*
* Generate the parts of H times right batch that is in a given 
* left batch 
*
* The complete set of blocks of the vector is specified by IBLOCK
*
* The right batch is specified by IROFF,NRBLOCK
* The left  batch is specified by ILOFF,NLBLOCK
*
* the input right batch is assumed to reside in CB 
*
* The output blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order 
*
* If Irestr.eq.1 then we are calculating :
* sigma(IBLK) = SUM(JBLK.LE.IBLK) (2-delta(iblk,jblk))/2 <IBLK!H!JBLK>C(JBLK)
*
* The blocks are scaled and the summation over JBLKS is restricted.
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      INCLUDE 'cands.inc'
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
*./ORBINP/ : NACOB used
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
*. Used : NSMOB
      INCLUDE 'lucinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
*
      INCLUDE 'csmprd.inc'
c      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
c      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
c     &              ADSXA(MXPOBS,2*MXPOBS),
c     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBFCB ')
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1              
      IBTP = 2             
*. alpha and beta strings with an electron removed
      IATPM1 = 3 
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5 
      IBTPM2 = 6
*
*. Number of supergroups 
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups 
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')  !done
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')  !done
      CALL STSTSM(dbl_mb(KSTSTS),dbl_mb(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR           
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
C?    WRITE(6,*) ' MXTSOB = ', MXTSOB
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSOOB_AS
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(IPRCIX.GE.1)
     &WRITE(6,*) ' MV7 : LCSBLK = ', LCSBLK
      LSCR1 = MAX(LCSBLK,LSCR1)
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
*.vectors able to hold strings of given sym and type
C     MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
C     LSCR3 = MAX(MAXIK*MXTSOB*MXTSOB,MXSTBL0)
C     CALL MEMMAN(KI1,LSCR3       ,'ADDL  ',1,'I1    ')
C     CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
*
C     CALL MEMMAN(KI2,LSCR3       ,'ADDL  ',1,'I2    ')
C     CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
*
C     CALL MEMMAN(KI3,MAXIK*MXTSOB,'ADDL  ',1,'I3    ')
C     CALL MEMMAN(KXI3S,MAXIK*MXTSOB,'ADDL  ',2,'XI3S  ')
*
C     CALL MEMMAN(KI4,MAXIK*MXTSOB,'ADDL  ',1,'I4    ')
C     CALL MEMMAN(KXI4S,MAXIK*MXTSOB,'ADDL  ',2,'XI4S  ')
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')   !done
*. Arrays giving allowed type combinations 
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'CIOIO ')
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'SIOIO ')
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
C          IAIBCM(ICISPC,IAIB)
      CALL IAIBCM(ISSPC,int_mb(KSIOIO))
      CALL IAIBCM(ICSPC,int_mb(KCIOIO))
*. Arrays giving block type
COLD  CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',1,'SVST  ')  !done
        CALL SIGVST(int_mb(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESC(int_mb(KSIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
      IF(IPRCIX.GE.2) THEN
        WRITE(6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
        WRITE(6,*) 'SBLOCK : MXADKBLK ', MXADKBLK 
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.2)  WRITE(6,*) ' LSCR2 = ', LSCR2
*
      LSCR12 = MAX(LSCR1,2*LSCR2)  
      CALL MEMMAN(KC2,LSCR12,'ADDL  ',2,'KC2   ')   !done
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')  !done
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')  !done
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')  !done
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')  !done
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')  !done
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')  !done
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')  !done
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')  !done
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,int_mb(KCBLTP),int_mb(KSVST))
*.Some TTS arrays 
      NOOS = NSMCI*NOCTPA*NOCTPB
      NTTS = MXNTTS
*
      CALL MEMMAN(KOOS1,NOOS,'ADDL  ',2,'OOS1  ')  !done
      CALL MEMMAN(KOOS2,NOOS,'ADDL  ',2,'OOS2  ')  !done
      CALL MEMMAN(KOOS3,NOOS,'ADDL  ',2,'OOS3  ')  !done
      CALL MEMMAN(KOOS4,NOOS,'ADDL  ',2,'OOS4  ')  !done
      CALL MEMMAN(KOOS5,NOOS,'ADDL  ',2,'OOS5  ')  !done
*. Space for four blocks of string occupations and arrays of 
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      DO I1234 = 1, 4
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')         
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*. Place perturbation integrals over one body integrals
      IF(IPERTOP.NE.0) then !CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        ktmp=kint1
        kint1 = kfi
        kfi = ktmp
      end if
c     IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
*. Scale and reformat
      NTEST = 0
      IF(IDC.EQ.2) THEN    
*. scale
        CALL SCDTTS(CB,IBLOCK(1,IROFF),NRBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,2,NTEST)
        CALL SCDTTS(HCB,IBLOCK(1,ILOFF),NLBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,2,NTEST)
*. reform 
        CALL RFTTS(CB,dbl_mb(KC2),IBLOCK(1,IROFF),NRBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,2,NTEST)
        CALL RFTTS(HCB,dbl_mb(KC2),IBLOCK(1,ILOFF),NLBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,2,NTEST)
      END IF
C?    WRITE(6,*) ' HCB before call to SBLOCKS2'
C?      CALL WRTTTS(HCB,IBLOCK(1,ILOFF),NLBLOCK,
C?   &              NSMST,NOCTPA,NOCTPB,
C?   &               WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),1)
C     RFTTS(BLOCKSI,BLOCKSO,IBLOCK,NBLOCK,
C    &                 ICOPY,NSMST,NOCTPA,NOCTPB,
C    &                 NSASO,NSBSO,IDC,PS,IWAY,IPRNT)
C     SCDTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
*
      CALL SBLOCKS2(NLBLOCK,IBLOCK(1,ILOFF),NRBLOCK,IBLOCK(1,IROFF),
     &            CB,HCB,dbl_mb(KC2),
     &            int_mb(KCIOIO),ISMOST(1,ICSM),int_mb(KCBLTP),
     &            NACOB,int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &            NAEL,IATP,NBEL,IBTP,
     &            IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &            ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
     &            LSCR1,dbl_mb(KINSCR),dbl_mb(KCSCR),dbl_mb(KSSCR),
     &            SXSTSM,dbl_mb(KSTSTS),dbl_mb(KSTSTD),SXDXSX,
     &            ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &            dbl_mb(KOOS1),dbl_mb(KOOS2),dbl_mb(KOOS3),
     &            dbl_mb(KOOS4),
     &            dbl_mb(KOOS5),int_mb(KI1),dbl_mb(KXI1S),
     &            int_mb(KI2),dbl_mb(KXI2S),IDOH2,MXPOBS,int_mb(KSVST),
     &            PSSIGN,IPRDIA,LUC,ICJKAIB,dbl_mb(KCJRES),
     &            dbl_mb(KSIRES),int_mb(KI3),dbl_mb(KXI3S),
     &            int_mb(KI4),dbl_mb(KXI4S),MXSXST,MXSXBL,
     &            MOCAA,MOCAB,IAPR,IRESTRICT)
*
      IF(IDC.EQ.2) THEN    
*. scale
C       CALL SCDTTS(HCB,IBLOCK,NBLOCK,NSMST,NOCTPA,NOCTPB,
C    &              WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
C    &              IDC,1,NTEST)
*. reform 
        CALL RFTTS(HCB,dbl_mb(KC2),IBLOCK(1,ILOFF),NLBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,ILOFF),NLBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF
C     RFTTS(BLOCKSI,BLOCKSO,IBLOCK,NBLOCK,
C    &                 ICOPY,NSMST,NOCTPA,NOCTPB,
C    &                 NSASO,NSBSO,IDC,PS,IWAY,IPRNT)
C     SCDTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,IDC,IWAY,IPRNT)

*. restore order
      IF(IPERTOP.NE.0) then !CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        ktmp=kint1
        kint1 = kfi
        kfi = ktmp
      end if
c     IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SBFCB ')
      RETURN
      END
      SUBROUTINE SBLOCKS2(NSBLOCK,ISBLOCK,NCBLOCK,ICBLOCK,CB,SB,C2,
     &                  ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA,NSSOB,
     &                  NAEL,IAGRP,NBEL,IBGRP,
     &                  IOCTPA,IOCTPB, NOCTPA,NOCTPB,
     &                  NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &                  ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &                  XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &                  SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC, 
     &                  ICOOSC,NCOOSC,ICOOSE,NCOOSE,IACOOS,
     &                  I1,XI1S,I2,XI2S,IDOH2,MXPOBS,ISTRFL,
     &                  PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &                  I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,IRESTR)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen , October 1995    
*
* Note : Sigma blocks are not zeroed !!
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block      
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type 
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If Irestr.ne. 0 then the sigma blocks to be 
* calculated are sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)
*. The C and sigma blocks are on input assume to be in determinant form 

      IMPLICIT REAL*8(A-H,O-Z)
*. Specific input : Blocks of sigma and C 
      INTEGER ISBLOCK(8,NSBLOCK)
      INTEGER ICBLOCK(8,NCBLOCK)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      DIMENSION ICOOSC(NOCTPA,NOCTPB,NSMST),NCOOSC(NOCTPA,NOCTPB,NSMST)
      DIMENSION ICOOSE(NOCTPA,NOCTPB,NSMST),NCOOSE(NOCTPA,NOCTPB,NSMST)
      DIMENSION IACOOS(NOCTPA,NOCTPB,NSMST)
      DIMENSION I1(*),I2(*),I3(*),XI1S(*),XI2S(*),XI3S(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
*.
      REAL * 8 INPROD
*
      CALL LUCIAQENTER('SBLO2')
      NTEST =  00
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' SBLOCKS2 speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)  
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS for each sigma block'
          DO IBLOCK = 1, NSBLOCK                                        
            WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          END DO
        WRITE(6,*) ' Number of C blocks included ', NCBLOCK
        WRITE(6,*) ' TTSS for each C block '
          DO IBLOCK = 1, NCBLOCK
            WRITE(6,'(10X,4I3,2I8)') (ICBLOCK(II,IBLOCK),II=1,4)
          END DO
          WRITE(6,*) ' IDC = ', IDC 
      END IF
*
C?    WRITE(6,*) ' LUC = ', LUC
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial blocks of C vector '
        CALL WRTTTS(CB,ICBLOCK,NCBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,1)
      END IF
* We are working with transposed sigma blocks , so 
      DO JSBLOCK = 1, NSBLOCK
        IATP = ISBLOCK(1,JSBLOCK)
        IBTP = ISBLOCK(2,JSBLOCK)
        IASM = ISBLOCK(3,JSBLOCK)
        IBSM = ISBLOCK(4,JSBLOCK)
        ISOFF = ISBLOCK(5,JSBLOCK)
        NASTR = NSSOA(IASM,IATP)
        NBSTR = NSSOA(IBSM,IBTP)
        CALL TRPMT3(SB(ISOFF),NASTR,NBSTR,C2)   
        CALL COPVEC(C2,SB(ISOFF),NASTR*NBSTR)
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP)THEN
          HALF = 0.5D0
          CALL SCALVE(SB(ISOFF),HALF,NASTR*NBSTR)
        END IF
      END DO
      IF(NTEST.GE.50) WRITE(6,*) ' S blocks zeroed '
*. Transpose C blocks ( for easier accesss in ab loop 
      IF(ICJKAIB.NE.0) THEN
        DO JCBLOCK = 1, NCBLOCK            
          JATP = ICBLOCK(1,JCBLOCK)
          JBTP = ICBLOCK(2,JCBLOCK)
          JASM = ICBLOCK(3,JCBLOCK)
          JBSM = ICBLOCK(4,JCBLOCK)
          JCOFF =ICBLOCK(5,JCBLOCK)
C?        write(6,*) ' JCBLOCK JCOFF ',JCBLOCK,JCOFF
*. Transposed C block is used so :
          CALL TRPMT3(CB(JCOFF),NSSOA(JASM,JATP),  
     &         NSSOB(JBSM,JBTP),C2)   
          CALL COPVEC(C2,CB(JCOFF),
     &         NSSOA(JASM,JATP)*NSSOB(JBSM,JBTP))    
        END DO
      END IF
      IF(NTEST.GE.50) WRITE(6,*) ' C blocks transposed '
*. Loop over sigma and C blocks in core and obtain  contribution from
* given C block to given S block
      DO 10000 ISBLK = 1, NSBLOCK
        IATP = ISBLOCK(1,ISBLK)
        IBTP = ISBLOCK(2,ISBLK)
        IASM = ISBLOCK(3,ISBLK)
        IBSM = ISBLOCK(4,ISBLK)
        ISOFF = ISBLOCK(5,ISBLK)
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Next s block in batch : '
          write(6,*) ' ISBLK IASM IBSM IATP IBTP'
          write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
        END IF
        IF(NIA.NE.0.AND.NIB.NE.0) THEN
*. Loop over C blocks
          DO 9000 ICBLK = 1, NCBLOCK                                     
            JATP = ICBLOCK(1,ICBLK)
            JBTP = ICBLOCK(2,ICBLK)
            JASM = ICBLOCK(3,ICBLK)
            JBSM = ICBLOCK(4,ICBLK)
            ICOFF = ICBLOCK(5,ICBLK)
            NJA = NSSOA(JASM,JATP)
            NJB = NSSOB(JBSM,JBTP)
            IF(IRESTR.EQ.1.AND.
     &         (JASM.GT.IASM.OR.
     &         JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &         JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &           GOTO 9000   
            XNORM2 = INPROD(CB(ICOFF),CB(ICOFF),NJA*NJB,NJA*NJB)
C!          XNORM2 = 1.0D0
            IF(XNORM2.NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
              PL = 1.0D0
              CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                    LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
              DO 8765 IPERM = 1, NPERM
                LLASM = LASM(IPERM)
                LLBSM = LBSM(IPERM)
                LLATP = LATP(IPERM)
                LLBTP = LBTP(IPERM)
                NLLA = NSSOA(LLASM,LLATP)
                NLLB = NSSOB(LLBSM,LLBTP)
                IF(LTRP(IPERM).EQ.1) THEN
                  LROW = NSSOA(LASM(IPERM-1),LATP(IPERM-1))
                  LCOL = NSSOB(LBSM(IPERM-1),LBTP(IPERM-1))
                  IF(ICJKAIB.EQ.0) THEN
                    CALL TRPMT3(CB(ICOFF),LROW,LCOL,C2)
                  ELSE
                    CALL TRPMT3(CB(ICOFF),LCOL,LROW,C2)
                  END IF
                  CALL COPVEC(C2,CB(ICOFF),LROW*LCOL)
                END IF
*Jeppe define LROW and LCOL
                IF(LSGN(IPERM).EQ.-1)
     &          CALL SCALVE(CB(ICOFF),-1.0D0,LROW*LCOL)
                IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &            (LLBSM.GT.LLASM.OR.
     &            (LLASM.EQ.LLBSM.AND.LLBTP.GT.LLATP))            )
     &            GOTO 8764
C               IF(IRESTR.EQ.1.AND.
C    &             (LLASM.GT.IASM.OR.
C    &              LLASM.EQ.IASM.AND.LLATP.GT.IATP.OR.
C    &              LLASM.EQ.IASM.AND.LLATP.EQ.IATP.AND.LLBTP.GT.IBTP))
C    &           GOTO 8764   
              
*
                IF(NTEST.GE.10) THEN
                  WRITE(6,*) ' RSSBCB will be called for '
                  WRITE(6,*) ' Sigma block : '
                  WRITE(6,*) ' ISOFF ', ISOFF
                  WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(6,*) ' C     block : '
                  WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                  WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                  WRITE(6,*) ' ICOFF ', ICOFF
                END IF
C?              write(6,*) ' Memory check before RSSBCB'
C?              call memchk
C?              write(6,*) ' checked passed '
C               IF(IRESTR.EQ.1.AND.IASM.EQ.LLASM.AND.
C    &                             IBSM.EQ.LLBSM.AND.
C    &                             IATP.EQ.LLATP.AND.
C    &                             IBTP.EQ.LLBTP) THEN
*
                IF(IRESTR.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     TWO = 2.0D0
                     CALL SCALVE(SB(ISOFF),TWO,NIA*NIB)
                END IF
*
                SCLFAC = 1.0D0
                CALL RSSBCB(IASM,IATP,IOCTPA,
     &               IBSM,IBTP,IOCTPB,
     &               LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &               NELFSPGP(1,IATP+IOCTPA-1),
     &               NELFSPGP(1,IBTP+IOCTPB-1),
     &               NELFSPGP(1,LLATP+IOCTPA-1),
     &               NELFSPGP(1,LLBTP+IOCTPB-1),
     &               NAEL,NBEL,
     &               IAGRP,IBGRP,
     &               SB(ISOFF),CB(ICOFF),IDOH2,
     &               ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &               NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &               NSMOB,NSMST,NSMSX,NSMDX,
     &               NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &               ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &               MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,SCLFAC,IPRNT)
*
C               IF(IRESTR.EQ.1.AND.IASM.EQ.LLASM.AND.
C    &                             IBSM.EQ.LLBSM.AND.
C    &                             IATP.EQ.LLATP.AND.
C    &                             IBTP.EQ.LLBTP) THEN
                IF(IRESTR.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     HALF= 0.5D0
                     CALL SCALVE(SB(ISOFF),HALF,NIA*NIB)
                END IF
C?              write(6,*) ' Memory check after RSSBCB'
C?              call memchk
C?              write(6,*) ' checked passed '
                
 8764         CONTINUE
 8765         CONTINUE
*. Transpose or scale to restore order in c 
              IF(LTRP(NPERM+1).EQ.1) THEN
                IF(ICJKAIB.EQ.0) THEN
                  CALL TRPMT3(CB(ICOFF),NJB,NJA,C2)  
                ELSE
                  CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)  
                END IF
                CALL COPVEC(C2,CB(ICOFF),NJA*NJB)  
              END IF
              IF(LSGN(NPERM+1).EQ.-1)
     &        CALL SCALVE(CB(ICOFF),-1.0D0,NJA*NJB)  
*
            END IF
 9000       CONTINUE
          END IF
*. End of loop over C blocks in Batch
10000     CONTINUE
*. End of loop over S blocks in batch
*
* Order 
      DO  ISBLK = 1 , NSBLOCK
        IATP = ISBLOCK(1,ISBLK)
        IBTP = ISBLOCK(2,ISBLK)
        IASM = ISBLOCK(3,ISBLK)
        IBSM = ISBLOCK(4,ISBLK)
        ISOFF  = ISBLOCK(5,ISBLK)
        ISOFFP = ISBLOCK(6,ISBLK)
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
        IF(ICJKAIB.NE.0) THEN
*. Tranpose sigma block was obtained, transpose to obtain correct block
           WRITE(6,*) ' Test: SB before TRPMT3'
           CALL WRTMAT(SB(ISOFF),NIB,NIA,NIB,NIA)
           CALL TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),
     &                 NSSOA(IASM,IATP),C2)
           CALL COPVEC(C2,SB(ISOFF),
     &                 NSSOA(IASM,IATP)* NSSOB(IBSM,IBTP))
        END IF
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) 
     &  CALL TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))
         WRITE(6,*) ' Test: SB as will be delivered'
           CALL WRTMAT(SB(ISOFF),NIA,NIB,NIA,NIB)
      END DO
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' output blocks from SBLOCKS2 '
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,1)
      END IF
*
C?    write(6,*) ' Memory check , exit from SBLOCKS2'
C?    call memchk
C?    write(6,*) ' checked passed '
      CALL LUCIAQEXIT('SBLO2')
      RETURN
      END
      FUNCTION SCALAR_PRODUCT_OF_B
     &(VEC1,VEC2,NBLOCK1,IBLOCK1,I2EQ1,NBLOCK2,I2TO1,IOFF)
*
* Two blocked vectors VEC1 and VEC2 are given.
* VEC1 contains NBLOCK1 blocks defined by IBLOCK1.
* If I2EQ1.ne.0  VEC2 has the same block structure as VEC1
* IF I2EQ1 .eq.0 VEC2 contains NBLOCK2 blocks, and I2TO1 gives 
* the mapping from blocks of VEC2 to blocks of vec1
*
* Find the scalar product between these two vectors 
*
* Jeppe Olsen, October 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION VEC1(*),VEC2(*),IBLOCK1(8,*)
      REAL*8 INPROD  
*. In use of I2EQ1 = 0
      DIMENSION I2TO1(*)
*
      X = 0.0D0
      IOFF2 = 1
      DO JBLOCK = 1, NBLOCK2
        IF(I2EQ1.NE.0) THEN
          JJBLOCK = JBLOCK-1+IOFF
          IOFF1   = IBLOCK1(6,JJBLOCK)
          IOFF2   = IOFF1
        ELSE
          JJBLOCK = I2TO1(JBLOCK)
          IOFF1 = IBLOCK1(6,JJBLOCK)
        END IF
        NELMNT = IBLOCK1(8,JJBLOCK)
        X = X + INPROD(VEC1(IOFF1),VEC2(IOFF2),NELMNT)
        IF(I2EQ1.EQ.0) IOFF2 = IOFF2 + NELMNT
      END DO
      SCALAR_PRODUCT_OF_B = X
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Output from INNER_PRODUCT_OF_BLOCKED_VECTORS'
         WRITE(6,*)
         WRITE(6,*) ' Number of blocks included ', NBLOCK2
         WRITE(6,*) ' Inner product ', X
      END IF
*
      RETURN
      END
      SUBROUTINE GET_BATCH_OF_FIRST
     &           (N1BLOCK,I1OFF,N1ELMNT,I1BLOCK,I1SM,I1SPC,
     &            E0,E1,ALPHA,LU0,IDC,
     &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &            N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,
     &            HAPR01,LU1,N1BLOCKT )
*
* Construct a part (a batch of blocks ) of the first order correction 
*
*  |1> = -(H(apr)-E0) -1 (H-E0-E1)|0> - ALPHA * (H(apr)-E0)-1|0>
*
* Obtain  also  overlap <0|H(apr)|1> (Only  LU1 .EQ.0 )
*
* Jeppe Olsen, Written : October 1995
*              Modified : December 1995
*                         February 1996 : LU1 added 
* Last Modification; Oct. 30, 2012; Jeppe Olsen; Aligned with current code
*
*
* =======
*  Input  
* =======
*         N1BLOCK  : Number of blocks to be calculated
*         I1OFF    : Absolute number of first block to be calculated 
*         N1ELMNT  : Number of ELEMENTS to be calculated
*         I1BLOCK  : List of all blocks in |1>
*         I1SM     : Symmetry of |1>
*         I1SPC    : Space of |1>
*         E0       : reference energy
*         SIN      : assumed overlap between |1> and |0>
*         LU0      : file containing 0 
*         IDC      : packing in use
*         IOIO0    : allowed combination of alpha and beta supergroups for |0>
*         NOCTPA,NOCTPB : Number of alpha and beta types 
*         NSMST   : Number of symmetries of strings 
*         NSASO   : Number of alpha strings per sym and occupation type
*         NSBSO   : Number of beta  strings per sym and occupation type
*         N0BLOCK : Total number of blocks in |0>
*         I0BLOCK : Blocks in |0>
*         LU1     : file containing |1> ( if constructed previously)
*         N1BLOCKT: Total number of blocks in |1>
*
* ======
* Output 
* ======
*
*
*     XB : Contains blocks of first order correction 
*     SB : Contains blocks of (H-E)|0>                '
*     SOUT : Overlap <0| (H(apr)-E0) -1 (H-E0)|0>  (if LU1 = 0)
*     HAPR01 : Overlap <0|H(apr)|1>                (if LU1 = 0)
* ========
*  Scratch
* ========
* IXBLOCK,IX1,IX2
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IOIO0(NOCTPA,NOCTPB)
*. Block structure of |1>
      DIMENSION I1BLOCK(8,*)
*. Block structure of |0>
      DIMENSION I0BLOCK(8,*)
*
      DIMENSION SB(*),CB(*),XB(*)
*
      INCLUDE 'oper.inc'
      INCLUDE 'cecore.inc'
      CHARACTER*6 CTASK 
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_FI')
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Information from  GET_FIRST_ORDER_CORRECTION '
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks to be constructed ', N1BLOCK
        WRITE(6,*) ' TTSS for Blocks to be constructed '
        DO IBLOCK = I1OFF, I1OFF-1+N1BLOCK
          WRITE(6,'(10X,4I3,2I8)') (I1BLOCK(II,IBLOCK),II=1,4)
        END DO
        IF(NTEST.GE.2) WRITE(6,*) ' First order space : ', I1SPC 
      END IF 
      CALL GFLUSH(6)
*
      IF(LU1.GT.0) THEN
*
*. Easy Living, first order correction is stored on LU1
* ======================================================
        
C             GET_BLOCKS_FROM_DISC(LU,NBLOCK,IOFF,IBLOCK,NBLOCKT,C)
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,N1BLOCK,I1OFF,I1BLOCK,N1BLOCKT,XB,1)
*. Blocks should be scaled if IDC.NE.1 ???????
         HAPR01 = 0.0D0
*. We are also interested in 
      ELSE
*
*. Well, not so easy living, I/WE ( maybe I am a parallel code by now)
*. have to construct the first order correction from SCRATCH
* ===================================================================
*. 
*
*
* ============================
*. Obtain sigma blocks: H |0> 
* ============================
*
      I12 = 2
      IPERTOP = 0
      IAPR = 0
C     write(6,*) ' call to SBLOCK : KVEC3 defined ?'
C     STOP'before call to SBLOCK'
      CTASK = 'SIGMA '
      CALL SBLOCK(N1BLOCK,I1BLOCK,I1OFF,CB,SB,LU0,0,0,0,0,0,0.0D0,
     &            ECORE,CTASK)
C     SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
C    &                  ECORE,ITASK) 
      IF(NTEST.GE.2000) THEN
        WRITE(6,*) ' Batch of H!0>: '
        CALL WRTTTS(SB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
C?    IF(NTEST.GE.2) write(6,*) ' Returning from Sblock '
*. Sigma Blocks that are in actual C space 
      CALL REDBLK(N1BLOCK,I1BLOCK,I1OFF,IOIO0,
     &            NOCTPA,NOCTPB,N10BLOCK,IXBLOCK,IX1)
C?    IF(NTEST.GE.10)
C?   &write(6,*) ' Returning from REDBLK '
*.Fetch the corresponding C blocks. 
      IREW = 1   
      ISCAL = 0 
      CALL GET_TTS_BATCH(CB,N10BLOCK,IXBLOCK,N0BLOCK,I0BLOCK,
     &                   NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                   IDC,LU0,IX2,IREW,ISCAL)
C?    IF(NTEST.GE.10)
C?   &write(6,*) ' Returning from GET_TTS'
* ====================================
*. Subtract E0+E1 |0> to give (H-(E0+E1)) |0>
* ====================================
      FACTOR = - (E0 + E1)
      CALL ADDBLKV(SB,CB,FACTOR,
     &            N10BLOCK,IX1,I1BLOCK,I1OFF)
C?    IF(NTEST.GE.10)
C?   &write(6,*) ' Returning from ADDBLKV'
      IF(NTEST.GE.2000) THEN
        WRITE(6,*) ' Batch of (H-(E0+E1))!0>'
        CALL WRTTTS(SB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
*
* ====================================
* obtain block of (H0-E0)-1**(H-E) !0>
* ====================================
*
      CALL COPVEC(SB,XB,N1ELMNT)
      FACTOR = - E0 
      ECOREX = 0.0D0
      JPERT = 1
      CALL DIATERM_GAS(FACTOR,1,XB,N1BLOCK,
     &            I1BLOCK,I1OFF,JPERT,0,0)
      IF(NTEST.GE.2000) THEN
        WRITE(6,*) ' Batch of (H0-E0)-1**(H-E) !0> '
        CALL WRTTTS(XB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
C?    IF(NTEST.GE.10)
C?   &write(6,*) ' Returning from DIATERM'
      ONEM = -1.0D0
      CALL SCALVE(XB,ONEM,N1ELMNT)
* =====================================
*. overlap <0|(H(apr)-E0) -1 (E0-H)|0>
* =====================================
      I2EQ1 = 0
      SOUT = SCALAR_PRODUCT_OF_B
     &       (XB,CB,N1BLOCK,I1BLOCK,I2EQ1,N10BLOCK,IX1,I1OFF)
*
* =====================================
* Add -ALPHA * (H0-E0)**-1 |0>
* =====================================
*
      IF(ALPHA.NE.0.0D0) THEN
*. (H(apr)-E0)**-1|0>
        FACTOR = -E0    
        ECOREX = 0
        JPERT = 1
        CALL DIATERM_GAS(FACTOR,1,CB,N10BLOCK,IXBLOCK,
     &                    1,JPERT,0,0)
        FACTOR = -ALPHA
        CALL ADDBLKV(XB,CB,FACTOR,N10BLOCK,IX1,I1BLOCK,I1OFF)
      END IF
*
*
* ================================
*. Contribution to <0| H(apr) |1>
* ================================
*
*. Fetch the C blocks again
      IREW = 1   
      ISCAL = 0
      CALL GET_TTS_BATCH(CB,N10BLOCK,IXBLOCK,N0BLOCK,I0BLOCK,
     &                   NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                   IDC,LU0,IX2,IREW,ISCAL)
*. H(apr)|0>
      FACTOR = 0
      ECOREX = 0
      JPERT =  1
      CALL DIATERM_GAS(FACTOR,2,CB,N10BLOCK,IXBLOCK,
     &                  1,JPERT,0,0)
      HAPR01 = SCALAR_PRODUCT_OF_B
     &       (XB,CB,N1BLOCK,I1BLOCK,I2EQ1,N10BLOCK,IX1,I1OFF)
*. 
*
      END IF
*. ( End of LU1=0/not 0 Branching)
*
*.
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' Batch of first order correction'
        CALL WRTTTS(XB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_FI')
*
      RETURN
      END 
      SUBROUTINE CLASS_PROD2(VEC1,VEC2,IBLOCK1,IMAP21,NBLOCK2,
     &                      NOCCLS,IOCCLS,CLSVEC)
*
* Two vectors in blocked form are given.
* Find contributions to each occupation class
*
* In this version VEC1 and VEC2 do not need to have identical
* block structure. Vec2 can contain only a subset of the 
* blocks in VEC1.
* IMAP21 gives mapping from the reduced vector to 
* the complete vector 
      IMPLICIT REAL*8(A-H,O-Z)
*
      REAL*8 INPROD
*. Input
      DIMENSION VEC1(*),VEC2(*)
      INTEGER IBLOCK1(8,*)
      INTEGER IMAP21(*)
*. Input/output
      DIMENSION CLSVEC(*)
      
*
      IOFF = 1
      DO JBLOCK2 = 1, NBLOCK2
        JBLOCK = IMAP21(JBLOCK2)
*
        JATP = IBLOCK1(1,JBLOCK)
        JBTP = IBLOCK1(2,JBLOCK)
        JOFF = IBLOCK1(6,JBLOCK)
        NELMNT = IBLOCK1(8,JBLOCK)      
C?    write(6,*) ' JBLOCK2 JBLOCK JATP JBTP JOFF IOFF NELMNT',
C?   &             JBLOCK2,JBLOCK,JATP,JBTP,JOFF, IOFF, NELMNT
*. Corresponding occupation class
        CALL SPGSPG_TO_CLASS(JATP,JBTP,JOCCLS,NOCCLS,IOCCLS)
        IF(JOCCLS.EQ.0) THEN
          WRITE(6,*) ' JOCCLS = 0 returned from  SPGSPG_TO_CLASS'
        ELSE 
          XTERM = INPROD(VEC1(JOFF),VEC2(IOFF),NELMNT)
C?        WRITE(6,*) ' CLASS_PROD : XTERM = ' , XTERM
          CLSVEC(JOCCLS) = CLSVEC(JOCCLS) + XTERM
        END IF
        IOFF = IOFF + NELMNT
      END DO
*
      NTEST = 0
      IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Updated CLSVEC '
         CALL WRTMAT(CLSVEC,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END  
      SUBROUTINE CLASS_PROD(VEC1,VEC2,NOCTPA,NOCTPB,
     &                      IBLOCK_OFF,NBLOCK,IBLOCK,NOCCLS,IOCCLS,
     &                      CLSVEC)
*
* Two vectors in blocked form are given.
* Find contributions to each occupation class
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      REAL*8 INPROD
*. Input
      DIMENSION VEC1(*),VEC2(*)
      INTEGER IBLOCK(8,*)
*. Input/output
      DIMENSION CLSVEC(*)
      
*
      IOFF = 1
      DO JBLOCK = IBLOCK_OFF,IBLOCK_OFF-1+NBLOCK
*
        JATP = IBLOCK(1,JBLOCK)
        JBTP = IBLOCK(2,JBLOCK)
        JASM = IBLOCK(3,JBLOCK)
        NELMNT = IBLOCK(8,JBLOCK)      
C       NELMNT = NOOS(JATP,JBTP,JASM)
*. Corresponding occupation class
        CALL SPGSPG_TO_CLASS(JATP,JBTP,JOCCLS,NOCCLS,IOCCLS)
C?      WRITE(6,*) 
C?   &  ' CLASS_PROD : JBLOCK, IOFF NELMNT ' ,JBLOCK,  IOFF, NELMNT 
        IF(JOCCLS.EQ.0) THEN
          WRITE(6,*) ' JOCCLS = 0 returned from  SPGSPG_TO_CLASS'
        ELSE 
          XTERM = INPROD(VEC1(IOFF),VEC2(IOFF),NELMNT)
C?        WRITE(6,*) ' CLASS_PROD : XTERM = ', XTERM
          CLSVEC(JOCCLS) = CLSVEC(JOCCLS) + XTERM
        END IF
        IOFF = IOFF + NELMNT
      END DO
*
      NTEST = 0
      IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Updated CLSVEC '
         CALL WRTMAT(CLSVEC,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END  
      SUBROUTINE SPGSPG_TO_CLASS(JOCTPA,JOCTPB,JOCCLS,NOCCLS,IOCCLS)
*
* Two supergroups of strings are given.
* Find the corresponding occupation class
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Absolute numbers 
      IATP = 1
      IBTP = 2
      JJOCTPA = IBSPGPFTP(IATP) - 1 + JOCTPA
      JJOCTPB = IBSPGPFTP(IBTP) - 1 + JOCTPB
*
      JOCCLS = 0
      DO KOCCLS= 1, NOCCLS
        I_AM_OKAY = 1
        DO IGAS = 1, NGAS
          IF(NELFSPGP(IGAS,JJOCTPA)+NELFSPGP(IGAS,JJOCTPB) 
     &       .NE. IOCCLS(IGAS,KOCCLS)                     )I_AM_OKAY=0
        END DO
        IF(I_AM_OKAY.EQ.1) JOCCLS = KOCCLS
      END DO
*
      IF(JOCCLS.EQ.0) THEN
        WRITE(6,*) ' Problem in SPGSPG_TO_CLASS : No match found ! '
      END IF
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info form SPGSPG_TO_CLASS'
        WRITE(6,*) ' input supergroups ',JJOCTPA,JJOCTPB
        WRITE(6,*) ' Output occupation class ', JOCCLS
      END IF
*
      RETURN
      END 
      SUBROUTINE DIA0TRM_GAS(ITASK,LUIN,LUOUT,VEC,VEC2,
     &                       FACTOR)
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 2)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 1)
*
* Note : Opposite to normal DIATRM_GAS !!!
*
*
* Jeppe Olsen, August 1995
*
* Last Modification; Nov. 1, 2012; Jeppe Olsen; Enabling inactive orbitals
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cecore.inc'
*
      INCLUDE 'cintfo.inc'
      INCLUDE 'cands.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
      DIMENSION VEC(*)
*
      CALL LUCIAQENTER('DIA0T')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIA0TR')
*
      ISM = ISSM
      ISPC = ISSPC
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRDIA)
      IF(NTEST.GE.1) WRITE(6,*) ' DIA0TRM : ISSM ISSPC :', ISSM,ISSPC
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' DIA0TRM speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        write(6,*) ' ISPC = ', ISPC
        WRITE(6,*) ' LUIN,LUOUT ', LUIN,LUOUT
        WRITE(6,*) ' IPART = ', IPART
        WRITE(6,*) ' FACTOR = ', FACTOR
      END IF
*
*
*. Perturbation operator
*
      IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
        I12 = 1
        IPERTOP = 1
      ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
       I12 = 2
       IPERTOP = 0
      END IF
*
      IF(NTEST.GE.10)WRITE(6,*) ' I12, IPERTOP',I12,IPERTOP

*. A bit of scracth 
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')  !done
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')  !done
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')  !done
      CALL MEMMAN(KLXA  ,NTOOB,   'ADDL  ',2,'KLXA  ')  !done
      CALL MEMMAN(KLXB  ,NTOOB,   'ADDL  ',2,'KLXB  ')  !done
      CALL MEMMAN(KLSCR ,2*NTOOB, 'ADDL  ',2,'KLSCR ')  !done
      CALL MEMMAN(KLH1D ,NTOOB,   'ADDL  ',2,'KLH1D ')  !done
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')  !done
      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR')  !done
      MAXA = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')  !done
*. Diagonal of one-body integrals and coulomb and exchange integrals
c     IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFIO),WORK(KINT1),NINT1)
      IF(IPERTOP.NE.0) then !CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        ktmp=kint1
        kint1 = kfi
        kfi = ktmp
      end if
      IF(IUSE_PH.EQ.1) THEN
c       CALL SWAPVE(WORK(KHINA),WORK(KINT1),NINT1)
        ktmp = khina
        khina = kint1
        kint1 = ktmp
        ECORE_SAVE = ECORE
        ECORE = ECORE_INA
      END IF
      CALL GT1DIA(dbl_mb(KLH1D))
      IF(IUSE_PH.EQ.1) THEN
c       CALL SWAPVE(WORK(KHINA),WORK(KINT1),NINT1)
        ktmp = khina
        khina = kint1
        kint1 = ktmp
      END IF
c     IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFIO),WORK(KINT1),NINT1)
      if(ipertop.ne.0) then
        ktmp = kfio
        kfio = kint1
        kint1 = ktmp
      end if
      WRITE(6,*) ' DIA0TRM_GAS : IPERTOP ', IPERTOP
*
      IF(I12.EQ.2) 
     &CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,dbl_mb(KLSCR2),IREOTS)
*. Interchange 
      IF(ITASK.EQ.1) THEN
        JTASK = 2
      ELSE 
        JTASK = 1
      END IF
*
*
*. Iblock driven, so just set up ....
*
      NTTS = MXNTTS                
      NOOS = NOCTPA*NOCTPB*NSMCI
      CALL MEMMAN(KLVIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'VIOIO ')  !done
      CALL MEMMAN(KLVBLTP,NSMST,'ADDL  ',2,'VBLTP ')  !done
*
      CALL IAIBCM(ISPC,int_mb(KLVIOIO))
      KSVST = 0
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,dbl_mb(KLVBLTP),int_mb(KSVST))
*
*. Space for partitioning of vectors
      NTTS = MXNTTS              
      CALL MEMMAN(KLVLBT ,NTTS  ,'ADDL  ',1,'VLBT  ')  !done
      CALL MEMMAN(KLVLEBT ,NTTS  ,'ADDL  ',1,'VLEBT ')  !done
      CALL MEMMAN(KLVI1BT,NTTS  ,'ADDL  ',1,'VI1BT ')  !done
      CALL MEMMAN(KLVIBT ,8*NTTS,'ADDL  ',1,'VIBT  ')  !done
*
      LSCR1 = 0
      CALL PART_CIV2(IDC,dbl_mb(KLVBLTP),int_mb(KNSTSO(IATP)),
     &               int_mb(KNSTSO(IBTP)),NOCTPA,NOCTPB,
     &               NSMST,LSCR1,int_mb(KLVIOIO),ISMOST(1,ISM),
     &               NXBATCH,int_mb(KLVLBT),int_mb(KLVLEBT),
     &               int_mb(KLVI1BT),int_mb(KLVIBT),1,
     &               ISIMSYM)
     &               
c     NBLOCKT = IFRMR(WORK(KLVLBT),1,1)
      NBLOCKT = int_mb(KLVLBT)
C     WRITE(6,*) ' NBLOCKT = ', NBLOCKT

      ECORES = ECORE
      WRITE(6,*) ' ECORE, FACTOR = ', ECORE, FACTOR
      CALL DIATERMS_GAS(NAEL,int_mb(KLASTR),NBEL,int_mb(KLBSTR),
     &             NTOOB,VEC,NSMST,dbl_mb(KLH1D),
     &             IDC,dbl_mb(KLXA),dbl_mb(KLXB),dbl_mb(KLSCR),
     &             dbl_mb(KLJ),dbl_mb(KLK),
     &             int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &             ECORES,LUIN,LUOUT,
     &             IPRDIA,NTOOB,
     &             dbl_mb(KLRJKA),
     &             I12,int_mb(KLVIBT),NBLOCKT,JTASK,FACTOR,0,0)
*. Clean up
      IF(IUSE_PH.EQ.1) THEN
        ECORE = ECORE_SAVE 
      END IF
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIA0TR')
      CALL LUCIAQEXIT('DIA0T')
*
      IF(NTEST.GE.100.AND.LUOUT.EQ.0) THEN
        WRITE(6,*)  ' output vector from DIA0TRM '
        CALL WRTTTS(VEC,int_mb(KLVIBT),NBLOCKT,
     &              NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),IDC)
      END IF
*
      RETURN
      END 
      SUBROUTINE DIATERM2_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
     &                       JPERT,J12,JDC)
* = DIATERM_GAS, just J12 added !
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
*
* If JPERT.NE.0, the perturbation operator as defined by IPART is used.
*
* Jeppe Olsen, August 1995
*
*. Last modification; Nov.1, 2012; Jeppe Olsen; Enabling inactive orbitals + PH
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'crun.inc'
*
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      INTEGER IBLOCK(8,*)
*
      DIMENSION VEC(*)
*
      CALL LUCIAQENTER('DIATR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIATRM')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
C     IF(JPERT.EQ.0) THEN
*. Use full Hamiltonian
C       I12 = 2
C       IPERTOP = 0
C     ELSE
*. Use perturbation operator 
C       IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
C         I12 = 1
C         IPERTOP = 1
C       ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
C         I12 = 2
C         IPERTOP = 0
C       END IF
C     END IF

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   DIATERM2_GAS speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        WRITE(6,*) ' JPERT,IPART,J12,IPERTOP',JPERT,IPART,J12,IPERTOP
      END IF
*. A bit of scracth 
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')  !done
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')  !done
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')  !done
      CALL MEMMAN(KLXA  ,NTOOB,   'ADDL  ',2,'KLXA  ')  !done
      CALL MEMMAN(KLXB  ,NTOOB,   'ADDL  ',2,'KLXB  ')  !done
      CALL MEMMAN(KLSCR ,2*NTOOB, 'ADDL  ',2,'KLSCR ')  !done
      CALL MEMMAN(KLH1D ,NTOOB,   'ADDL  ',2,'KLH1D ')  !done
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR') !done
      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR') !done
      MAXA = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')  !done
*. Diagonal of one-body integrals and coulomb and exchange integrals
      IF(IUSE_PH.EQ.1) THEN
c       CALL SWAPVE(WORK(KHINA),WORK(KINT1),NINT1)
        ktmp = khina
        khina = kint1
        kint1 = ktmp
        ECORE_SAVE = ECORE
        ECORE = ECORE_INA
      END IF
      CALL GT1DIA(dbl_mb(KLH1D))
      IF(IUSE_PH.EQ.1) THEN
c       CALL SWAPVE(WORK(KHINA),WORK(KINT1),NINT1)
        ktmp = khina
        khina = kint1
        kint1 = ktmp
      END IF
      IF(J12.EQ.2)
     &CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,dbl_mb(KLSCR2),IREOTS)
*. Core energy not included 
      ECOREP = 0.0D0
      CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,dbl_mb(KLSCR2),IREOTS)
*
      SHIFT = ECORE_ORIG-ECORE
      FACTORX = FACTOR + SHIFT
*
      CALL DIATERMS_GAS(NAEL,int_mb(KLASTR),NBEL,int_mb(KLBSTR),
     &             NTOOB,VEC,NSMST,dbl_mb(KLH1D),
     &             JDC,dbl_mb(KLXA),dbl_mb(KLXB),dbl_mb(KLSCR),
     &             dbl_mb(KLJ),dbl_mb(KLK),
     &             int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &             ECOREP,0,0,
     &             IPRDIA,NTOOB,
     &             dbl_mb(KLRJKA),
     &             J12,IBLOCK(1,IOFF),NBLOCK,ITASK,FACTORX,0,0)
C    &                  IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)    
*. Clean up
      IF(IUSE_PH.EQ.1) THEN
        ECORE = ECORE_SAVE 
      END IF
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIATRM')
      CALL LUCIAQEXIT('DIATR')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' output vector from DIATRM '
        CALL WRTTTS(VEC,IBLOCK(1,IOFF),NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),IDC)
      END IF
*
      RETURN
      END 
      SUBROUTINE DIATERM_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
     &                       JPERT,I0CHK,I0BLK)     
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
*
* If JPERT.NE.0, the perturbation operator as defined by IPART is used.
*
* IF ICHBLKS = 1, entries in IZBLKS are checked for zero blocks
*
* Jeppe Olsen, August 1995
*
* Last modification; Nov. 1, 2012; Jeppe Olsen; Allowing use of inactive orbitals + ph
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'crun.inc'
*
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      INTEGER IBLOCK(8,*)
*
      DIMENSION VEC(*)
*
      CALL LUCIAQENTER('DIATR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIATRM')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      IF(JPERT.EQ.0) THEN
*. Use full Hamiltonian
        I12 = 2
        IPERTOP = 0
      ELSE
*. Use perturbation operator 
        IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
          I12 = 1
          IPERTOP = 1
        ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
          I12 = 2
          IPERTOP = 0
        END IF
      END IF

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' DIATERM speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        WRITE(6,*) ' JPERT,IPART,I12,IPERTOP',JPERT,IPART,I12,IPERTOP
      END IF
*. A bit of scracth 
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')  !done
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')  !done
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ') !done
      CALL MEMMAN(KLXA  ,NTOOB,   'ADDL  ',2,'KLXA  ')  !done
      CALL MEMMAN(KLXB  ,NTOOB,   'ADDL  ',2,'KLXB  ')  !done
      CALL MEMMAN(KLSCR ,2*NTOOB, 'ADDL  ',2,'KLSCR ')  !done
      CALL MEMMAN(KLH1D ,NTOOB,   'ADDL  ',2,'KLH1D ')  !done
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')  !done
      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR')  !done
      MAXA = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')  !done
*. Diagonal of one-body integrals and coulomb and exchange integrals
c     IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1O),NINT1)
      IF(IPERTOP.NE.0) THEN
        ktmp = kfi
        kfi = kint1o
        kint1o = ktmp
      ENDIF
      IF(IUSE_PH.EQ.1) THEN
c       CALL SWAPVE(WORK(KHINA),WORK(KINT1),NINT1)
        ktmp = khina
        khina = kint1
        kint1 = ktmp
        ECORE_SAVE = ECORE
        ECORE = ECORE_INA
      END IF
      CALL GT1DIA(dbl_mb(KLH1D))
      IF(IUSE_PH.EQ.1) THEN
c       CALL SWAPVE(WORK(KHINA),WORK(KINT1),NINT1)
        ktmp = khina
        khina = kint1
        kint1 = ktmp
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Diagonal of 1-el ints '
        CALL WRTMAT(dbl_mb(KLH1D),1,NTOOB,1,NTOOB)
      END IF
c     IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1O),NINT1)
      if(ipertop.ne.0) then
        ktmp = kfio
        kfio = kint1
        kint1 = ktmp
      end if
      IF(I12.EQ.2)
     &CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,dbl_mb(KLSCR2),IREOTS)
      IF(NTEST.GE.10 ) THEN
        WRITE(6,*) ' ECORE and FACTOR', ECORE,FACTOR
C       WRITE(6,*) ' ECORE_ORIG, ECORE', ECORE_ORIG,ECORE
      END IF
*
      CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,dbl_mb(KLSCR2),IREOTS)
      CALL DIATERMS_GAS(NAEL,int_mb(KLASTR),NBEL,int_mb(KLBSTR),
     &             NTOOB,VEC,NSMST,dbl_mb(KLH1D),
     &             IDC,dbl_mb(KLXA),dbl_mb(KLXB),dbl_mb(KLSCR),
     &             dbl_mb(KLJ),dbl_mb(KLK),
     &             int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &             ECORE,0,0,
     &             IPRDIA,NTOOB,
     &             dbl_mb(KLRJKA),
     &             I12,IBLOCK(1,IOFF),NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)    
      IF(IUSE_PH.EQ.1) THEN
        ECORE = ECORE_SAVE
      END IF
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIATRM')
      CALL LUCIAQEXIT('DIATR')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' output vector from DIATRM '
        CALL WRTTTSC(VEC,IBLOCK(1,IOFF),NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),IDC,
     &              I0CHK,I0BLK)
      END IF
*
      RETURN
      END 
      SUBROUTINE DIATERMS_GAS(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,VEC,NSMST,H,
     &                  IDC,XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,
     &                  ECORE,LUIN,LUOUT,
     &                  IPRNT,NTOOB,RJKAA,I12,
     &                  IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)    
*
* Terms from diagonal to specific blocks
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* Calculate determinant diagonal
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, July 1995, GAS version                
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
      IMPLICIT REAL*8           (A-H,O-Z)
*.General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION H(NORB)
      DIMENSION IBLOCK(8,*)
*. 
      INTEGER I0BLK(*) 
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION VEC (*)
*
      NTEST =  000
      NTEST = MAX(NTEST,IPRNT)
C?    WRITE(6,*) ' NTEST = ',NTEST
*
      IF(LUIN.GT.0) REWIND LUIN
      IF(LUOUT.GT.0) REWIND LUOUT
 
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' DIATERMS_GAS in action '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' LUIN,LUOUT = ', LUIN,LUOUT
        WRITE(6,*) ' NBLOCK =', NBLOCK 
        WRITE(6,*) ' I0CHK = ', I0CHK   
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NORB,1,NORB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*) ' I12 and ITASK = ', I12,ITASK
        END IF
      WRITE(6,*) ' FACTOR = ',FACTOR
      END IF
*
**3 Diagonal elements according to Handys formulae
*   (corrected for error)
*
*   DIAG(IDET) = HII*(NIA+NIB)
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
*              +         J(I,J) * NIA*NJB
*
*. K goes to J - K
      IF(I12.EQ.2) 
     &CALL VECSUM(RK,RK,RJ,-1.0D0,+1.0D0,NTOOB **2)
*
      ITDET = 0
      IDET = 0
      DO JBLOCK = 1, NBLOCK
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
        IOFF = IBLOCK(6,JBLOCK)
        IF(NTEST.GE.20) THEN
         WRITE(6,*) ' Block in action : IATP IBTP IASM IBSM ',
     &               IATP,IBTP,IASM,IBSM
        END IF
*
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
*
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' After GETSTR for A strings '
          WRITE(6,*) ' alpha strings obtained '
          NAST = NSSOA(IASM,IATP)
          CALL IWRTMA(IASTR,NAEL,NAST,NAEL,NAST) 
        END IF
*
        IOFF =  1                 
        NIA = NSSOA(IASM,IATP)
        DO IA = 1 ,NSSOA(IASM,IATP)
          EAA = 0.0D0
          DO IEL = 1, NAEL
            IAEL = IASTR(IEL,IA)
            EAA = EAA + H(IAEL)
            IF(I12.EQ.2) THEN
              DO JEL = 1, NAEL
                EAA =   EAA + 0.5D0*RK(IASTR(JEL,IA),IAEL )
              END DO   
            END IF
          END DO
          RJKAA(IA-IOFF+1) = EAA 
        END DO
*. Obtain alpha strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
        NIB =  NSSOB(IBSM,IBTP)
*
        IMZERO=0
        IF(LUIN.GT.0) THEN
          CALL IFRMDS(LDET,1,-1,LUIN)
          IDET = 0
          CALL FRMDSC(VEC(1),LDET,-1,LUIN,IMZERO,IAMPACK)
        END IF
*
        IF(I0CHK.EQ.1) THEN
          IMZERO = I0BLK(JBLOCK)
          IF(IMZERO.EQ.1) THEN
*.Update offset to next block
            IF(IPACK.EQ.1.AND.IATP.EQ.IBTP) THEN
              IDET = IDET + NIA*(NIA+1)/2
            ELSE
              IDET = IDET + NIA*NIB      
            END IF
          END IF
        END IF
C?      WRITE(6,*) ' DIATERMS_GAS : I0CHK,JBLOCK IMZERO',
C?   &  I0CHK,JBLOCK,IMZERO
*
        IF(IMZERO.NE.1) THEN
*. Calculate ...
*
        DO IB = 1 ,NIB    
*
*. Terms depending only on IB
*
          HB = 0.0D0
          RJBB = 0.0D0
          CALL SETVEC(XB,0.0D0,NORB)

          DO IEL = 1, NBEL
            IBEL = IBSTR(IEL,IB)
            HB = HB + H(IBEL )
*
            IF(I12.EQ.2) THEN
              DO JEL = 1, NBEL
                RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
              END DO
*
              DO IORB = 1, NORB
                XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
              END DO
            END IF
          END DO
          EB = HB + 0.5D0*RJBB + ECORE
*
          IF(IPACK.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1                 
          END IF
*
          IASTOP = NSSOA(IASM,IATP) 
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            ITDET = ITDET + 1
            X = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NAEL
              X = X +XB(IASTR(IEL,IA)) 
            END DO
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
            IF(ITASK.EQ.1) THEN 
              IF(ABS(X+FACTOR) .GT. 1.0D-10) THEN
                VEC(IDET) = VEC(IDET)/(X+FACTOR)
              ELSE
                VEC(IDET) = 0.0D0
              END IF
            ELSE
              VEC(IDET) = VEC(IDET)*(X+FACTOR)
            END IF
           IF(NTEST.GE.1000)
     &     WRITE(6,*) ' IDET,DIAG,VEC(IDET) ', IDET,X,VEC(IDET)
          END DO
        END DO
        END IF
        
*
        IF(LUOUT.GT.0) THEN
          CALL ITODS(LDET,1,-1,LUOUT)
          CALL TODSC(VEC,LDET,-1,LUOUT)
C?        WRITE(6,*) ' Number of elements transferred to DISC ',
C?   &    LDET
          IDET = 0
        END IF
*
      END IF
      END DO
*
      IF(LUOUT.GT.0) THEN
       IONEM = -1
       CALL ITODS(IONEM,1,-1,LUOUT)
      END IF
*
C?    WRITE(6,*) ' Mission DIATERMS finished '
*
      RETURN
      END
      SUBROUTINE ADDBLKV(VEC,VECA,FACTOR,NBLKA,IBLKA,IBLOCK,IOFF)
*
* A blocked vector VEC is given 
* add factor * VECADD to this
* mapping of VECADD blocks to VEC is given by IBLKA
*
* Jeppe Olsen, August 95
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*),VECA(*)
*. Base and length of blocks in complete vector
      INTEGER IBLOCK(8,*)           
*. Blocks to be added 
      INTEGER IBLKA(NBLKA)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to ADDBLKV '
        WRITE(6,*) ' Factor ', FACTOR
        WRITE(6,*) ' blocks to be updated '
        CALL IWRTMA(IBLKA,1,NBLKA,1,NBLKA)
      END IF
*
      IOFFI = 1
      MAXE = 0
      DO JBLK = 1, NBLKA
C       IOFFO  = IBLOCK(6,IBLKA(JBLK)-1+IOFF)
C       NELMNT = IBLOCK(8,IBLKA(JBLK)-1+IOFF)
        IOFFO  = IBLOCK(6,IBLKA(JBLK))
        NELMNT = IBLOCK(8,IBLKA(JBLK))
        MAXE = MAX(MAXE,IOFFO+NELMNT-1)
        ONE = 1.0D0
        CALL VECSUM(VEC(IOFFO),VEC(IOFFO),VECA(IOFFI),
     &              ONE,FACTOR,NELMNT)
        IOFFI = IOFFI + NELMNT
       END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' output from ADDBLKV '
       CALL WRTMAT(VEC,1,MAXE,1,MAXE)
      END IF
*
      RETURN
      END 
      SUBROUTINE ADD_TTS_BLK_TO_VECTOR
     &           (NBLOCKI,IBLOCKI,IOFFI,VECI,
     &            NBLOCKO,IBLOCKO,IOFFO,VECO,FACTOR)
* A vector VECI containing NBLOCKI blocks defined by IBLOCKI
* are given. Add those blocks to a vector VECO,
* defined by NBLOCKO,IBLOCKO,IOFFO,VECO as 
*
*  VECO = VECO + FACTOR*VECI 
*
* Jeppe Olsen,  November 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLOCKI(8,*)
      INTEGER IBLOCKO(8,*)
      DIMENSION VECI(*)
*. output
      DIMENSION VECO(*)
*
      DO JBLOCKI = IOFFI,IOFFI+NBLOCKI-1
        DO JBLOCKO = IOFFO, IOFFO + NBLOCKO-1
          IF(IBLOCKI(1,JBLOCKI).EQ. IBLOCKO(1,JBLOCKO).AND.
     &    IBLOCKI(2,JBLOCKI).EQ. IBLOCKO(2,JBLOCKO)   .AND. 
     &    IBLOCKI(3,JBLOCKI).EQ. IBLOCKO(3,JBLOCKO)   .AND. 
     &    IBLOCKI(4,JBLOCKI).EQ. IBLOCKO(4,JBLOCKO)          ) THEN 
            IFROM = IBLOCKI(6,JBLOCKI)-IBLOCKI(6,IOFFI)+1
            ITO   = IBLOCKO(6,JBLOCKO)-IBLOCKO(6,IOFFO)+1
            NELMNT = IBLOCKI(8,JBLOCKI)
            ONE = 1.0D0
            CALL VECSUM(VECO(ITO),VECO(ITO),VECI(IFROM),
     &                  ONE,FACTOR,NELMNT)
          END IF
        END DO
      END DO
*
      RETURN
      END 
      SUBROUTINE GET_TTS_BLK_IN_VECTOR
     &           (NBLOCKI,IBLOCKI,IOFFI,VECI,
     &            NBLOCKO,IBLOCKO,IOFFO,
     &            NBLOCKIO,LBLOCKIO,IBLOCKIO,VECIO )
* A vector VECI containing NBLOCKI blocks defined by IBLOCKI
* are given. Obtain those blocks in VECI that corresponds to
* IBLOCKO(1,IOFFO) - IBLOCKO(1,IOFFO+NBLOCKO-1)
*
*. The number of common blocks is NBLOCKIO
*. The number of elements in these blocks is LBLOCKIO
*. Save the common blocks in VECIO and the corresponding block info in
*  IBLOCKIO
*  
*
* Jeppe Olsen,  November 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLOCKI(8,*)
      INTEGER IBLOCKO(8,*)
      DIMENSION VECI(*)
*. output
      INTEGER IBLOCKIO(8,*)
      DIMENSION VECIO(*)
*. Output
*
      IOFFP = 1
      IOFF = 1
      NBLOCKIO = 0
      LBLOCKIO = 0
      DO JBLOCKI = IOFFI,IOFFI+NBLOCKI-1
        DO JBLOCKO = IOFFO, IOFFO + NBLOCKO-1
          IF(IBLOCKI(1,JBLOCKI).EQ. IBLOCKO(1,JBLOCKO).AND.
     &    IBLOCKI(2,JBLOCKI).EQ. IBLOCKO(2,JBLOCKO)   .AND. 
     &    IBLOCKI(3,JBLOCKI).EQ. IBLOCKO(3,JBLOCKO)   .AND. 
     &    IBLOCKI(4,JBLOCKI).EQ. IBLOCKO(4,JBLOCKO)          ) THEN 
            NBLOCKIO = NBLOCKIO + 1
            CALL ICOPVE(IBLOCKI(1,JBLOCKI),IBLOCKIO(1,NBLOCKIO),8)
            IBLOCKIO(5,NBLOCKIO) = IOFF
            IBLOCKIO(6,NBLOCKIO) = IOFFP
            CALL COPVEC(VECI(IBLOCKI(6,JBLOCKI)),VECIO(IOFFP),
     &                  IBLOCKI(8,JBLOCKI))
            IOFF = IOFF + IBLOCKI(7,JBLOCKI)
            IOFFP= IOFFP+ IBLOCKI(8,JBLOCKI)
          END IF
        END DO
      END DO
      LBLOCKIO = IOFFP - 1
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
     &   ' Number of blocks obtained in GET_TTS_FROM_VECTOR ',NBLOCKIO
        WRITE(6,*)
     &   ' Number of elements obtained ',LBLOCKIO
        WRITE(6,*) ' Blocks obtained : '
        CALL IWRTMA(IBLOCKIO,8,NBLOCKIO,8,NBLOCKIO)
        WRITE(6,*) ' corresponding vector '
        NELMNTIO = IOFFP - 1
        CALL WRTMAT(VECIO,NELMNTIO,1,NELMNTIO,1)
      END IF
*
      RETURN
      END 
      SUBROUTINE REDBLK(NBLOCKI,IBLOCKI,IBBLOCKI,IOCOC,
     &                  NOCTPA,NOCTPB,
     &                  NBLOCKO,IBLOCKO,IMAP)
*
* a set of TTS blocks are given in IBLOCKI starting from 
* block IBBLOCKI,
*
* Obtain those blocks that are allowed in space defined by IOCOC 
*
* Jeppe Olsen , August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLOCKI(8,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Output
      INTEGER IBLOCKO(8,*)
      INTEGER IMAP(*)
*
      IOFF = 1
      IOFFP = 1
      NBLOCKO = 0
      DO JBLOCKI = IBBLOCKI,IBBLOCKI+NBLOCKI-1
        IATP = IBLOCKI(1,JBLOCKI)
        IBTP = IBLOCKI(2,JBLOCKI)
        IF(IOCOC(IATP,IBTP).GT.0) THEN
          NBLOCKO = NBLOCKO + 1
          IMAP(NBLOCKO) = JBLOCKI
          CALL ICOPVE(IBLOCKI(1,JBLOCKI),IBLOCKO(1,NBLOCKO),8)
          IBLOCKO(5,NBLOCKO) = IOFF
          IBLOCKO(6,NBLOCKO) = IOFFP
          IOFF = IOFF + IBLOCKO(7,NBLOCKO)
          IOFFP= IOFFP+ IBLOCKO(8,NBLOCKO)
        END IF
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of blocks obtained in REDBLK ', NBLOCKO
        WRITE(6,*) ' output => input block map '
        CALL IWRTMA(IMAP,1,NBLOCKO,1,NBLOCKO)
        WRITE(6,*) ' Blocks obtained : '
        CALL IWRTMA(IBLOCKO,8,NBLOCKO,8,NBLOCKO)
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_TTS_BATCHN(CTTS,NBLOCK,IBLOCK,NBLOCKC,IBLOCKC,
     &                  NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                  IDC,LUC,IREW,ISCALE)
*
* Read in a batch of C blocks from file LUC.
*
* The complete file is defined by NBLOCKC,IBLOCKC,
* and the blocks of the actual batch is defined by NBLOCK,IBLOCK.
* Vector packing is defined by IDC
*
* Should be initialized with rewind on LUC i.e. IREW = 1
*
*. Feb 96 : Modified to accomodate packed files
*. March 97 : Shaved, assumes IBLOCK and IBLOCKC have identical ordering
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CTTS(*),NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IBLOCKC(8,NBLOCKC),IBLOCK(8,NBLOCK)
*
      COMMON/HIDLUC/IBLK,IATP,IBTP,IASM,IBSM
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to GET_TTS_BATCH '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' Number of blocks in batch ', NBLOCK
        WRITE(6,*)
        WRITE(6,*) ' IATP  IBTP  IASM  IBSM '
        WRITE(6,*) ' ====================== '
        DO JBLOCK = 1, NBLOCK
          WRITE(6,'(4I4)') (IBLOCK(II,JBLOCK),II = 1, 4 )
        END DO
      END IF
*

*
      IF(IREW.EQ.1) THEN
        REWIND (LUC)
        IBLK = 1
      END IF
*
      IAMPACK = 1
      IMZERO= 0
*. Ready to loop over blocks
*
*. Loop over blocks to be read
      DO JBLOCK = 1, NBLOCK
          IF(IBLOCK(1,JBLOCK).GT.0) THEN
            ISKIP = 0
          ELSE 
            ISKIP = 1
          END IF
*
          IF(IBLK.GT.NBLOCKC) THEN
            REWIND LUC
            WRITE(6,*) ' Notice : LUC rewinded in GET_TTS'
            WRITE(6,*) ' Less than optimal programming ? '
            IBLK = 1
          END IF
*. Loop over blocks on file
  999   CONTINUE
*. is this block identical to next block on disc
          IF(ABS(IBLOCK(1,JBLOCK)).EQ.ABS(IBLOCKC(1,IBLK)) .AND. 
     &       IBLOCK(2,JBLOCK).EQ.IBLOCKC(2,IBLK) .AND. 
     &       IBLOCK(3,JBLOCK).EQ.IBLOCKC(3,IBLK) .AND. 
     &       IBLOCK(4,JBLOCK).EQ.IBLOCKC(4,IBLK) ) THEN
             IF(ISKIP.EQ.0) THEN
*. Read me !
              CALL IFRMDS(LBL,1,-1,LUC)
              IOFFO   = IBLOCK(6,JBLOCK)
              CALL FRMDSC(CTTS(IOFFO),LBL,-1,LUC,IMZERO,IAMPACK)
              IBLK = IBLK + 1
            ELSE
*. Skip me, you are not interested in me at all
              CALL IFRMDS(LBL,1,-1,LUC)
              CALL SKPRCD2(LBL,-1,LUC)
              IBLK = IBLK + 1
            END IF
          ELSE 
*. Skip me
            CALL IFRMDS(LBL,1,-1,LUC)
            CALL SKPRCD2(LBL,-1,LUC)
            IBLK = IBLK + 1
            GOTO 999
          END IF
C     END IF
*. End of loop over blocks on file
      END DO

*
*. Rescale from combination form to determinant
*
      IF(IDC.EQ.2.AND.ISCALE.NE.0) THEN
        CALL SCDTTS(CTTS,IBLOCK,NBLOCK,
     &  NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,
     &  IDC,2,NTEST)
C     SCDTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Batch of blocks read in '
         WRITE(6,*) ' ========================'
         WRITE(6,*)
         CALL WRTTTS(CTTS,IBLOCK,NBLOCK,NSMST,NOCTPA,NOCTPB,
     &               NSASO,NSBSO,IDC)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_TTS_BATCH(CTTS,NBLOCK,IBLOCK,NBLOCKC,IBLOCKC,
     &                  NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                  IDC,LUC,ISCR,IREW,ISCALE)
*
* Read in a batch of C blocks from file LUC.
*
* The complete file is defined by NBLOCKC,IBLOCKC,
* and the blocks of the actual batch is defined by NBLOCK,IBLOCK.
* Vector packing is defined by IDC
*
* Should be initialized with rewind on LUC i.e. IREW = 1
*
*. Feb 96 : Modified to accomodate packed files
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CTTS(*),NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IBLOCKC(8,NBLOCKC),IBLOCK(8,NBLOCK)
*
      COMMON/HIDLUC/IBLK,IATP,IBTP,IASM,IBSM
*. Local scratch : should atleast be 3*NBLOCK
      DIMENSION ISCR(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to GET_TTS_BATCH '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' Number of blocks in batch ', NBLOCK
        WRITE(6,*)
        WRITE(6,*) ' IATP  IBTP  IASM  IBSM '
        WRITE(6,*) ' ====================== '
        DO JBLOCK = 1, NBLOCK
          WRITE(6,'(4I4)') (IBLOCK(II,JBLOCK),II = 1, 4 )
        END DO
      END IF
*

*
      IF(IREW.EQ.1) THEN
        REWIND (LUC)
        IBLK = 1
        IATP = IBLOCKC(1,1)
        IBTP = IBLOCKC(2,1)
        IASM = IBLOCKC(3,1)
        IBSM = IBLOCKC(4,1)
      END IF
*
*. Map batch blocks to global blocks
*
      DO JBLOCK = 1, NBLOCK
        DO JBLOCKC = 1, NBLOCKC
          IF(IBLOCK(1,JBLOCK).EQ.IBLOCKC(1,JBLOCKC) .AND. 
     &       IBLOCK(2,JBLOCK).EQ.IBLOCKC(2,JBLOCKC) .AND. 
     &       IBLOCK(3,JBLOCK).EQ.IBLOCKC(3,JBLOCKC) .AND. 
     &       IBLOCK(4,JBLOCK).EQ.IBLOCKC(4,JBLOCKC)       ) THEN  
           ISCR(JBLOCK) = JBLOCKC
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Blocks mapped to global order '
        CALL IWRTMA(ISCR,1,NBLOCK,1,NBLOCK)
      END IF
*. Order
C     ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
*
* ORDER A STRING OF INTEGERS TO ASCENDING ORDER
*
* IINST : INPUT STRING
* IOUTST : OUTPUT STRING
* NELMNT : NUMBER OF INTEGERS 
* INO : Mapping array from new to old order             
      KLORIG = 1
      KLORDER = KLORIG + NBLOCK
      KLNTOO = KLORDER + NBLOCK
      CALL ORDINT(ISCR(KLORIG),ISCR(KLORDER),NBLOCK,ISCR(KLNTOO),NTEST)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Blocks ordered and original place '
        CALL IWRTMA(ISCR(KLORDER),1,NBLOCK,1,NBLOCK)
        CALL IWRTMA(ISCR(KLNTOO),1,NBLOCK,1,NBLOCK)
      END IF
*
*. Ready to loop over blocks
*
C     IOFFO = 1 
      DO JJBLK = 1, NBLOCK
        JBLK = ISCR(KLORDER-1+JJBLK)
        IOFFO   = IBLOCK(6,ISCR(KLNTOO-1+JJBLK))
        IBSM = IBLOCKC(4,1)
*. Loop to this block
        IPACK = 1
        IMZERO= 0
        DO KBLK = IBLK,JBLK-1
          CALL IFRMDS(LBL,1,-1,LUC)
          CALL SKPRCD2(LBL,-1,LUC)
C         CALL IFRMDS(LBL,1,-1,LUC)
C         IF(IPACK.EQ.1) CALL IFRMDS(IMZERO,1,-1,LUC)
C         IF(IPACK.EQ.0.OR.IMZERO.EQ.0) THEN
C           READ(LUC)
C         END IF 
        END DO
*. and : READIN !!!
        CALL IFRMDS(LBL,1,-1,LUC)
        CALL FRMDSC(CTTS(IOFFO),LBL,-1,LUC,IMZERO,IAMPACK)
C?      WRITE(6,*) ' GET_TTS, JJBLK,IOFFO LBL ', JJBLK,IOFFO,LBL
C       IOFFO = IOFFO + LBL 
        IBLK = JBLK+1
      END DO 
*
*. Rescale from combination form to determinant
*
C     ISCALE = 0
      IF(IDC.EQ.2.AND.ISCALE.NE.0) THEN
        CALL SCDTTS(CTTS,IBLOCK,NBLOCK,
     &  NSMST,NICTOA,NOCTPB,NSASO,NSBSO,
     &  IDC,2,NTEST)
C     SCDTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Batch of blocks read in '
         WRITE(6,*) ' ========================'
         WRITE(6,*)
         CALL WRTTTS(CTTS,IBLOCK,NBLOCK,NSMST,NOCTPA,NOCTPB,
     &               NSASO,NSBSO,IDC)
      END IF
C     WRTTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,ISC)

      RETURN
      END
      SUBROUTINE ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
*
* ORDER A STRING OF INTEGERS TO ASCENDING ORDER
*
* IINST : INPUT STRING
* IOUTST : OUTPUT STRING
* NELMNT : NUMBER OF INTEGERS 
* INO : Mapping array from new to old order             
*
* THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
* ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IINST(NELMNT),IOUTST(NELMNT),INO(NELMNT)
*
      IF(NELMNT.EQ.0) GOTO 1001
      CALL ICOPVE(IINST,IOUTST,NELMNT)
      DO  5 I = 1, NELMNT
        INO(I) = I
    5 CONTINUE
C
C       BEGIN TO ORDER
C
        JOE = 1
  10    I = JOE
  20    CONTINUE
        IF(I.EQ.NELMNT) GO TO 50
        IF(IOUTST(I).LE.IOUTST(I+1)) GO TO 40
        JOE = I + 1
  30    SWAP = IOUTST(I)
        IOUTST(I) = IOUTST(I+1)
        IOUTST(I+1) = SWAP
        ISWAP = INO(I)
        INO(I) = INO(I+1)
        INO(I+1) = ISWAP
        IF(I.EQ.1) GO TO 10
        I = I - 1
        IF(IOUTST(I).GT.IOUTST(I+1)) GO TO 30
        GO TO 10
 40     I = I + 1
      GO TO 20
C
C     END ORDER
C
 50   CONTINUE
*
 1001 CONTINUE
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GE.200) THEN
        WRITE(6,*) ' Result from ORDINT '
        WRITE(6,*) 
        WRITE(6,*)  ' Input string '
        CALL IWRTMA(IINST,1,NELMNT,1,NELMNT)
        WRITE(6,*)  ' Ordered string '
        CALL IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' New to old order '
        CALL IWRTMA(INO,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
      SUBROUTINE WRTTTS2(BLOCKS,IBLOCK,NBLOCK,IOFF,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,ISC)
*
* Print a batch of TTS blocks as given by IBLOCK.
* The blocks starts from block IOFF
*
*
* ISC = 1 : In slater determinant form 
* ISC = 2 : In Combination        form 
*
*. Jeppe Olsen, August 1995           
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input 
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,*)
*
*
      WRITE(6,*) '  Batch of blocks '
      WRITE(6,*) ' ================= '
      WRITE(6,*)
      WRITE(6,'(A,I4)') ' Number of blocks in batch ', NBLOCK
*
      DO JBLOCK = IOFF,IOFF - 1 + NBLOCK
C?      WRITE(6,*) ' TESTY, JBLOCK = ', JBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
C?      WRITE(6,*) ' IATP, IBTP, IASM, IBSM = ',
C?   &               IATP, IBTP, IASM, IBSM
*
        IF (ISC.EQ.1 ) THEN
          JOFF = IBLOCK(5,JBLOCK)
        ELSE
          JOFF = IBLOCK(6,JBLOCK)
        END IF
*
*. Is this block diagonal 
        IF(ISC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE 
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
C?      write(6,*) ' iatp ibtp iasm ibsm nia nib ',
C?   &  iatp,ibtp,iasm,ibsm,nia,nib
*
        IF(IPACK.EQ.1) THEN 
          NELMNT = NIA*(NIA+1)/2
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL PRSM2(BLOCKS(JOFF) ,NIA)
          END IF
        ELSE
          NELMNT = NIA*NIB
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL WRTMAT(BLOCKS(JOFF) ,NIA,NIB,NIA,NIB)
          END IF
        END IF
      END DO
*
      RETURN
      END 
      SUBROUTINE WRTTTSC(BLOCKS,IBLOCK,NBLOCK,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,ISC,I0CHK,I0BLK)
*
* Print a batch of TTS blocks as given by IBLOCK
*
*
* ISC = 1 : In slater determinant form 
* ISC = 2 : In Combination        form 
*
*. Jeppe Olsen, from WRTTTS, July 97 
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input 
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,NBLOCK)
      INTEGER I0BLK(NBLOCK)
*
*
      WRITE(6,*) '  Batch of blocks '
      WRITE(6,*) ' ================= '
      WRITE(6,*)
      WRITE(6,'(A,I4)') ' Number of blocks in batch ', NBLOCK
*
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF (ISC.EQ.1 ) THEN
          IOFF = IBLOCK(5,JBLOCK)
        ELSE
          IOFF = IBLOCK(6,JBLOCK)
        END IF
*
*. Is this block diagonal 
        IF(ISC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE 
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
C?      write(6,*) ' iatp ibtp iasm ibsm nia nib ',
C?   &  iatp,ibtp,iasm,ibsm,nia,nib
*
        IF(I0CHK.NE.0) THEN
          IF(I0BLK(JBLOCK).NE.0) THEN
            I_AM_ZERO = 1
          ELSE  
            I_AM_ZERO = 0
          END IF
        ELSE
          I_AM_ZERO = 0
        END IF
        IF(IPACK.EQ.1) THEN 
          NELMNT = NIA*(NIA+1)/2
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            IF(I_AM_ZERO.EQ.0) THEN
              CALL PRSM2(BLOCKS(IOFF) ,NIA)
            ELSE
              WRITE(6,*) ' Vanishing block'
            END IF
          END IF
        ELSE
          NELMNT = NIA*NIB
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            IF(I_AM_ZERO.EQ.0) THEN
              CALL WRTMAT(BLOCKS(IOFF) ,NIA,NIB,NIA,NIB)
            ELSE
              WRITE(6,*) ' Vanishing block'
            END IF
          END IF
        END IF
*
        END IF
      END DO
*
      RETURN
      END 
      SUBROUTINE WRTTTS(BLOCKS,IBLOCK,NBLOCK,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,ISC)
*
* Print a batch of TTS blocks as given by IBLOCK
*
*
* ISC = 1 : In slater determinant form 
* ISC = 2 : In Combination        form 
*
*. Jeppe Olsen, August 1995           
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input 
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,NBLOCK)
*
*
      WRITE(6,*) '  Batch of blocks '
      WRITE(6,*) ' ================= '
      WRITE(6,*)
      WRITE(6,'(A,I4)') ' Number of blocks in batch ', NBLOCK
*
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF (ISC.EQ.1 ) THEN
          IOFF = IBLOCK(5,JBLOCK)
        ELSE
          IOFF = IBLOCK(6,JBLOCK)
        END IF
*
*. Is this block diagonal 
        IF(ISC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE 
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
C?      write(6,*) ' iatp ibtp iasm ibsm nia nib ',
C?   &  iatp,ibtp,iasm,ibsm,nia,nib
*
        IF(IPACK.EQ.1) THEN 
          NELMNT = NIA*(NIA+1)/2
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL PRSM2(BLOCKS(IOFF) ,NIA)
          END IF
        ELSE
          NELMNT = NIA*NIB
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL WRTMAT(BLOCKS(IOFF) ,NIA,NIB,NIA,NIB)
          END IF
        END IF
*
        END IF
      END DO
*
      RETURN
      END 
      SUBROUTINE SCDTTS(BLOCKS,IBLOCK,NBLOCK,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
*
* Scale batch of
* blocks between determinant and combination form 
*
*
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
*
* The blocks are assumed to be in packed form !! 
*
*. Jeppe Olsen, August 1995           
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input 
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,NBLOCK)
*
C?    LOGICAL DIAGBL
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Information from SCDTTS '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
C     SUBROUTINE WRTTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,ISC,IPRNT)
        CALL WRTTTS(BLOCKS,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,2)
      END IF
*
      SQ2 = SQRT(2.0D0)
      SQ2I = 1.0D0/SQ2
*
      IBASE = 1
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IOFFP= IBLOCK(6, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*. Is this block diagonal in packed form
        IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE 
          IPACK = 0
        END IF
        NIA   = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
        IF(IPACK .EQ. 1 ) THEN                       
          NELMNT =  NIA*(NIA+1)/2
        ELSE
          NELMNT =  NIA*NIB
        END IF
*Ms combinations
        IF(IDC.EQ.2) THEN
          IF(IWAY.EQ.1) THEN
            FACTOR = SQ2
          ELSE
            FACTOR = SQ2I
          END IF
          CALL SCALVE(BLOCKS(IOFFP),FACTOR,NELMNT)
          IF(IPACK.EQ.1 ) THEN                   
            FACTOR = 1.0D0/FACTOR
            CALL SCLDIA(BLOCKS(IOFFP),FACTOR,NIA,1)     
          END IF
        END IF
*
        END IF
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output vector '  
        CALL WRTTTS(BLOCKS,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,2)
      END IF
*
      RETURN
      END
      SUBROUTINE RFTTS(BLOCKSI,BLOCKSO,IBLOCK,NBLOCK,
     &                 ICOPY,NSMST,NOCTPA,NOCTPB,
     &                 NSASO,NSBSO,IDC,PS,IWAY,IPRNT)
*
* Reformat between determinant and combination form of 
* matrices. No scaling is performed .
*
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
*
* Combination storage mode is defined BY IDC
*
*. Jeppe Olsen, August 1995           
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input 
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKSI(*),BLOCKSO(*)
      INTEGER IBLOCK(8,NBLOCK)
*

      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
      LENGTH = 0
      IF(IWAY.EQ.1) THEN
        ISCI = 1 
        ISCO = 2
      ELSE
        ISCI = 2
        ISCO = 1
      END IF
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from RFTTS  '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
        CALL WRTTTS(BLOCKSI,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,ISCI)
      END IF
*
      SQ2 = SQRT(2.0D0)
      SQ2I = 1.0D0/SQ2
*
      IBASE = 1
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF(IWAY.EQ.1) THEN
          IOFFI = IBLOCK(5,JBLOCK)
          IOFFO = IBLOCK(6,JBLOCK)
        ELSE 
          IOFFO = IBLOCK(5,JBLOCK)
          IOFFI = IBLOCK(6,JBLOCK)
        END IF
*. Is this block diagonal in packed form
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE 
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
*. Number of elements in output block
        IF(IPACK .EQ. 1 .AND. ISCO.EQ.2 ) THEN                       
          NELMNT =  NIA*(NIA+1)/2
        ELSE
          NELMNT =  NIA*NIB
        END IF
C?     WRITE(6,*) ' JBLOCK, NELMNT = ', JBLOCK,NELMNT
C?     write(6,*) 
C?   & ' RFTTS : IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
C?     WRITE(6,*)
C?   & ' RFTTS : NIA NIB IOFFI,IOFFO',NIA,NIB,IOFFI,IOFFO
*
        IF(IPACK.EQ.0) THEN
*. Just copy
          CALL COPVEC(BLOCKSI(IOFFI),BLOCKSO(IOFFO),NELMNT)
        ELSE
          IF(IWAY.EQ.1) THEN
*. unpacked => packed
C TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
            CALL TRIPK3(BLOCKSI(IOFFI),BLOCKSO(IOFFO),1,NIA,NIA,PS)
          ELSE
*. Packed => unpacked
            CALL TRIPK3(BLOCKSO(IOFFO),BLOCKSI(IOFFI),2,NIA,NIA,PS)
          END IF
        END IF
        LENGTH = LENGTH + NELMNT
        END IF
      END DO
*
      IF(ICOPY.NE.0) THEN
        CALL COPVEC(BLOCKSO,BLOCKSI,LENGTH)
      END IF
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from RFTTS  '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Output vector '
        CALL WRTTTS(BLOCKSO,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,ISCO)
      END IF
*
      RETURN
      END
      SUBROUTINE SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
     &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
     &                  ECORE,ITASK,ILEN,CNORM,SCLFAC,CONVOG,nnod) 
*
* Generate a set of sigma blocks, 
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order 
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
* ICONSPA,ICONSPB  added October 1996
* ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
*
* If ICBAT_RES .eq.1 then it as assumed that only
* Cbatches ICBAT_INI to ICBAT_END are stored on  LUC
*
* CV added to allow for incore calc, May 2010
* ECORE added, June 2010
*
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      CHARACTER*6 ITASK
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      INCLUDE 'cands.inc'
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*. Input (if ICISTR = 1) 
CNW   DIMENSION CV(*)
      integer CV,ILEN,CONVOG(nnod,2)
*
*./ORBINP/ : NACOB used
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
*
      INCLUDE 'csmprd.inc'
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      INCLUDE 'cintfo.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBLOCK')
*
      CALL LUCIAQENTER('SBLOK')
      NTEST = 000
      NTEST = MAX(NTEST,IPRCIX)
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Output from SBLOCK '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IPERTOP = ', IPERTOP
        WRITE(6,*) ' SBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
        WRITE(6,*) ' LUC in SBLOCK ', LUC
        WRITE(6,*) ' I12 in SBLOCK = ', I12
        WRITE(6,*) ' Ecore = ', ECORE
        WRITE(6,*) ' NBLOCK = ', NBLOCK
        WRITE(6,'(A,A6)') ' ITASK = ', ITASK
*
        IF(ICBAT_RES.EQ.1) THEN
          WRITE(6,*) ' C-batches restricted to interval ',
     &    ICBAT_INI,ICBAT_END 
        END IF
*
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL REWINO(LUCBLK)
      END IF
*
* Info for this internal space
*. type of alpha and beta strings - as H does not change 
*. the number of electrons, I do not distinguish between spaces for C and S
      IF(ICSPC.LE.NCMBSPC) THEN
       IATP = 1              
       IBTP = 2             
      ELSE
       IATP = IALTP_FOR_GAS(ICSPC)
       IBTP = IBETP_FOR_GAS(ICSPC)
       IF(NTEST.GE.100) WRITE(6,*) ' SBLOCK : IATP, IBTP = ', IATP, IBTP
      END IF
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. alpha and beta strings with an electron removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-1,IATPM1)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-1,IBTPM1)
*. alpha and beta strings with two electrons removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-2,IATPM2)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-2,IBTPM2)
*. Number of supergroups 
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups 
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
*. connection matrices for supergroups
*
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')  !done
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')  !done
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPA),IPRCIX,ga_nodeid())
cVOG &              int_mb(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPB),IPRCIX,ga_nodeid())
cVOG &              int_mb(KCONSPB),IPRCIX)
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')  !done
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')  !done
      CALL STSTSM(dbl_mb(KSTSTS),dbl_mb(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR           
*. Largest number of strings of given symmetry and type
      MAXA = 0
      MAXA0 = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXA = MAX(MAXA,MAXA0)
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
*
      MAXB = 0
      MAXB0 = IMNMX(int_mb(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MAXB = MAX(MAXB,MAXB0)
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.3 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
C
C BERT Make MXSTBL dependent on number of processors, with minimum of 1
C
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
C?    WRITE(6,*) ' MXTSOB = ', MXTSOB
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(ISIMSYM.EQ.0) THEN
        LSCR1 = MXSOOB
      ELSE
        LSCR1 = MXSOOB_AS
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MAX(MXTSOB ** 4, NTOOB**2)
      IF(IPRCIX.GE.5)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')  !done
      CALL MEMMAN(KINSCR2,INTSCR,'ADDL  ',2,'INSCR2')  !done
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IF( ISSPC.GE.ICSPC) THEN
         KSCIOIO = KSIOIO
      ELSE
         KSCIOIO = KCIOIO
      END IF
      IPRCIXL = 0
      CALL MXRESCPH(int_mb(KSCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIXL,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &            IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,
     &            MX_NSPII)
      IF(IPRCIX.GE.3) THEN
        WRITE(6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM
         WRITE(6,*) 'SBLOCK : MXADKBLK ', MXADKBLK
         WRITE(6,*) ' MX_NSPII = ', MX_NSPII
      END IF
      IF(ISIMSYM.EQ.1) THEN 
        MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
        MXADKBLK = MXADKBLK_AS
      END IF
*. For hardwired routines MXCIJAB is also used
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.3)  WRITE(6,*) ' LSCR2 = ', LSCR2
C  I assume memory was allocated for blocks, so
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      IF(ISIMSYM.EQ.1) LSCR3 = MAX(LSCR3,NSMST*MXADKBLK)
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')  !done
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')  !done
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')  !done
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')  !done
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')  !done
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')  !done
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')  !done
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')  !done
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,int_mb(KCBLTP),int_mb(KSVST))
*.Some TTS arrays 
      NOOS = NOCTPA*NOCTPB*NSMCI 
      NTTS = MXNTTS
*. For scaling for each TTS block
*. Moved to RASSG3
C     CALL MEMMAN(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')  !done

*. Space for four blocks of string occupations and arrays of 
*. reordering arrays
*. Also used to hold an NORB*NORB matrix  
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
*. Set up to two blocks for orbital conserving operator 
C     DO I1234 = 1, 2
      DO I1234 = 1, 1
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')  
      END DO
   
      DO I1234 = 1, 2
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
* 4 arrays containing all strings of given sym. Dimension can  be 
*   reduced to largest number of strings in alpha or beta.
C?    WRITE(6,*) ' SBLOCKS : MAX_STR_SPGP = ', MAX_STR_SPGP
c..dongxia: the following 4 arrays ar not used anywhere currently
c     IF(ISIMSYM.NE.0) THEN
c       CALL MEMMAN(KLREOJA,MAX_STR_SPGP,'ADDL  ',1,'REOJA ')
c       CALL MEMMAN(KLREOJB,MAX_STR_SPGP,'ADDL  ',1,'REOJB ')
c       CALL MEMMAN(KLREOIA,MAX_STR_SPGP,'ADDL  ',1,'REOIA ')
c       CALL MEMMAN(KLREOIB,MAX_STR_SPGP,'ADDL  ',1,'REOIB ')
c     END IF
c..
*
      I_PRINT_MEM_MAP = 0
      IF(ga_nodeid().eq.0.and.I_PRINT_MEM_MAP.EQ.1) THEN
       WRITE(6,*) ' Memory map in SBLOCK '
       IDUM = 0
       CALL MEMMAN(IDUM,IDUM,'PRINT', IDUM, 'SBLOCK')
      END IF
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*. Place perturbation integrals over one body integrals
CINSERT_START
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
*. Prepare for perturbation calculation  
*
C     IF(IPERTOP.NE.0) THEN
*. Matrix specifying partiotioned spaces
        CALL MEMMAN(KLH0SPC,NOCTPA*NOCTPB,'ADDL  ',2,'H0SPC ')  !done
        CALL H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                NGAS,MXPNGAS,dbl_mb(KLH0SPC),NELFGP)
C       IF(IH0SPC.EQ.0) THEN
*. Form of perturbation in subspace has not been defined,
*. Use current IPART
          IH0INSPC(1) = IPART
C       END IF
C     END IF
*
C?    WRITE(6,*) ' KCI1BT, WORK(KCI1BT)(3) ', KCI1BT
C?    CALL IWRTMA(WORK(KCI1BT),1,1,1,1)
*
C     IF(ISIMSYM.EQ.0) THEN
*. Active-passive division not in use, so
       KSIPA = 0
       KCJPA = 0
*. June 30 : I will try to use go through the normal route with ISIMSYM :
c      cpu_vog  = util_cpusec()
c      wall_vog = util_wallsec()
      CALL SBLOCKS(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,dbl_mb(KC2),
     &     int_mb(KCIOIO),ISMOST(1,ICSM),int_mb(KCBLTP),
     &     NACOB,int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &     NAEL,IATP,NBEL,IBTP,
     &     IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &     ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
     &     LSCR1,dbl_mb(KINSCR),dbl_mb(KCSCR),dbl_mb(KSSCR),
     &     SXSTSM,dbl_mb(KSTSTS),dbl_mb(KSTSTD),SXDXSX,
     &     ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &     int_mb(KI1),dbl_mb(KXI1S),
     &     int_mb(KI2),dbl_mb(KXI2S),IDOH2,MXPOBS,int_mb(KSVST),
     &     PSSIGN,IPRDIA,LUC,ICJKAIB,dbl_mb(KCJRES),
     &     dbl_mb(KSIRES),int_mb(KI3),dbl_mb(KXI3S),
     &     int_mb(KI4),dbl_mb(KXI4S),MXSXST,MXSXBL,
     &     MOCAA,MOCAB,IAPR,
     &     int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCI1BT),int_mb(KCIBT),
     &     IRESTRICT,int_mb(KCONSPA),int_mb(KCONSPB),CNORM,SCLFAC,
     &     LUCBLK,IPERTOP,IH0INSPC,dbl_mb(KLH0SPC),
     &     ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,I_RES_AB,
     &     IUSE_PA,int_mb(KCJPA),int_mb(KSIPA),ISIMSYM,dbl_mb(KINSCR2),
     &     MXADKBLK,ICISTR,CV,ECORE,NCBATCH,ITASK,ILEN,CONVOG,nnod)
c      cpu_vog  = util_cpusec() - cpu_vog
c      wall_vog = util_wallsec()- wall_vog
c      write(*,00002) ga_nodeid(), cpu_vog , wall_vog !VOG timings
c00002     format(/'SBLK Proc',i6,' Task  times  cpu:',f11.1,
c     &           's     wall:',f11.1,'s')
C     ELSE
C     CALL SBLOCKSN(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
C    &     WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
C    &     NACOB,WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &     NAEL,IATP,NBEL,IBTP,
C    &     IOCTPA,IOCTPB,NOCTPA,NOCTPB,
C    &     NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
C    &     ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
C    &     LSCR1,WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
C    &     SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
C    &     ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
C    &     WORK(KI1),WORK(KXI1S),
C    &     WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,WORK(KSVST),
C    &     PSSIGN,IPRDIA,LUC,ICJKAIB,WORK(KCJRES),
C    &     WORK(KSIRES),WORK(KI3),WORK(KXI3S),
C    &     WORK(KI4),WORK(KXI4S),MXSXST,MXSXBL,
C    &     MOCAA,MOCAB,IAPR,
C    &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),int_mb(KLIBT),
C    &     IRESTRICT,WORK(KCONSPA),WORK(KCONSPB),WORK(KLSCLFAC),
C    &     LUCBLK,IPERTOP,IH0INSPC,WORK(KLH0SPC),
C    &     ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,I_RES_AB,
C    &     IUSE_PA,WORK(KCJPA),WORK(KSIPA),ISIMSYM,
C    &     WORK(KLREOIA),WORK(KLREOIB),WORK(KLREOJA),WORK(KLREOJB),
C    &     IUSE_HW,WORK(KINSCR2))
C     END IF
*
      IF(NTEST.GE.200) 
     &WRITE(6,*) ' In SBLOCK after SBLOCKS, IDC =  ', IDC
    
      IF(IDC.EQ.2) THEN    
*. reform 
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL ITODS(-1,1,-1,LUCBLK)
      END IF
*. Eliminate local memory
      IDUM = 0
      IF(NTEST.GE.100) WRITE(6,*) ' Leaving SBLOCK '
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SBLOCK')
      CALL LUCIAQEXIT('SBLOK')
      RETURN
      END
      SUBROUTINE SBLOCKS(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA,NSSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB, NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &           ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &           XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC, 
     &           I1,XI1S,I2,XI2S,IDOH2,MXPOBS,ISTRFL,
     &           PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &           I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &           ICONSPA,ICONSPB,CNORM,SCLFAC,
     &           LUCBLK,IPERTOP,IH0INSPC,IH0SPC,
     &           ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &           I_RES_AB,IUSE_PA,CJPA,SIPA,ISIMSYM,XINT2,MXADKBLK,
     &           ICISTR,CV,ECORE,NCBATCH,ITASK,ILEN,CONVOG,nnod)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Contribution to Sigma-blocks from all C-blocks
*
* Jeppe Olsen , Winter of 1991
*               Last modification : April 99
*               ISIMSYM added again june 2003
*               CV added, May 2010
*               ECORE added, June 2010
*               ITASK added, March 2012
*
* Last modification; Oct. 31, 2012; Jeppe Olsen; Aligning IRESTRICT
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type 
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If IRESTRICT.NE. 0 THEN we are after :
* sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)


      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      INTEGER ISBLOCK(8,*)
* (if ICISTR = 1) 
CNW   DIMENSION CV(*)
      integer CV
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
      integer CONVOG(nnod,2)
*
      CHARACTER*6 ITASK
      LOGICAL LOCAL
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),XINT2(*),CSCR(*),SSCR(*)
*
      INTEGER   I1(MXADKBLK,*),I2(MXADKBLK,*)
      INTEGER   I3(MXADKBLK,*),I4(MXADKBLK,*) 
      DIMENSION XI1S(MXADKBLK,*),XI2S(MXADKBLK,*)
      DIMENSION XI3S(MXADKBLK,*),XI4S(MXADKBLK,*)
*. Space for storing info on batches
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*. Zero order Hamiltonian
      INTEGER IH0SPC(NOCTPA,NOCTPB)
      INTEGER IH0INSPC(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION SCLFAC(*),CNORM(*)
*
      COMMON/H_OCC_CONS/IH_OCC_CONS
* IH_OCC_CONS =1 implies that we should employ occupation conserving 
* part of Hamiltonian
      INCLUDE 'comjep.inc'
      INCLUDE 'chktyp.inc'
*.
      REAL * 8 INPROD
*
      CALL LUCIAQENTER('SBLOS')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBLOCS')
*
C?    WRITE(6,*) ' IPERTOP in SBLOCKS = ', IPERTOP
      IF(IH_OCC_CONS.EQ.1.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Occupation conserving part of Hamiltonian '
      END IF
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' SBLOCKS speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)  
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,'(A,A6)') ' ITASK = ', ITASK
      END IF
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' TTSS for each ACTIVE sigma block'
        DO IBLOCK = 1, NSBLOCK
          IF(ISBLOCK(1,IBLOCK).GT.0) 
     &    WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' Number of blocks per batch in C '
        CALL IWRTMA(LCBLOCK,1,NCBATCH,1,NCBATCH)
        NCBLOCK_T = IELSUM(LCBLOCK, NCBATCH)
        WRITE(6,'(A,I3,I8)') 
     &  ' Number of batches and active blocks in C ', NCBATCH, NCBLOCK_T
        WRITE(6,*) ' TTSS for each ACTIVE C block'
        DO IBLOCK = 1, NCBLOCK_T
          IF(ICBLOCK(1,IBLOCK).GT.0) 
     &    WRITE(6,'(10X,4I3,2I8)') (ICBLOCK(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' IDC PS IPERTOP', IDC,PS,IPERTOP
        WRITE(6,*) ' IDOH2 = ',IDOH2
        WRITE(6,*) ' I_RES_AB=',I_RES_AB
        WRITE(6,*) ' ECORE = ', ECORE
        WRITE(6,*) ' I_AM_DOING_BK = ', I_AM_DOING_BK
        WRITE(6,*) ' LUC, ICISTR = ', LUC, ICISTR
*
        IF(ICBAT_RES.EQ.1) THEN
         WRITE(6,*) ' Restricted set of C batches '
         WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        ELSE
         WRITE(6,*) ' Unrestricted C batches, NCBATCH = ', NCBATCH
        END IF
      END IF
*
      IF(NTEST.GE.500.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Initial C vector '
        IF(ICISTR.NE.1) THEN 
          CALL WRTVCD(CB,LUC,1,-1)
        ELSE
          call ga_get(CV,1,1,1,1,cfirst,1)
          WRITE(6,*) ' First element of C ',cfirst
CNW       CALL WRTMAT(CV,1,1,1,1)
*. Print using ICLBLOCK
        END IF
      END IF
*. Number of batches, blocks, and offsets
      IF(ICBAT_RES.EQ.1) THEN
        NCBATCH_L = ICBAT_END - ICBAT_INI + 1
        ICBLK_INI = I1CBLOCK(ICBAT_INI)
        NCBLOCK_T = IELSUM(LCBLOCK(ICBAT_INI), NCBATCH_L)
      ELSE 
        NCBATCH_L = NCBATCH
        ICBLK_INI = 1
        NCBLOCK_T = IELSUM(LCBLOCK(1), NCBATCH_L)
      END IF
      IF(NTEST.GE.500) WRITE(6,*) ' ICBLK_INI = ', ICBLK_INI
*
* ===========================
* 1 : Arrays for accessing C 
* ============================
*. Find the active blocks on LUC, store info in SCLFAC
C     IF(ICISTR.NE.1) THEN
*. Active blocks on LUC
C       CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC(ICBLK_INI),CB)
C     ELSE
*. Active blocks in vectors CV 
cVOG Old bug fixed
C BERT NEW NEED TO SERIALIZE FIND_ACTIVE_BLOCKS_IN_VECTOR
C ONLY NEED TO DO MY OWN BLOCK, ISSUE = FIRST ITERATION ONLY ONE BLOCK ACTIVE
C INSTEAD OF ASSUME C LOCAL ASSUME SIGMA LOCAL, REQUIRES GETS INSTEAD OF ACC
C WHICH IS MORE EFFICIENT, IN THAT CASE ROUTINE BELOW IS FINE MOSTLY
C
C PULL OUTSIDE BEFORE RASSG3 OVER BATCHES OF SIGMA, SO WE DO IT ONLY ONCE
C
C       CALL FIND_ACTIVE_BLOCKS_IN_VECTOR(CV,SCLFAC(ICBLK_INI),
C    &       NCBLOCK_T,ICBLOCK(1,ICBLK_INI))
cVOG &       NCBLOCK_T,ICBLOCK)
C     END IF
C     IF(NTEST.GE.1000) THEN
C       WRITE(6,*) ' Array giving active blocks '
C       CALL WRTMAT(CNORM(ICBLK_INI),1,NCBLOCK_T,1,NCBLOCK_T)
C     END IF
*
*. If ISIMSYM is active, enforce the same activity of all blocks 
*. belonging to a TT block
*. Number of BLOCKS in  C
        NBLOCK_C =  I1CBLOCK(NCBATCH)-1+ LCBLOCK(NCBATCH)
C       WRITE(6,*) ' Number of blocks in C', NBLOCK_C
*
       CALL ENFORCE_SIMSYM_ACT_BLK(NBLOCK_C,ICBLOCK,SCLFAC)
C           ENFORCE_SIMSYM_ACT_BLK(NBLK,IBLK,XSCA)
*
* Initialize sigma blocks
C?     CALL MEMCHK2('CHECAD')
C
      DO JSBLOCK = 1, NSBLOCK
        IATP = ISBLOCK(1,JSBLOCK)
        IBTP = ISBLOCK(2,JSBLOCK)
        IASM = ISBLOCK(3,JSBLOCK)
        IBSM = ISBLOCK(4,JSBLOCK)
        IOFF = ABS(ISBLOCK(5,JSBLOCK))
        LOCAL= ISBLOCK(5,JSBLOCK).GT.0
        NASTR = NSSOA(IASM,IATP)
        NBSTR = NSSOB(IBSM,IBTP)
        ZERO = 0.0D0
        IF((ISBLOCK(1,JSBLOCK).GT.0).AND.LOCAL)
     &  CALL SETVEC(SB(IOFF),ZERO,NASTR*NBSTR)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Initializing sigma-block:'
          WRITE(6,*) 
     &    ' JSBLOCK, IATP, IASM, IBSM =', JSBLOCK, IATP,IASM,IBSM
          WRITE(6,*) ' IOFF, NASTR, NBSTR = ', IOFF, NASTR,NBSTR
        END IF
        
      END DO
C?     CALL MEMCHK2('CHECAE')
* Loop over batches over C blocks      
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE 
        MXEXC = 1
      END IF
      IF(ICISTR.EQ.1) THEN
        IOFF_C = 1
      ELSE 
        REWIND LUC
      END IF
      IF(ICBAT_RES.EQ.1) THEN
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Restricted set of C batches '
          WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        END IF
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
*
      IOFF_C_BATCH = 0
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
        IF(NTEST.GE.100.and.ga_nodeid().eq.0)
     &  WRITE(6,*) ' 20000, JCBATCH = ', JCBATCH
        IF(LUC.LE.0) THEN
*. Address of start of batch of C
          IF(JCBATCH.EQ.JCBAT_INI) THEN
           IOFF_C_BATCH = 1
          ELSE
           IOFF_C_BATCH = IOFF_C_BATCH+LECBLOCK(JCBATCH-1)
          END IF
         END IF
*
*. Read C blocks into core
*
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
        DO JJCBLOCK = 1, NJBLOCK            
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
          IF(NTEST.GE.100) THEN
           WRITE(6,*) ' JJCBLOCK, JBLOCK  NJBLOCK = ',
     &                  JJCBLOCK, JBLOCK, NJBLOCK
           WRITE(6,*) ' I1CBLOCK(1) = ',  I1CBLOCK(1)
          END IF
       CALL MEMCHK2('CHECKA')
*. Offset for this block in complete vector
          IF(LUC.LE.0) THEN
            IOFF_C_BLOCK = IOFF_C_BATCH-1+ICBLOCK(6,JBLOCK)
          ELSE
            IOFF_C_BLOCK = 1
          END IF
*. Will this block be needed ??
          INTERACT = 0
          INTERACTL = 0
          IF(SCLFAC(JBLOCK).EQ. 1.0D0) THEN
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          PL = 1.0D0
          CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
*.Loop over Sigma blocks in batch
            DO JSBLOCK = 1, NSBLOCK
            IDENT = 0
            LOCAL = ISBLOCK(5,JSBLOCK).GT.0
            IF(ISBLOCK(1,JSBLOCK).GT.0.AND.LOCAL) THEN
              IATP = ISBLOCK(1,JSBLOCK)
              IBTP = ISBLOCK(2,JSBLOCK)
              IASM = ISBLOCK(3,JSBLOCK)
              IBSM = ISBLOCK(4,JSBLOCK)
*. Are the two blocks connected by allowed excitation
              CALL CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,
     &             LLASM,LLBSM,IOCTPA,IOCTPB,ICONSPA,ICONSPB,
     &             NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)
*
              IF(ITASK.EQ.'TRACID') THEN
*. Traci with only diagonal excitations 
                IF(.NOT.(IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &                   LLATP.EQ.IATP.AND.LLBTP.EQ.IBTP))
     &          INTERACT = 0
              END IF
*
              IF(INTERACT.EQ.1) INTERACTL = 1
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
            END IF
            END DO
          END DO
*.        ^ End of checking whether C-block is needed
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
CNW       IF(NTEST.GE.10000) WRITE(6,*) ' Interactl = ', INTERACTL
          IF(INTERACTL.EQ.1) THEN
CNW         IF(NTEST.GE.1000) THEN
CNW         WRITE(6,*) 
CNW  &      ' GSTTBL will be called for JATP, JBTP, JASM, JBSM = ',
CNW  &                                  JATP, JBTP, JASM, JBSM 
CNW         WRITE(6,*) ' IOFF_C_BLOCK = ', IOFF_C_BLOCK
CNW         END IF
CNW         CALL GSTTBL(CV(IOFF_C_BLOCK),CB(JOFF),JATP,JASM,JBTP,JBSM,
CNW  &                  ICOCOC,
CNW  &                  NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
CNW  &                  PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
C
            CALL GSTTBL(CV,CNORM(JBLOCK),IOFF_C_BLOCK,CB(JOFF),JATP,
     &                  JASM,JBTP,JBSM,ICOCOC,
     &                  NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                  PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
             IF(NTEST.GE.10000)  WRITE(6,*) ' GSTTBL was called '
*. Note in GSTTBL : ICOOSC only used for CI vectors in core, 
          ELSE
*. not relevant
            IF(LUC.GT.0) THEN
              CALL IFRMDS(LBL,-1,1,LUC)
              CALL SKPRCD2(LBL,-1,LUC)
              SCLFAC(JBLOCK) = 0.0D0
            END IF 
          END IF
*
          IF(NTEST.GE.1000) THEN
            IF(INTERACTL.EQ.1) THEN
              WRITE(6,*) ' TTSS for C block read in  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            ELSE
              WRITE(6,*) ' TTSS for C block skipped  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            END IF
          END IF
*
        END DO
*       ^ End of loop over Blocks
*
*. Loop over blocks of sigma and C in core and obtain  contribution from
*  given C block to given S block
*. Loop over C blocks
        ICOFF_TT = -0
        DO 9000 ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK 
          IF(NTEST.GE.1000) WRITE(6,*) ' 9000, ICBLK = ', ICBLK
*. Is this a new TT block -or is it identical to the previous 
*. TT block
          INEW_CTT = 0
          IF(ICBLK.EQ. I1CBLOCK(JCBATCH)) THEN
             INEW_CTT = 1
          ELSE 
             JATP_PR = JATP
             JBTP_PR = JBTP 
             IF(ICBLOCK(1,ICBLK).NE.JATP_PR.OR.
     &          ICBLOCK(2,ICBLK).NE.JBTP_PR) INEW_CTT = 1
          END IF
*. last block
          IEND_CTT = 0
          IF(ICBLK.EQ.I1CBLOCK(JCBATCH)-1+NJBLOCK) THEN
            IEND_CTT = 1
          ELSE 
            IF(ICBLOCK(1,ICBLK).NE.ICBLOCK(1,ICBLK+1).OR.
     &         ICBLOCK(2,ICBLK).NE.ICBLOCK(2,ICBLK+1)   ) IEND_CTT = 1
          END IF
C?        WRITE(6,*) ' Loop 9000, ICBLK, INEW_CTT, SCLFAC = ', 
C?   &                            ICBLK, INEW_CTT,SCLFAC(ICBLK)
*
          JATP = ICBLOCK(1,ICBLK)
          JBTP = ICBLOCK(2,ICBLK)
          JASM = ICBLOCK(3,ICBLK)
          JBSM = ICBLOCK(4,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
          IF(INEW_CTT.EQ.1) ICOFF_TT = ICOFF
          NJA = NSSOA(JASM,JATP)
          NJB = NSSOB(JBSM,JBTP)
C?        WRITE(6,*) ' ICBLK, JATP, JBTP =', ICBLK,JATP,JBTP
*
          IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
            PL = 1.0D0
            CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                  LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
            DO 8765 IPERM = NPERM,1, -1
              LLASM = LASM(IPERM)
              LLBSM = LBSM(IPERM)
              LLATP = LATP(IPERM)
              LLBTP = LBTP(IPERM)
              NLLA = NSSOA(LLASM,LLATP)
              NLLB = NSSOB(LLBSM,LLBTP)
C?            WRITE(6,*) ' LLASM, LLATP = ', LLASM, LLATP
C?            WRITE(6,*) ' LLBSM, LLBTP = ', LLBSM, LLBTP
C?            WRITE(6,*) ' SBLOCKS: NLLA, NLLB = ', NLLA, NLLB
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
              IF(IPERM.EQ.1) THEN       
*. Revert initial (IPERM = 2) negation 
                IF(IDC.EQ.2.AND.PS.EQ.-1.0D0.AND.NPERM.EQ.2) THEN
                  CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                END IF
                IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM) THEN
*. Diagonal blocks, Transposing corresponds to scaling
                  IF(PS.EQ.-1.0D0) THEN
                    CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                  END IF
                ELSE
*. ofdiagonal blocks, explicit transposing
                  CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)   
                  CALL COPVEC(C2,CB(ICOFF),NJA*NJB)   
                END IF
              ELSE
                IF(IDC.EQ.2.AND.PS.EQ.-1.0D0) THEN
                  CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
C?                CALL SCALVE(CB(ICOFF),PSSIGN,NJA*NJB)
                END IF
              END IF
*
C?     CALL MEMCHK2('CHECKB')
              DO 10000 ISBLK = 1, NSBLOCK
*. Start of a new TT block of Sigma ?
              INEW_STT = 0
              IF(ISBLK.EQ. 1) THEN
                 INEW_STT = 1
              ELSE 
                 IATP_PR = IATP
                 IBTP_PR = IBTP 
                 IF(ISBLOCK(1,ISBLK).NE.IATP_PR.OR.
     &              ISBLOCK(2,ISBLK).NE.IBTP_PR) INEW_STT = 1
              END IF

              LOCAL = ISBLOCK(5,ISBLK).GT.0
              IF(ISBLOCK(1,ISBLK) .GT. 0 .AND. LOCAL) THEN
                IATP = ISBLOCK(1,ISBLK)
                IBTP = ISBLOCK(2,ISBLK)
                IASM = ISBLOCK(3,ISBLK)
                IBSM = ISBLOCK(4,ISBLK)
                ISOFF = ABS(ISBLOCK(5,ISBLK))
                NIA = NSSOA(IASM,IATP)
                NIB = NSSOB(IBSM,IBTP)
*
                IF(ISIMSYM.EQ.0.AND.NIA*NIB.EQ.0) GOTO 10000
C ISM,IB,IA
                IF(IRESTRICT.EQ.1.AND.
     &             (JATP.GT.IATP.OR.
     &              JATP.EQ.IATP.AND.JBTP.GT.IBTP.OR.
     &              JATP.EQ.IATP.AND.JBTP.EQ.IBTP.AND.JASM.GT.IASM))
     &           GOTO 10000
COLD            IF(IRESTRICT.EQ.1.AND.
COLD &             (JASM.GT.IASM.OR.
COLD &             JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
COLD &             JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
COLD &          GOTO 10000
*. Are the two blocks connected by allowed excitation
                CALL CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,
     &               LLASM,LLBSM,IOCTPA,IOCTPB,ICONSPA,ICONSPB,
     &               NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)
*
                IF(ITASK.EQ.'TRACID') THEN
*. Traci between two FCI blocks with only a single active active space..
*. Only interactions if blocks are identical
                  IF(.NOT.(IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &                     IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP ))
     &            INTERACT = 0
                END IF
*. IF BK approximation is active, check whether block should 
* be calculated exactly, by diagonal or is set to zero
                I_DO_EXACT_BLK = 1
                IF(I_AM_DOING_BK.EQ.1) THEN
                 CALL CHECK_BLOCKS_FOR_BK_APPROX(
     &           IATP,IBTP,LLATP,LLBTP,IASM,IBSM,
     &           LLASM,LLBSM,IOCTPA,IOCTPB,I_DO_EXACT_BLK)
C?               WRITE(6,*) ' I_DO_EXACT_BLK = ', I_DO_EXACT_BLK
C     CHECK_BLOCKS_FOR_BK_APPROX(
C    &          IATP,IBTP,JATP,JBTP,IASM,IBSM,JASM,JBSM,
C    &          IOCTPA,IOCTPB,I_DO_EXACT_BLOCK)
                END IF
 
                IF(INTERACT.EQ.0.OR.I_DO_EXACT_BLK.EQ.0) GOTO 10000

                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' PROCESSOR ', ga_nodeid(), 'Before RSSB'
                  WRITE(6,*) ' Next s block in batch : '
                  write(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(6,*) ' I_DO_EXACT_BLK = ',
     &                         I_DO_EXACT_BLK
                END IF
*
                IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &            ((LLBSM.GT.LLASM).OR.
     &            (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &            GOTO 8764
*
                IF(NTEST.GE.200) THEN
                  WRITE(6,*) ' PROCESSOR ', ga_nodeid(), 'IS SPEAKING'
                  WRITE(6,*) ' RSSBCB will be called for '
                  WRITE(6,*) ' Sigma block : '
                  WRITE(6,*) ' ISOFF ', ISOFF
                  WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(6,*) ' C     block : '
                  WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                  WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                  WRITE(6,*) ' ICOFF ', ICOFF
                  WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
                END IF
*
                IF(IRESTRICT.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     XFAC = 0.5D0*SCLFAC(ICBLK)
                ELSE
                     XFAC = SCLFAC(ICBLK)
                END IF
*. Form of operator in action
C               IF(IPERTOP.NE.0) THEN
*. Not exact Hamiltonian in use 
C?                 WRITE(6,*) ' IATP, IBTP, JATP, JBTP ', 
C?   &                          IATP, IBTP, JATP, JBTP
                   IPTSPC = IH0SPC(IATP,IBTP)
                   JPTSPC = IH0SPC(JATP,JBTP)
                   IJOP   = IH0INSPC(IPTSPC)
C               ELSE
C                  IPTSPC = -1
C                  JPTSPC = -1
C                  IJOP   = -1
C               END IF
                IF(ICHKTYP.EQ.1.AND.IPTSPC.NE.JPTSPC) GOTO 8764
*
                IDUM = 0
*
*. Add core-energy to diagonal terms
*
                IF(IASM.EQ.LLASM.AND.IATP.EQ.LLATP.AND.
     &             IBSM.EQ.LLBSM.AND.IBTP.EQ.LLBTP.AND.
     &             ITASK(1:5).NE.'TRACI') THEN
                   IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
                     FACTOR = ECORE/2.0D0/SQRT(2.0D0)
                   ELSE IF(IDC.EQ.2) THEN
                     FACTOR = ECORE/SQRT(2.0D0)
                   ELSE
                     FACTOR = ECORE
                   END IF
                   CALL VECSUM(SB(ISOFF),SB(ISOFF),CB(ICOFF),
     &             1.0D0,FACTOR,NIA*NIB)
c?                 WRITE(6,*) ' Core energy added = ', ECORE
                END IF
     &             
       CALL MEMCHK2('CHECKC')
                IF(I_DO_EXACT_BLK.EQ.1) THEN
!      cpu_vogr  = util_cpusec()
!      wall_vogr = util_wallsec()
c                 if (ga_nodeid().eq.0) then
c                  CONVOG(nnod,1) = CONVOG(nnod,1) + 1
c                 else
c                  CONVOG(ga_nodeid(),1) = CONVOG(ga_nodeid(),1) + 1
c                 end if
                 CALL RSSBCB2(IASM,IATP,IOCTPA,
     &               IBSM,IBTP,IOCTPB,
     &               LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &               NELFSPGP(1,IATP+IOCTPA-1),
     &               NELFSPGP(1,IBTP+IOCTPB-1),
     &               NELFSPGP(1,LLATP+IOCTPA-1),
     &               NELFSPGP(1,LLBTP+IOCTPB-1),
     &               NAEL,NBEL,
     &               IAGRP,IBGRP,
     &               SB(ISOFF),CB(ICOFF),IDOH2,
     &               ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &               NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &               NSMOB,NSMST,NSMSX,NSMDX,
     &               NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &               ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &               MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
     &               IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
     &               IUSE_PH,IPHGAS,I_RES_AB,IUSE_PA,CJPA,SIPA,
     &               IH_OCC_CONS,ISIMSYM,XINT2,
     &               MXADKBLK,INEW_CTT,INEW_STT,IEND_CTT,CB(ICOFF_TT),
     &               ITASK,ILEN)
!      cpu_vogr  = util_cpusec() - cpu_vogr
!      wall_vogr = util_wallsec()- wall_vogr
!      write(*,00002) ga_nodeid(), cpu_vogr , wall_vogr, ISBLK, ICBLK !VOG timings
!00002     format(/'RSSB Proc',i6,' Task  times  cpu:',f11.1,
!     &           's     wall:',f11.1,'s',2i4)
                ELSE IF(I_DO_EXACT_BLK.EQ.-1) THEN
*. Multiply with diagonal
C                   ADDDIA_TERM(FACTOR,CVEC,SVEC,IASPGP,IBSPGP,IASM,IBSM,
C    &                       JPERT)
                      FACTOR = 0.0D0
                      CALL TRPMT3(SB(ISOFF),NIB,NIA,C2)
                      CALL COPVEC(C2,SB(ISOFF),NIA*NIB)
                      CALL TRPMT3(CB(ICOFF),NLLB,NLLA,C2)
                      CALL COPVEC(C2,CB(ICOFF),NLLB*NLLA)
*
                      CALL ADDDIA_TERM(FACTOR,CB(ICOFF),SB(ISOFF),
     &                      IATP,IBTP,IASM,IBSM)
* 
                      CALL TRPMT3(SB(ISOFF),NIA,NIB,C2)
                      CALL COPVEC(C2,SB(ISOFF),NIA*NIB)
                      CALL TRPMT3(CB(ICOFF),NLLA,NLLB,C2)
                      CALL COPVEC(C2,CB(ICOFF),NLLB*NLLA)
                END IF! Exact_block
*
 8764           CONTINUE
                END IF
*               ^ End if S-block should be calculated
10000         CONTINUE
*.            ^  End of loop over sigma blocks 
 8765       CONTINUE
          END IF
*         ^ End of C-block is nonvanishing
 9000   CONTINUE
*.      ^ End of loop over C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
*
* Order 
      DO  ISBLK = 1 , NSBLOCK
      LOCAL = ISBLOCK(5,ISBLK).GT.0
      IF(ISBLOCK(1,ISBLK).GT.0.AND.LOCAL) THEN
        IATP = ISBLOCK(1,ISBLK)
        IBTP = ISBLOCK(2,ISBLK)
        IASM = ISBLOCK(3,ISBLK)
        IBSM = ISBLOCK(4,ISBLK)
        ISOFF  = ISBLOCK(5,ISBLK)
        ISOFFP = ISBLOCK(6,ISBLK)
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
C BERT NEW ONLY DO IF UPDATED!!!!
        IF(ICJKAIB.NE.0) THEN
*. Tranpose sigma block was obtained, transpose to obtain correct block
           CALL TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),
     &                 NSSOA(IASM,IATP),C2)
           CALL COPVEC(C2,SB(ISOFF),
     &                 NSSOA(IASM,IATP)* NSSOB(IBSM,IBTP))
        END IF
        IF(IDC.EQ.2) THEN
         IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
           CALL TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))
         ELSE 
C?         CALL SCALVE(SB(ISOFF),PS,
C?   &          NSSOA(IASM,IATP)* NSSOB(IBSM,IBTP))
         END IF
        END IF
*
      END IF
      END DO
*
      IF(NTEST.GE.500.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' output blocks from SBLOCKS '
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,1)
      END IF
      IF(NTEST.GE.20) WRITE(6,*) ' Leaving SBLOCKS'
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SBLOCS')
      CALL LUCIAQEXIT('SBLOS')
      RETURN
      END
      SUBROUTINE IAIBCM(ICISPC,IAIB)
*
* obtain allowed combinbation of alpha- and beta- supergroups 
* for CI space ICISPC
*
* Master for IAIBCM_GAS
*
*      Jeppe Olsen, august 1995
*                   I_RE_MS2 added, May 99
*      Small changes added to treat extra space = NCISPC + 1, march 2007
*      Constraint on orbitals in  ensemble gas space added, Febr. 2012
*      (Geneva, just before darkness....)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'crun.inc'
*. Output
      INTEGER IAIB(*)
*
      IF(ICISPC.LE.NCMBSPC) THEN
        IATP = 1
        IBTP = 2
      ELSE
        IATP = IALTP_FOR_GAS(ICISPC)
        IBTP = IBETP_FOR_GAS(ICISPC)
      END IF
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        write(6,*) ' IAIBCM speaking '
        WRITE(6,*) ' =============== '
        write(6,*) ' ICISPC, LCMBSPC, ICISPC ' 
        WRITE(6,*) ICISPC,  LCMBSPC(ICISPC)
        WRITE(6,*) (ICMBSPC(II,ICISPC),II=1, LCMBSPC(ICISPC))
        WRITE(6,*) ' IATP, IBTP, NOCTPA, NOCTPB ', 
     &               IATP, IBTP, NOCTPA, NOCTPB
      END IF

      CALL IAIBCM_GAS(LCMBSPC(ICISPC),ICMBSPC(1,ICISPC),
     &            IGSOCCX,NOCTPA,
     &            NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,MXPNGAS,NGAS,IAIB,IPRDIA,I_RE_MS2_SPACE,
     &            I_RE_MS2_VALUE,I_CHECK_ENSGS)
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) 
     &  ' Matrix giving allowed combinations of al-, be-supergroups'
        CALL IWRTMA(IAIB,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END 
      SUBROUTINE OCCLS_IN_CI(NOCCLS,IOCCLS,ICISPC,NINCCLS,INCCLS)
*
* A set of occupation classes are given.
* FInd out the classes that are allowed for CI space ICISPC 
* 
*     Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc' 
*. Specific input
      INTEGER IOCCLS(NGAS,NOCCLS)
*. Scratch
C     INTEGER IACOC(MXPNGAS)
*. Output
      INTEGER INCCLS(*)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' OCCLS_IN_CI :  Input classes '
        CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
      END IF
*
      NINCCLS = 0
      DO JOCCLS = 1, NOCCLS

*
        INCLUDE = 0
        IEL = -2303
        DO JJCMBSPC = 1, LCMBSPC(ICISPC)
          JCMBSPC = ICMBSPC(JJCMBSPC,ICISPC)
          I_AM_OK = 1
          DO IGAS = 1, NGAS
            IF(IGAS.EQ.1) THEN
              IEL = IOCCLS(1,JOCCLS)
            ELSE 
              IEL = IEL + IOCCLS(IGAS,JOCCLS)
            END IF
            IF(IEL.LT.IGSOCCX(IGAS,1,JCMBSPC).OR.
     &      IEL.GT.IGSOCCX(IGAS,2,JCMBSPC)) I_AM_OK = 0
          END DO
          IF(I_AM_OK .EQ. 1 ) INCLUDE = 1
        END DO
*
        IF(INCLUDE.EQ.1) THEN
          NINCCLS = NINCCLS + 1
          INCCLS(JOCCLS) = 1      
        ELSE
          INCCLS(JOCCLS) = 0
        END IF
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from OCCLS_IN_CI '
        WRITE(6,*) ' ========================'
        WRITE(6,*) ' CI space under study ', ICISPC
        WRITE(6,*) ' Number of occupation classes included ',NINCCLS
        WRITE(6,*) ' Active classes : 1 => active 0=> inactive '
        CALL IWRTMA(INCCLS,1,NOCCLS,1,NOCCLS)  
      END IF
*
      RETURN
      END 
      SUBROUTINE DIRDIR(INSPC,IUTSPC,ISM,LUIN,LEVEL,E0,LUSCR,LUSCR2,
     &                  LU1,EOUT) 


* A vector is given in space INSPC . Calculate approximations 
* to energies in larger space UTSPC
*
* Input vector is on file LUIN, but is assumed expanded to
* IUTSPC ( Suitable empty blocks have been inserted )
*
* LEVEL = 1 : Calculate norm of (H-E0)!0>
*                          !1> = (H0-E0)**(-1) (H-E0) !0>
*             Second order energy = <0!(H-E0) (H0-E0)**(-1) (H-E0)!0>
*             <1!1>,<0!1>
* The second order energy is only meaningful for spaces not included in !0>

*
* LEVEL = 2 : Calculate also <1!H!1>, 
*      
* Jeppe Olsen, August 1995 : Trying to get my multibillion debt paid
*                            before going to Cambridge 
*
*. Last modification; Oct. 2012; Jeppe Olsen; Bringing it back alive
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      IDUM  = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DIRDI')
*
      WRITE(6,'(20X,A)') '***********************************'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '*       Welcome to DIRDIR         *'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '*   vector free calculations      *'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '* Jeppe Olsen, written   oct.1995 *'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '* Version of Oct. 2012            *'
      WRITE(6,'(20X,A)') '***********************************'
      WRITE(6,'(20X,A)') 
*
      WRITE(6,*) ' Input space ', INSPC
      WRITE(6,*) ' Output space ', IUTSPC
      WRITE(6,*) ' Level of calculation ', LEVEL
      WRITE(6,*) ' Form of perturbation operator ', MPORENP
      WRITE(6,*) 
*
      IF(NOCSF.EQ.0) THEN
        WRITE(6,*) ' Stop as you have entered DIRDIR using CSFs'
        STOP ' Stop as you have entered DIRDIR using CSFs'
      END IF
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL = NAEL + NBEL
*. Number of occupation classes in largest possible space ( compund space )
*. The occupation classes in compound space have already been calculated and stored in 
*. KIOCCLS (dimension: NOCCLS_MAX)
*. What classes of compound space are active in input and output space
      CALL MEMMAN(KLACCLI,NOCCLS_MAX,'ADDL  ',1,'ACCLI ')  !done
      CALL MEMMAN(KLACCLO,NOCCLS_MAX,'ADDL  ',1,'ACCLO ')  !done
*          OCCLS_IN_CI(NOCCLS_MAX,IOCCLS,ICISPC,NINCCLS,INCCLS)
      CALL OCCLS_IN_CI(NOCCLS_MAX,int_mb(KIOCCLS),INSPC,NACTI,
     &     int_mb(KLACCLI))
      CALL OCCLS_IN_CI(NOCCLS_MAX,int_mb(KIOCCLS),IUTSPC,NACTO,
     &     int_mb(KLACCLO))
*
*
*. Space for saving results : S11 : <1!1> in different classes
*.                            S01 : <0!1> in different classes
*.                            RES0 : norm of (H-E)!0> in different classes
*                             E2    : Second order energy contributions 
*                                     from different classes
*
      CALL  MEMMAN(KLS11,NOCCLS_MAX,'ADDL  ',2,'KLS11 ')  !done
      CALL  MEMMAN(KLS01,NOCCLS_MAX,'ADDL  ',2,'KLS01 ')  !done
      CALL  MEMMAN(KLRES0,NOCCLS_MAX,'ADDL  ',2,'KLRES0')  !done
      CALL  MEMMAN(KLE2 ,NOCCLS_MAX,'ADDL  ',2,'KLE2  ')  !done
*. Some terms defining perturbation
      WRITE(6,*) ' DIRDIR: MPORENP = ', MPORENP
      IF(MPORENP.EQ.1) THEN
* E0 is expectation value of H0
        IE0AVEX = 1
      ELSE
*. EN => Current energy is used
        IE0EVEX = 2
      END IF
*. Calculate FI+FA if required
      IF(MPORENP.EQ.1) THEN
        WRITE(6,*) ' M-P Hamiltonian will be constructed '
c       CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
        ktmp = kfi
        kfi = kint1o
        kint1o = ktmp
        CALL FIFAM(KFI) 
c       CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
        ktmp = kfi
        kfi = kfio
        kfio = ktmp
        ECORE_H = 0.0D0
        IF(IUSE_PH.EQ.1) THEN
          CALL FI(KFI,ECORE_H,0)
C         ECORE = ECORE_ORIG + ECCP 
        END IF
      END IF

      
        CALL DIRDIR1N(E0,ECORE,INSPC,IUTSPC,ISM,LUIN,
     &       dbl_mb(KLS11),dbl_mb(KLS01),dbl_mb(KLRES0),dbl_mb(KLE2),
     &       NOCCLS_MAX,int_mb(KIOCCLS),NACTI,int_mb(KLACCLI),
     &       NACTO,int_mb(KLACCLO),LUSCR,LUSCR2,LEVEL,LU1,EOUT)
*
      RETURN
      END
      SUBROUTINE DIRDIR1N(EREF,ECORE,INSPC,IIUTSPC,ISM,LUC,
     &             S11,S01,RES0,E2,
     &             NOCCLS,IOCCLS,NACTI,IACTI,
     &             NACTO,IACTO,LUSCR,LUSCR2,LEVEL,LU1,EOUT)
*
*
* First order correction to  the CI vector   
* !1> = -(H0-E0) ** -1 ( H-E0 ) ! 
*
* is calculated in pieces and the information is processed to give
*
* S11 : Overlap matrix S11(IOCCLS) = <1(IOCCLS)!1(IOCCLS)>
* S01 : Overlap matrix S01(IOCCLS) = <0(IOCCLS)!1(IOCCLS)>
* RES0 : Norm of residual =  (H-E) !0> in different classes
* E2   : <0!(H-E0)(H0-E0)**-1(H-E0)!0> in different classes
*
* NOCCLS is number of occupation classes in compound space, NACTI, IACTI, NACTO, IACTO
* are the active occupation classes in input and output
*
* Jeppe Olsen, August 1995
*
* Last modification; Oct. 2012; Jeppe Olsen; Bringing it back alive
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD, INPRDD
*
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'cprnt.inc'
*. Common block for communicating with sigma
      INCLUDE 'cands.inc'
      INCLUDE 'oper.inc'
C     COMMON/OPER/I12,IPERTOP,IAPR,MNRS1E,MXRS3E,IPART
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*. 
      DIMENSION IOCCLS(NGAS,*)
      DIMENSION S11(NOCCLS),RES0(NOCCLS),E2(NOCCLS)
      DIMENSION S01(NOCCLS)
      CHARACTER*6 CTASK
*
      IUTSPC = IIUTSPC
      CTASK = 'SIGMA '

*
      NTEST =   10
*. LBLK should be defined outside for future use
      LBLK = -1
      
      
*
C     IF(NTEST.GE.3) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ================================= '
        WRITE(6,*) '        Welcome to DIRDIR1N '
        WRITE(6,*) ' ================================= '
        WRITE(6,*) 
        WRITE(6,*) '     Input  CI space ', INSPC
        WRITE(6,*) '     Output CI space ', IUTSPC
        WRITE(6,*) '     Level parameter ', LEVEL
        WRITE(6,*)
        WRITE(6,*) ' Initial CI vector is in space ', INSPC
        WRITE(6,*)
        WRITE(6,*) ' Calculations to be performed : '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        IF(LEVEL.GE.1) 
     &  WRITE(6,*) '     PT2 calculation in space .......', IUTSPC
        IF(LEVEL.GE.2) 
     &  WRITE(6,*) '   + 1 CI iteration + PT3 in space ..', IUTSPC
        IF(LEVEL.EQ.3)
     &  WRITE(6,*) '   + PT4 in space ...................', IUTSPC
        IF(LEVEL.EQ.4)
     &  WRITE(6,*) '   + PT4 in space ...................', IUTSPC+1
        WRITE(6,*)
        IF(MPORENP.EQ.1) THEN
          WRITE(6,*) 
     &    '     Perturbation operator: diagonal Moller-Plesset '
        ELSE IF(MPORENP.EQ.2) THEN
          WRITE(6,*) 
     &    '     Perturbation operator: diagonal Epstein-Nesbet'
        ELSE
          WRITE(6,*) ' Unknown perturbation parameter = ', MPORENP
          STOP ' Unknown perturbation parameter '
        END IF
        WRITE(6,*)
    
       
*
        IF(IC1DSC.EQ.0) THEN 
          WRITE(6,*) ' First order correction is not stored '
        ELSE
          WRITE(6,*) ' First order correction stored on file ', LU1
        END IF
        WRITE(6,*) ' Scratch files ', LUSCR,LUSCR2
        WRITE(6,*) ' LU1 = ', LU1
*
C     END IF
*. Zero - before I forget
      ZERO = 0.0D0
      CALL SETVEC(S11,ZERO,NOCCLS)
      CALL SETVEC(S01,ZERO,NOCCLS)
      CALL SETVEC(RES0,ZERO,NOCCLS)
      CALL SETVEC(E2,ZERO,NOCCLS)
*.Tired of compiler warnings
      KLB = -2803
*. Transfer to /CANDS
      ICSPC = INSPC
      ISSPC = IUTSPC
*. zero and first order space
      ISPC0 = INSPC
      ISPC1 = IUTSPC
*
      ISSM = ISM
      ICSM = ISM
*. I will not work with complete vectors in a single batch
      ICOMP = 0
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIRDR1')
* The calculation of (H-E0) will be processed in blocks,
* allocate blocks for (H-E0)!0> and !0>
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSOOB_AS
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      WRITE(6,*) ' Size of BATCH ', LSCR1
      CALL MEMMAN(KCB,LSCR1,'ADDL  ',2,'KCB   ')
      CALL MEMMAN(KSB,LSCR1,'ADDL  ',2,'KSB   ')
CM    CALL MEMMAN(KLB,LSCR1,'ADDL  ',2,'KLB   ')
*
*. Structure and dimensions of output CI vector.
*
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
C START
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR           
*. type of alpha and beta strings
      IATP = 1              
      IBTP = 2             
*. alpha and beta strings with an electron removed
      IATPM1 = 3 
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5 
      IBTPM2 = 6
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. Block for storing complete or partial CI-vector
      IF(ISIMSYM.EQ.1.OR.ICISTR.EQ.2) THEN
        LBLOCK = MXSOOB_AS
      ELSE
        LBLOCK = MXSOOB
      END IF
      IF(NOCSF.EQ.0) THEN
        LBLOCK  = MAX(NSD_FOR_OCCLS_MAX,MXSOOB)
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' TEST, DIRDIR LCSBLK, LBLOCK, MXSOOB  = ',
     &                          LCSBLK, LBLOCK, MXSOOB
      END IF
      ICOMP = 0
      ILTEST = -3006

      ISFIRST = 1
      IF(ISFIRST.EQ.1) THEN
*. Actually, as the code is running  only for SD's the occupation class info below is 
*. redundant..
        CALL Z_BLKFO_FOR_CISPACE(ISSPC,ISSM,LBLOCK,ICOMP,
     &       NTEST,NSBLOCK,NSBATCH,
     &       int_mb(KSIOIO),int_mb(KSBLTP),NSOCCLS_ACT,
     &       int_mb(KSIOCCLS_ACT),
     &       int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSLBLK),int_mb(KSI1BT),
     &       int_mb(KSIBT),IDUM,IDUM,ILTEST)
      END IF
* KSLBT: Number of blocks per batch
* KSLEBT: Number of elements per batch
* KSLBLK: Length of each block (packed)
* KSI1BT: Is first blovck of each batch
* KSIBT is IBLKFO (the eight column blockfo)
*
      ICFIRST = 1
      IF(ICFIRST.EQ.1) THEN
        CALL Z_BLKFO_FOR_CISPACE(ICSPC,ICSM,LBLOCK,ICOMP,
     &       NTEST,NCBLOCK,NCBATCH,int_mb(KCIOIO),int_mb(KCBLTP),
     &       NCOCCLS_ACT,int_mb(KCIOCCLS_ACT),
     &       int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &       int_mb(KCIBT),
     &       int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
      END IF
      NMBLOCK = MAX(NCBLOCK,NSBLOCK)
      CALL MEMMAN(KLISCR,NMBLOCK,'ADDL  ',2,'ISCR  ')

*. a bit of local scratch
      NTTS = MXNTTS                
      CALL MEMMAN(KLXIBT ,8*NTTS,'ADDL  ',1,'XIBT  ')
      CALL MEMMAN(KLXISCR,3*NTTS,'ADDL  ',1,'KLXISC')
      CALL MEMMAN(KLXISC2,3*NTTS,'ADDL  ',1,'KLXIS2')
C KLXIBT, KLXISCR, KLXISC2
*
*. Blocks and scratch for sigma
      IF(LEVEL.LE.2) THEN
        IPRCIXL = 0
        CALL MXRESCPH(int_mb(KSIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &              NSMST,NSTFSMSPGP,MXPNSMST,
     &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIXL,MAXK,
     &              NELFSPGP,
     &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &              IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,
     &              MXADKBLK_AS,MX_NSPII)
        IF(IPRCIX.GE.2) THEN
          WRITE(6,*) 'DIRDIR : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                         MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
           WRITE(6,*) 'DIRDIR : MXADKBLK,MX_NSPII ', MXADKBLK,MX_NSPII
           WRITE(6,*) ' DIRDIR: MXCJ_ALLSYM = ', MXCJ_ALLSYM
        END IF
        IF(ISIMSYM.EQ.1) THEN
          MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
          MXADKBLK_AS = MXADKBLK
        END IF
*
        LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
        IF(IPRCIX.GE.3) WRITE(6,'(A,3(2X, I9))')
     &  'MXCJ,MXCIJA,MXCIJB,MXCIJAB = ',
     &   MXCJ,MXCIJA,MXCIJB,MXCIJAB
        IF(IPRCIX.GE.2)
     &  WRITE(6,*) ' Space for resolution matrices ',LSCR2
        LSCR12 = MAX(LSCR1,2*LSCR2)  
        WRITE(6,*) ' LSCR12 = ', LSCR12
        CALL MEMMAN(KVEC3,LSCR12,'ADDL  ',2,'KLB   ')
        KLB = KVEC3
      END IF
*
* Information on second order correction
*
     
      ISPC2 =  0
CT    IF(LEVEL.GE.3) THEN
*
CT      IF(LEVEL.EQ.3) THEN
CT       ISPC2 = IUTSPC    
CT      ELSE
CT       ISPC2 = IUTSPC + 1
CT      END IF 
*
CT      CALL MEMMAN(KLIOIO2,NOCTPA*NOCTPB,'ADDL  ',2,'IOIO2 ')
CT      CALL MEMMAN(KLBLTP2,NSMST,'ADDL  ',2,'BLTP2 ')
*. Partitioning of second order correction 
CT      NTTS = MXNTTS
CT      CALL MEMMAN(KLLBT2 ,NTTS  ,'ADDL  ',1,'LBT2  ')
CT      CALL MEMMAN(KLLEBT2 ,NTTS  ,'ADDL  ',1,'LEBT2 ')
CT      CALL MEMMAN(KLI1BT2,NTTS  ,'ADDL  ',1,'I1BT2 ')
CT      CALL MEMMAN(KLIBT2 ,8*NTTS,'ADDL  ',1,'IBT2  ')
*
CT      CALL Z_BLKFO_FOR_CISPACE(ISPC2,ISSM,LBLOCK,ICOMP,
CT   &       NTEST,NC2BLOCK,NC2BATCH,
CT   &       WORK(KLIOIO2),WORK(KLBLTP2),NC2OCCLS_ACT,
CT   &       WORK(KC2IOCCLS_ACT),
CT   &       WORK(KLC2LBT),WORK(KLLEBT2),int_mb(KSLBLK),int_mb(KSI1BT),
CT   &       int_mb(KSIBT),
CT   &       WORK(KSNOCCLS_BAT),WORK(KSIBOCCLS_BAT),0,ILTEST)
*
CT      CALL IAIBCM(ISPC2,WORK(KLIOIO2))
CT      KSVST = 1
CT      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLTP2),WORK(KSVST))
CT      WRITE(6,*) ' Number of batches for |2>   ', NBATCH2
CT      CALL MXRESC(WORK(KLIOIO2),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
CT   &              NSMST,NSTFSMSPGP,MXPNSMST,
CT   &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
CT   &              NELFSPGP,
CT   &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
CT      IF(IPRCIX.GE.2) THEN
CT        WRITE(6,*) 'DIRDIR : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
CT   &                         MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
CT         WRITE(6,*) 'DIRDIR : MXADKBLK ', MXADKBLK
CT      END IF
CT      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
CT      IF(IPRCIX.GE.2)
CT   &  WRITE(6,*) ' Space for resolution matrices ',LSCR2
CT      LSCR12 = MAX(LSCR1,2*LSCR2)  
CT      CALL MEMMAN(KVEC3,LSCR12,'ADDL  ',2,'KLB   ')
CT      KLB = KVEC3
CT     END IF! Level = 3 
*. Total number of blocks in C
      NTBLOCKC = NCBLOCK
*. Total number of blocks in In Sigma ( First order correction )
      NTBLOCKS = NSBLOCK
      IF(NTEST.GE.1) THEN
        WRITE(6,'(A,2I5)')
     &  ' Number of batches for !0> and !1> ', NCBATCH, NSBATCH
        WRITE(6,'(A,2I8)')
     &  ' Number of blocks for !0> and !1> ', NTBLOCKC, NTBLOCKS
       END IF
       
*
* ======================================================
* Calculate first order correction to wave function and 
* second order correction to energy 
* ======================================================
*
       WRITE(6,*) 
       WRITE(6,'(20X,A)') ' **********************************'
       WRITE(6,'(20X,A)') '  Second order energy calculation : '
       WRITE(6,'(20X,A)') ' **********************************'
       WRITE(6,*) 
       WRITE(6,*) 
*
       IF(IC1DSC.EQ.1) THEN
         CALL REWINO(LU1)
       END IF
*
*. Start by studying whether |0> is an eigenfunction for H(apr).
*  This has some bearings upon how to evaluate (H0-E0)**-1
*. Calculate H(apr)|0> - <0!H(apr)!0> |0>
      CALL REWINO(LUC)
      CALL REWINO(LUSCR)
      ICSPC = ISPC0
      ISSPC = ISPC0
*. H(apr) |0>
      ZERO = 0.0D0
      IPART = MPORENP
      CALL DIA0TRM_GAS(1,LUC,LUSCR,dbl_mb(KCB),dbl_mb(KSB),ZERO)
*. <0|H(apr)|0>
      ONE = 1.0D0
      HAPR00 = INPRDD(dbl_mb(KCB),dbl_mb(KSB),LUC,LUSCR,1,LBLK)
      CALL VECSMD(dbl_mb(KCB),dbl_mb(KSB),ONE,-HAPR00,
     &            LUSCR,LUC,LUSCR2,1,LBLK)
      XNORM = INPRDD(dbl_mb(KCB),dbl_mb(KSB),LUSCR2,LUSCR2,1,LBLK)
*
      IF(ABS(XNORM/HAPR00) .LE. 1.0D-12) THEN
       IHAPREIG = 1
      ELSE
       IHAPREIG = 0
      END IF
*. Choice of zero order energy 
      IF(IE0AVEX.EQ.1) THEN
        E0 = HAPR00
      ELSE
        E0 = EREF
      END IF
      WRITE(6,*) ' E0 = ', E0
*
      HAPRM100 = 0.0D0
      IF(IHAPREIG.EQ.0) THEN
*<0!(H(apr)-E0)**-1 |0>
        CALL DIA0TRM_GAS(2,LUC,LUSCR,dbl_mb(KCB),dbl_mb(KSB),-E0)
        HAPRM100  = INPRDD(dbl_mb(KCB),dbl_mb(KSB),LUC,LUSCR,1,LBLK)
        WRITE(6,*) ' HAPRM100', HAPRM100
      ENDIF
*. First order energy 
      E1T = EREF - E0
*
      WRITE(6,*) ' Zero order energy  ', E0
      WRITE(6,*) ' First order energy ', E1T


      WRITE(6,*) '  HAPR00,  XNORM, IHAPREIG, HAPRM100 : ',
     &              HAPR00,  XNORM, IHAPREIG, HAPRM100
      CALL GFLUSH(6)
C     stop ' Enforced stop in DIRDIR1N'

*
*. Loop over partitionings of (H-E)!0> vector
      IFIRST = 1
      S01T = 0.0D0
      ICSPC = ISPC0
      ISSPC = ISPC1
      DO IBATCH = 1, NSBATCH
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Information about sigma BATCH ', IBATCH
          CALL GFLUSH(6)
        END IF
*. Information about this batch 
        NBLOCK = IFRMR(int_mb(KSLBT),1,IBATCH)
        NELMNT = IFRMR(int_mb(KSLEBT),1,IBATCH)
        IBLOCK_OFF = IFRMR(int_mb(KSI1BT),1,IBATCH)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' NBLOCK, NELMNT = ', NBLOCK, NELMNT
        END IF
*. Obtain blocks of first order corrections 
        SIN = 0.0D0
        CALL GET_BATCH_OF_FIRST
     &       (NBLOCK,IBLOCK_OFF,NELMNT,int_mb(KSIBT),ISM,IUTSPC,
     &        E0,E1T,SIN,LUC,IDC,
     &        int_mb(KCIOIO),NOCTPA,NOCTPB,NSMST,
     &        int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &        NTBLOCKC,int_mb(KCIBT),int_mb(KLXIBT),int_mb(KLXISCR),
     &        int_mb(KLXISC2),dbl_mb(KSB),dbl_mb(KCB),dbl_mb(KLB),SOUT,
     &        HAPR01,0,NTBLOCKS)
C KLXIBT, KLXISCR, KLXISC2
C     GET_BATCH_OF_FIRST
C    &           (N1BLOCK,I1OFF,N1ELMNT,I1BLOCK,I1SM,I1SPC,
C    &            E0,E1,ALPHA,LU0,IDC,
C    &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
C    &            N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,
C    &            HAPR01,LU1,N1BLOCKT )
*. On output : LB : Batch of first order correction
*              SB : Batch of (H-E0)|0>
*
        IF(IC1DSC.EQ.1) THEN
*. Save batch of first order corrections on DISC '
C              EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
          CALL EXTRROW2(int_mb(KSIBT),8,IBLOCK_OFF,8,
     &                  NBLOCK,dbl_mb(KLISCR))
          CALL TODSCN(dbl_mb(KLB),NBLOCK,dbl_mb(KLISCR),LBLK,LU1)
        END IF
*
        S01T = S01T + SOUT
*. Find contributions to 
*. 1 : norm of nonorthogonalized first order correction vector
        CALL CLASS_PROD(dbl_mb(KLB),dbl_mb(KLB),NOCTPA,NOCTPB,
     &                  IBLOCK_OFF,NBLOCK,
     &                  int_mb(KSIBT),NOCCLS,IOCCLS,S11)
*. 2 : norm of residual vector
        CALL CLASS_PROD(dbl_mb(KSB),dbl_mb(KSB),NOCTPA,NOCTPB,
     &                  IBLOCK_OFF,NBLOCK,
     &                  int_mb(KSIBT),NOCCLS,IOCCLS,RES0)
*. 3 :  contributions to second order energy correction
        CALL CLASS_PROD(dbl_mb(KSB),dbl_mb(KLB),NOCTPA,NOCTPB,
     &                  IBLOCK_OFF,NBLOCK,
     &                  int_mb(KSIBT),NOCCLS,IOCCLS,E2)
      END DO
*. To put an end of to first order correction file  
      IF(IC1DSC.EQ.1) THEN
        CALL ITODS(-1,1,LBLK,LU1)
      END IF
*. We wanted the norm of the residual so
      XNORMR = 0.0D0 
      XNORM1 = 0.0D0 
      DO JOCCLS = 1, NOCCLS
       XNORMR = XNORMR + RES0(JOCCLS) 
       XNORM1 = XNORM1 + S11(JOCCLS) 
       RES0(JOCCLS) = SQRT(RES0(JOCCLS))
      END DO 
      XNORMR = SQRT(XNORMR)
      XNORM1 = SQRT(XNORM1)
*
      IF(IPART.EQ.1) THEN
        IHYLLE2 = 1
      ELSE
        IHYLLE2 = 0
      END IF
      WRITE(6,*) ' TESTY, IHYLLE2, IPART = ', IHYLLE2, IPART
*. (Not invoked for EN hamiltonian )
* 
*. Calculate <1!H0-E0!1> to Obtain Hylleraas expression for second order 
*  energy : <1!H0-E0!1> + 2 <0!H!1>
* 
* |1> = |1'> - alpha1*(H apr -E )**-1 |0>
*
* <1!H0-E0!1> = <1!H apr -E0!1> = <1'!H apr -E0!1'> 
*                               - <0|1'>**2 / <0!(H0-E0)-1!0>
      IF(IHYLLE2.EQ.1.AND.IPART.EQ.1) THEN
*
*
* Calculate <1'|H apr |'1> = sum(l,r) <0(l)|H apr |0(r)> 
*                   = sum(l)    <0(l)|H apr |0(l)>
*                   + 2sum(r.lt.l) <0(l)|H apr |0(r)>
*
*. Loop over batches of <1'| and |1'> 
        WRITE(6,*)
        WRITE(6,*) ' ==============================='
        WRITE(6,*) '  Hylleraas second order energy  '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL REWINO(LUSCR)
        H011 = 0.0D0
*. Only one-electron part in sigma 
        I12 = 1
*. Loop over left batches, i.e.  batches of H(apr)!1>
        DO ILBATCH = 1, NSBATCH
          WRITE(6,*) ' Starting outer loop for left  batch ', ILBATCH
          NLBLOCK = IFRMR(int_mb(KSLBT),1,ILBATCH)
          NLELMNT = IFRMR(int_mb(KSLEBT),1,ILBATCH)
          ILOFF = IFRMR(int_mb(KSI1BT),1,ILBATCH)
*. Loop over right batches, i.e.  batches of !1>
          DO IRBATCH  = 1, ILBATCH 
            WRITE(6,*) 
     &      '     Information about batches (left and right)',
     &      ILBATCH,IRBATCH
            NRBLOCK = IFRMR(int_mb(KSLBT),1,IRBATCH)
            NRELMNT = IFRMR(int_mb(KSLEBT),1,IRBATCH)
            IROFF = IFRMR(int_mb(KSI1BT),1,IRBATCH)
*
* Obtain first order correction |1'> in IRBATCH 
*
            IF(IC1DSC.EQ.0) THEN
              SIN = 0.0D0
              CALL GET_BATCH_OF_FIRST
     &             (NRBLOCK,IROFF,NRELMNT,int_mb(KSIBT),ISM,IUTSPC,
     &              E0,E1T,SIN   ,LUC,IDC,
     &              int_mb(KCIOIO),NOCTPA,NOCTPB,NSMST,
     &              int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &              NTBLOCKC,int_mb(KCIBT),int_mb(KLXIBT),
     &              int_mb(KLXISCR),int_mb(KLXISC2),
     &              dbl_mb(KSB),dbl_mb(KCB),dbl_mb(KLB),SOUT,
     &              HAPR01,0,NTBLOCKS)
C KLXIBT, KLXISCR, KLXISC2
C     GET_BATCH_OF_FIRST
C    &           (N1BLOCK,I1OFF,N1ELMNT,I1BLOCK,I1SM,I1SPC,
C    &            E0,E1,ALPHA,LU0,IDC,
C    &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
C    &            N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,
C    &            HAPR01,LU1,N1BLOCKT )
            ELSE
*. Read in from DISC
              CALL GET_BLOCKS_FROM_DISC
     &        (LU1,NRBLOCK,IROFF,int_mb(KSIBT),NTBLOCKS,dbl_mb(KLB),1)
            END IF 
*
*. Transfer batch of first order correction to  LUSCR2 
*  ( Sblock/Sblocks assumes that blocks are  on disc)
C              EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
          CALL EXTRROW2(int_mb(KSIBT),8,IROFF,8,
     &                  NRBLOCK,int_mb(KLXISCR))
          CALL REWINO(LUSCR2)
          CALL TODSCNP(dbl_mb(KLB),NRBLOCK,int_mb(KLXISCR),-1,LUSCR2)
          CALL ITODS(-1,1,-1,LUSCR2)
C         CALL ITODS(LREC(IREC),1,LBLK,LU)
*
*. Obtain (left) batch of H apr  |1(right_batch>  
*
            I12 = 1
*. MP one-electron integrals
C           CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
            IRESTR = 1

            ICSPC = IUTSPC
            ISSPC = IUTSPC
            CALL SBLOCK(NLBLOCK,int_mb(KSIBT),ILOFF,dbl_mb(KCB),
     &                  dbl_mb(KSB),
     &                  LUSCR2,IRESTR,0,1,IRBATCH,IRBATCH,0.0D0,ECORE,
     &                  CTASK)
C     SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
C    &                  ECORE,ITASK) 
COLD             SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END) 
*. Add to previous obtained contributions to left batch of H apr !1>
            IF(IRBATCH.NE.1) THEN
              CALL REWINO(LUSCR)
              CALL FRMDSC(dbl_mb(KLB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
              ONE = 1.0D0
              CALL VECSUM(dbl_mb(KSB),dbl_mb(KSB),dbl_mb(KLB),ONE,ONE,
     &                    NLELMNT)
            END IF
*. And transfer to disc to save for future generations
            CALL REWINO(LUSCR)
            CALL TODSC(dbl_mb(KSB),NLELMNT,-1,LUSCR)
*. Restore usual one-electron integrals
C           CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
          END DO
*         ^ End of loop over batches of !1'> (right batches)
*. Left Batch of H apr |1'> is now completed 
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Batch of H0 |1> '
            CALL WRTTTS2(dbl_mb(KSB),int_mb(KSIBT),NLBLOCK,ILOFF,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),IDC)
          END IF
*
*. Obtain left batch of first order corrections
*
          ICSPC = INSPC
          ISSPC = IUTSPC
          IF(IC1DSC.EQ.0) THEN
            SIN = 0.0D0
            CALL GET_BATCH_OF_FIRST
     &           (NLBLOCK,ILOFF,NLELMNT,int_mb(KSIBT),ISM,IUTSPC,
     &            E0,E1T,SIN,LUC,IDC,
     &            int_mb(KCIOIO),NOCTPA,NOCTPB,NSMST,
     &            int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &            NTBLOCKC,int_mb(KCIBT),int_mb(KLXIBT),int_mb(KLXISCR),
     &            int_mb(KLXISC2),dbl_mb(KSB),dbl_mb(KCB),dbl_mb(KLB),
     &            SOUT,
     &            HAPR01,0,NTBLOCKS)
C KLXIBT, KLXISCR, KLXISC2
C     GET_BATCH_OF_FIRST
C    &           (N1BLOCK,I1OFF,N1ELMNT,I1BLOCK,I1SM,I1SPC,
C    &            E0,E1,ALPHA,LU0,IDC,
C    &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
C    &            N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,
C    &            HAPR01,LU1,N1BLOCKT )

          ELSE
*. Read in from DISC
              CALL GET_BLOCKS_FROM_DISC
     &        (LU1,NLBLOCK,ILOFF,int_mb(KSIBT),NTBLOCKS,dbl_mb(KLB),1)
          END IF 
*. Obtain saved sigma block 
          CALL REWINO(LUSCR)
          CALL FRMDSC(dbl_mb(KSB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
          H011 = H011 + INPROD(dbl_mb(KLB),dbl_mb(KSB),NLELMNT)
        END DO
*       ^ End of loop over left batches
*
        IF(IRESTR.EQ.1) THEN
*. Only half of H011 was constructed, 
          H011 = 2.0D0*H011
        END IF
*
        WRITE(6,*) ' H011 ', H011
        H0ME011P = H011 - E0*S11P
        WRITE(6,*) ' <1''!H0-E0!1''> ', H0ME011P
        IF(IHAPREIG.EQ.0) THEN
          H0ME011 = H0ME011P - S01T**2/HAPRM100
        ELSE
          H0ME011 = H0ME011P
        END IF
        WRITE(6,*)  ' <1!H0-E0!1> ',  H0ME011
      END IF
*     ^ End if for Hylleraas calculation
*
*
* Collect second order energy and print
*
      S11P = XNORM1 ** 2
* 
      IF(IHAPREIG.EQ.0) THEN
        ALPHA1 = S01T/HAPRM100
      ELSE
        ALPHA1 = 0.0D0
      END IF
*
*. Orthogonalizetion term to second order energy 
      E2ORT = ALPHA1*S01T
*. Total energy through second order
      E2T = 0.0D0
      DO JOCCLS = 1, NOCCLS
       E2T = E2T + E2(JOCCLS)
      END DO
      E2T = E2T + E2ORT
*. Hylleraas form of E2
      IF(IHYLLE2.EQ.1) E2TH = 2*E2T +  H0ME011
*. Total energy through second order 
      E012T =  E0+E1T+E2T
      EOUT = E012T
*. Total energy through second order , Hylleraas form
      IF(IHYLLE2.EQ.1) E012TH =  E0+E1T+E2TH
*, Normalization factor 
      FACN = 1.0D0/XNORM1
      WRITE(6,*) ' Occupation classes '
      CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
*
C?    IF(NTEST.GE.1 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Contributions to second order energy :  '
        WRITE(6,*) ' ======================================= '
        WRITE(6,*)
        WRITE(6,'(6X,A)') 
     &  'Class  Norm of !1(u)> Norm of (H-E)!0>   Second order energy'
        WRITE(6,'(6X,A)') 
     &  ' ============================================================'
        DO JOCCLS = 1, NOCCLS
         WRITE(6,'(6x,I3,2X,E13.6,2X,E13.6,7X,E16.9)')
     &   JOCCLS,S11(JOCCLS),RES0(JOCCLS),E2(JOCCLS)
        END DO

        WRITE(6,*) 
        WRITE(6,'(A,F18.9)') 
     &  ' Second order energy contribution ................... ',E2T
        WRITE(6,*) 
        WRITE(6,'(A,F18.9)') 
     &  ' Second order energy Orthogonalization term ......... ',E2ORT
        WRITE(6,*)
        WRITE(6,'(A,F18.9)') 
     &  ' Energy through second order ........................ ',E012T
        WRITE(6,*)
        WRITE(6,'(A,F18.9)') 
     &  ' Norm of first order residual ....................... ',XNORMR
        WRITE(6,'(A,F18.9)') 
     &  ' Norm of unorthogonalized first order vector ........ ',XNORM1
        WRITE(6,'(A,F18.9)') 
     &  ' Overlap <0|(H0-E0)-1(H-E)|0> ....................... ',S01T
*
        IF(IHYLLE2.EQ.1) THEN
        WRITE(6,'(A,F18.9)') 
     &  ' Hylleraas second order energy ...................... ',E2TH
        WRITE(6,'(A,F18.9)') 
     &  ' Energy through second order ( Hylleraas )  ......... ',E012TH
        END IF
C     END IF
      CALL GFLUSH(6)
*
      IF(LEVEL.GE.2 .AND. IC1DSC.GT.0 .AND. ALPHA1 .NE. 0.0D0 ) THEN
*. Construct complete first order correction and store on disc
        WRITE(6,*)
        WRITE(6,*) ' ==============================================='
        WRITE(6,*) ' Construction of complete first order correction'
        WRITE(6,*) ' ==============================================='
        WRITE(6,*)
*
        HAPR01T = 0.0D0
        S11T = 0.0D0
        CALL REWINO(LU1)
        DO IBATCH = 1, NSBATCH
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' Information about sigma BATCH ', IBATCH
          END IF
*. Information about this batch
          NBLOCK = IFRMR(int_mb(KSLBT),1,IBATCH)
          NELMNT = IFRMR(int_mb(KSLEBT),1,IBATCH)
          IBLOCK_OFF = IFRMR(int_mb(KSI1BT),1,IBATCH)
*. Obtain blocks of first order corrections
          CALL GET_BATCH_OF_FIRST
     &         (NBLOCK,IBLOCK_OFF,NELMNT,int_mb(KSIBT),ISM,IUTSPC,
     &          E0,E1T,ALPHA1,LUC,IDC,
     &          int_mb(KCIOIO),NOCTPA,NOCTPB,NSMST,
     &          int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &          NTBLOCKC,int_mb(KCIBT),int_mb(KLXIBT),int_mb(KLXISCR),
     &          int_mb(KLXISC2),dbl_mb(KSB),dbl_mb(KCB),dbl_mb(KLB),
     &          SOUT,
     &          HAPR01,0,NTBLOCKS)
C KLXIBT, KLXISCR, KLXISC2
C     GET_BATCH_OF_FIRST
C    &           (N1BLOCK,I1OFF,N1ELMNT,I1BLOCK,I1SM,I1SPC,
C    &            E0,E1,ALPHA,LU0,IDC,
C    &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
C    &            N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,
C    &            HAPR01,LU1,N1BLOCKT )
*. On output : LB : Batch of first order correction
          S11T = S11T + INPROD(dbl_mb(KLB),dbl_mb(KLB),NELMNT)
          HAPR01T = HAPR01T + HAPR01
*. Save batch of first order corrections on DISC '
          CALL EXTRROW2(int_mb(KSIBT),8,IBLOCK_OFF,8,
     &                  NBLOCK,int_mb(KLXISCR))
          CALL TODSCN(dbl_mb(KLB),NBLOCK,int_mb(KLXISCR),LBLK,LU1)
        END DO
      END IF


      IF(LEVEL.GE.2) THEN 
        WRITE(6,*) 
        WRITE(6,'(20X,A)') ' *****************************'
        WRITE(6,'(20X,A)') '    One direct CI iteration : '
        WRITE(6,'(20X,A)') ' *****************************'
        WRITE(6,*) 
        WRITE(6,*) 
*
*
* Calculate <1|H|1> = sum(l,r) <0(l)|H|0(r)> 
*                   = sum(l)    <0(l)|H|0(l)>
*                   + 2sum(r.lt.l) <0(l)|H|0(r)>
*
*. and <1|H0|1>
*. Loop over batches of <1| and |1> 
        CALL REWINO(LUSCR)
        H11 = 0.0D0
        H011 = 0.0D0
        IF(IC1DSC.EQ.0) THEN
          HAPR01T = 0.0D0
        END IF
        S11T = 0.0D0
*
        IRESTRICT =  1
CTEST   IRESTRICT =  0
CTEST   WRITE(6,*) ' Restriction of batches removed for test '
*
        DO ILBATCH = 1, NSBATCH
          WRITE(6,*) ' Starting outer loop for left  batch ', ILBATCH
*. Construct this batch of first order corrections and store in LUSCR
          NLBLOCK = IFRMR(int_mb(KSLBT),1,ILBATCH)
          NLELMNT = IFRMR(int_mb(KSLEBT),1,ILBATCH)
          ILOFF = IFRMR(int_mb(KSI1BT),1,ILBATCH)
*
          IF(IRESTRICT.EQ.1) THEN
            IRMAX = ILBATCH
          ELSE
            IRMAX =  NSBATCH
          END IF
C         IRMAX = NSBATCH
          DO IRBATCH  = 1, IRMAX   
*
            WRITE(6,*) 
     &      '     Information about batches (left and right)',
     &      ILBATCH,IRBATCH
            CALL GFLUSH(6)
*
            NRBLOCK = IFRMR(int_mb(KSLBT),1,IRBATCH)
            NRELMNT = IFRMR(int_mb(KSLEBT),1,IRBATCH)
            IROFF = IFRMR(int_mb(KSI1BT),1,IRBATCH)
* calculate first order corrections in R batch 
            ICSPC = INSPC
            ISSPC = IUTSPC
            CALL GET_BATCH_OF_FIRST
     &           (NRBLOCK,IROFF,NRELMNT,int_mb(KSIBT),ISM,IUTSPC,
     &            E0,E1T,ALPHA1,LUC,IDC,
     &            int_mb(KCIOIO),NOCTPA,NOCTPB,NSMST,
     &            int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &            NTBLOCKC,int_mb(KCIBT),int_mb(KLXIBT),int_mb(KLXISCR),
     &            int_mb(KLXISC2),dbl_mb(KSB),dbl_mb(KCB),dbl_mb(KLB),
     &            SOUT,
     &            HAPR01,LU1,NTBLOCKS)
C KLXIBT, KLXISCR, KLXISC2
*. calculate those blocks of H |1(right_batch)> that resides in 
*  left batch 
*
*. Transfer batch of first order correction to  LUSCR2 
*  ( Sblock/Sblocks assumes that blocks are  on disc)
C                EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
            CALL EXTRROW2(int_mb(KSIBT),8,IROFF,8,
     &                    NRBLOCK,int_mb(KLXISCR))
            CALL REWINO(LUSCR2)
            CALL TODSCNP(dbl_mb(KLB),NRBLOCK,int_mb(KLXISCR),-1,LUSCR2)
            CALL ITODS(-1,1,-1,LUSCR2)
*
*. Obtain (left) batch of H  |1(right_batch>  
*
            I12 = 2
            IPERTOP = 0
            IAPR = 0
*
            ICSPC = IUTSPC
            ISSPC = IUTSPC
*. Change the pointers to C so they are to S
*
            KCLBT_SAVE =  KCLBT
            KCLBLK_SAVE =  KCLBLK
            KCLEBT_SAVE =  KCLEBT
            KCI1BT_SAVE = KCI1BT
            KCIBT_SAVE =  KCIBT
            KC2B_SAVE = KC2B
            NCBATCH_SAVE = NCBATCH
*
            KCLBT =  KSLBT
            KCLBLK =  KSLBLK
            KCLEBT =  KSLEBT
            KCI1BT = KSI1BT
            KCIBT =  KSIBT
            KC2B = KS2B
            NCBATCH = NSBATCH
*
C?          WRITE(6,*) ' < Left batch !H! Right batch > '
*. We are only constructing half of the sigma, add Ecore later to avoid problems
            ECOREX = 0.0D0
            CALL SBLOCK(NLBLOCK,int_mb(KSIBT),ILOFF,dbl_mb(KCB),
     &                  dbl_mb(KSB),
     &                  LUSCR2,IRESTRICT,0,1,IRBATCH,IRBATCH,0.0D0,
     &                  ECOREX, CTASK )
C     SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
C    &                  ECORE,ITASK) 
*. And restore order
            KCLBT =  KCLBT_SAVE
            KCLBLK =  KCLBLK_SAVE
            KCLEBT =  KCLEBT_SAVE
            KCI1BT = KCI1BT_SAVE
            KCIBT =  KCIBT_SAVE
            KC2B = KC2B_SAVE
            NCBATCH = NCBATCH_SAVE
*
*. Add to previous obtained contributions to left batch of H !1>
            IF(IRBATCH.NE.1) THEN
              CALL REWINO(LUSCR)
              CALL FRMDSC(dbl_mb(KLB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
              ONE = 1.0D0
              CALL VECSUM(dbl_mb(KSB),dbl_mb(KSB),dbl_mb(KLB),ONE,ONE,
     &                    NLELMNT)
            END IF
*. In last final batch, multiply with two to obtain full term
            IF(IRESTRICT.EQ.1.AND.IRBATCH.EQ.IRMAX) THEN
              TWO = 2.0D0
              CALL SCALVE(dbl_mb(KSB),TWO,NLELMNT)
            END IF
*. And transfer to disc to save for future generations
            CALL REWINO(LUSCR)
            CALL TODSCP(dbl_mb(KSB),NLELMNT,-1,LUSCR)
          END DO
*         ^ End of loop over right batches
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Batch of H |1> '
            CALL WRTTTS2(dbl_mb(KSB),int_mb(KSIBT),NLBLOCK,ILOFF,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),IDC)
          END IF
*. construct left batch of first order corrections
          ICSPC = INSPC
          ISSPC = IUTSPC
          CALL GET_BATCH_OF_FIRST
     &         (NLBLOCK,ILOFF,NLELMNT,int_mb(KSIBT),ISM,IUTSPC,
     &          E0,E1T,ALPHA1,LUC,IDC,
     &          int_mb(KCIOIO),NOCTPA,NOCTPB,NSMST,
     &          int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &          NTBLOCKC,int_mb(KCIBT),int_mb(KLXIBT),int_mb(KLXISCR),
     &          int_mb(KLXISC2),dbl_mb(KSB),dbl_mb(KCB),dbl_mb(KLB),
     &          SOUT,
     &          HAPR01,LU1,NTBLOCKS)
C KLXIBT, KLXISCR, KLXISC2
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Left Batch of first order correction'
            CALL WRTTTS2(dbl_mb(KLB),int_mb(KSIBT),NLBLOCK,ILOFF,
     &                   NSMST,NOCTPA,NOCTPB,
     &                   int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),IDC)
          END IF

          IF(IC1DSC.EQ.0)  HAPR01T = HAPR01T + HAPR01
          S11T = S11T + INPROD(dbl_mb(KLB),dbl_mb(KLB),NLELMNT)
* <1|H0!1>
          CALL COPVEC(dbl_mb(KLB),dbl_mb(KSB),NLELMNT)
          FACTOR = 0.0D0
          ITASK = 2
          JPERT = 1
          CALL DIATERM_GAS(FACTOR,ITASK,dbl_mb(KSB),NLBLOCK,
     &         int_mb(KSIBT),ILOFF,JPERT,0,0)
          H011 = H011 + INPROD(dbl_mb(KLB),dbl_mb(KSB),NLELMNT)
*. <1|H|1>
          CALL REWINO(LUSCR)
          CALL FRMDSC(dbl_mb(KSB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
*. Update H11
          H11 = H11 + INPROD(dbl_mb(KLB),dbl_mb(KSB),NLELMNT)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' LB and SB vectors to H11 '
            CALL WRTMAT(dbl_mb(KLB),1,NLELMNT,1,NLELMNT)
            CALL WRTMAT(dbl_mb(KSB),1,NLELMNT,1,NLELMNT)
          END IF
        END DO
*. Normalize to unit norm of first order correction 
        XNORM = SQRT(S11T)
*
COLD    FACN = 1.0D0/XNORM1
        FACN = 1.0D0/XNORM
        H00 = E0+ E1T 
        H11 = H11/S11T + ECORE
        H01 = E2T/XNORM
C       H011 = H011 + ECORE*S11T
        WRITE(6,*) ' 2 by 2 matrix (orthonormal vectors) '
        WRITE(6,*) '======================================'
        WRITE(6,*) 
        WRITE(6,*) '          <0|H|0> '
        WRITE(6,*) '          <0|H|1>             <1|H|1> '
        WRITE(6,*) 
        WRITE(6,'(4X,F22.12)') H00
        WRITE(6,'(4X,2F22.12)') H01,H11
*. In honour of B. Bechman :
        A =  1.0D0
        B = -(H00+H11)
        C = H00*H11-H01**2
*. 
        EA = -B/(2*A) - SQRT(B**2 - 4*A*C)/(2*A)      
        EB = -B/(2*A) + SQRT(B**2 - 4*A*C)/(2*A)
        ELOW = MIN(EA,EB)
        EHIGH = MAX(EA,EB)
        EOUT = ELOW
        WRITE(6,'(A,2F22.12)') '   Eigenvalues ', ELOW,EHIGH
        WRITE(6,'(A)') '   Eigenvector for lowest eigenvalue '
        IF(H01.NE.0.0D0) THEN
          XLOW = (ELOW-H00)/H01
          XLNORM = SQRT(1.0D0+XLOW**2)
          FACTOR = 1.0D0/XLNORM
          WRITE(6,'(3X,2F22.12)')  FACTOR, XLOW*FACTOR 
        ELSE
          IF(H00.LE.H11) THEN
            WRITE(6,'(A)') '1 ,  0 '
          ELSE
            WRITE(6,'(A)') '0 ,  1 '
          END IF
        END IF
        CALL GFLUSH(6)
*. Third order energy 
        WRITE(6,*) 
        WRITE(6,'(20X,A)') 
     &  ' ****************************************'
        WRITE(6,'(20X,A)')
     &   '    Third order perturbation theory : '
        WRITE(6,'(20X,A)')
     &   ' ****************************************'
        WRITE(6,*) 
        WRITE(6,*) 
        H00 = E0 
*. We scaled H11 above, remove this
C       H11 = (H11 - ECORE)*S11T
        H11 = H11*S11T
        WRITE(6,*) ' <1| H0 |1> ', H011 
        WRITE(6,*) '  <1| H0 -E0 |1> ',H011-E0*S11T
C       E3T = (H11-E1T)*S11T - H011
        E3T = H11-E1T*S11T - H011
        WRITE(6,*) ' Third order energy correction ', E3T
        E0123T = E0 + E1T +  E2T + E3T
        WRITE(6,*) ' Energy through third order ',
     &              E0+E2T+E3T
          
        WRITE(6,*) '  <0| H apr | 1 > :' , HAPR01T
        WRITE(6,*) '   S11T = ', S11T 
*
        XNORM1 = SQRT(S11T)
*
        WRITE(6,*) ' Contributions to third order energy :  '
        WRITE(6,*) ' ===================================== '
        WRITE(6,*) 
        WRITE(6,'(A,F18.9)') 
     &  ' Third order energy contribution :................... ',E3T
        WRITE(6,*)
        WRITE(6,'(A,F18.9)') 
     &  ' Energy through third order :........................ ',E0123T
        WRITE(6,'(A,F18.9)') 
     &  ' Norm of orthogonalized first order vector ..::...... ',XNORM1
      END IF
*
      IF(LEVEL.GE.3) THEN 
        WRITE(6,*) 
        WRITE(6,'(20X,A)') 
     &  ' *****************************************'
        WRITE(6,'(20X,A)')
     &   '    Fourth order perturbation theory : '
        WRITE(6,'(20X,A)')
     &   ' *****************************************'
        WRITE(6,*) 
        WRITE(6,*) 
*. The fourth energy looks - in its most stable form as
*  E4 = <2|H0-E0|2> + 2<1|V-E1|2> - E2<1|1>
*. But the most easy - at hte moment is 
*. E4 = <1|V-E1|2> - E2<1|1>
        IWAY = 3
        WRITE(6,*) ' expression used for fourth order energy '
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' E4 = <2|H0-E0|2> + 2<1|V-E1|2> - E2<1|1> '
        ELSE IF(IWAY.EQ.2) THEN
          WRITE(6,*) ' E4 = <1|V-E1|2> - E2<1|1> '
        ELSE IF(IWAY.EQ.3) THEN
          WRITE(6,*) ' E4 =-<2|H0-E0|2> - E2<1|1> '
        END IF
*. Loop over batches of second order correction
        CALL REWINO(LUSCR)
        H022 = 0.0D0
        V12 = 0.0D0
        S22U = 0.0D0
        S02T= 0.0D0
*
        DO IBATCH2 = 1, NBATCH2
          WRITE(6,*) ' Starting outer loop for batch ', IBATCH2
*. Construct this batch of first order corrections and store in LUSCR
          NBLOCK2 = IFRMR(int_mb(KLLBT2) ,1,IBATCH2)
          NELMNT2 = IFRMR(int_mb(KLLEBT2),1,IBATCH2)
          IOFF2 = IFRMR(int_mb(KLI1BT2)  ,1,IBATCH2)
          write(6,*) ' NBLOCK2 NELMNT2 IOFF2 ',
     &                 NBLOCK2,NELMNT2,IOFF2
*
* calculate second order corrections in  batch
*
C         HAPR01TX = 0.0D0
          CALL GET_BATCH_OF_SECOND
     &         (NBLOCK2,IOFF2,NELMNT2,int_mb(KLIBT2),ISM,ISPC2,
     &         NSBATCH,int_mb(KLSLBT),int_mb(KLSLEBT),int_mb(KLSIBT),
     &         ISM,ISPC1,
     &         NBATCHC,int_mb(KLCLBT),int_mb(KLCIBT),ISPC0,
     &         E0,E1T,E2T,ALPHA1,LUC,LUSCR2,IDC,
     &         int_mb(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &         int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &         int_mb(KLIBT),dbl_mb(KLISCR),dbl_mb(KLISCR2),
     &         dbl_mb(KLB),dbl_mb(KCB),dbl_mb(KSB),SOUT,HAPR01T,LU1,
     &         NTBLOCKS)
          S22U = S22U + INPROD(dbl_mb(KSB),dbl_mb(KSB),NELMNT2)
          S02T = S02T + SOUT
*. Well now we have the second order correction in SB , save it !
          CALL REWINO(LUSCR)
          CALL TODSC(dbl_mb(KSB),NELMNT2,-1,LUSCR)
*. The term <2|H0-E0|2>
* (H0-E0)|2>
          CALL COPVEC(dbl_mb(KSB),dbl_mb(KLB),NELMNT2)
          FACTOR = - E0   
          ITASK = 2
          ECOREX = 0.0D0
          IUTSPC = ISPC2
          JPERT = 1
          CALL DIATERM_GAS(FACTOR,ITASK,dbl_mb(KLB),NBLOCK2,
     &         int_mb(KLIBT2),IOFF2,JPERT,0,0)
C         DIATERM_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
C    &                     ISM,ISPC,ECORE,I12)
          H022 = H022 + INPROD(dbl_mb(KLB),dbl_mb(KSB),NELMNT2)
*. The terms E2<1|2> and <1|V|2> : Loop over batches of first correction 
          IF(IWAY.NE.3) THEN
          DO IBATCH1 = 1, NSBATCH
*
            WRITE(6,*)
     &      '     Information about batches (First - and second order)',
     &      IBATCH1,IBATCH2
*
            NBLOCK1 = IFRMR(int_mb(KLSLBT),1,IBATCH1)
            NELMNT1 = IFRMR(int_mb(KLSLEBT),1,IBATCH1)
            IOFF1 = IFRMR(int_mb(KLSI1BT),1,IBATCH1)
* calculate first order corrections in batch
            ICSPC = ISPC0
            ISSPC = ISPC1
            CALL GET_BATCH_OF_FIRST
     &           (NBLOCK1,IOFF1,NELMNT1,int_mb(KLSIBT),ISM,ISPC1,
     &            E0,E1T,ALPHA1,LUC,IDC,
     &            int_mb(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &            int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &            NTBLOCKC,int_mb(KLCIBT),int_mb(KLIBT),int_mb(KLI1BT),
     &            dbl_mb(KLISCR),dbl_mb(KSB),dbl_mb(KCB),dbl_mb(KLB),
     &            SOUT,
     &            HAPR01,LU1,NTBLOCKS)
*. First order correction is now in WORK(KLB)
*. save on LUSCR after second order batch
*. Position
            CALL REWINO(LUSCR)
            CALL FRMDSC(dbl_mb(KCB),NELMNT2,-1,LUSCR,IMZERO,IAMPACK)
*
            CALL TODSC(dbl_mb(KLB),NELMNT1,-1,LUSCR)
*<2|H|1>
*
            ICSPC = ISPC1
            ISSPC = ISPC2
*
            ZERO = 0.0D0
            CALL SETVEC(dbl_mb(KCB),ZERO,NELMNT2)
            I12 = 2
            IREST = 0
            CALL  SBATCH2_FROM_CBATCH(NBLOCK2,IOFF2,NBLOCK1,IOFF1,
     &                                 int_mb(KLIBT2),int_mb(KLSIBT),
     &                                 dbl_mb(KLB),dbl_mb(KCB),IREST)
*. Result is in in KCB
*. read second order batch into KSB
            CALL REWINO(LUSCR)
            CALL FRMDSC(dbl_mb(KSB),NELMNT2,-1,LUSCR,IMZERO,IAMPACK)
*
            V12 = V12 + INPROD(dbl_mb(KSB),dbl_mb(KCB),NELMNT2)
*. Obtain <1|H0-E1|2> and subtract
*. Second order blocks that are in first order  batch :
            CALL GET_TTS_BLK_IN_VECTOR(NBLOCK2,int_mb(KLIBT2),IOFF2,
     &           dbl_mb(KSB),NBLOCK1,int_mb(KLSIBT),IOFF1,
     &           NBLOCK12,NELMNT12,int_mb(KLIBT),dbl_mb(KLB))
C           GET_TTS_BLK_IN_VECTOR
C    &             (NBLOCKI,IBLOCKI,IOFFI,VECI,
C    &              NBLOCKO,IBLOCKO,IOFFO,
C    &              NBLOCKIO,IBLCKIO,VECIO )
*. Extracted Second order blocks are now in KLB
*  H0+E1|2>
            FACTOR = E1T
            ITASK = 2
            ECOREX = 0.0D0
            JPERT = 1
            CALL DIATERM_GAS(FACTOR,ITASK,dbl_mb(KLB),NBLOCK12,
     &           int_mb(KLIBT),1,JPERT,0,0)
*. extract the corresponding first order blocks
*. Read in first order corrections in KCB
            CALL FRMDSC(dbl_mb(KCB),NELMNT1,-1,LUSCR,IMZERO,IAMPACK)
*. And extract the common blocks
            CALL GET_TTS_BLK_IN_VECTOR(NBLOCK1,int_mb(KLSIBT),IOFF1,
     &           dbl_mb(KCB),NBLOCK12,int_mb(KLIBT),1,
     &           NBLOCK12A,NELMNT12A,int_mb(KLIBT),dbl_mb(KSB))
*. Common   first order blocks are now in KSB
            V12 = V12 - INPROD(dbl_mb(KLB),dbl_mb(KSB),NELMNT12)
          END DO
          END IF
*.End of loop over batches of first order corrections
        END DO 
*.End of loop over batches of second order corrections
*. The complete second order correction was not obtained,
*. It reads : |2> = |2'> - alpha2(H apr -E0 )-1 |0>
*. Only |2'> was obtained above
*. correct the V12 term
*. And <2|H0-E0|2>
        IF(IHAPREIG.NE.0) THEN
         V12 = V12
         ALPHA2 = 0.0D0
        ELSE
         ALPHA2 = S02T/HAPRM100
         V12 = V12 - ALPHA2*(-S02T+(E2T-HAPR01T)*HAPRM100)
         H022 = H022 -2*ALPHA2*S02T + ALPHA2**2 * HAPRM100
        END IF
*. and the norm
        S22 = S22U - S02T**2
        WRITE(6,'(A,F18.9)') 
     &  ' Norm of orthogonalized second order correction vector ',
     &  SQRT( S22 ) 
*. The time has come to present the conclusions :
        WRITE(6,*) 
        WRITE(6,*) ' Contributions to fourth order energy :  '
        WRITE(6,*) ' ======================================= '
        WRITE(6,*)
        IF(IWAY.EQ.1) THEN
          WRITE(6,'(A,F18.9)') 
     &    '    <2|H0-E0|2> ..................... ', H022
          WRITE(6,'(A,F18.9)') 
     &    '    2 * <1|V-E1|2> .................. ', 2.0D0*V12
          WRITE(6,'(A,F18.9)') 
     &    '    -E2*<1|1> ....................... ', -E2T*S11T
          E4T = H022 + 2.0D0*V12 - E2T*S11T
        ELSE IF(IWAY.EQ.2) THEN
          WRITE(6,'(A,F18.9)')
     &    '    <1|V-E1|2> ...................... ', V12
          WRITE(6,'(A,F18.9)')
     &    '    -E2*<1|1> ....................... ', -E2T*S11T
          E4T = V12 - E2T*S11T
        ELSE IF (IWAY.EQ.3) THEN
          WRITE(6,'(A,F18.9)')
     &    '   -<2|H0-E0|2> ..................... ',-H022
          WRITE(6,'(A,F18.9)')
     &    '    -E2*<1|1> ....................... ', -E2T*S11T
          E4T = -H022 - E2T*S11T
        END IF
*
        WRITE(6,*)
        WRITE(6,'(A,F18.9)') 
     &  ' Fourth order energy contribution ... ',E4T
        WRITE(6,*)
        E01234T =  E0+ECORE +E1T+E2T+E3T+E4T
        EOUT = E01234T
        WRITE(6,'(A,F18.9)') 
     &  ' Energy through fourth order ........ ',E01234T            
        WRITE(6,*)

        WRITE(6,'(A,E18.9)') 
     &  ' overlap between 2(unortho) and 0 ........ ',S02T
         
      END IF
*. End of fourth order part 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DIRDR1')
      RETURN
      END 
      SUBROUTINE HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                   IABSPC,JABSPC,IABOP,JABOP,IIF,JDOH2,
     &                   IDOH2,IMZERO,IDIAG)
*
* Decide upon the treatment of matrix element
*
* <IASM IATP IBSM IBTP | H(apr) | JASM, JATP JBSM, JBTP>
*
* and do preparations (IIF = 1 )
* or counteract preparations (IIF = 2)
*
* Jeppe Olsen, The last day of January 1996
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      IF((IABSPC.NE.JABSPC).OR.
     &   (IABSPC.EQ.JABSPC.AND.(IABOP.EQ.1.OR.IABOP.EQ.2).AND.
     &    (JASM.NE.IASM.OR.IATP.NE.JATP.OR.IBTP.NE.JBTP))) THEN
*. Zero
        IMZERO = 1
      ELSE
*. Not zero !
        IMZERO = 0
      END IF
*. Diagonal approximation?
        IF(IABOP.EQ.1.OR.IABOP.EQ.2) THEN
          IDIAG = 1
        ELSE
          IDIAG = 0
        END IF
*. Moller Plesset or normal operator ?
        IF(IABOP.EQ.1.OR.IABOP.EQ.3.OR.IABOP.EQ.5)THEN
          IMP = 1
        ELSE
          IMP = 0
        END IF
C     END IF
*. Two - or one- electron operator
      IF(IABOP.EQ.1.OR.IABOP.EQ.3) THEN
        IDOH2 = 0
      ELSE
C       IDOH2 = JDOH2
        IDOH2 = 1
      END IF
*. Put MP integrals in place ( Or put good old one-electron integrals 
*. back where they belong
      IF(IMP.EQ.1) THEN
c       CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        ktmp = kint1
        kint1 = kfi
        kfi = ktmp
c       CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
        ktmp = kint1o
        kint1o = kfio
        kfio = ktmp
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
         WRITE(6,*) ' HMATAPR speaking '
         WRITE(6,*) ' ================='
         WRITE(6,*) ' | IASM IATP IBSM IBTP > :',
     &              '|',IASM,IATP,IBSM,IBTP,'>'
         WRITE(6,*) ' | JASM JATP JBSM JBTP > :',
     &              '|',JASM,JATP,JBSM,JBTP,'>'
         WRITE(6,*) ' IABSPC,JABSPC :', IABSPC,JABSPC
         WRITE(6,*) ' IABOP ', IABOP
         WRITE(6,*) ' Results : IMP IDIAG IDOH2 IMZERO : ',
     &                IMP,IDIAG,IDOH2, IMZERO 
      END IF
*
      RETURN
      END
      SUBROUTINE RSSBCBN(IASM,IATP,IOCPTA,
     &                  IBSM,IBTP,IOCTPB,
     &                  JASM,JATP,JBSM,JBTP,NGAS,
     &                  IAOC,IBOC,JAOC,JBOC, 
     &                  NAEL,NBEL,
     &                  IJAGRP,IJBGRP,
     &                  SB,CB,JDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IDC,PS,
     &                  ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &                  MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,IPRNT,IHAPR,
     &                  IPTSPC,JPTSPC,IJOP,NNSEL2E,ISEL2E,SCLFAC)
*
* Contributions to sigma block (iasm iatp, ibsm ibtp ) from
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS      : Number of active spaces in calculation
* IAOC,IBOC : Number of electrons in each AS for sigma supergroups
* JAOC,JBOC : Number of electrons in each AS for C     supergroups
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* IDOH2 : = 0 => no two electron operator
* IDOH2 : = 1 =>    two electron operator
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ICJKAIB =1 =>  construct C(Ka,Jb,j) and S(Ka,Ib,i) as intermediate 
*                 matrices in order to reduce overhead
*
* IHAPR : .ne. 0 implies thatt the exact Hamiltonian shoulf not be uses
* In the case IPTSPC and JPTSPC defined the perturbation spaces
* a nonvanishing perturbation is allowed inside each subspace.
* The actual type of approximate Hamiltonian in each subspace is defined by
* IHFORM
* NNSEL2E : Only selected 2e terms will be included
* ISEL2E : orbital spaces in which 2e terms are included
*          (Currently : all indeces identical )
*
* ======
* Output
* ======
* SB : fresh sigma block
*
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*. Output
      DIMENSION CB(*),SB(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*)
      DIMENSION C2(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION IBLOCK(8)
*. For H(apr)
      DIMENSION ISEL2E(*)
COLD  DIMENSION IHSPC(NOCPTA,NOCTPB)
COLD  DIMENSION IHFORM(*)
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
*
C?    WRITE(6,*) ' Memcheck entering RSSBCB '
C?    CALL MEMCHK
      IF(NTEST.GE.000) THEN
        WRITE(6,*) ' =================='
        WRITE(6,*) ' RSSBCBN :  C block '
        WRITE(6,*) ' ==================='
        IF(ICJKAIB.EQ.0) THEN
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
        ELSE
        CALL WRTMAT(CB,NJB,NJA,NJB,NJA)
        END IF
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' RSSBCBN : Initial  S block '
        WRITE(6,*) ' =========================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
        WRITE(6,*) ' Overall scalefactor ',SCLFAC
      END IF
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' IAOC and IBOC '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' JAOC and JBOC  : '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' IHAPR and IDIAG ', IHAPR,IDIAG
        WRITE(6,*) ' IASM IATP JASM JATP ', IASM,IATP,JASM,JATP
        WRITE(6,*) ' IBSM IBTP JBSM JBTP ', IBSM,IBTP,JBSM,JBTP
        WRITE(6,*) ' NAEL NBEL ', NAEL, NBEL
      END IF
* Should the corresponding Hamiltonian matrix block be 
* calculated exactly or approximately
      IF(IHAPR.NE.0) THEN
C       HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
C    &                   IABSPC,JABSPC,IABOP,JABOP,IIF,JDOH2,
C    &                   IDOH2,IMZERO,IDIAG)
        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &               IMZERO,IDIAG)
COLD    WRITE(6,*) ' RSSBCBN : ', NNSEL2E, ISEL2E(1)
        NSEL2E = NNSEL2E
        IF(IMZERO.NE.0) GOTO 9999
      ELSE
*. Operator specified by input 
        IAPRLEV =-1    
        IDOH2 = JDOH2
        IDIAG = 0
        NSEL2E = 0
      END IF
*
*
      IF(IDC.EQ.2.AND.IATP.EQ.IBTP.AND.IASM.EQ.IBSM .AND.
     &            JASM.EQ.JBSM.AND.JATP.EQ.JBTP) THEN
*. Diagonal sigma block, use alpha-beta symmetry to reduce computations.
        IUSEAB = 1
      ELSE
        IUSEAB = 0
      END IF
*
      IF(IDIAG.EQ.0) THEN
*
* Nondiagonal operator    
*
      IF(IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
* Sigma beta beta contribution
* =============================
*
* Sigma aa(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*                 * ((ij!kl)-(il!kj)) C(IA,JB)
*                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*.One electron part
*. If ICJKAIB is active matrices are transposed, so back transpose
        IF(ICJKAIB.NE.0) THEN
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
   
        IF(NBEL.GE.1) THEN
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB before RSBB1E'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
          IF(NTEST.GE.200.and.ga_nodeid().eq.0)
     &    WRITE(6,*) ' I am going to call RSBB1E'
          CALL RSBB1E(IBSM,IBTP,IOCTPB,JBSM,JBTP,IOCTPB,
     &         IJBGRP,NIA,
     &         NGAS,IBOC,JBOC,
     &         SB,CB,
     &         ADSXA,SXSTST,STSTSX,
     &         MXPNGAS,NOBPTS,IOBPTS,
     &         ITSOB,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &         NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
     &         NIB,MXSXST,CJRES(1),CJRES(1+MXSXBL),MOCAA,SCLFAC,NTEST)
C?    WRITE(6,*) ' Memcheck after RSBB1E    '
C?    CALL MEMCHK
           IF(NTEST.GE.101.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' SB after RSBB1E, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' SB after RSBB1E'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
        END IF
        IF(IDOH2.NE.0.AND.NBEL.GE.2) THEN
*. Two electron part
          IF(NTEST.GE.100.and.ga_nodeid().eq.0)
     &    WRITE(6,*) ' I am going to call RSBB2A'
          CALL RSBB2A(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NIB,
     &                NGAS,IBOC,JBOC,                
     &                SB,CB,
     &                ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &                NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &                CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC,NTEST,
     &                NSEL2E,ISEL2E)
C?    WRITE(6,*) ' Memcheck after RSBB2A    '
C?    CALL MEMCHK
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB after RSBB2A, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' SB after RSBB2a'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
        END IF
*. If ICJKAIB is active matrices are transposed, so back transpose
        IF(ICJKAIB.NE.0) THEN
          CALL TRPMT3(SB,NIA,NIB,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
      END IF
*
* =============================
* Sigma alpha beta contribution
* =============================
*
      IF(IDOH2.NE.0.AND.NAEL.GE.1.AND.NBEL.GE.1) THEN
        IF(NTEST.GE.101)
     &  WRITE(6,*) ' I am going to call RSBB2B'
        IIITRNS = 1
        IF(IIITRNS.EQ.1.AND.NIB.GT.NIA.AND.NJB.GE.NJA) THEN
           JJJTRNS = 1
        ELSE
           JJJTRNS = 0
        END IF
*
        IF (JJJTRNS.EQ.0) THEN
          CALL RSBB2B(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                IJAGRP,IJBGRP,NGAS,
     &                IAOC,IBOC,JAOC,JBOC,
     &                SB,CB,
     &                ADSXA,STSTSX,MXPNGAS,
     &                NOBPTS,IOBPTS,ITSOB,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &                IUSEAB,ICJKAIB,CJRES,SIRES,C2,SCLFAC,NTEST,
     &                NSEL2E,ISEL2E)
*
         ELSE IF ( JJJTRNS.EQ.1) THEN
*. well lets give the transpose routine some more practice : Transpose back
          CALL TRPMAT(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
*
          CALL TRPMAT(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
*
          CALL RSBB2B(IBSM,IBTP,IASM,IATP,NIB,NIA,
     &                JBSM,JBTP,JASM,JATP,NJB,NJA,
     &                IJBGRP,IJAGRP,NGAS,
     &                IBOC,IAOC,JBOC,JAOC,
     &                SB,CB,
     &                ADSXA,STSTSX,MXPNGAS,
     &                NOBPTS,IOBPTS,ITSOB,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &                IUSEAB,ICJKAIB,CJRES,SIRES,C2,SCLFAC,NTEST,
     &                NSEL2E,ISEL2E)

*. Transpose ( To compensate later transposition )
          CALL TRPMAT(SB,NIA,NIB,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMAT(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
C?    WRITE(6,*) ' Memcheck after RSBB2B    '
C?    CALL MEMCHK
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB after RSBB2B, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' SB after RSBB2b'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
      END IF
*
* =============================
* Sigma alpha alpha contribution
* =============================
*
      IF(IUSEAB.EQ.0) THEN
*. Transpose for alpha excitations
      IF(NAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
        IF(ICJKAIB.EQ.0) THEN
          CALL TRPMT3(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
          CALL TRPMT3(SB,NIA,NIB,C2)
          CALL COPVEC(C2,NIA*NIB,SB)
        END IF
*
* alpha single excitation
*
        IF(NTEST.GE.101)
     &  WRITE(6,*) ' I am going to call RSBB1E (last time )'
        CALL RSBB1E(IASM,IATP,IOCTPA,JASM,JATP,IOCTPA,
     &                   IJAGRP,NIB,
     &                   NGAS,IAOC,JAOC,
     &                   SB,CB,
     &                   ADSXA,SXSTST,STSTSX,
     &                   MXPNGAS,NOBPTS,IOBPTS,
     &                   ITSOB,MAXI,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                   NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
     &                   NIA,MXSXST,CJRES(1),CJRES(1+MXSXBL),
     &                   MOCAA,SCLFAC,NTEST)
C?    WRITE(6,*) ' Memcheck after RSBB1E, 2 '
C?    CALL MEMCHK
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB transposed after RSBB1, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' SB transposed  after RSBB1E'
             call wrtmat(SB,nib,nia,nib,nia)
           END IF
*
* alpha double excitation
*
        IF(IDOH2.NE.0.AND.NAEL.GE.2) THEN
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB2A (last time )'
          CALL RSBB2A(IASM,IATP,JASM,JATP,IJAGRP,NIB,NIA,
     &         NGAS,IAOC,JAOC,  
     &         SB,CB,
     &         ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &         NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC,NTEST,
     &         NSEL2E,ISEL2E)
C?    WRITE(6,*) ' Memcheck after RSBB2A, 2 '
C?    CALL MEMCHK
        END IF
        IF(ICJKAIB.EQ.0) THEN
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NJA*NJB)
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
*
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB transposed after RSBB2A, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' SB after RSBB2A'
          call wrtmat(sb,nia,nib,nia,nib)
        END IF
      END IF
      END IF
*
      ELSE IF (IDIAG.EQ.1) THEN
*. Diagonal approxiation           
       IBLOCK(1) = IATP
       IBLOCK(2) = IBTP
       IBLOCK(3) = IASM
       IBLOCK(4) = IBSM
       IBLOCK(5) = 1
       IBLOCK(6) = 1
       IF(IDOH2.EQ.0) THEN
         I12 = 1
       ELSE
         I12 = 2
       END IF
       WRITE(6,*) ' IDOH2, I12 ', IDOH2,I12
       ITASK = 2
       FACTOR = 0.0D0
*. Input is in det basis 
       IIDC = 1
*. Well, we are working with transposed matrices so
       CALL TRPMT3(CB,NIB,NIA,C2)
       WRITE(6,*) ' DIATERM2_GAS will be called '
C     DIATERM2_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
C    &                       JPERT,J12,JDC)
       CALL DIATERM2_GAS(FACTOR,ITASK,C2,1,IBLOCK,1,0,I12,IIDC)
       CALL TRPMT3(C2,NIA,NIB,CB)
C Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
C Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)

C      CALL VVTOV(CB,C2,C2,NIA*NIB)
       ONE = 1.0D0
*. what to do when IUSEAB = 1 is in use ???
       IF(IUSEAB.EQ.0) THEN
         FACTOR = 1.0D0
       ELSE
         FACTOR = 0.5D0
       END IF
       CALL VECSUM(SB,SB,CB,ONE,FACTOR,NIA*NIB)
      END IF
*
 9999 CONTINUE
      IF(IHAPR.NE.0) THEN
*. Clean up
C       HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
C    &                   IABSPC,JABSPC,IABOP,JABOP,IIF,JDOH2,
C    &                   IDOH2,IMZERO,IDIAG)
        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &               IMZERO,IDIAG)
      END IF
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' RSSBCB : Final S block '
        WRITE(6,*) ' =========================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
      END IF
      NTESTO = NTEST
C?    STOP ' Jeppe forced me to stop in RSSBCB '
      RETURN
      END
      SUBROUTINE WRTBLKN_EP(VEC,NBLOCK,LBLOCK)
*
* Write the NBLOCK blocks of VEC
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*)
      DIMENSION LBLOCK(NBLOCK)
      REAL*8 INPROD
      INTEGER, PARAMETER :: MAXLEN = 500
*
      IOFF = 1
      DO IBLOCK = 1, NBLOCK
        IF( LBLOCK(IBLOCK).GT.0) THEN
          XNRM2 = INPROD(VEC(IOFF),VEC(IOFF),LBLOCK(IBLOCK))
          WRITE(6,*) ' Block : ', IBLOCK
          WRITE(6,*) ' ==================='
C?        WRITE(6,*)
C?        WRITE(6,*) ' Length : ', LBLOCK(IBLOCK)
C?        WRITE(6,*) ' Norm   : ', SQRT(XNRM2)
C?        WRITE(6,*) ' Norm^2 : ', XNRM2
          WRITE(6,*)
          IF (LBLOCK(IBLOCK).GT.MAXLEN) THEN
            WRITE(6,*) 'printing only the first ',MAXLEN,' elements'
            CALL WRTMAT_EP(VEC(IOFF),1,MAXLEN,1,MAXLEN)
          ELSE
            LL = LBLOCK(IBLOCK)
            CALL WRTMAT_EP(VEC(IOFF),1,LL,1,LL)
          END IF
          IOFF = IOFF + LBLOCK(IBLOCK)
          WRITE(6,*)
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE WRTBLKN(VEC,NBLOCK,LBLOCK)
*
* Write the NBLOCK blocks of VEC
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*)
      DIMENSION LBLOCK(NBLOCK)
      REAL*8 INPROD
      INTEGER, PARAMETER :: MAXLEN = 500
*
      IOFF = 1
      DO IBLOCK = 1, NBLOCK
        IF( LBLOCK(IBLOCK).GT.0) THEN
          XNRM2 = INPROD(VEC(IOFF),VEC(IOFF),LBLOCK(IBLOCK))
          WRITE(6,*) ' Block : ', IBLOCK
          WRITE(6,*) ' ==================='
C?        WRITE(6,*)
C?        WRITE(6,*) ' Length : ', LBLOCK(IBLOCK)
C?        WRITE(6,*) ' Norm   : ', SQRT(XNRM2)
C?        WRITE(6,*) ' Norm^2 : ', XNRM2
          WRITE(6,*)
          IF (LBLOCK(IBLOCK).GT.MAXLEN) THEN
            WRITE(6,*) 'printing only the first ',MAXLEN,' elements'
            CALL WRITVE(VEC(IOFF),MAXLEN)
          ELSE
            CALL WRITVE(VEC(IOFF),LBLOCK(IBLOCK))
          END IF
          IOFF = IOFF + LBLOCK(IBLOCK)
          WRITE(6,*)
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
*
* Extract row IROW from integer matrix INMAT, starting from column
* ICOLOFF
*
* Jeppe Olsen, Winter 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION INMAT(NROW,ICOLOFF-1+NCOL)
      DIMENSION IOUTVEC(NCOL)
*
      DO ICOL = ICOLOFF, ICOLOFF - 1 +  NCOL
        IOUTVEC(ICOL-ICOLOFF+1) = INMAT(IROW,ICOL)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Output vector from EXTRROW '  
        WRITE(6,*) ' Extracted ROW ', IROW
        CALL IWRTMA(IOUTVEC,1,NCOL,1,NCOL)
      END IF
*
      RETURN 
      END 
      SUBROUTINE EXTRROW(INMAT,IROW,NROW,NCOL,IOUTVEC)
*
* Extract row IROW from integer matrix INMAT
*
* Jeppe Olsen, Winter 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION INMAT(NROW,NCOL)
      DIMENSION IOUTVEC(NCOL)
*
      DO ICOL = 1, NCOL
        IOUTVEC(ICOL) = INMAT(IROW,ICOL)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Output vector from EXTRROW '  
        WRITE(6,*) ' Extracted ROW ', IROW
        CALL IWRTMA(IOUTVEC,1,NCOL,1,NCOL)
      END IF
*
      RETURN 
      END 
      SUBROUTINE FRMDSCN(VEC,NREC,LBLK,LU)
*
* Read  VEC as multiple record file, NREC records read. 
*. If NREC .lt. 0, then records are read untill end of 
* vector mark ( aka -1) is obtained
*
* Modified Aug. 2005 ( NREC = -1 option added )
*
      IMPLICIT REAL*8(A-H,O-Z)
*. OUtput
      DIMENSION VEC(*)
*
      IOFF = 1
      IF(NREC.GE.0) THEN 
        DO IREC = 1, NREC
          CALL IFRMDS(LREC,1,LBLK,LU)
          CALL FRMDSC(VEC(IOFF),LREC,LBLK,LU,IMZERO,IAMPACK)
          IOFF = IOFF + LREC
        END DO
      ELSE
*. Loop untill -1 is obtained 
  999   CONTINUE
          CALL IFRMDS(LREC,1,LBLK,LU)
          IF(LREC.GE.0) THEN
            CALL FRMDSC(VEC(IOFF),LREC,LBLK,LU,IMZERO,IAMPACK)
            IOFF = IOFF + LREC
        GOTO 999
          END IF
*         ^ End of recordlength  is larger than 0
      END IF
*     ^ End of NREC switch
*
      RETURN
      END
      SUBROUTINE TODSCNP(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
* Only nonzero elements are written to DISC  
*
* Statistics are handled in TODSCP
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
        IF(LREC(IREC).GE.0) THEN
*. Normal complete record
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSCP(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
*. zero record of length -LREC(IREC)
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,1)
        END IF
         
*
      END DO
*
      RETURN
      END
      SUBROUTINE TODSCN(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
*
* Statistcs are handled in TODSC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSC(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,0)
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE GET_BLOCKS_FROM_DISC
     &           (LU,NBLOCK,IOFF,IBLOCK,NBLOCKT,C,IREW)
*
* Obtain blocks IOFF - IOFF + NBLOCK from file LU
*
* Jeppe Olsen, January 1996
*
c     IMPLICIT REAL*8(A-H,O-Z)
c     INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      DIMENSION C(*)
      DIMENSION IBLOCK(8,*)
*  
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'stinf.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*. Local Scratch
*. Should be uncommented if NEW_OR_OLD is set to 2
*
C     PARAMETER (MXPNBLK = 33000)
C     DIMENSION ISCR(3*MXPNBLK)
*
C?    WRITE(6,*) ' ENTERING GET_BLOCKS_FROM_DISC'
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
CM    IREW = 1
      ISCAL = 0
      NEW_OR_OLD=1
      IF(NEW_OR_OLD.EQ.2) THEN
        MXPNBLK = 0
        IF(NBLOCKT.GT.MXPNBLK) THEN
          WRITE(6,*) ' Increase parameter MXPNBLK in'
          WRITE(6,*) ' GET_BLOCKS_FROM_DISC '
          WRITE(6,*) ' Current and required values ',
     &     MXPNBLK,NBLOCKT
          STOP 'GET_BLOCKS_FROM_DISC'
        END IF
        CALL GET_TTS_BATCH(C,NBLOCK,IBLOCK(1,IOFF),NBLOCKT,IBLOCK,
     &                     NOCTPA,NOCTPB,NSMST,
     &                     int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &                     IDC,LU,ISCR,IREW,ISCAL)
      ELSE
        CALL GET_TTS_BATCHN(C,NBLOCK,IBLOCK(1,IOFF),NBLOCKT,IBLOCK,
     &                     NOCTPA,NOCTPB,NSMST,
     &                     int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &                     IDC,LU,IREW,ISCAL)
      END IF
      
*
      RETURN
      END

*


      FUNCTION LBATCH_F(LBLOCK,NBLOCK,IOFF,LBATCHMX,LBATCHA)
*
* obtain number of blocks in batch so length of batch is less of equal 
* to LBATCH
*
* Jeppe Olsen, January 1996
*
* Input 
* =====
* LBLOCK : Length of each block
* NBLOCK : Total number of BLOCKS
* IOFF   : First block in current batch
* LBATCHMX : Allowed number of elements in Batch
* LBATCH : Number of elements in BATCH
*
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION LBLOCK(NBLOCK)
*
      LENGTH = 0
      LENGTHB = 0

      DO IBLOCK = IOFF, NBLOCK
        IF(LENGTH + LBLOCK(IBLOCK).LE.LBATCHMX) THEN
          LENGTHB = LENGTHB + 1
          LENGTH = LENGTH +  LBLOCK(IBLOCK)
        ELSE
          GOTO 999
        END IF
      END DO
  999 CONTINUE
*
      LBATCH_F = LENGTHB
      LBATCHA = LENGTH
*
      IF(LENGTHB.EQ.0) THEN
        WRITE(6,*) ' Sorry mate '
        WRITE(6,*) ' Not a single block can be stored in batch'
        WRITE(6,*) ' Change batch size - or find BUG '
        STOP ' Problem in LBATCH_F'
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of blocks in BATCH ', LBATCH_F
      END IF
*
      RETURN
      END
      SUBROUTINE RSSBCB2(IASM,IATP,IOCTPA,IBSM,IBTP,IOCTPB,
     &                  JASM,JATP,JBSM,JBTP,NGAS,
     &                  IAOC,IBOC,JAOC,JBOC,NAEL,NBEL,
     &                  IJAGRP,IJBGRP,SB,CB,JDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,MXPNGASX,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBSX,IDC,PS,
     &                  ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &                  MXSXBL,MXSXST,MOCAA,MOCAB,IAPR_,IPRNT,IHAPR,
     &                  IPTSPC,JPTSPC,IJOP,NNSEL2E,ISEL2E,SCLFAC,
     &                  IUSE_PH,IPHGAS,I_RES_AB_,IUSE_PA,CJPA,SIPA,
     &                  IH_OCC_CONS,ISIMSYM,XINT2,MXADKBLK,INEW_CTT,
     &                  INEW_STT,IEND_CTT,CBTT,ITASK,ILEN)
*
* Contributions to sigma block (iasm iatp, ibsm ibtp ) from
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS      : Number of active spaces in calculation
* IAOC,IBOC : Number of electrons in each AS for sigma supergroups
* JAOC,JBOC : Number of electrons in each AS for C     supergroups
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* IDOH2 : = 0 => no two electron operator
* IDOH2 : = 1 =>    two electron operator
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* IHAPR : .ne. 0 implies thatt the exact Hamiltonian shoulf not be uses
* In the case IPTSPC and JPTSPC defined the perturbation spaces
* a nonvanishing perturbation is allowed inside each subspace.
* The actual type of approximate Hamiltonian in each subspace is defined by
* IHFORM
* NNSEL2E : Only selected 2e terms will be included
* ISEL2E : orbital spaces in which 2e terms are included
*          (Currently : all indeces identical )
* CBTT   : CB at start of TT block - used for SIMSYM 
*
* ======
* Output
* ======
* SB : fresh sigma block
*
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*               Modified July 2003, to allow use of several symmetry-
*               blocks in alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
      CHARACTER*6 ITASK
      INCLUDE 'multd2h.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'oper.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION CB(*),SB(*),CBTT(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*)
      DIMENSION C2(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION IBLOCK(8)
      DIMENSION IPHGAS(*)
      DIMENSION XINT(*),XINT2(*)
*. For H(apr)
      DIMENSION ISEL2E(*)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
*
C?    WRITE(6,*) ' Entering RSSBCB2, NTEST = ', NTEST
      CALL MEMCHK2('RSSBEN')
      CALL LUCIAQENTER('RSSBE')
*
CT    IF(ITASK(1:5).EQ.'SIGMA') THEN
CT    DO I = 1, 100
CT      WRITE(6,*) ' Fusk is going on '
CT    END DO
CT    CB(1) = 0.0D0
CT    CB(4) = 0.0D0
*
CT    END IF
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' RSSBCB2 speaking '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' ITASK = ', ITASK
      END IF
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ==============================='
        WRITE(6,*) ' RSSBCB2 :  C block (transposed)'
        WRITE(6,*) ' ================================'
        CALL WRTMAT(CB,NJB,NJA,NJB,NJA)
        WRITE(6,*) ' ======================================='
        WRITE(6,*) ' RSSBCB2 : Initial  S block(transposed) '
        WRITE(6,*) ' ======================================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
        WRITE(6,*) ' Overall scalefactor ',SCLFAC
        WRITE(6,*) ' IHAPR,JDOH2 = ', IHAPR,JDOH2
        WRITE(6,*) ' IUSE_PH,I_RES_AB_ = ', IUSE_PH,I_RES_AB_
        WRITE(6,*) ' INEW_CTT,INEW_STT, ISIMSYM = ',
     &               INEW_CTT,INEW_STT, ISIMSYM
      END IF
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' IAOC and IBOC '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' JAOC and JBOC  : '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' IASM IATP JASM JATP ', IASM,IATP,JASM,JATP
        WRITE(6,*) ' IBSM IBTP JBSM JBTP ', IBSM,IBTP,JBSM,JBTP
        WRITE(6,*) ' NAEL NBEL ', NAEL, NBEL
        WRITE(6,'(A,4(2X,I4))') ' NJA, NJB, NIA, NIB = ',
     &                            NJA, NJB, NIA, NIB
      END IF
* Should the corresponding Hamiltonian matrix block be 
* calculated exactly or approximately
      IF(IHAPR.NE.0) THEN
        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &               IMZERO,IDIAG)
        IF(NTEST.GE. 20) 
     &  WRITE(6,*) ' RSSBCBN : ', NNSEL2E, ISEL2E(1)
        NSEL2E = NNSEL2E
        IF(IMZERO.NE.0) GOTO 9999
      ELSE
*. Operator specified by input 
        IAPRLEV =-1
        IDOH2 = JDOH2
        IDIAG = 0
        NSEL2E = 0
      END IF
      IF(NTEST.GE. 20) 
     &WRITE(6,*) ' IHAPR, IDIAG IDOH2 ' , IHAPR,IDIAG, IDOH2
*
      IF(ITASK(1:5).EQ.'SIGMA') THEN
*
* Perform a standard sigma calculation
        IF(IDC.EQ.2.AND.IATP.EQ.IBTP.AND.IASM.EQ.IBSM .AND.
     &     I_RES_AB_.EQ.0.AND.JASM.EQ.JBSM.AND.JATP.EQ.JBTP) THEN
*. Diagonal sigma block, use alpha-beta symmetry to reduce computations.
          IUSEAB = 1
        ELSE
          IUSEAB = 0
          END IF
*
        IF(IDIAG.EQ.0.AND.IH_OCC_CONS.EQ.0) THEN
*
* =======================
* Calculate block exactly 
* =======================
*
        IF(I_RES_AB_.NE.1.AND.IUSEAB.EQ.0.
     &     AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
* Sigma beta beta contribution
* =============================
*
      CALL MEMCHK2('SIGBI ')
* Sigma aa(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*                 * ((ij!kl)-(il!kj)) C(IA,JB)
*                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*.One electron part
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
*. set common oper (just in case we have unrestricted orbitals)
          ISPCAS=2 ! beta beta
          IF(NBEL.GE.0) THEN
            IF(NTEST.GE.200) THEN
              WRITE(6,*) ' SB before RSBB1E'
              call wrtmat(sb,nia,nib,nia,nib)
            END IF
            IF(NTEST.GE.101.and.ga_nodeid().eq.0)
     &      WRITE(6,*) ' I am going to call RSBB1E(beta)'
            CALL RSBB1E(IBSM,IBTP,IOCTPB,JBSM,JBTP,IOCTPB,
     &           IJBGRP,NIA,NGAS,IBOC,JBOC,
     &           SB,CB,ADSXA,SXSTST,STSTSX,
     &           MXPNGAS,NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &           NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
     &           NIB,MXSXST,CJRES(1),CJRES(1+MXSXBL),MOCAA,SCLFAC,
     &           IUSE_PH,IPHGAS,NTEST)
            IF(NTEST.GE.101) THEN
              WRITE(6,*) ' SB after RSBB1E(beta)'
              call wrtmat(sb,nib,nia,nib,nia)
            END IF
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' first element of SB after RSBB1E',
     &        SB(1) 
            END IF
          END IF ! NBEL.GE.0
          IF(IDOH2.NE.0.AND.NBEL.GE.0) THEN
*. Two electron part
            IF(NTEST.GE.101)
     &      WRITE(6,*) ' I am going to call RSBB2A(beta-beta)'
      CALL MEMCHK2('SIGBBI')
            CALL RSBB2A(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NIB,
     &           NGAS,IBOC,JBOC,SB,CB,
     &           ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &           NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &           NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &           CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC,NTEST,0,0,
     &           IUSE_PH,IPHGAS,XINT2)
             IF(NTEST.GE.200) THEN
               WRITE(6,*) ' SB after RSBB2a(beta-beta)'
               call wrtmat(sb,nib,nia,nib,nia)
             END IF
          END IF ! IDOH2.NE.0.AND.NBEL.GE.0
          CALL TRPMT3(SB,NIA,NIB,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF ! I_RES_AB_.NE.1.AND.IUSEAB.EQ.0
               !  .AND.IATP.EQ.JATP.AND.JASM.EQ.IASM
*
* =============================
* Sigma alpha beta contribution
* =============================
*
      CALL MEMCHK2('SIGABI')
*. set common oper (just in case we have unrestricted orbitals)
        ISPCAS=3 ! alpha beta
        IF(IDOH2.NE.0.AND.NAEL.GE.0.AND.NBEL.GE.0) THEN
* 
          IF(NTEST.GE.101)  
     &        WRITE(6,*) ' I am going to call RSBB2B'
          IIITRNS = 1
          IF(IIITRNS.EQ.1) THEN
*. Call advice routine
C     ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
C            CALL ADVICE_SIGMA(IAOC,IBOC,JAOC,JBOC,1,LADVICE)
             CALL ADVICE_SIGMA3(IAOC,IBOC,JAOC,JBOC,1,LADVICE,
     &                          NIA,NIB,NJA,NJB)
*. LADVICE = 2 => implies transpose
             IF(LADVICE.EQ.2) THEN  
               JJJTRNS = 1
             ELSE
               JJJTRNS = 0
             END IF
          ELSE
             JJJTRNS = 0
          END IF
*
C!        JJJTRNS = 0
          IVNEW = 1
          ISIMSYM = 0
          IF(ISIMSYM.EQ.0) THEN
*. Route with one TTS block at a time 
            IF (JJJTRNS.EQ.0) THEN
              IF(IVNEW.EQ.0) THEN
                CALL RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &             JASM,JATP,JBSM,JBTP,NJA,NJB,
     &             IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &             SB,CB,ADSXA,STSTSX,MXPNGAS,
     &             NOBPTS,MAXK,
     &             SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &             XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &             IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &             XINT2)
              ELSE 
                CALL RSBB2BN2(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &             JASM,JATP,JBSM,JBTP,NJA,NJB,
     &             IJAGRP,IJBGRP,IOCTPA,IOCTPB,NGAS,IAOC,IBOC,JAOC,JBOC,
     &             SB,CB,ADSXA,STSTSX,MXPNGAS,
     &             NOBPTS,MAXK,
     &             SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &             XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &             IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &             XINT2,NSTFSMSPGP,
     &             ISPGPFTP(1,IATP+IOCTPA-1),ISPGPFTP(1,IBTP+IOCTPB-1),
     &             ISPGPFTP(1,JATP+IOCTPA-1),ISPGPFTP(1,JBTP+IOCTPB-1)) 
              END IF ! IVNEW switch
*
            ELSE IF ( JJJTRNS.EQ.1) THEN
*. well lets give the transpose routine some more practice : Transpose back
              CALL TRPMT3(SB,NIB,NIA,C2)
              CALL COPVEC(C2,SB,NIA*NIB)
              CALL TRPMT3(CB,NJB,NJA,C2)
              CALL COPVEC(C2,CB,NJA*NJB)
C             WRITE(6,*) ' RSSBCB2 : Transpose path choosen'
              IF(IVNEW.EQ.0) THEN
                CALL RSBB2BN(IBSM,IBTP,IASM,IATP,NIB,NIA,
     &             JBSM,JBTP,JASM,JATP,NJB,NJA,
     &             IJBGRP,IJAGRP,NGAS,IBOC,IAOC,JBOC,JAOC,
     &             SB,CB,ADSXA,STSTSX,MXPNGAS,
     &             NOBPTS,MAXK,
     &             SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &             XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &             IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &             XINT2)

              ELSE 
                CALL RSBB2BN2(IBSM,IBTP,IASM,IATP,NIB,NIA,
     &             JBSM,JBTP,JASM,JATP,NJB,NJA,
     &             IJBGRP,IJAGRP,IOCTPB,IOCTPA,NGAS,IBOC,IAOC,JBOC,JAOC,
     &             SB,CB,ADSXA,STSTSX,MXPNGAS,
     &             NOBPTS,MAXK,
     &             SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &             XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &             IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &             XINT2,NSTFSMSPGP,
     &             ISPGPFTP(1,IBTP+IOCTPB-1),ISPGPFTP(1,IATP+IOCTPA-1),
     &             ISPGPFTP(1,JBTP+IOCTPB-1),ISPGPFTP(1,JATP+IOCTPA-1))
              END IF
*. Transpose ( To compensate later transposition )
              CALL TRPMT3(SB,NIA,NIB,C2)
              CALL COPVEC(C2,SB,NIA*NIB)
              CALL TRPMT3(CB,NJA,NJB,C2)
              CALL COPVEC(C2,CB,NJA*NJB)
            END IF
*           ^ End if transposing should be used
          ELSE IF(ISIMSYM.EQ.1.AND.IEND_CTT.EQ.1.AND.INEW_STT.EQ.1) THEN
*
* Route where all symmetryblocks belonging to given TT blocks of C 
* and S are treated simultaneously
*
* If time-reversal symmetry is used and the total symmetry of 
* C or S vector is not total symmetric, and two types are identical, 
* I have not completed the coding so ...
          
            IRES_AB = -0803
            IF(IRES_AB.EQ.1.AND.(IASM.NE.IBSM.OR.JASM.NE.JBSM)) THEN
             WRITE(6,*) 
     &       ' ISIMSYM = 1 may not be working for nonsym cases'
             WRITE(6,*) ' So, I will stop '
             STOP       
     &       ' ISIMSYM = 1 may not be working for nonsym cases'
            END IF 
*
            ICSM = MULTD2H(JASM,JBSM)
            ISSM = MULTD2H(IASM,IBSM)
*
            IA_ABS = IATP + IOCTPA -1
            IB_ABS = IBTP + IOCTPB -1
            JA_ABS = JATP + IOCTPA -1
            JB_ABS = JBTP + IOCTPB -1
*
*. transpose is temp eliminated 
            JJJTRNS = 0
            IF(JJJTRNS.EQ.0) THEN
*. No need for transposing so 
              CALL RSBB2BN_MS(ICSM,ISSM,IATP,IBTP,
     &             NSTFSMSPGP(1,IA_ABS),NSTFSMSPGP(1,IB_ABS),
     &             JATP,JBTP,
     &             NSTFSMSPGP(1,JA_ABS),NSTFSMSPGP(1,JB_ABS),
     &             IJAGRP,IJBGRP,IOCPTPA,IOCTPB,NGAS,IAOC,IBOC,JAOC,JBOC,
     &             SB,CBTT,ADSXA,STSTSX,MXPNGAS,
     &             NOBPTS,MAXK,
     &             SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &             XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &             IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &             XINT2,MXADKBLK)
            ELSE IF (JJJTRNS.EQ.1) THEN
*. Start by transposing all symmetry-blocks in C and S blocks 
C  TRP_TT_BLK(C,NA,NB,ISM,NSMST,IFB,SCR)
              CALL TRP_TT_BLK(CBTT,NSTFSMSPGP(1,JA_ABS),
     &             NSTFSMSPGP(1,JB_ABS),ICSM,NSMST,2,C2)
              CALL TRP_TT_BLK(SB,NSTFSMSPGP(1,IA_ABS),
     &             NSTFSMSPGP(1,IB_ABS),ISSM,NSMST,2,C2)
              CALL RSBB2BN_MS(ICSM,ISSM,IBTP,IATP,
     &             NSTFSMSPGP(1,IB_ABS),NSTFSMSPGP(1,IA_ABS),
     &             JBTP,JATP,
     &             NSTFSMSPGP(1,JB_ABS),NSTFSMSPGP(1,JA_ABS),
     &             IJBGRP,IJAGRP,IOCTPB,IOCTPA,NGAS,IBOC,IAOC,JBOC,JAOC,
     &             SB,CBTT,ADSXA,STSTSX,MXPNGAS,
     &             NOBPTS,MAXK,
     &             SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &             XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &             IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &             XINT2,MXADKBLK)
              CALL TRP_TT_BLK(CBTT,NSTFSMSPGP(1,JA_ABS),
     &             NSTFSMSPGP(1,JB_ABS),ICSM,NSMST,1,C2)
              CALL TRP_TT_BLK(SB,NSTFSMSPGP(1,IA_ABS),
     &             NSTFSMSPGP(1,IB_ABS),ISSM,NSMST,1,C2)
            END IF
*           ^ End of transpose switch 
          END IF
*         ^ End of switch for call to RSBB2BN_MS 
          IF(NTEST.GE.200) THEN
            WRITE(6,*) ' SB after RSBB2b'
            call wrtmat(sb,nib,nia,nib,nia)
          END IF
        END IF
*       ^ End if alpha-beta loop should be called
*
* =============================
* Sigma alpha alpha contribution
* =============================
*
      CALL MEMCHK2('SIGAI ')
        IF(I_RES_AB_.NE.-1.AND.
     &     NAEL.GE.0.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* alpha single excitation
*
*. set common oper (just in case we have unrestricted orbitals)
          ISPCAS=1 ! alpha alpha
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB1E(alpha)'
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' SB transposed before RSBB1, first element '
            call wrtmat(sb,1,1    ,nia,nib)
          END IF
          CALL MEMCHK2('ASI1BE')
          CALL RSBB1E(IASM,IATP,IOCTPA,JASM,JATP,IOCTPA,
     &                     IJAGRP,NIB,
     &                     NGAS,IAOC,JAOC,
     &                     SB,CB,
     &                     ADSXA,SXSTST,STSTSX,
     &                     MXPNGAS,NOBPTS,IOBPTS,
     &                     ITSOB,MAXI,MAXK,
     &                     SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                     NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
     &                     NIA,MXSXST,CJRES(1),CJRES(1+MXSXBL),
     &                     MOCAA,SCLFAC,IUSE_PH,IPHGAS,NTEST)
          CALL MEMCHK2('ASI1AF')
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' SB transposed after RSBB1, first element '
            call wrtmat(sb,1,1    ,nia,nib)
          END IF
          IF(NTEST.GE.200) THEN
            WRITE(6,*) ' SB transposed  after RSBB1E(alpha)'
            call wrtmat(SB,nib,nia,nib,nia)
          END IF
*
* alpha double excitation
*
          IF(IDOH2.NE.0.AND.NAEL.GE.0) THEN
            IF(NTEST.GE.101)
     &      WRITE(6,*) ' I am going to call RSBB2A (alpha-alpha )'
      CALL MEMCHK2('SIGAAI')
            CALL RSBB2A(IASM,IATP,JASM,JATP,IJAGRP,NIB,NIA,
     &           NGAS,IAOC,JAOC,  
     &           SB,CB,
     &           ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &           NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &           NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &           CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC,NTEST,0,0,
     &           IUSE_PH,IPHGAS,XINT2)
          CALL MEMCHK2('SIGAAA')
          END IF
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' SB transposed after RSBB2A, first element '
            call wrtmat(sb,1,1    ,nia,nib)
          END IF
          IF(NTEST.GE.200) THEN
            WRITE(6,*) ' SB after RSBB2A(alpha-alpha)'
            call wrtmat(sb,nia,nib,nia,nib)
          END IF
        END IF ! I_RES_AB_.NE.-1 etc.

*. reset common oper (to track uninitialized calls)
        ISPCAS=0
*
        ELSE IF (IDIAG.EQ.1.OR.IH_OCC_CONS.EQ.1) THEN
*
*. Diagonal approxiation (IDIAG = 1) 
*  or complete orbital conserving part of Ham (IH_OCC_CONS = 1)
*
          IBLOCK(1) = IATP
          IBLOCK(2) = IBTP
          IBLOCK(3) = IASM
          IBLOCK(4) = IBSM
          IBLOCK(5) = 1
          IBLOCK(6) = 1
          IF(IDOH2.EQ.0) THEN
            I12 = 1
          ELSE
            I12 = 2
          END IF
C?        WRITE(6,*) ' IDOH2, I12 ', IDOH2,I12
          ITASKL = 2
          FACTOR = 0.0D0
*. Well, we are not using transposed matrices here so
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
          IF(IATP.EQ.JATP.AND.IBTP.EQ.JBTP.AND.
     &       IASM.EQ.JASM.AND.IBSM.EQ.JBSM) THEN
C?          WRITE(6,*) ' DIATERM2_GAS will be called '
            CALL COPVEC(CB,C2,NJA*NJB)
*. Input is in det basis 
            IIDC = 1
            CALL DIATERM2_GAS(FACTOR,ITASKL,C2,1,IBLOCK,1,0,I12,IIDC)
          ELSE
            ZERO = 0.0D0
            CALL SETVEC(C2,ZERO,NIA*NIB)
          END IF
*. Remaining occupation conserving operator
          IF(IH_OCC_CONS.EQ.1) THEN
            CALL HCONFDIA_BBM(NAEL,NBEL,IJAGRP,IJBGRP,
     &           IASM,IATP,IAOC,NIA,IBSM,IBTP,IBOC,NIB,
     &           JASM,JATP,JAOC,NJA,JBSM,JBTP,JBOC,NJB,XINT,CB,C2)
          END IF
          IF(IUSEAB.EQ.0) THEN
            FACTOR = 1.0D0*SCLFAC
          ELSE
            FACTOR = 0.5D0*SCLFAC
          END IF
C             MAT_P_MATT(A,B,NR,NC,COEF)
          CALL MAT_P_MATT(SB,C2,NIB,NIA,FACTOR)
          CALL TRPMT3(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
*
c       call luciaqenter('gdgop')
c       call ga_dgop(1, SB, nia*nib, '+')
c       call luciaqexit('gdgop')
 9999 CONTINUE
        IF(IHAPR.NE.0) THEN
*. Clean up
          CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                 IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &                 IMZERO,IDIAG)
        END IF
      ELSE IF(ITASK(1:6).EQ.'TRACID') THEN
C?      WRITE(6,*) ' Routine for blockdiagonal TRACI will be called '
        CALL GET_TPAM(XINT)
*. (Remember that blocks are transposed at entrance!)
*. alpha excitations
*. 
C     TRACI_BLKDIA(ICSM,ICTP,ICEL,
C    &           IGRP,NROW,NCOL,NGAS,SB,CB,NOBPTS,
C    &           I1,XI1S,T,NSMOB,NTESTG)
*. TRACI changes C, so C2 will be used
        CALL COPVEC(CB,C2,NJA*NJB)
        CALL TRACI_BLKDIA(JASM,JATP,JAOC,IJAGRP,NJB,NJA,
     &       NGAS,SB,CB,NOBPTS,I1,XI1S,
     &       XINT,NSMOB,NTEST)
C?      WRITE(6,*) ' Home from TRACI_FCIBLK for alpha'
*. and beta-excitations
        CALL TRPMT3(SB,NIB,NIA,C2)
        CALL TRACI_BLKDIA(JBSM,JBTP,JBOC,IJBGRP,NJA,NJB, 
     &       NGAS,SB,C2,NOBPTS,I1,XI1S,XINT,NSMOB,NTEST)
*. We now have the result in SB, in untransposed form. As program
*. expects transposed form, transpose
        CALL TRPMAT(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ==================================='
        WRITE(6,*) ' RSSBCB2 : Final S block (transposed)'
        WRITE(6,*) ' ==================================='
        CALL WRTMAT(SB,NIB,NIA,NIB,NIA)
      END IF
      CALL LUCIAQEXIT('RSSBE')
*
CT    IF(ITASK(1:5).EQ.'SIGMA')
CT   &STOP ' Enforced stop at end of RSSBCB2' 
*
      NTESTO = NTEST
      RETURN
      END
      SUBROUTINE SBLOCKS_BEFORE_JUNE03(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA,NSSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB, NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &           ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &           XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC, 
     &           I1,XI1S,I2,XI2S,IDOH2,MXPOBS,ISTRFL,
     &           PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &           I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &           ICONSPA,ICONSPB,SCLFAC,
     &           LUCBLK,IPERTOP,IH0INSPC,IH0SPC,
     &           ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &           I_RES_AB,IUSE_PA,CJPA,SIPA,ISIMSYM,XINT2)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen , Winter of 1991
*               Last modification : April 99
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block      
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type 
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If IRESTRICT.NE. 0 THEN we are after :
* sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)


      IMPLICIT REAL*8(A-H,O-Z)
*. Specific input
      INTEGER ISBLOCK(8,*)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),XINT2(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),XI1S(*),XI2S(*),XI3S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*. Zero order Hamiltonian
      INTEGER IH0SPC(NOCTPA,NOCTPB)
      INTEGER IH0INSPC(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION SCLFAC(*)
*
      COMMON/H_OCC_CONS/IH_OCC_CONS
* IH_OCC_CONS =1 implies that we should employ occupation conserving 
* part of Hamiltonian
      INCLUDE 'comjep.inc'
      INCLUDE 'chktyp.inc'
*.
      REAL * 8 INPROD
*
      CALL LUCIAQENTER('SBLOC')
*
C?    WRITE(6,*) ' IPERTOP in SBLOCKS = ', IPERTOP
      IF(IH_OCC_CONS.EQ.1) THEN
        WRITE(6,*) ' Occupation conserving part of Hamiltonian '
      END IF
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' SBLOCKS speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)  
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0) 
     &      WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          END DO
          WRITE(6,*) ' IDC PS IPERTOP', IDC,PS,IPERTOP
          WRITE(6,*) ' IDOH2 = ',IDOH2
          WRITE(6,*) ' I_RES_AB=',I_RES_AB
      END IF
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
* ===========================
* 1 : Arrays for accessing C 
* ============================
*. Find batches of C - strings
      CALL PART_CIV2(IDC,ICBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,
     &              LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &              I1CBLOCK,ICBLOCK,0,ISIMSYM)
*. Find the active blocks on LUC, store info in SCLFAC
      CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)
*
      JOFF = -2803
* Initialize sigma blocks
      DO JSBLOCK = 1, NSBLOCK
        IATP = ISBLOCK(1,JSBLOCK)
        IBTP = ISBLOCK(2,JSBLOCK)
        IASM = ISBLOCK(3,JSBLOCK)
        IBSM = ISBLOCK(4,JSBLOCK)
        IOFF = ISBLOCK(5,JSBLOCK)
        NASTR = NSSOA(IASM,IATP)
        NBSTR = NSSOB(IBSM,IBTP)
        ZERO = 0.0D0
        IF(ISBLOCK(1,JSBLOCK).GT.0) 
     &  CALL SETVEC(SB(IOFF),ZERO,NASTR*NBSTR)
      END DO
* Loop over batches over C blocks      
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE 
        MXEXC = 1
      END IF
      REWIND LUC
      IF(ICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
*
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
*
*. Read C blocks into core
*
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
        DO JJCBLOCK = 1, NJBLOCK            
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
*. Will this block be needed ??
          INTERACT = 0
          IF(SCLFAC(JBLOCK).EQ. 1.0D0) THEN
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          PL = 1.0D0
          CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
*.Loop over Sigma blocks in batch
            DO JSBLOCK = 1, NSBLOCK
            IDENT = 0
            IF(ISBLOCK(1,JSBLOCK).GT.0) THEN
              IATP = ISBLOCK(1,JSBLOCK)
              IBTP = ISBLOCK(2,JSBLOCK)
              IASM = ISBLOCK(3,JSBLOCK)
              IBSM = ISBLOCK(4,JSBLOCK)
*. Are the two blocks connected by allowed excitation
              CALL CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,
     &             LLASM,LLBSM,IOCTPA,IOCTPB,ICONSPA,ICONSPB,
     &             NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
            END IF
            END DO
          END DO
*.        ^ End of checking whether C-block is needed
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
          IF(INTERACT.EQ.1) THEN
            CALL GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                   NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                   PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
*. Note in GSTTBL : ICOOSC only used for CI vectors in core, 
          ELSE
*. not relevant
            CALL IFRMDS(LBL,-1,1,LUC)
            CALL SKPRCD2(LBL,-1,LUC)
            SCLFAC(JBLOCK) = 0.0D0
          END IF
*
          IF(NTEST.GE.100) THEN
            IF(INTERACT.EQ.1) THEN
              WRITE(6,*) ' TTSS for C block read in  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            ELSE
              WRITE(6,*) ' TTSS for C block skipped  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            END IF
          END IF
*
        END DO
*       ^ End of loop over Blocks
*
*. Loop over blocks of sigma and C in core and obtain  contribution from
*  given C block to given S block
*. Loop over C blocks
        DO 9000 ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK 
          JATP = ICBLOCK(1,ICBLK)
          JBTP = ICBLOCK(2,ICBLK)
          JASM = ICBLOCK(3,ICBLK)
          JBSM = ICBLOCK(4,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
          NJA = NSSOA(JASM,JATP)
          NJB = NSSOB(JBSM,JBTP)
C?        WRITE(6,*) ' ICBLK, JATP, JBTP =', ICBLK,JATP,JBTP
*
          IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
            PL = 1.0D0
            CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                  LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
            DO 8765 IPERM = NPERM,1, -1
              LLASM = LASM(IPERM)
              LLBSM = LBSM(IPERM)
              LLATP = LATP(IPERM)
              LLBTP = LBTP(IPERM)
              NLLA = NSSOA(LLASM,LLATP)
              NLLB = NSSOB(LLBSM,LLBTP)
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
              IF(IPERM.EQ.1) THEN       
                IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM) THEN
*. Diagonal blocks, Transposing corresponds to scaling
                  IF(PS.EQ.-1.0D0) THEN
                    CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                  END IF
                ELSE
*. ofdiagonal blocks, explicit transposing
                  CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)   
                  CALL COPVEC(C2,CB(ICOFF),NJA*NJB)   
                END IF
              END IF
*
              DO 10000 ISBLK = 1, NSBLOCK
              IF(ISBLOCK(1,ISBLK) .GT. 0 ) THEN
                IATP = ISBLOCK(1,ISBLK)
                IBTP = ISBLOCK(2,ISBLK)
                IASM = ISBLOCK(3,ISBLK)
                IBSM = ISBLOCK(4,ISBLK)
                ISOFF = ISBLOCK(5,ISBLK)
                NIA = NSSOA(IASM,IATP)
                NIB = NSSOB(IBSM,IBTP)
*
                IF(NIA*NIB.EQ.0) GOTO 10000
                IF(IRESTRICT.EQ.1.AND.
     &             (JASM.GT.IASM.OR.
     &             JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &             JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &          GOTO 10000
*. Are the two blocks connected by allowed excitation
                CALL CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,
     &               LLASM,LLBSM,IOCTPA,IOCTPB,ICONSPA,ICONSPB,
     &               NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)
C
                IF(INTERACT.EQ.0) GOTO 10000

                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' Next s block in batch : '
                  write(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                END IF
*
                IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &            ((LLBSM.GT.LLASM).OR.
     &            (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &            GOTO 8764
*
                IF(NTEST.GE.60) THEN
                  WRITE(6,*) ' RSSBCB will be called for '
                  WRITE(6,*) ' Sigma block : '
                  WRITE(6,*) ' ISOFF ', ISOFF
                  WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(6,*) ' C     block : '
                  WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                  WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                  WRITE(6,*) ' ICOFF ', ICOFF
                  WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
                END IF
*
                IF(IRESTRICT.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     XFAC = 0.5D0*SCLFAC(ICBLK)
                ELSE
                     XFAC = SCLFAC(ICBLK)
                END IF
*. Form of operator in action
C               IF(IPERTOP.NE.0) THEN
*. Not exact Hamiltonian in use 
C?                 WRITE(6,*) ' IATP, IBTP, JATP, JBTP ', 
C?   &                          IATP, IBTP, JATP, JBTP
                   IPTSPC = IH0SPC(IATP,IBTP)
                   JPTSPC = IH0SPC(JATP,JBTP)
                   IJOP   = IH0INSPC(IPTSPC)
C               ELSE
C                  IPTSPC = -1
C                  JPTSPC = -1
C                  IJOP   = -1
C               END IF
                IF(ICHKTYP.EQ.1.AND.IPTSPC.NE.JPTSPC) GOTO 8764
*
C?              WRITE(6,*) ' Memcheck before RSSBCB2 '
C?              CALL MEMCHK
C?              WRITE(6,*) ' Memory checked passed '
* 
                IDUM = 0
                CALL RSSBCB2(IASM,IATP,IOCTPA,
     &               IBSM,IBTP,IOCTPB,
     &               LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &               NELFSPGP(1,IATP+IOCTPA-1),
     &               NELFSPGP(1,IBTP+IOCTPB-1),
     &               NELFSPGP(1,LLATP+IOCTPA-1),
     &               NELFSPGP(1,LLBTP+IOCTPB-1),
     &               NAEL,NBEL,
     &               IAGRP,IBGRP,
     &               SB(ISOFF),CB(ICOFF),IDOH2,
     &               ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &               NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &               NSMOB,NSMST,NSMSX,NSMDX,
     &               NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &               ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &               MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
     &               IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
     &               IUSE_PH,IPHGAS,I_RES_AB,IUSE_PA,CJPA,SIPA,
     &               IH_OCC_CONS,ISIMSYM,XINT2)
*
C?              WRITE(6,*) ' Memcheck after RSSBCB2 '
C?              CALL MEMCHK
C?              WRITE(6,*) ' Memory checked passed '
C?                WRITE(6,*) ' JATP, JBTP after RSSBCB2 = ',
C?   &            JATP, JBTP
C?                WRITE(6,*) ' LLATP, LLBTP after RSSBCB2 = ',
C?   &            LLATP, LLBTP
 8764           CONTINUE
                END IF
*               ^ End if S-block should be calculated
10000         CONTINUE
*.            ^  End of loop over sigma blocks 
 8765       CONTINUE
          END IF
*         ^ End of C-block is nonvanishing
 9000   CONTINUE
*.      ^ End of loop over C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
*
* Order 
      DO  ISBLK = 1 , NSBLOCK
      IF(ISBLOCK(1,ISBLK).GT.0) THEN
        IATP = ISBLOCK(1,ISBLK)
        IBTP = ISBLOCK(2,ISBLK)
        IASM = ISBLOCK(3,ISBLK)
        IBSM = ISBLOCK(4,ISBLK)
        ISOFF  = ISBLOCK(5,ISBLK)
        ISOFFP = ISBLOCK(6,ISBLK)
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
        IF(ICJKAIB.NE.0) THEN
*. Tranpose sigma block was obtained, transpose to obtain correct block
           CALL TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),
     &                 NSSOA(IASM,IATP),C2)
           CALL COPVEC(C2,SB(ISOFF),
     &                 NSSOA(IASM,IATP)* NSSOB(IBSM,IBTP))
        END IF
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) 
     &  CALL TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))
*
      END IF
      END DO
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' output blocks from SBLOCKS '
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,1)
      END IF
*
      CALL LUCIAQEXIT('SBLOC')
      RETURN
      END
      SUBROUTINE SIGDEN3_CI(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
     &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
     &                  ECORE,ISIGDEN) 
*
* Generate a set of sigma blocks or contributions to density,  
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB (if sigma-generation)
*
* The blocks are scaled and reformed to combination order 
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
*
* If ICBAT_RES .eq.1 then it as assumed that only
* Cbatches ICBAT_INI to ICBAT_END are stored on  LUC
*
*
* Jeppe Olsen, April 2011, from SBLOCK
*
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      INCLUDE 'cands.inc'
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*. Input (if ICISTR = 1) 
CNW   DIMENSION CV(*)
      integer CV
*
*./ORBINP/ : NACOB used
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
*
      INCLUDE 'csmprd.inc'
c      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
c      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
c     &              ADSXA(MXPOBS,2*MXPOBS),
c     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIDE3 ')
*
C?    WRITE(6,*) ' IPERTOP in SBLOCK = ', IPERTOP
      NTEST = 0
      NTEST = MAX(NTEST,IPRCIX)
      IF(NTEST.GE.5)
     &WRITE(6,*) ' SBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
C?    WRITE(6,*) ' I12 in SBLOCK = ', I12
      IF(LUCBLK.GT.0) THEN
        CALL REWINO(LUCBLK)
      END IF
*
* Info for this internal space
*. type of alpha and beta strings - as H does not change 
*. the number of electrons, I do not distinguish between spaces for C and S
      IF(ICSPC.LE.NCMBSPC) THEN
       IATP = 1              
       IBTP = 2             
      ELSE
       IATP = IALTP_FOR_GAS(ICSPC)
       IBTP = IBETP_FOR_GAS(ICSPC)
       IF(NTEST.GE.100) WRITE(6,*) ' SBLOCK : IATP, IBTP = ', IATP, IBTP
      END IF
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. alpha and beta strings with an electron removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-1,IATPM1)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-1,IBTPM1)
*. alpha and beta strings with two electrons removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-2,IATPM2)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-2,IBTPM2)
*. Number of supergroups 
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups 
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
*. connection matrices for supergroups
*
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPB),IPRCIX)
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(dbl_mb(KSTSTS),dbl_mb(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR           
*. Largest number of strings of given symmetry and type
      MAXA = 0
      MAXA0 = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXA = MAX(MAXA,MAXA0)
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
*
      MAXB = 0
      MAXB0 = IMNMX(int_mb(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MAXB = MAX(MAXB,MAXB0)
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.3 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
C?    WRITE(6,*) ' MXTSOB = ', MXTSOB
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(ISIMSYM.EQ.0) THEN
        LSCR1 = MXSOOB
      ELSE
        LSCR1 = MXSOOB_AS
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MAX(MXTSOB ** 4, NTOOB**2)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
      CALL MEMMAN(KINSCR2,INTSCR,'ADDL  ',2,'INSCR2')
*. Arrays giving allowed type combinations 
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'CIOIO ')
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',1,'SIOIO ')
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
      CALL IAIBCM(ISSPC,int_mb(KSIOIO))
      CALL IAIBCM(ICSPC,int_mb(KCIOIO))
*. Arrays giving block type
COLD  CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',1,'SVST  ')
        CALL SIGVST(int_mb(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IF( ISSPC.GE.ICSPC) THEN
         KSCIOIO = KSIOIO
      ELSE
         KSCIOIO = KCIOIO
      END IF
      CALL MXRESCPH(int_mb(KSCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &            IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,
     &            MX_NSPII)
      IF(IPRCIX.GE.3) THEN
        WRITE(6,*) 'SIGDEN3_CI:MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM',
     &                         MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM
         WRITE(6,*) 'SIGDEN3_CI : MXADKBLK ', MXADKBLK
         WRITE(6,*) ' MX_NSPII = ', MX_NSPII
      END IF
      IF(ISIMSYM.EQ.1) THEN 
        MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
        MXADKBLK = MXADKBLK_AS
      END IF
*. For hardwired routines MXCIJAB is also used
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.3)  WRITE(6,*) ' LSCR2 = ', LSCR2
C  I assume memory was allocated for blocks, so
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      IF(ISIMSYM.EQ.1) LSCR3 = MAX(LSCR3,NSMST*MXADKBLK)
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,int_mb(KCBLTP),int_mb(KSVST))
*.Some TTS arrays 
      NOOS = NOCTPA*NOCTPB*NSMCI 
      NTTS = MXNTTS
*
*. for partitioning of vector 
      CALL MEMMAN(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      CALL MEMMAN(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      CALL MEMMAN(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      CALL MEMMAN(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      CALL MEMMAN(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')

*. Space for four blocks of string occupations and arrays of 
*. reordering arrays
*. Also used to hold an NORB*NORB matrix  
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
*. Set up to two blocks for orbital conserving operator 
C     DO I1234 = 1, 2
      DO I1234 = 1, 1
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')  
      END DO
   
      DO I1234 = 1, 2
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
* 4 arrays containing all strings of given sym. Dimension can  be 
*   reduced to largest number of strings in alpha or beta.
C?    WRITE(6,*) ' SBLOCKS : MAX_STR_SPGP = ', MAX_STR_SPGP
      IF(ISIMSYM.NE.0) THEN
        CALL MEMMAN(KLREOJA,MAX_STR_SPGP,'ADDL  ',1,'REOJA ')
        CALL MEMMAN(KLREOJB,MAX_STR_SPGP,'ADDL  ',1,'REOJB ')
        CALL MEMMAN(KLREOIA,MAX_STR_SPGP,'ADDL  ',1,'REOIA ')
        CALL MEMMAN(KLREOIB,MAX_STR_SPGP,'ADDL  ',1,'REOIB ')
      END IF
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*. Place perturbation integrals over one body integrals
CINSERT_START
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
*. Prepare for perturbation calculation  
*
C     IF(IPERTOP.NE.0) THEN
*. Matrix specifying partiotioned spaces
        CALL MEMMAN(KLH0SPC,NOCTPA*NOCTPB,'ADDL  ',2,'H0SPC ')
        CALL H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                NGAS,MXPNGAS,dbl_mb(KLH0SPC),NELFGP)
C       IF(IH0SPC.EQ.0) THEN
*. Form of perturbation in subspace has not been defined,
*. Use current IPART
          IH0INSPC(1) = IPART
C       END IF
C     END IF
*
      IF(ISIGDEN.EQ.2.AND.IDC.EQ.2) THEN    
*. reform 
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,
     &             2,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,2,NTEST)
      END IF
*
C     IF(ISIMSYM.EQ.0) THEN
*. June 30 : I will try to use go through the normal route with ISIMSYM :
      KSIPA = 0
      KCJPA = 0
      CALL SIGDEN4_CI(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,dbl_mb(KC2),
     &     int_mb(KCIOIO),ISMOST(1,ICSM),int_mb(KCBLTP),
     &     NACOB,int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &     NAEL,IATP,NBEL,IBTP,
     &     IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &     ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
     &     LSCR1,dbl_mb(KINSCR),dbl_mb(KCSCR),dbl_mb(KSSCR),
     &     SXSTSM,dbl_mb(KSTSTS),dbl_mb(KSTSTD),SXDXSX,
     &     ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &     int_mb(KI1),dbl_mb(KXI1S),
     &     int_mb(KI2),dbl_mb(KXI2S),IDOH2,MXPOBS,int_mb(KSVST),
     &     PSSIGN,IPRDIA,LUC,ICJKAIB,dbl_mb(KCJRES),
     &     dbl_mb(KSIRES),int_mb(KI3),dbl_mb(KXI3S),
     &     int_mb(KI4),dbl_mb(KXI4S),MXSXST,XSXBL,
     &     MOCAA,MOCAB,IAPR,
     &     int_mb(KLLBT),int_mb(KLLEBT),int_mb(KLI1BT),int_mb(KLIBT),
     &     IRESTRICT,int_mb(KCONSPA),int_mb(KCONSPB),dbl_mb(KLSCLFAC),
     &     LUCBLK,IPERTOP,IH0INSPC,dbl_mb(KLH0SPC),
     &     ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,I_RES_AB,
     &     IUSE_PA,int_mb(KCJPA),int_mb(KSIPA),ISIMSYM,dbl_mb(KINSCR2),
     &     MXADKBLK,ICISTR,CV,ECORE)
C     ELSE
C     CALL SBLOCKSN(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
C    &     WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
C    &     NACOB,WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &     NAEL,IATP,NBEL,IBTP,
C    &     IOCTPA,IOCTPB,NOCTPA,NOCTPB,
C    &     NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
C    &     ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
C    &     LSCR1,WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
C    &     SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
C    &     ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
C    &     WORK(KI1),WORK(KXI1S),
C    &     WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,WORK(KSVST),
C    &     PSSIGN,IPRDIA,LUC,ICJKAIB,WORK(KCJRES),
C    &     WORK(KSIRES),WORK(KI3),WORK(KXI3S),
C    &     WORK(KI4),WORK(KXI4S),MXSXST,MXSXBL,
C    &     MOCAA,MOCAB,IAPR,
C    &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),int_mb(KLIBT),
C    &     IRESTRICT,WORK(KCONSPA),WORK(KCONSPB),WORK(KLSCLFAC),
C    &     LUCBLK,IPERTOP,IH0INSPC,WORK(KLH0SPC),
C    &     ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,I_RES_AB,
C    &     IUSE_PA,WORK(KCJPA),WORK(KSIPA),ISIMSYM,
C    &     WORK(KLREOIA),WORK(KLREOIB),WORK(KLREOJA),WORK(KLREOJB),
C    &     IUSE_HW,WORK(KINSCR2))
C     END IF
*
      IF(NTEST.GE.100) WRITE(6,*) ' In SIGDEN3 after SIGDEN4 '
      IF(ISIGDEN.EQ.1.AND.IDC.EQ.2) THEN    
*. reform 
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NOCTPA,NOCTPB,
     &              int_mb(KNSTSO(IATP)), int_mb(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL ITODS(-1,1,-1,LUCBLK)
      END IF
*. Eliminate local memory
      IDUM = 0
      IF(NTEST.GE.100) WRITE(6,*) ' Leaving SIGDEN3_CI '
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SIDE3 ')
      RETURN
      END
      SUBROUTINE SIGDEN4_CI(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA,NSSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB, NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS_X,
     &           ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &           XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC, 
     &           I1,XI1S,I2,XI2S,IDOH2,MXPOBS_X,ISTRFL,
     &           PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &           I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &           ICONSPA,ICONSPB,SCLFAC,
     &           LUCBLK,IPERTOP,IH0INSPC,IH0SPC,
     &           ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &           I_RES_AB,IUSE_PA,CJPA,SIPA,ISIMSYM,XINT2,MXADKBLK,
     &           ICISTR,CV,ECORE,ISIGDEN)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
*
* Jeppe Olsen, April 2011, from SBLOCKS
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block      
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type 
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If IRESTRICT.NE. 0 THEN we are after :
* sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
*. For location of density matrices
      INCLUDE 'glbbas.inc'
*
*. Specific input
      INTEGER ISBLOCK(8,*)
* (if ICISTR = 1) 
CNW   DIMENSION CV(*)
      integer CV
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),XINT2(*),CSCR(*),SSCR(*)
*
      INTEGER   I1(MXADKBLK,*),I2(MXADKBLK,*)
      INTEGER   I3(MXADKBLK,*),I4(MXADKBLK,*) 
      DIMENSION XI1S(MXADKBLK,*),XI2S(MXADKBLK,*)
      DIMENSION XI3S(MXADKBLK,*),XI4S(MXADKBLK,*)
*
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*. Zero order Hamiltonian
      INTEGER IH0SPC(NOCTPA,NOCTPB)
      INTEGER IH0INSPC(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION SCLFAC(*)
*
      COMMON/H_OCC_CONS/IH_OCC_CONS
* IH_OCC_CONS =1 implies that we should employ occupation conserving 
* part of Hamiltonian
      INCLUDE 'comjep.inc'
      INCLUDE 'chktyp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*.
      REAL * 8 INPROD
*
      CALL LUCIAQENTER('SIDE4')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIDE4 ')
*
C?    WRITE(6,*) ' IPERTOP in SBLOCKS = ', IPERTOP
      IF(IH_OCC_CONS.EQ.1) THEN
        WRITE(6,*) ' Occupation conserving part of Hamiltonian '
      END IF
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' SIGDEN4_CI speaking :'
        WRITE(6,*) ' ======================'
        WRITE(6,*)  
        WRITE(6,*) ' ISIGDEN = ', ISIGDEN
        WRITE(6,*) ' Number of sigma blocks active ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0) 
     &      WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          END DO
          WRITE(6,*) ' IDC PS IPERTOP', IDC,PS,IPERTOP
          WRITE(6,*) ' IDOH2 = ',IDOH2
          WRITE(6,*) ' I_RES_AB=',I_RES_AB
          WRITE(6,*) ' ECORE = ', ECORE
      END IF
*
      IF(ISIGDEN.EQ.2) THEN
        KRHO1_A = KRHO1
        KRHO2_A = KRHO2
      ELSE
        KRHO1_A = 0
        KRHO2_A = 0
      END IF
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' Initial C vector '
        IF(ICISTR.NE.1) THEN 
          CALL WRTVCD(CB,LUC,1,-1)
        ELSE
*. Print using ICLBLOCK
        END IF
      END IF
*
* ===========================
* 1 : Arrays for accessing C 
* ============================
*. Find batches of C - strings
      IF(ICISTR.EQ.1) THEN
        ICOMP = 1
      ELSE
        ICOMP = 0
      END IF
*
      CALL PART_CIV2(IDC,ICBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,
     &              LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &              I1CBLOCK,ICBLOCK,ICOMP,ISIMSYM)
*. Find the active blocks on LUC, store info in SCLFAC
      IF(ICISTR.NE.1) THEN
*. Active blocks on LUC
        CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)
      ELSE
*. Active blocks in vectors CV = number of blocks in BATCH 1
        NCBLOCK_T = LCBLOCK(1)
        CALL FIND_ACTIVE_BLOCKS_IN_VECTOR(CV,SCLFAC,NCBLOCK_T,ICBLOCK)
      END IF
*. If ISIMSYM is active, enforce the same activity of all blocks 
*. belonging to a TT block
*. Number of BLOCKS in  C
        NBLOCK_C =  I1CBLOCK(NCBATCH)-1+ LCBLOCK(NCBATCH)
C       WRITE(6,*) ' Number of blocks in C', NBLOCK_C
*
       CALL ENFORCE_SIMSYM_ACT_BLK(NBLOCK_C,ICBLOCK,SCLFAC)
C           ENFORCE_SIMSYM_ACT_BLK(NBLK,IBLK,XSCA)
*
* Initialize sigma blocks
      DO JSBLOCK = 1, NSBLOCK
        IATP = ISBLOCK(1,JSBLOCK)
        IBTP = ISBLOCK(2,JSBLOCK)
        IASM = ISBLOCK(3,JSBLOCK)
        IBSM = ISBLOCK(4,JSBLOCK)
        IOFF = ISBLOCK(5,JSBLOCK)
        NASTR = NSSOA(IASM,IATP)
        NBSTR = NSSOB(IBSM,IBTP)
        ZERO = 0.0D0
        IF(ISBLOCK(1,JSBLOCK).GT.0.AND.ISIGDEN.EQ.1) 
     &  CALL SETVEC(SB(IOFF),ZERO,NASTR*NBSTR)
      END DO
* Loop over batches over C blocks      
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE 
        MXEXC = 1
      END IF
      IF(ICISTR.EQ.1) THEN
        IOFF_C = 1
      ELSE 
        REWIND LUC
      END IF
      IF(ICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
*
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
*
*. Read C blocks into core
*
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
        DO JJCBLOCK = 1, NJBLOCK            
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
*. Will this block be needed ??
          INTERACT = 0
          INTERACTL = 0
          IF(SCLFAC(JBLOCK).EQ. 1.0D0) THEN
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          PL = 1.0D0
          CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
*.Loop over Sigma blocks in batch
            DO JSBLOCK = 1, NSBLOCK
            IDENT = 0
            IF(ISBLOCK(1,JSBLOCK).GT.0) THEN
              IATP = ISBLOCK(1,JSBLOCK)
              IBTP = ISBLOCK(2,JSBLOCK)
              IASM = ISBLOCK(3,JSBLOCK)
              IBSM = ISBLOCK(4,JSBLOCK)
*. Are the two blocks connected by allowed excitation
              CALL CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,
     &             LLASM,LLBSM,IOCTPA,IOCTPB,ICONSPA,ICONSPB,
     &             NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)
              IF(INTERACT.EQ.1) INTERACTL = 1
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
            END IF
            END DO
          END DO
*.        ^ End of checking whether C-block is needed
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
          IF(INTERACTL.EQ.1) THEN
            CALL GSTTBL(C,CB(ICOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                   NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                   PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK),CV)
*. Note in GSTTBL : ICOOSC only used for CI vectors in core, 
          ELSE
*. not relevant
            CALL IFRMDS(LBL,-1,1,LUC)
            CALL SKPRCD2(LBL,-1,LUC)
            SCLFAC(JBLOCK) = 0.0D0
          END IF
*
          IF(NTEST.GE.100) THEN
            IF(INTERACT.EQ.1) THEN
              WRITE(6,*) ' TTSS for C block read in  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            ELSE
              WRITE(6,*) ' TTSS for C block skipped  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            END IF
          END IF
*
        END DO
*       ^ End of loop over Blocks
*
*. Loop over blocks of sigma and C in core and obtain  contribution from
*  given C block to given S block
*. Loop over C blocks
        ICOFF_TT = -2303
        DO 9000 ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK 
*. Is this a new TT block -or is it identical to the previous 
*. TT block
          INEW_CTT = 0
          IF(ICBLK.EQ. I1CBLOCK(JCBATCH)) THEN
             INEW_CTT = 1
          ELSE 
             JATP_PR = JATP
             JBTP_PR = JBTP 
             IF(ICBLOCK(1,ICBLK).NE.JATP_PR.OR.
     &          ICBLOCK(2,ICBLK).NE.JBTP_PR) INEW_CTT = 1
          END IF
*. last block
          IEND_CTT = 0
          IF(ICBLK.EQ.I1CBLOCK(JCBATCH)-1+NJBLOCK) THEN
            IEND_CTT = 1
          ELSE 
            IF(ICBLOCK(1,ICBLK).NE.ICBLOCK(1,ICBLK+1).OR.
     &         ICBLOCK(2,ICBLK).NE.ICBLOCK(2,ICBLK+1)   ) IEND_CTT = 1
          END IF
C?        WRITE(6,*) ' Loop 9000, ICBLK, INEW_CTT, SCLFAC = ', 
C?   &                            ICBLK, INEW_CTT,SCLFAC(ICBLK)
*
          JATP = ICBLOCK(1,ICBLK)
          JBTP = ICBLOCK(2,ICBLK)
          JASM = ICBLOCK(3,ICBLK)
          JBSM = ICBLOCK(4,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
          IF(INEW_CTT.EQ.1) ICOFF_TT = ICOFF
          NJA = NSSOA(JASM,JATP)
          NJB = NSSOB(JBSM,JBTP)
C?        WRITE(6,*) ' ICBLK, JATP, JBTP =', ICBLK,JATP,JBTP
*
          IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
            PL = 1.0D0
            CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                  LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
            DO 8765 IPERM = NPERM,1, -1
              LLASM = LASM(IPERM)
              LLBSM = LBSM(IPERM)
              LLATP = LATP(IPERM)
              LLBTP = LBTP(IPERM)
              NLLA = NSSOA(LLASM,LLATP)
              NLLB = NSSOB(LLBSM,LLBTP)
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
              IF(IPERM.EQ.1) THEN       
                IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM) THEN
*. Diagonal blocks, Transposing corresponds to scaling
                  IF(PS.EQ.-1.0D0) THEN
                    CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                  END IF
                ELSE
*. ofdiagonal blocks, explicit transposing
                  CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)   
                  CALL COPVEC(C2,CB(ICOFF),NJA*NJB)   
                END IF
              END IF
*
              DO 10000 ISBLK = 1, NSBLOCK
*. Start of a new TT block of Sigma ?
              INEW_STT = 0
              IF(ISBLK.EQ. 1) THEN
                 INEW_STT = 1
              ELSE 
                 IATP_PR = IATP
                 IBTP_PR = IBTP 
                 IF(ISBLOCK(1,ISBLK).NE.IATP_PR.OR.
     &              ISBLOCK(2,ISBLK).NE.IBTP_PR) INEW_STT = 1
              END IF

              IF(ISBLOCK(1,ISBLK) .GT. 0 ) THEN
                IATP = ISBLOCK(1,ISBLK)
                IBTP = ISBLOCK(2,ISBLK)
                IASM = ISBLOCK(3,ISBLK)
                IBSM = ISBLOCK(4,ISBLK)
                ISOFF = ISBLOCK(5,ISBLK)
                NIA = NSSOA(IASM,IATP)
                NIB = NSSOB(IBSM,IBTP)
*
                IF(ISIMSYM.EQ.0.AND.NIA*NIB.EQ.0) GOTO 10000
                IF(IRESTRICT.EQ.1.AND.
     &             (JASM.GT.IASM.OR.
     &             JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &             JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &          GOTO 10000
*. Are the two blocks connected by allowed excitation
                CALL CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,
     &               LLASM,LLBSM,IOCTPA,IOCTPB,ICONSPA,ICONSPB,
     &               NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)
C
                IF(INTERACT.EQ.0) GOTO 10000

                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' Next s block in batch : '
                  write(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                END IF
*
                IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &            ((LLBSM.GT.LLASM).OR.
     &            (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &            GOTO 8764
*
                IF(NTEST.GE.60) THEN
                  WRITE(6,*) ' RSSBCB will be called for '
                  WRITE(6,*) ' Sigma block : '
                  WRITE(6,*) ' ISOFF ', ISOFF
                  WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                  WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(6,*) ' C     block : '
                  WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                  WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                  WRITE(6,*) ' ICOFF ', ICOFF
                  WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
                END IF
*
                IF(IRESTRICT.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     XFAC = 0.5D0*SCLFAC(ICBLK)
                ELSE
                     XFAC = SCLFAC(ICBLK)
                END IF
*. Form of operator in action
C               IF(IPERTOP.NE.0) THEN
*. Not exact Hamiltonian in use 
C?                 WRITE(6,*) ' IATP, IBTP, JATP, JBTP ', 
C?   &                          IATP, IBTP, JATP, JBTP
                   IPTSPC = IH0SPC(IATP,IBTP)
                   JPTSPC = IH0SPC(JATP,JBTP)
                   IJOP   = IH0INSPC(IPTSPC)
C               ELSE
C                  IPTSPC = -1
C                  JPTSPC = -1
C                  IJOP   = -1
C               END IF
                IF(ICHKTYP.EQ.1.AND.IPTSPC.NE.JPTSPC) GOTO 8764
*
                IDUM = 0
*
*. Add core-energy to diagonal terms
*
                IF(ISIGDEN.EQ.1.AND.
     &             IASM.EQ.LLASM.AND.IATP.EQ.LLATP.AND.
     &             IBSM.EQ.LLBSM.AND.IBTP.EQ.LLBTP) THEN
                   IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
                     FACTOR = ECORE/2.0D0/SQRT(2.0D0)
                   ELSE
                     FACTOR = ECORE
                   END IF
                   CALL VECSUM(SB(ISOFF),SB(ISOFF),CB(ICOFF),
     &             1.0D0,FACTOR,NIA*NIB)
C                  WRITE(6,*) ' Core energy added = ', ECORE
                END IF
     &             
                IF(ISIGDEN.EQ.1) THEN
*
*. Sigma generation
*
                CALL RSSBCB2(IASM,IATP,IOCTPA,
     &               IBSM,IBTP,IOCTPB,
     &               LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &               NELFSPGP(1,IATP+IOCTPA-1),
     &               NELFSPGP(1,IBTP+IOCTPB-1),
     &               NELFSPGP(1,LLATP+IOCTPA-1),
     &               NELFSPGP(1,LLBTP+IOCTPB-1),
     &               NAEL,NBEL,
     &               IAGRP,IBGRP,
     &               SB(ISOFF),CB(ICOFF),IDOH2,
     &               ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &               NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &               NSMOB,NSMST,NSMSX,NSMDX,
     &               NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &               ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &               MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
     &               IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
     &               IUSE_PH,IPHGAS,I_RES_AB,IUSE_PA,CJPA,SIPA,
     &               IH_OCC_CONS,ISIMSYM,XINT2,
     &               MXADKBLK,INEW_CTT,INEW_STT,IEND_CTT,CB(ICOFF_TT))
                ELSE
*
* Density matrix construction
*
                     IF(IDOH2.EQ.1) THEN
                       I12_L = 2
                     ELSE
                       I12_L = 1
                     END IF
*. At the moment, no spindensities
                     IDOSRHO1 = 0
                     IDOSRHO2 = 0
                     NDACTOB = NACOB
                   
                     CALL GSDNBB2(I12_L,dbl_mb(KRHO1_A),dbl_mb(KRHO2_A),
     &                    IIASM,IIATP,IIBSM,IIBTP,
     &                    JASM,JATP,JBSM,JBTP,NGAS,
     &                    NELFSPGP(1,IOCTPA-1+IIATP),
     &                    NELFSPGP(1,IOCTPB-1+IIBTP),
     &                    NELFSPGP(1,IOCTPA-1+JATP),
     &                    NELFSPGP(1,IOCTPB-1+JBTP),
     &                    NAEL,NBEL,IAGRP,IBGRP,
     &                    SB(ISOFF),CB(ICOFF),C2,
     &                    ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                    MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,
     &                    SSCR,CSCR,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                    X,NSMOB,NSMST,NSMSX,NSMDX,
     &                    NIIA,NIIB,NJA,NJB,MXPOBS,
     &                    IPRNT,NACOB,RHO1S,SCLFAC,
     &                    S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &                    IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &                    NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,
     &                    IOBPTS_SEL,NINOB)
                END IF
*
 8764           CONTINUE
                END IF
*               ^ End if S-block should be activated
10000         CONTINUE
*.            ^  End of loop over sigma blocks 
 8765       CONTINUE
          END IF
*         ^ End of C-block is nonvanishing
 9000   CONTINUE
*.      ^ End of loop over C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
*
* Order 
      DO  ISBLK = 1 , NSBLOCK
      IF(ISBLOCK(1,ISBLK).GT.0) THEN
        IATP = ISBLOCK(1,ISBLK)
        IBTP = ISBLOCK(2,ISBLK)
        IASM = ISBLOCK(3,ISBLK)
        IBSM = ISBLOCK(4,ISBLK)
        ISOFF  = ISBLOCK(5,ISBLK)
        ISOFFP = ISBLOCK(6,ISBLK)
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
        IF(ISIGDEN.EQ.1.AND.ICJKAIB.NE.0) THEN
*. Tranpose sigma block was obtained, transpose to obtain correct block
           CALL TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),
     &                 NSSOA(IASM,IATP),C2)
           CALL COPVEC(C2,SB(ISOFF),
     &                 NSSOA(IASM,IATP)* NSSOB(IBSM,IBTP))
        END IF
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) 
     &  CALL TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))
*
      END IF
      END DO
*
      IF(ISIGDEN.EQ.1.AND.NTEST.GE.50) THEN
        WRITE(6,*) ' output blocks from SIGDEN4_CI '
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,1)
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' Leaving SIGDEN4_CI'
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SIDE4 ')
      CALL LUCIAQEXIT('SIDE4')
      RETURN
      END
      SUBROUTINE CHECK_BLOCKS_FOR_BK_APPROX(
     &          IATP,IBTP,JATP,JBTP,IASM,IBSM,JASM,JBSM,
     &          IOCTPA,IOCTPB,I_DO_EXACT_BLOCK)

*. Check whether block <IATP, IBTP! H! JATP, JBTP> should
* be calculated exactly or by BK approx
* 
*
*
*. Input
* ======
* IATP IBTP JATP JBTP : Supergroups, relative numbers
* IOCTPA, IOBTPB : Offset for type
*
*. Output
*
* I_DO_EXACT_BLOCK = 1 => Do exact block
*                  = 0 => Set block to zero
*                  =-1 => Use diagonal aproximation
* Giovanni +Jeppe Olsen, Sept 2011, on a bench at Torre Normanna, Sicily
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
*. Local 
      INTEGER IOCC(MXPNGAS), JOCC(MXPNGAS)
*
      INTERACT = 0
*. Obtain Occupations of IASTR*IBSTR, JASTR*JBSTR
      IONE = 1
      CALL IVCSUM(IOCC,
     &     NELFSPGP(1,IOCTPA-1+IATP),NELFSPGP(1,IOCTPB-1+IBTP),
     &     IONE,IONE,NGAS)
      CALL IVCSUM(JOCC,
     &     NELFSPGP(1,IOCTPA-1+JATP),NELFSPGP(1,IOCTPB-1+JBTP),
     &     IONE,IONE,NGAS)
*
      IOCC_IN = ICHECK_OCC_IN_ACCSPC(IOCC,IOCCPSPC,NGAS,MXPNGAS)
      JOCC_IN = ICHECK_OCC_IN_ACCSPC(JOCC,IOCCPSPC,NGAS,MXPNGAS)
*. If both occupation are outside of IOCCPSPC, we make approximations
      IF(IOCC_IN.EQ.0.AND.JOCC_IN.EQ.0) THEN
*. If the blocks are identical use diagonal approximation
        IF(IATP.EQ.JATP.AND.IASM.EQ.JASM.AND.
     &     IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*. Compute as diagonal
           I_DO_EXACT_BLOCK = -1
        ELSE
*. Neglect
           I_DO_EXACT_BLOCK =  0
        END IF
      ELSE
*. Atleast one block is in PSPC, so calculate exactly
        I_DO_EXACT_BLOCK = 1
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*)  ' CHECK_BLOCKS_FOR_BK_APPROX speaking'
       WRITE(6,'(A, 4I4)') ' Input blocks IA, IB, JA, JB = ',
     &                     IATP,IBTP, JATP, JBTP
       WRITE(6,'(A,I4)') ' I_DO_EXACT_BLOCK = ', I_DO_EXACT_BLOCK
      END IF
*
      RETURN
      END
      FUNCTION ICHECK_OCC_IN_ACCSPC(IOCC,IMINMAX,NGAS,MXPNGAS)
*
* Is Occupation of GAS Spaces defined by IOCC in 
* accumulated space defined by IMINMAX
*
* Jeppe Olsen, for BK impmentation, take II
*
      INCLUDE 'implicit.inc' 
*. Input
      INTEGER IOCC(NGAS)
      INTEGER IMINMAX(MXPNGAS,2)
*
      IACC = 0
      I_AM_IN = 1
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          IACC = IOCC(1) 
        ELSE 
          IACC = IACC + IOCC(IGAS)
        END IF
        IF(IMINMAX(IGAS,1).GT.IACC.OR.IMINMAX(IGAS,2).LT.IACC) 
     &  I_AM_IN = 0
      END DO
* 
       ICHECK_OCC_IN_ACCSPC = I_AM_IN
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Input to ICHECK_OCC_IN_ACCSPC, IOCC, IMINMAX'
       CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
       CALL IWRTMA(IMINMAX,NGAS,2,MXPNGAS,2)
       WRITE(6,*) ' And the verdict is ', I_AM_IN
      END IF
*
      RETURN
      END
      SUBROUTINE ADDDIA_TERM(FACTOR,CVEC,SVEC,IASPGP,IBSPGP,IASM,IBSM,
     &                       JPERT)
*. Update Sigma vector with diagonal terms for a given block
*     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP) 
*                         + (FACTOR+DIAG(IASPGP,IBSPGP))CVEC(IASPGP,IBSPGP)
*
* If JPERT.NE.0, the perturbation operator as defined by IPART is used.
*
* Jeppe Olsen, September 2011
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cecore.inc'
*
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION CVEC(*)
*. Output
      DIMENSION SVEC(*)
*
      CALL LUCIAQENTER('ADDDI')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ADDDIA')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
C     IF(JPERT.EQ.0) THEN
*. Use full Hamiltonian
C       I12 = 2
C       IPERTOP = 0
C     ELSE
*. Use perturbation operator 
C       IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
C         I12 = 1
C         IPERTOP = 1
C       ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
C         I12 = 2
C         IPERTOP = 0
C       END IF
C     END IF

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   ADDDIA_TERM speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        WRITE(6,*) ' IASPGP, IBSPGP = ', IASPGP, IBSPGP
        WRITE(6,*) ' JPERT,IPART,J12,IPERTOP',JPERT,IPART,J12,IPERTOP
      END IF
*. A bit of scracth 
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')  !done
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')  !done
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')  !done
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')  !done
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')  !done
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')  !done
      CALL MEMMAN(KLH1D ,NTOOB,   'ADDL  ',2,'KLH1D ')
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')  !done
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')  !done
      IATP = 1
      NOCTPA = NOCTYP(IATP)
      MAXA = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')  !done
*. Diagonal of one-body integrals and coulomb and exchange integrals
*. Integrals assumed in place so :
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      IF(IUSE_PH.EQ.1) THEN
c       CALL SWAPVE(WORK(KH),WORK(KINT1),NINT1)
        ktmp = kh
        kh = kint1
        kint1 = ktmp
        ECORE_SAVE = ECORE
        ECORE = ECORE_ORIG
      END IF
      CALL GT1DIA(dbl_mb(KLH1D))
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,dbl_mb(KLSCR2),IREOTS)
      CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,dbl_mb(KLSCR2),IREOTS)
*

      ECOREP = ECORE
      FACTORX = 0.0D0
C?    WRITE(6,*) ' ECOREP = ',  ECOREP
*
      CALL ADDDIA_TERMS(NAEL,int_mb(KLASTR),NBEL,int_mb(KLBSTR),
     &             NACOB,CVEC,SVEC,NSMST,dbl_mb(KLH1D),
     &             IDC,dbl_mb(KLXA),dbl_mb(KLXB),dbl_mb(KLSCR),
     &             dbl_mb(KLJ),dbl_mb(KLK),
     &             int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &             ECOREP,
     &             IPRDIA,NTOOB,
     &             dbl_mb(KLRJKA),
     &             IASPGP,IASM,IBSPGP,IBSM,FACTORX,I12)
C     SUBROUTINE ADDDIA_TERMS(NAEL,IASTR,NBEL,IBSTR,
C    &                  NORB,CVEC,SVEC,NSMST,H,
C    &                  IDC,XA,XB,SCR,RJ,RK,
C    &                  NSSOA,NSSOB,
C    &                  ECORE,
C    &                  IPRNT,NTOOB,RJKAA,
C    &                  IASPGP,IASM,IBSPGP,IBSM,FACTOR,I12)    
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'ADDDIA')
      CALL LUCIAQEXIT('ADDDI')
*
      RETURN
      END 
      SUBROUTINE ADDDIA_TERMS(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,CVEC,SVEC,NSMST,H,
     &                  IDC,XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,
     &                  ECORE,
     &                  IPRNT,NTOOB,RJKAA,
     &                  IASPGP,IASM,IBSPGP,IBSM,FACTOR,I12)    
*
*. Update Sigma vector with diagonal terms for a given block
*     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP) 
*                         + FACTOR*DIAG(IASPGP,IBSPGP)CVEC(IASPGP,IBSPGP)
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, September 2011
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
      IMPLICIT REAL*8           (A-H,O-Z)
*.General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION CVEC(*)
*. 
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION SVEC(*)
*
      NTEST =  00
      NTEST = MAX(NTEST,IPRNT)
C?    WRITE(6,*) ' NTEST = ',NTEST
*
COLD  IF(LUIN.GT.0) REWIND LUIN
COLD  IF(LUOUT.GT.0) REWIND LUOUT
 
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' ADDDIA_TERMS in action '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' IASM, IASPGP, IBSM, IBSPGP = ',
     &               IASM, IASPGP, IBSM, IBSPGP
        WRITE(6,*) ' I12 = ', I12
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NORB,1,NORB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*) ' I12 and ITASK = ', I12,ITASK
        END IF
      WRITE(6,*) ' FACTOR = ',FACTOR
      END IF
*
**3 Diagonal elements according to Handys formulae
*   (corrected for error)
*
*   DIAG(IDET) = HII*(NIA+NIB)
*              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,A)*N(J,A)
*              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,B)*N(J,B)
*              +         J(I,J) * N(I,A)*N(J,B)
* N(X) are occupation numbers
*
*. K goes to J - K
      IF(I12.EQ.2) 
     &CALL VECSUM(RK,RK,RJ,-1.0D0,+1.0D0,NTOOB **2)
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IASPGP
      IDUM = 0
      CALL GETSTR_TOTSM_SPGP(1,IASPGP,IASM,NAEL,NASTR1,IASTR,
     &                         NORB,0,IDUM,IDUM)
      IF(NTEST.GE.1000) THEN
        write(6,*) ' After GETSTR for A strings '
        WRITE(6,*) ' alpha strings obtained '
        NIA = NSSOA(IASM,IASPGP)
        CALL IWRTMA(IASTR,NAEL,NIA,NAEL,NIA) 
      END IF
*
      NIA = NSSOA(IASM,IASPGP)
      DO IA = 1 ,NIA
        EAA = 0.0D0
        DO IEL = 1, NAEL
          IAEL = IASTR(IEL,IA)
          EAA = EAA + H(IAEL)
          IF(I12.EQ.2) THEN
            DO JEL = 1, NAEL
              EAA =   EAA + 0.5D0*RK(IASTR(JEL,IA),IAEL )
            END DO   
          END IF
        END DO
        RJKAA(IA) = EAA 
      END DO
*. Obtain alpha strings of sym IBSM and type IBTP
      CALL GETSTR_TOTSM_SPGP(2,IBSPGP,IBSM,NBEL,NBSTR1,IBSTR,
     &                       NORB,0,IDUM,IDUM)
C?    WRITE(6,*) ' IBSM, IBSPGP = ', IBSM, IBSPGP
      NIB =  NSSOB(IBSM,IBSPGP)
*
      IDET = 0
      DO IB = 1 ,NIB    
C?      WRITE(6,*) ' TEST, IB, BETA STRING ', IB,
C?   &  (IBSTR(JEL,IB),JEL=1, NBEL)
*
*. Terms depending only on IB
*
        HB = 0.0D0
        RJBB = 0.0D0
        CALL SETVEC(XB,0.0D0,NORB)

        DO IEL = 1, NBEL
          IBEL = IBSTR(IEL,IB)
          HB = HB + H(IBEL )
*
          IF(I12.EQ.2) THEN
            DO JEL = 1, NBEL
              RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
            END DO
*
            DO IORB = 1, NORB
              XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
            END DO
          END IF
        END DO
*. Core energy is added another place (outside, so)
        EB = HB + 0.5D0*RJBB + ECORE - ECORE
C?      WRITE(6,*) ' EB, HB, RJBB = ', EB, HB, RJBB
*
        DO IA = 1,NSSOA(IASM,IASPGP)
C?      WRITE(6,*) ' TEST, IA, Alpha STRING ', IA,
C?   &  (IASTR(JEL,IA),JEL=1, NAEL)
          IDET = IDET + 1
          X = EB + RJKAA(IA)
          DO IEL = 1, NAEL
            X = X +XB(IASTR(IEL,IA)) 
          END DO
*
C?        WRITE(6,*) ' IDET, SVEC(IDET) before add ', IDET, SVEC(IDET)
          SVEC(IDET) = SVEC(IDET) + CVEC(IDET)*(X+FACTOR)
C?       write(6,*) ' IDET,X,SVEC(IDET) ', IDET,X,SVEC(IDET)
        END DO ! IA
      END DO ! IB
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input and output vectors, ADDDIA_TERMS '
        CALL WRTMAT(CVEC,1,IDET,1,IDET)
        CALL WRTMAT(SVEC,1,IDET,1,IDET)
      END IF
*
      RETURN
      END
      SUBROUTINE IACCU(IAR,NELMNT)
*
* Change array to accumulated form:
* IAR(I) = SUM(J=1,I) IAR(J)
*
*
      INCLUDE 'implicit.inc'
*. input and output
      INTEGER IAR(NELMNT)
*
      DO J = 2, NELMNT
       IAR(J) = IAR(J) + IAR(J-1)
      END DO
*
      RETURN
      END
      SUBROUTINE TRACI_BLKDIA(ICSM,ICTP,ICEL,
     &           IGRP,NROW,NCOL,NGASX,SB,CB,NOBPTSX,
     &           I1,XI1S,T,NSMOBX,NTESTG)
*
* Peer-Aake biotransformation:
* One-electron excitations on column strings
* Special version for block diagonal transformations, 
*
*. Note that CB is overwritten in the process
*
* =====
* Input
* =====
*
* ICSM,ICTP : Symmetry and type of sigma columns
* IBCTP : Base for sigma column types 
* ICEL : Number of electrons per gas space in column strings
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NCOL : Number of columns in S and C block
* NGAS : Number of active sets 
* CB   : Input C block
* NSMOB: Number of symmetries of orbitals
* T : T-coefficients
*
* ======
* Output
* ======
* SB : updated sigma block
* =======
* Scratch
* =======
*
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
*
* Jeppe Olsen, Indian summer of 2011 for the VB project
*
      IMPLICIT REAL*8(A-H,O-Z)
*. MAX dimensions 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
*.Specific Input
      DIMENSION CB(*), T(*)
      INTEGER  ICEL(NGAS)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION I1(*),XI1S(*)
*
      REAL*8 INPROD
*
      CALL LUCIAQENTER('TRADBL')
*
      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
*
*
      IF(NTEST.GE.10)THEN
        WRITE(6,*) 
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Information from TRACI_BLKDIA '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ICEL : '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' The T- matrix '
        CALL APRBLM2(T,NTOOBS,NTOOBS,NSMOB,0)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input block of transposed C'
        CALL WRTMAT(CB,NROW,NCOL,NROW,NCOL)
      END IF
*
C?    X2C = INPROD(CB,CB,NROW*NROW)
*. Obtain groups
C     GET_SPGP_INF(ISPGP,ITP,ICGRP)
      CALL GET_SPGP_INF(ICTP,IGRP,ICGRP)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ICTP,IGRP,ICGRP = ', ICTP,IGRP,ICGRP
      END IF
*. Loop over symmetry and type of orbitals - identical for I and J
      DO IJSM = 1, NSMOB
        IF(NTEST.GE.1000) WRITE(6,*) ' IJSM = ', IJSM
*. Offset for start of block with symmetry IJSM in T
        IF(IJSM.EQ.1) THEN
          IJ_SM_OFF = 1
        ELSE 
*. Update pointer to start of matrix block using NOB_IJSM from previous IJSM
          IJ_SM_OFF = IJ_SM_OFF + NOB_IJSM**2
        END IF
*. Number of orbitals with this symmetry
        NOB_IJSM = NTOOBS(IJSM)
        IJ_GAS_OFF = 1
        DO IJTP = 1, NGAS
         IF(NTEST.GE.1000) write(6,*) ' IJSM, IJTP ', IJSM, IJTP
         IF(IJTP.GT.1) IJ_GAS_OFF = IJ_GAS_OFF + NOBPTS(IJTP-1,IJSM)
         NIJORB = NOBPTS(IJTP,IJSM)
         KFRST = 1
         IFRST = 1
         IF(NIJORB.NE.0) THEN
*
*. Obtain annihilation/creation maps for all K strings
*
          ONE = 1.0D0
          CALL ADAST_GAS(IJSM,IJTP,NGAS,ICGRP,ICSM,
     &         I1,XI1S,NKASTR,IEND,1,1,KACT,ONE,2)
C              ADAST_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
C    &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC,
C    &                    IAC)
          IF(NTEST.GE.100) WRITE(6,*) ' NKSTR = ', NKSTR
          DO JORB = 1, NIJORB
            ZERO = 0.0D0
            CALL SETVEC(SB,ZERO,NROW*NCOL)
            DO KSTR = 1, NKASTR
*. Is <Kstr!a_j .ne. 0
              JCSTR =   I1((JORB-1)*NKASTR+KSTR)
              SIGNJ = XI1S((JORB-1)*NKASTR+KSTR)
              IF(NTEST.GE.10000) WRITE(6,*) ' KSTR, JORB, JCSTR  = ',
     &        KSTR, JORB, JCSTR
              IF(JCSTR.NE.0) THEN
                DO IORB = 1, NIJORB
                  ICSTR =   I1((IORB-1)*NKASTR+KSTR)
                  SIGNI = XI1S((IORB-1)*NKASTR+KSTR)
                  IF(ICSTR.NE.0) THEN
                    SIGN = SIGNI*SIGNJ
                    ICOFF = (JCSTR-1)*NROW + 1
                    ISOFF = (ICSTR-1)*NROW + 1
C?                  IIORB = IORB + IJ_GAS_OFF - 1
C?                  JJORB = JORB + IJ_GAS_OFF - 1
*. Orbitals IIORB and JJORB are type ordered GAS orbitals. Obtain
*  the corresponding type orbitals
                    IABS = IORB+IOBPTS_GN(IJTP,IJSM)-1
                    IREO = IREOTS(IABS) - IBSO(IJSM) + 1
                    JABS = JORB+IOBPTS_GN(IJTP,IJSM)-1
                    JREO = IREOTS(JABS) - IBSO(IJSM) + 1
                    IF(NTEST.GE.10000) THEN
                      WRITE(6,*) ' IJ_SM_OFF, IJSM = ',
     &                             IJ_SM_OFF, IJSM
                      WRITE(6,*) ' IREO, JREO = ',IREO,JREO
                      WRITE(6,*) ' IJ = ',
     &                IJ_SM_OFF-1 + (JREO-1)*NOB_IJSM + IREO
                      WRITE(6,*) ' ISOFF, ICOFF = ', ISOFF, ICOFF
                    END IF
*
                    TIJ =T(IJ_SM_OFF-1 + (JREO-1)*NOB_IJSM + IREO)
*
 
                    IF(IORB.EQ.JORB) TIJ = TIJ - 1.0D0
                    FACTOR = TIJ*SIGN
                    IF(NTEST.GE.10000) WRITE(6,*) ' FACTOR = ', FACTOR
                    DO IRSTR = 1, NROW
                      SB(ISOFF - 1+ IRSTR) = SB(ISOFF - 1+IRSTR) + 
     &                FACTOR*CB(ICOFF-1+IRSTR)
                    END DO
                  END IF
                END DO! IORB
              END IF !JCSCTR.NE.0
            END DO! KSTR
            ONE = 1.0D0
            CALL VECSUM(CB,CB,SB,ONE,ONE,NROW*NCOL)
            IF(NTEST.GE.200) THEN
              WRITE(6,*) ' SB after JORB = ', JORB
              CALL WRTMAT(SB,NROW,NCOL,NROW,NCOL)
              WRITE(6,*) ' CB after JORB = ', JORB
              CALL WRTMAT(CB,NROW,NCOL,NROW,NCOL)
            END IF
          END DO ! JORB
         END IF
        END DO! IJTP
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' CB after IJSM = ', IJSM
          CALL WRTMAT(CB,NROW,NCOL,NROW,NCOL)
        END IF
*.(end of loop over symmetries)
      END DO !IJSM
*. And copy to SB- where we normally have the results
      CALL COPVEC(CB,SB,NROW*NCOL)
*
C?    X2S = INPROD(SB,SB,NROW*NROW)
C?    WRITE(6,*) ' Square norm of input and output block ',
C?   &           X2C, X2S
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated SB '
        CALL WRTMAT(SB,NROW,NCOL,NROW,NCOL)
      END IF
*
      CALL LUCIAQEXIT('TRADBL')
*
      RETURN
      END
      SUBROUTINE TODSCN2(VEC,NREC,LREC,LBLK,LU,IWR_EOV)
*
* Write VEC as multiple record file accordin to NREC and LREC
* IWR_EOV = 1, then end of vector mark, -1 is written
*
* Statistcs are handled in TODSC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSC(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,0)
        END IF
      END DO
*
      IF(IWR_EOV.EQ.1) THEN
        CALL ITODS(-1,1,LBLK,LU)
      END IF
*
      RETURN
      END
c $Id$
