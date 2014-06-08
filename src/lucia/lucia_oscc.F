      SUBROUTINE GEN_CONFEX(NTCONF,ITCONF,NOCCLS,IBSPOX_FOR_CLS,
     &                      ISPOX_FOR_CLS)
*
* Generate all T-configuration 
* A configuration excitation is by the orbital configurations of 
* of T | ref>
* A T configuration is defined by 
* 1) A set of orbital numbers in ascending order 
*    giving the occupied orbital of T |ref>
* 2) Occupation (1 or 2) or each of the above orbitals
*
* The T-configurations are obtained by generating 
* the spinororbital excitations that are the first 
* members of the various classes of class excitations
*
* Jeppe Olsen, July 2001
*
*. The configuration excitations 
      RETURN
      END
      SUBROUTINE SPOBEX_FOR_OCCLS(
     &           IEXTP_TO_OCCLS,NOCCLS,ISOX_TO_OX,NSOX,
     &           NSOX_FOR_OCCLS,ISOX_FOR_OCCLS,IBSOX_FOR_OCCLS)
*
* Obtain spin-orbital excitation corresponding to given occupation class
*
* Jeppe Olsen, April 24 in Kerkrade The Netherlands
*
* Note : in input NSOX should be the extended number including 
*        the unit operator
      INCLUDE 'implicit.inc'
*.  Input 
      INTEGER IEXTP_TO_OCCLS(NOCCLS),ISOX_TO_OX(NSOX)
*. Output
      INTEGER NSOX_FOR_OCCLS(NOCCLS),IBSOX_FOR_OCCLS(NOCCLS)
      INTEGER ISOX_FOR_OCCLS(NSOX)
*
      IB = 1
      DO IOCCLS = 1, NOCCLS
C?      WRITE(6,*) ' IOCCLS ', IOCCLS
        IBSOX_FOR_OCCLS(IOCCLS) = IB
        N = 0
        DO ISOX = 1, NSOX
C?        WRITE(6,*) ' ISOX, ISOX_TO... ', ISOX,
C?   &                 ISOX_TO_OX(ISOX)
          IF(IEXTP_TO_OCCLS(ISOX_TO_OX(ISOX)).EQ.IOCCLS) THEN
             N = N + 1
             L =  IB - 1 + N
             ISOX_FOR_OCCLS(L) = ISOX
          END IF
        END DO
        IB = IB + N
        NSOX_FOR_OCCLS(IOCCLS) = N
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' spin-orbital excitations to occupations '
        DO IOCCLS = 1, NOCCLS
          WRITE(6,*) ' Occupation class : ', IOCCLS            
          N = NSOX_FOR_OCCLS(IOCCLS)
          IB =  IBSOX_FOR_OCCLS(IOCCLS)
          WRITE(6,*) ' Number of spin orbital orbital types : ', N
          WRITE(6,*) ' The spin-orbital excitations'
          CALL IWRTMA(ISOX_FOR_OCCLS(IB),1,N,1,N)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE DECODE2_LINE(LINE,NCHAR,NENTRY,IENTRY,MXENTRY)
*
* A CHAR line is given.
* Find number of separate items, with each item
* being separated by a , or a space
*
* Jeppe Olsen, Someday in 97 where I really should be doing more
* important things
*
*. Entry
      CHARACTER*(*) LINE
*. Output
      CHARACTER*102 IENTRY(MXENTRY)
*. Local scratch
      CHARACTER*102 CSCR
*
*  
      JITEM=0
      JEFF = 0
*. Last item with a character being nonblank
      ILAST = 0
      DO ICHAR = 1, NCHAR
        IF(LINE(ICHAR:ICHAR).NE.' ')ILAST = ICHAR
      END DO
      IF(ILAST.EQ.0) GOTO 999
*
      DO ICHAR = 0, ILAST
        IF(ICHAR.EQ.0.OR.LINE(ICHAR:ICHAR).EQ.','.OR.
     &     (JITEM.GE.1.AND.LINE(ICHAR:ICHAR).EQ.' ')   ) THEN
*Start of new item, make sure there is space and clean up
          JITEM = JITEM + 1
          IF(JITEM .GT.MXENTRY) THEN
            WRITE(6,*) 'DECODE_LINE:MXENTRY too small'
            WRITE(6,*) ' Number of entries larger than MXENTRY'
            WRITE(6,*) ' JITEM, MXENTRY', JITEM, MXENTRY
            WRITE(6,*)
            WRITE(6,*)
            DO JENTRY = 1, JITEM  
              WRITE(6,'(A,I3,102A)') 'Entry ',JENTRY,IENTRY(JENTRY)
            END DO
            STOP'DECODE_LINE:MXENTRY too small'
          END IF
*. Copy previous entry to permanent
          IF(JITEM.NE.1) THEN
            IENTRY(JITEM-1) = CSCR
          END IF
*. and clean
          DO JCHAR = 1, NCHAR
            CSCR(JCHAR:JCHAR) = ' '
          END DO
          JEFF = 0
        ELSE
*. Continuation of previous item
          JEFF = JEFF + 1
          CSCR(JEFF:JEFF) = LINE(ICHAR:ICHAR)
        END IF
*
      END DO
  999 CONTINUE
*. Transfer last item to permanant residence
      IF(JEFF.NE.0) IENTRY(JITEM) = CSCR
*
      NENTRY = JITEM
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from Decode line '
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        WRITE(6,*) ' Number of separate entries', NENTRY
        WRITE(6,*)
        DO JENTRY = 1, NENTRY
          WRITE(6,'(A,I3,102A)') 'Entry ',JENTRY,IENTRY(JENTRY)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE MASK_CCVEC(ITSS_TP,NTSS_TP,CCVEC,ISM,
     &           MASK_SD,MSK_AEL,MSK_BEL,VALUE,
     &           MX_ST_TSOSO_BLK_MX)
*
* Use determinant MASK_SD to mask CC vector CCVEC
*
*  Operatorstrings containing  creation of spinorbitals occupied in MASK_SD
*  or containing annihilation of spinorbitals occupied in MASK_SD are set to 
*  VALUE
*
*
* Jeppe Olsen, April 23 in Kerkrade, The Netherlands
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'orbinp.inc'
*
      INTEGER MASK_SD(100,2)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MASK_V')
*. Four blocks of string occupations
      CALL MEMMAN(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
*. Expand MASK determinant to full form 
      CALL MEMMAN(KLMSKA,MXPORB,'ADDL  ',1,'MSKA  ')
      CALL MEMMAN(KLMSKB,MXPORB,'ADDL  ',1,'MSKB  ')
      IZERO = 0
      CALL ISETVC(WORK(KLMSKA),IZERO,NTOOB)
      CALL ISETVC(WORK(KLMSKB),IZERO,NTOOB)
      DO IAB = 1, 2
        IF(IAB.EQ.1) THEN
          MEL = MSK_AEL
          KLMSK = KLMSKA
        ELSE
          MEL = MSK_BEL
          KLMSK = KLMSKB
        END IF
        DO IEL = 1, MEL
          IORB = MASK_SD(IEL,IAB)
          CALL ISTVC3(WORK(KLMSK),IORB,1,1)
        END DO
      END DO
*. And then the work
      CALL MASK_CCVECS(ITSS_TP,NTSS_TP,CCVEC,ISM,VALUE,
     &     WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &     WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &     WORK(KLMSKA),WORK(KLMSKB))
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'MASK_V')
*
      RETURN
      END 
      SUBROUTINE MASK_CCVECS(ITSS_TP,NTSS_TP,CCVEC,ISM,VALUE,
     &           IOCC_CA,IOCC_CB,IOCC_AA,IOCC_AB,
     &           MASKA, MASKB)
*
*
* Use determinant MASK_SD=MASKA,MASKB to mask CC vector CCVEC, inner routine
*
*  Operatorstrings containing  creation of spinorbitals occupied in MASK_SD
*  or containing annihilation of spinorbitals occupied in MASK_SD are set to 
*  VALUE
*
*
* Jeppe Olsen, April 23 in Kerkrade, The Netherlands

*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER ITSS_TP(4*NGAS,NTSS_TP)
      INTEGER MASKA(*), MASKB(*)
*. Input and Output
      DIMENSION CCVEC(*)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*
      NTEST = 00
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' MASKA and MASKB vectors '
        CALL IWRTMA(MASKA,1,NTOOB,1,NTOOB)
        CALL IWRTMA(MASKB,1,NTOOB,1,NTOOB)
      END IF
*
      IT = 0
      DO ITSS = 1, NTSS_TP
*. Transform from occupations to groups
       CALL OCC_TO_GRP(ITSS_TP(1+0*NGAS,ITSS),IGRP_CA,1      )
       CALL OCC_TO_GRP(ITSS_TP(1+1*NGAS,ITSS),IGRP_CB,1      )
       CALL OCC_TO_GRP(ITSS_TP(1+2*NGAS,ITSS),IGRP_AA,1      )
       CALL OCC_TO_GRP(ITSS_TP(1+3*NGAS,ITSS),IGRP_AB,1      )
*
       NEL_CA = IELSUM(ITSS_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ITSS_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ITSS_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ITSS_TP(1+3*NGAS,ITSS),NGAS)
*. Diagonal block ?
       CALL DIAG_EXC_CC(ITSS_TP(1+0*NGAS,ITSS),
     &                  ITSS_TP(1+1*NGAS,ITSS),
     &                  ITSS_TP(1+2*NGAS,ITSS),
     &                  ITSS_TP(1+3*NGAS,ITSS),NGAS,IDIAG)
       IF(MSCOMB_CC.EQ.0.OR.IDIAG.EQ.0) THEN
         IRESTRICT = 0
       ELSE 
         IRESTRICT = 1
       END IF
       DO ISM_C = 1, NSMST
         ISM_A = MULTD2H(ISM,ISM_C) 
         DO ISM_CA = 1, NSMST
           ISM_CB = MULTD2H(ISM_C,ISM_CA)
           DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*
            ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
            ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
            IF(IRESTRICT.EQ.1.AND.ISM_BETA.GT.ISM_ALPHA) GOTO 777
            IF(IRESTRICT.EQ.0.OR.ISM_ALPHA.GT.ISM_BETA) THEN
             IRESTRICT_LOOP = 0
            ELSE
             IRESTRICT_LOOP = 1
            END IF
*. obtain strings
            CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,NSTR_CA,
     &           IOCC_CA, NORBT,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,NSTR_CB,
     &           IOCC_CB, NORBT,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,NSTR_AA,
     &           IOCC_AA, NORBT,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,NSTR_AB,
     &           IOCC_AB, NORBT,0,IDUM,IDUM)
*. Loop over T elements as  matric T(I_CA, I_CB, IAA, I_AB)
             DO I_AB = 1, NSTR_AB
              IF(IRESTRICT_LOOP.EQ.1) THEN
                I_AA_MIN = I_AB
              ELSE
                I_AA_MIN = 1
              END IF
              DO I_AA = I_AA_MIN, NSTR_AA
               DO I_CB = 1, NSTR_CB
                IF(IRESTRICT_LOOP.EQ.1.AND.I_AB.EQ.I_AA) THEN
                 ICA_MIN = I_CB
                ELSE
                 ICA_MIN = 1
                END IF
                DO I_CA = ICA_MIN, NSTR_CA
*. Test 
                 IMASK = 0
                 DO IOP = 1, NEL_CA
                   JORB = IOCC_CA((I_CA-1)*NEL_CA+IOP)
                   IF(MASKA(JORB).EQ.1) IMASK = 1
                 END DO
                 DO IOP = 1, NEL_CB
                   JORB = IOCC_CB((I_CB-1)*NEL_CB+IOP)
                   IF(MASKB(JORB).EQ.1) IMASK = 1
                 END DO
                 DO IOP = 1, NEL_AA
                   JORB = IOCC_AA((I_AA-1)*NEL_AA+IOP)
                   IF(MASKA(JORB).EQ.0) IMASK = 1
                 END DO
                 DO IOP = 1, NEL_AB
                   JORB = IOCC_AB((I_AB-1)*NEL_AB+IOP)
                   IF(MASKB(JORB).EQ.0) IMASK = 1
                 END DO
               
                 IT = IT + 1
                 IF(IMASK.EQ.1) CCVEC(IT) = VALUE
                END DO
               END DO
              END DO
             END DO
*            ^ End of loop over elements of block
  777       CONTINUE
            END DO
*           ^ End of loop over ISM_AA
         END DO
*        ^ End of loop over ISM_CA
       END DO
*      ^ End of loop over ISM_C
      END DO
*     ^ End of loop over ITSS
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' Masked vector  '
        CALL WRTMAT(CCVEC,1,IT,1,IT)
      END IF
*
      RETURN
      END
c $Id$
