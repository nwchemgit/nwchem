      SUBROUTINE ADD_SKAIIB(SB,NI,NIA,SKAIIB,NKA,NIB,
     &                          I,ISCA,SSCA)
*
* Update Transposed sigma block with contributions for given orbital index j 
* from the matrix S(Ka,i,Ib)
*
* S(Ib,Isca(Ka)) =  S(Ib,Isca(Ka)) + Ssca(Ka)*S(Ka,I,Ib)
*
*
* For efficient processing of alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
       DIMENSION SKAIIB(*),SSCA(*),ISCA(*)
*. Input and Output
       DIMENSION SB(NIB,NIA)
*. For statistics 
       INCLUDE 'rou_stat.inc'
*
       CALL LUCIAQENTER('SKAII')
       NCALL_ADD_SKAIIB = NCALL_ADD_SKAIIB  + 1
*. To get rid of annoying and incorrect compiler warnings
      ICOFF = 0
*
C     LBLK = 100
      LBLK = 40
      NBLK = NIB/LBLK
      IF(LBLK*NBLK.LT.NIB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NIB)
        ICONST = NKA*NI
        IADR0 =  (I-1)*NKA+(ICOFF-1-1)*NKA*NI
        IF(ICEND.GT.ICOFF) THEN
*. Use form with Inner loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              XOP_ADD_SKAIIB = XOP_ADD_SKAIIB  + ICEND - ICOFF + 1 
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (I-1)*NKA+(ICOFF-1-1)*NKA*NI
              IADR = IADR0 + KA
              DO IB = ICOFF,ICEND
*. Adress of S(Ka,i,Ib)
                IADR = IADR + ICONST
                SB(Ib,IROW) = SB(Ib,IROW)+S*SKAIIB(IADR)
              END DO
            END IF
          END DO
        ELSE
*. Form with no loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              XOP_ADD_SKAIIB = XOP_ADD_SKAIIB  + 1 
              S = SSCA(KA)
              IROW = ISCA(KA)
              IADR = IADR0 + KA + ICONST
              SB(ICOFF,IROW) = SB(ICOFF,IROW)+S*SKAIIB(IADR)
            END IF
          END DO
        END IF
*       ^ End of test of ICOFF=ICEND
      END DO
      CALL LUCIAQEXIT('SKAII')
*
      RETURN
      END
C               GET_CKAJJB(CB,NJ,NJA,CJRES,NKABTC,NJB,
C    &                          JJ,I1(1,JJ),XI1S(1,JJ)
      SUBROUTINE GET_CKAJJB(CB,NJ,NJA,CKAJJB,NKA,NJB,
     &                          J,ISCA,SSCA)
*
* Obtain for given orbital index j the gathered matrix
*
* C(Ka,j,Jb) = SSCA(Ka)C(Jb,ISCA(Ka))
*
* For efficient processing of alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CB(NJB,NJA), SSCA(*),ISCA(*)
*. Output
      DIMENSION CKAJJB(*)
*. For statistics 
      INCLUDE 'rou_stat.inc'
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_CKAJJB for J = ', J
        WRITE(6,*) ' NKA, NJB = ', NKA, NJB
      END IF
*
      NCALL_GET_CKAJJB = NCALL_GET_CKAJJB + 1
*. To get rid of annoying and incorrect compiler warnings
      ICOFF = 0
*
C     LBLK = 100
      LBLK = 40
      NBLK = NJB/LBLK
      IF(LBLK*NBLK.LT.NJB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NJB)
        ICONST = NKA*NJ 
        IADR0 =  (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
*
        IF(ICEND.GT.ICOFF) THEN
*. Inner loop over JB
          DO KA  = 1, NKA
C?          WRITE(6,*) ' KA = ', KA
            XOP_GET_CKAJJB = XOP_GET_CKAJJB + ICEND-ICOFF+1
*
            IF(ISCA(KA).LT.0) THEn
              WRITE(6,*) ' Problem in GET_CK... '
              WRITE(6,*) ' KA, ISCA(KA) = ', KA, ISCA(KA)
              STOP ' Problem IN GET_CKAJJB' 
            END IF
*
C?          WRITE(6,*) ' KA, ISCA(KA) = ', KA, ISCA(KA)
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
*. Adress of C(Ka,j,Jb)
                IADR = IADR + ICONST
C?              WRITE(6,*) 'JB, IADR, IROW, S = ', JB, IADR,IROW,S
                CKAJJB(IADR) = S*CB(JB,IROW)
              END DO
            ELSE  
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
                IADR = IADR + ICONST
C?              WRITE(6,*) ' IADR = ', IADR
                CKAJJB(IADR) = 0.0D0          
              END DO
            END IF
          END DO
        ELSE
*. No inner loop over JB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
              IADR = IADR0 + KA
*. Adress of C(Ka,j,Jb)
              IADR = IADR + ICONST
C?            WRITE(6,*) ' 
C?   &        JB, IADR, ICOFF, IROW, S = ', JB, IADR,ICOFF,IROW,S
              CKAJJB(IADR) = S*CB(ICOFF,IROW)
            ELSE  
              IADR = IADR0 + KA
              IADR = IADR + ICONST
C?            WRITE(6,*) ' IADR = ', IADR
              CKAJJB(IADR) = 0.0D0          
            END IF
          END DO
        END IF
*       ^ End of test ICEND,ICOFF
      END DO
*
      RETURN
      END
      SUBROUTINE NXTNUM(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL(I),MAXVAL(I). 
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen Oct 1994
*
*. Input
      DIMENSION MINVAL(*),MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
            INUM(JPLACE) = MINVAL(JPLACE)
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
      SUBROUTINE NUMST4(NEL,NORB1,NEL1MN,NEL1MX,NORB2,
     &                  NORB3,NEL3MN,NEL3MX,NSTTP)
*
* Number of strings per type for class of strings with 
*
* Between NEL1MN AND NEL1MX electrons in the first NORB1 orbitals
* Between NEL3MN AND NEL3MX electrons in the last  NORB3 orbitals
*
      INTEGER NSTTP(*)
      NTEST = 0
      NSTRIN = 0
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' NUMST4 '
        WRITE(6,*) ' NEL NEL1MN NEL1MX NEL3MN NEL3MX '
        WRITE(6,*) NEL,NEL1MN,NEL1MX,NEL3MN,NEL3MX
      END IF
*
      ITPMAX = 0
      DO 100 IEL1 = NEL1MN,MIN(NEL1MX,NORB1,NEL)
        NSTIN1 = IBION(NORB1,IEL1)
        IEL3MN = MAX ( NEL3MN,NEL-(IEL1+NORB2) )
        IEL3MX = MIN ( NEL3MX,NEL-IEL1)
        DO 80 IEL3 = IEL3MN, IEL3MX
         IEL2 = NEL - IEL1-IEL3
         NSTINT = NSTIN1*IBION(NORB2,IEL2)*IBION(NORB3,IEL3)
         NSTRIN = NSTRIN + NSTINT
         ITP = (NEL1MX-IEL1)
     &         * (NEL3MX-NEL3MN+1)          
     &         + IEL3-NEL3MN + 1        
         NSTTP(ITP) = NSTINT
         ITPMAX= MAX(ITPMAX,ITP)
  80   CONTINUE
 100  CONTINUE
*
      IF( NTEST .GE.1 ) THEN
      WRITE(6,'(/A,I6)') '  Number of strings generated ... ', NSTRIN
         WRITE(6,*)  
         WRITE(6,*) ' Largest type number ', ITPMAX
         WRITE(6,*) ' Number of strings per type '
         CALL IWRTMA(NSTTP,1,ITPMAX,1,ITPMAX)
      END IF
*
      RETURN
      END
      FUNCTION NCASTR(IAC,NSTTPI,NTPSTI,ICLSI,NOBATP,NOBTP,IELPTP)
*
* Number of allowed annihilation/creations from a given group 
* of stings 
*
*     Jeppe Olsen, June 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  NSTTPI(NTPSTI)
      INTEGER  NOBATP(NOBTP)
      INTEGEr IELPTP(NOBTP,*)
*
      LCA = 0
      DO IOBTP = 1, NOBTP
        DO ISTTP = 1, NTPSTI
*. Type of resulting string
          CALL NEWTYP(ISTTP,IAC,IOBTP,1,ITPO)
C?        WRITE(6,*) ' IOBTP ISTTP => ITPO, ICLSO '
C?        WRITE(6,*)   IOBTP,ISTTP,ITPO,ICLSO
          IF(IAC.EQ.1) THEN
            NENTRY = IELPTP(IOBTP,ISTTP)
          ELSE
            NENTRY = NOBATP(IOBTP)-IELPTP(IOBTP,ISTTP)
          END IF
          IF(ITPO.GT.0) THEN
            LCA = LCA + NENTRY*NSTTPI(ISTTP)
          END IF
*
        END DO
      END DO
*
      WRITE(6,*) ' Number of generated strings ', LCA
      NCASTR = LCA
*
      RETURN
      END
         
      SUBROUTINE CIDIA5(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,DIAG,NSMST,H,
     &                  ISMOST,IBLTP,XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,IOCOC,NOCTPA,NOCTPB,
     &                  ISSOA,ISSOB,LUDIA,ECORE,
     &                  PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,I12)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, February 1994 , obtained from CIDIA4
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
      IMPLICIT REAL*8           (A-H,O-Z)
*.General input
      DIMENSION NSSOA(NOCTPA,*),NSSOB(NOCTPB,* )
      DIMENSION ISSOA(NOCTPA,*),ISSOB(NOCTPB,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ISMOST(*),IBLTP(*)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION DIAG(*)
*
      NTEST =  0
      NTEST = MAX(NTEST,IPRNT)
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
*
*. To get rid of annoying and incorrect compiler warnings
      IOFF = 0
*
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NORB,1,NORB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
        END IF
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
      IDET = 0
      ITDET = 0
      IF(LUDIA.NE.0) CALL REWINO(LUDIA)
      DO 1000 IASM = 1, NSMST
        IBSM = ISMOST(IASM)
        IF(IBSM.EQ.0.OR.IBLTP(IASM).EQ.0) GOTO 1000
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
        DO 999  IATP = 1,NOCTPA
          IF(IREST1.EQ.1) THEN
            MXBTP = IATP
          ELSE
            MXBTP = NOCTPB
          END IF
*
*. Will strings of this type be used ?
*
          IUSED = 0
          DO IBTP = 1, MXBTP
            IF( NSSOB(IBTP,IBSM).NE.0.AND.
     &      IOCOC(IATP,IBTP).NE.0) IUSED = 1
          END DO 
          IF( IUSED .EQ. 0 ) GOTO 987
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
          IOFF =  ISSOA(IATP,IASM)
          DO IA = IOFF,IOFF+NSSOA(IATP,IASM)-1
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
  987     CONTINUE
*
          DO 900 IBTP = 1,MXBTP
          IF(IOCOC(IATP,IBTP) .EQ. 0 ) GOTO 900
          IBSTRT = ISSOB(IBTP,IBSM)
          IBSTOP = IBSTRT + NSSOB(IBTP,IBSM)-1
          DO 899 IB = IBSTRT,IBSTOP
            IBREL = IB - IBSTRT + 1
*
*. Terms depending only on IB
*
            HB = 0.0D0
            RJBB = 0.0D0
            CALL SETVEC(XB,0.0D0,NORB)
*
            DO 990 IEL = 1, NBEL
              IBEL = IBSTR(IEL,IB)
              HB = HB + H(IBEL )
*
              IF(I12.EQ.2) THEN
                DO 980 JEL = 1, NBEL
                  RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
  980           CONTINUE
*
                DO 970 IORB = 1, NORB
                  XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
  970           CONTINUE
              END IF
  990       CONTINUE
            EB = HB + 0.5D0*RJBB + ECORE
*
            IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
              IASTRT = ISSOA(IATP,IASM) - 1 + IBREL
            ELSE
              IASTRT = ISSOA(IATP,IASM)
            END IF
            IASTOP = ISSOA(IATP,IASM) + NSSOA(IATP,IASM) - 1
            DO 800 IA = IASTRT,IASTOP
              IDET = IDET + 1
              ITDET = ITDET + 1
              X = EB + RJKAA(IA-IOFF+1)
              DO 890 IEL = 1, NAEL
                X = X +XB(IASTR(IEL,IA)) 
  890         CONTINUE
              DIAG(IDET) = X
              IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
  800       CONTINUE
  899     CONTINUE
*. Yet a RAS block of the diagonal has been constructed
          IF(ICISTR.GE.2) THEN
            IF(NTEST.GE.100) THEN
              write(6,*) ' number of diagonal elements to disc ',IDET
              CALL WRTMAT(DIAG,1,IDET,1,IDET)
            END IF
            CALL ITODS(IDET,1,-1,LUDIA)
            CALL TODSC(DIAG,IDET,-1,LUDIA)
            IDET = 0
          END IF
  900   CONTINUE
  999   CONTINUE
*. Yet a symmetry block of the diagonal has been constructed
*
 1000 CONTINUE
 
      IF(NTEST.GE.10) WRITE(6,*)
     &' Number of diagonal elements generated ',ITDET
*
      IF(NTEST .GE.100 .AND.ICISTR.LE.1 ) THEN
        WRITE(6,*) ' CIDIAGONAL '
        CALL WRTMAT(DIAG(1),1,IDET,1,IDET)
      END IF
*
      IF ( ICISTR.GE.2 ) CALL ITODS(-1,1,-1,LUDIA)
*
      RETURN
      END


      SUBROUTINE DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
     &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
*
* Multiply/divide elements of two vectors residing on disc
* Elements corresponding to absolute adresses in ISCAT
* are set to the elements of XSCAT
*
* For elements not in ISCAT the operation is thus :
*
* For INV.NE. 0 :
*
*    V3(I) = (V1(I)+FAC)-1 * V2(I)
*    LU3      LU1            LU2
* 
* For INV.NE. 0 :
* 
*    V3(I) = (V1(I)+FAC) * V2(I)
*    LU3         LU1        LU2
*
* LBLK defines structure of files
*
* Overlap between input and output vector is also calculated
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  VEC1(*),VEC2(*)
      DIMENSION  XSCAT(NSCAT),ISCAT(NSCAT)
*
      REAL*8 INPROD
*
      XINOUT = 0.0D0
*
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          REWIND LU2
          REWIND LU3
        ELSE
          CALL REWINE( LU1,LBLK)
          CALL REWINE( LU2,LBLK)
          CALL REWINE( LU3,LBLK)
         END IF
      END IF
*
* Loop over blocks
*
        IOFF = 1
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LU1) LBL1
          READ(LU2) LBL2
          WRITE(LU3) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LU1)
          CALL IFRMDS(LBL2,1,-1,LU2)
          CALL ITODS (LBL1,1,-1,LU3)
        END IF
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVSD : '
     &                     , LBL1,LBL2
          STOP ' DIFFERENT BLOCKSIZES IN DMTVSD '
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          LENGTH = LBL1
          CALL FRMDSC(VEC1,LENGTH,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC(VEC2,LENGTH,KBLK,LU2,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
            IF(INV .NE. 0 ) THEN
              CALL DIAVC2(VEC1,VEC2,VEC1,FAC,LENGTH)
            ELSE
              CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
              CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FAC,LENGTH)
            END IF
          END IF
*
          IF(NSCAT.NE.0) THEN
            IFIRST = IOFF
            ILAST = IOFF + LENGTH - 1
            DO 100 I = 1, NSCAT
              IF(IFIRST .LE. ISCAT(I) .AND. ISCAT(I) .LE. ILAST ) 
     &        VEC1(ISCAT(I)-IOFF+1) = XSCAT(I)
  100       CONTINUE
          END IF
*
          XINOUT = XINOUT + INPROD(VEC1,VEC2,LENGTH)
          CALL TODSC(VEC1,LENGTH,KBLK,LU3)
          IOFF = IOFF + LENGTH 
        END IF
       
*
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      RETURN
      END
      SUBROUTINE H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
     &                  H0,SHIFT,XWORK,XH0PSX,
     &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG,IPRECOND,
     &                  ISBSPPR_ACT)
*
* Calculate inverted general preconditioner matrix times vector
*
* Disc version
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  LUOUT       LUDIA        LUIN
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* H0: Form depends on IPRECOND
*
* IPRECOND = 1
*==============
*
* H0 is the diagonal on LUDIA and a matrix of dimension NP1 
* Full diagonal form of matrix is stored/hidden in H0
*
* IPRECOND = 2:
* ==============
*
* H0 is the diagonal on file LUDIA
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements on LUDIA corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* If ISBSPPR_ACT ne 0, then this space is considered as the P1 space, no P2 space 
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , September 1993
*               May 99        : IPRECOND added 
*               Jan 13        : ISBSPPR_ACT added, includes extensive changes
*               July 13       : Going back to simple one-blocks precond.
*
*
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File Containing diagonal of H0 
* LUIN  : File Containing input vector   
* LBLK : Defines format of files
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* XWORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
* ISBSPPR_ACT: CI space defining P1 space, included through direct procedures
*
* Notice: NP1, NP2, NQ are MAX, not actual dimensions!! NP1.. defines
*         therefore offsets
* NPQDM is actual dimension
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector 
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
      INCLUDE 'implicit.inc'
COLD  REAL * 8  INPROD
*
      DIMENSION VEC1(*),VEC2(*)          
      DIMENSION IPNTR(*),H0(*)
      DIMENSION XWORK(*)
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from H0M1TD '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' IPRECOND, NPQDM = ', IPRECOND, NPQDM
      END IF
*
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
* 
*. Explicit Hamiltonian in subspace
*
*. Obtain subspace components of input vector
        IZERO = 0
C            GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,NTESTG)
        CALL GATVCD(LUIN,LBLK,NPQDM,IPNTR,XWORK(KLV1),VEC1,
     &              IZERO)
*
        IF(IPRECOND.EQ.1) THEN
*. Simple one space preconditioner
*. Eigen-vectors and -values are stored in H0, adresses...
          KLH0_SUBDT = 1
          KLH0_MAT = KLH0_SUBDT + NP1
          KLH0_EIGVAL = KLH0_MAT +  NP1*(NP1+1)/2
          KLH0_EIGVEC = KLH0_EIGVAL + NP1
*
C         XDXTV(VECUT,VECIN,X,DIA,NDIM,SCR,SHIFT,IINV)
          CALL XDXTV(XWORK(KLV2),XWORK(KLV1),
     &         H0(KLH0_EIGVEC),H0(KLH0_EIGVAL),
     &         NPQDM,XWORK(KLSCR),SHIFT,1)
        ELSE IF (IPRECOND.EQ.2) THEN
*. Solve linear equations in subspace
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               XWORK(KLV2),XWORK(KLV1),SHIFT,XWORK(KLSCR),
     &               NTEST )
        END IF ! PRECOND switch
      END IF
*
*. Calculate inverse diagonal and scatter results from subspace,
*. Write to file LUOUT
C     DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
C    &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
      CALL DMTVDS(VEC1,VEC2,LUDIA,LUIN,LUOUT,SHIFT,1,1,
     &            IPNTR,XWORK(KLV2),NPQDM,LBLK,XH0PSX)
      IF(IPRECOND.EQ.3) THEN
C     HCONFINTV(LURHS,LUX,SHIFTG,SHIFT_DIAG,VECIN,VECOUT,
C    &                LBLK,LUPROJ,LUPROJ2,LLUDIA)
        CALL HCONFINTV(LUIN,LUOUT,SHIFT,SHIFT,VEC1,VEC2,
     &                LBLK,0,0,LUDIA)
      END IF
*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TD '
        WRITE(6,*) ' ========================= '
*. Note : works only if result vector is first file on LUOUT 
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)   
        WRITE(6,*) ' Overlap between input and output vector',
     &               XH0PSX
      END IF
*
      RETURN
      END
*
      SUBROUTINE MICDV4(MV7,VEC1,VEC2,LU1,LU2,RNRM,EIG,
     &                  FINEIG,MAXIT,NVAR,
     &                  LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,CONVER,RNRM_CNV,ISBSPPR_ACT)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Input :
* =======
*        MV7: Name of sigma routine
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENLVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENLVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
C      REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
       DIMENSION RNRM_CNV(NROOT)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2                               
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(10)
       REAL*8 INPRDD
*
       EXTERNAL MV7
*
       IF(IPRT.GT.5) WRITE(6,*) ' MICDV4 in action '
C?     WRITE(6,*) ' TESTY: ISBSPPR_ACT = ', ISBSPPR_ACT
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
 
       IOLSTM = 0
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MICDV4 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop on MICDV4 '
         STOP 20
       END IF
       IF(IPRT.GE.5) WRITE(6,'(A,F13.7)') 'EIGSHF = ', EIGSHF
*
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       DO 10 IVEC = 1,NINVEC
         CALL REWINO(LU3)
         CALL REWINO(LU4)
         CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
         CALL MV7(VEC1,VEC2,LU3,LU4,0,0)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINO(LU4)
         CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINO(LU2)
         DO 8 JVEC = 1, IVEC
           CALL REWINO(LU3)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU3,0,LBLK)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.3 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
C        CALL PRSYM(APROJ,NINVEC)
         CALL PRSYM_EP(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGENL(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF(IPRT .GE. 10 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(6,'(5F21.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 10 ) THEN
         WRITE(6,*) ' Initial set of eigen values (no shift) '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        IF(IPRT  .GE. 10 ) THEN
         WRITE(6,*) ' Info from iteration .... ', ITER
        END IF
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
         EIGAPR = EIG(ITER-1,IROOT)
*
         CALL REWINO(LU1)
         CALL REWINO(LU2)
         EIGAPR = EIG(ITER-1,IROOT)
         DO 60 IVEC = 1, NVEC
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           IF(IVEC.EQ.1) THEN
             CALL REWINO( LU3 )
*                 SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
             CALL SCLVCD(LU2,LU3,FACTOR,VEC1,0,LBLK)
           ELSE
             CALL REWINO(LU3)
             CALL REWINO(LU4)
C                 VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU2,LU3,0,LBLK)
           END IF
C
           FACTOR = (-EIGAPR)*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
           CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU1,LU4,0,LBLK)
   60    CONTINUE
         IF ( IPRT  .GE. 10000 ) THEN
           WRITE(6,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU4,1,LBLK)
         END IF
*  Strange place to put convergence but ....
C                      INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
         IF(RNORM.LT. TEST .OR. 
     &      (ITER.GT.2.AND.
     &      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.E_CONV)) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU3
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.10) THEN
             WRITE(6,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINO(LUDIA)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
*. Assuming diagonal preconditioner
           IPRECOND = 1
           CALL H0M1TD(LU3,LUDIA,LU4,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND,ISBSPPR_ACT)
C               H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG,IPRECOND,ISBSPPR_ACT)
           IF ( IPRT  .GE.200) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
           X2NORM = INPRDD(VEC1,VEC2,LU3,LU3,1,-1)
C?         WRITE(6,*) 'TEST in MICDC4,Norm of resid and dia**-1*resid',
C?   &     RNORM**2, X2NORM
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* Current eigen-vector on LU5
             CALL REWINO(LU1)
             DO 66 IVEC = 1, NVEC
               FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
               IF(IVEC.EQ.1) THEN
                 IF(NVEC.EQ.1) THEN
                   CALL REWINO( LU5 )
                   CALL SCLVCD(LU1,LU5,FACTOR,VEC1,0,LBLK)
                 ELSE
                   CALL REWINO( LU4 )
                   CALL SCLVCD(LU1,LU4,FACTOR,VEC1,0,LBLK)
                 END IF
               ELSE
                 CALL REWINO(LU5)
                 CALL REWINO(LU4)
                 CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU1,LU5,0,LBLK)
                 CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
               END IF
   66        CONTINUE
             IF ( IPRT  .GE. 200) THEN
               WRITE(6,*) '  (current  X ) '
               CALL WRTVCD(VEC1,LU5,1,LBLK)
             END IF
* (H0 - E )-1  * X on LU4
             CALL REWINO(LU5)
             CALL REWINO(LU4)
             CALL REWINO(LUDIA)
*
             CALL H0M1TD(LU4,LUDIA,LU5,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND,ISBSPPR_ACT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPRDD(VEC1,VEC2,LU5,LU4,1,LBLK)
* is X an eigen vector for (H0 - 1 ) - 1
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU5,LU4,1,LBLK))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU5,LU3,1,LBLK)
                FACTOR = (-DELTA)/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
                CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*
                IF(IPRT.GE.200) THEN
                  WRITE(6,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU3,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
           CALL REWINE( LU1 ,LBLK)
           DO 80 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             WORK(IVEC) = INPRDD(VEC1,VEC2,LU1,LU3,0,LBLK)
C?       WRITE(6,*) ' MICDV4 : Overlap ', WORK(IVEC)
   80      CONTINUE
*
           CALL REWINE(LU1,LBLK)
           DO 82 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             CALL REWINE(LU4,LBLK)
             CALL VECSMD(VEC1,VEC2,-WORK(IVEC),1.0D0,LU1,LU3,
     &                   LU4,0,LBLK)
             CALL COPVCD(LU4,LU3,VEC1,1,LBLK)
   82      CONTINUE
           IF ( IPRT  .GE. 200 ) THEN
             WRITE(6,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU3,LBLK)
           CALL SCLVCD(LU3,LU1,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.200) THEN
             CALL SCLVCD(LU3,LU4,FACTOR,VEC1,1,LBLK)
             WRITE(6,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU4,1,LBLK)
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,NVEC,VEC1,1,LBLK)
      DO 150 IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
        CALL MV7(VEC1,VEC2,LU3,LU4,0,0)
        CALL REWINE(LU4,LBLK)
        CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO 140 JVEC = 1, NVEC+IVEC
          CALL REWINE(LU4,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU4,0,LBLK)
  140   CONTINUE
  150 CONTINUE
      NVEC = NVEC + IADD
      IF(IPRT.GE.100) THEN
         WRITE(6,*) ' Updated subspace matrix '
         CALL  PRSYM(APROJ,NVEC)
      END IF
*. Diagonalize projected matrix
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGENL(WORK(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = (-C1NRM)/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F21.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(6,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(6,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 10 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset      
*
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINE( LU5,LBLK)
        DO 320 IROOT = 1, NROOT
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
          CALL REWINE(LU3,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
  320   CONTINUE
*. Transfer C vectors to LU1
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 411 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU1,VEC1,0,LBLK)
  411   CONTINUE
*. corresponding sigma vectors
        CALL REWINE (LU5,LBLK)
        CALL REWINE (LU2,LBLK)
        DO 329 IROOT = 1, NROOT
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,WORK(IROOT),VEC1,0,LBLK)
  329   CONTINUE
*
* Transfer HC's to LU2
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 400 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU2,VEC1,0,LBLK)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
*
      END IF
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
 
* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      DO IROOT = 1, NROOT
        FINEIG(IROOT) = EIG(ITER,IROOT) + EIGSHF
        RNRM_CNV(IROOT) = RNRM(ITER,IROOT)
      END DO
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F21.10)
          IF(IPRT.GE.1000) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F21.13,2X,E12.5)
 1600   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
*
      SUBROUTINE TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
C
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C              LOWER HALF OF AUTPAK IS STORED IN APAK
C   IWAY = 2 : PACKED TO FULL FORM
C              APAK STORED IN LOWER HALF
C               SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
C.. NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
*
* Some considerations on cache minimization used for IMET = 2 Loop
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
*. To get rid of annoying and incorrect compiler warnings
      IOFF = 0
      JOFF = 0
*
*. Packing : No problem with cache misses    
*
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO J = 1,NDIM
          DO I = J , NDIM
            APAK(IJ+I) = AUTPAK(I,J)
          END DO
          IJ = IJ +NDIM-J
        END DO   
      END IF
*
* Unpacking : cache misses can occur so two routes
*
      IF( IWAY .EQ. 2 ) THEN
*. Use blocked algorithm
      IMET = 2
      IF(IMET.EQ.1) THEN
*. No blocking
        IJ = 0
        DO J = 1,NDIM
          DO I = J,NDIM
           AUTPAK(J,I) = SIGN*APAK(IJ+I)
           AUTPAK(I,J) = APAK(IJ+I)
          END DO
          IJ = IJ + NDIM-J
        END DO
      ELSE IF (IMET .EQ. 2 ) THEN
*. Blocking
        LBLK = 40
        NBLK = MATDIM/LBLK
        IF(LBLK*NBLK.LT.MATDIM) NBLK = NBLK + 1
        DO JBLK = 1, NBLK
          IF(JBLK.EQ.1) THEN
            JOFF = 1
          ELSE 
            JOFF = JOFF + LBLK
          END IF
          JEND = MIN(JOFF+LBLK-1,MATDIM)
          DO IBLK = JBLK, NBLK
            IF(IBLK.EQ.JBLK) THEN
              IOFF = JOFF
            ELSE
              IOFF = IOFF + LBLK
            END IF
            IEND = MIN(IOFF+LBLK-1,MATDIM)
              DO J = JOFF,JEND
                IF(IBLK.EQ.JBLK) THEN
                  IOFF2 = J
                ELSE 
                  IOFF2 = IOFF
                END IF
                IJOFF = (J-1)*MATDIM-J*(J-1)/2
                DO I = IOFF2,IEND
                  AUTPAK(J,I) = SIGN*APAK(IJOFF+I)
                  AUTPAK(I,J) = APAK(IJOFF+I)
                END DO
              END DO
*. End of loop over I and J
            END DO
          END DO
*. End of loop over blocks of I and J
        END IF
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPK3 '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSM2(APAK,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE TRIPK3G(AUTPAK,APAK,ICOFF,IWAY,MATDIM,NDIM,SIGN)
C
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C              LOWER HALF OF AUTPAK IS STORED IN APAK
C   IWAY = 2 : PACKED TO FULL FORM
C              APAK STORED IN LOWER HALF
C               SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
C.. NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
*
* Some considerations on cache minimization used for IMET = 2 Loop
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      DIMENSION AUTPAK(MATDIM,MATDIM)
      integer APAK,ICOFF,mybyte,icast,ihandle,KLMAT1
*. To get rid of annoying and incorrect compiler warnings
      IOFF = 0
      JOFF = 0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRIPK3')
      CALL MEMMAN(ITEMP,NDIM*(NDIM+1)/2,'ADDL  ',2,'ITEMP ')
*
*. Packing : No problem with cache misses    
*
      IF( IWAY .EQ. 1 ) THEN
        call luciaqenter('tr3g1')
        IJ = 0
        DO J = 1,NDIM
          DO I = J , NDIM
            dbl_mb(ITEMP+IJ) = AUTPAK(I,J)
            IJ = IJ + 1
          END DO
        END DO   
        call ga_put(APAK,ICOFF,ICOFF+NDIM*(NDIM+1)/2-1,1,1,
     &              dbl_mb(ITEMP),1)
        call luciaqexit('tr3g1')
      END IF
*
* Unpacking : cache misses can occur so two routes
*
      IF( IWAY .EQ. 2 ) THEN
*. No blocking
        call luciaqenter('tr3g2')
        call ga_get(APAK,ICOFF,ICOFF+NDIM*(NDIM+1)/2-1,1,1,
     &             dbl_mb(ITEMP),1)
        IJ = 0
        DO J = 1,NDIM
          DO I = J,NDIM
           AUTPAK(J,I) = SIGN*dbl_mb(ITEMP+IJ)
           AUTPAK(I,J) = dbl_mb(ITEMP+IJ)
           IJ = IJ + 1
          END DO
        END DO
        call luciaqexit('tr3g2')
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPK3 '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
CNW     CALL PRSM2(APAK,NDIM)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRIPK3')
      RETURN
      END
      SUBROUTINE GATRMT(MATIN,NROWIN,NCOLIN,MATUT,NROWUT,NCOLUT,
     &                  ISCA,SSCA)
*
* Gather rows of transposed matrix MATIN  to  MATUT
*
* MATUT(I,J) = SSCA(I)*MATIN(J,ISCA(I)),(ISCA(I) .ne. 0 )
*
*      
* Jeppe Olsen, Getting LUCIA in shape , Feb1994
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MATIN,MATUT
*.Input
      DIMENSION ISCA(*),SSCA(*),MATIN(NCOLIN,NROWIN)
*. ( MATIN Transposed )
*.Output
      DIMENSION MATUT(NROWUT,NCOLUT)
*
*. To get rid of annoying and incorrect compiler warnings
      ICOFF = 0
C     LBLK = 100
      LBLK = 40
      NBLK = NCOLUT/LBLK
      IF(LBLK*NBLK.LT.NCOLUT) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NCOLUT)
        DO I = 1, NROWUT
          IF(ISCA(I).NE.0) THEN
            S = SSCA(I)
            IROW = ISCA(I)
            DO J = ICOFF,ICEND
              MATUT(I,J) = S*MATIN(J,IROW)
            END DO
          ELSE IF (ISCA(I).EQ.0) THEN
            DO J = ICOFF,ICEND
              MATUT(I,J) = 0.0D0
            END DO
          END IF
        END DO
      END DO
*
      RETURN
      END 
      SUBROUTINE SCARMT(MATIN,NROWIN,NCOLIN,MATUT,NROWUT,NCOLUT,
     &                  ISCA,SSCA)
*
* Scatter-add  rows of MATIN to transposed matrix MATUT

*  MATUT(J,ISCA(I)) = MATUT(J,ISCA(I)) + SSCA(I)*MATIN(I,J) 
*  ( if INDEX(I).ne.0 )
*      
* Jeppe Olsen, Getting LUCIA in shape , Feb1994
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MATIN,MATUT
*.Input
      DIMENSION ISCA(*),SSCA(*),MATIN(NROWIN,NCOLIN)
*.Input and Output
      DIMENSION MATUT(NCOLUT,NROWUT)
*.                 (MATUT transposed !)
*. To get rid of annoying and incorrect compiler warnings
      ICINOF = 0
*
C     LBLK = 100
      LBLK = 40
      NBLK = NCOLIN/LBLK
      IF(LBLK*NBLK.LT.NCOLIN) NBLK = NBLK + 1
      DO ICINBL = 1, NBLK
        IF(ICINBL.EQ.1) THEN
          ICINOF = 1
        ELSE
          ICINOF = ICINOF + LBLK
        END IF
        ICINEN = MIN(ICINOF+LBLK-1,NCOLIN)
        DO I = 1, NROWIN
          IF(ISCA(I).NE.0) THEN
            S = SSCA(I)
            IROW = ISCA(I)
            DO ICOL = ICINOF,ICINEN
              MATUT(ICOL,IROW) = MATUT(ICOL,IROW)+S*MATIN(I,ICOL)
            END DO
          END IF
        END DO
      END DO
*
      RETURN
      END 
      SUBROUTINE TRPAD3(MAT,FACTOR,NDIM)
*
*  MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
*
*. With some considerations of effective cache use for large 
*  matrices
*
      IMPLICIT REAL*8           (A-H,O-Z)
      REAL*8            MAT(NDIM,NDIM)
              FAC2 = 1.0D0 - FACTOR**2
*
C     IWAY = 1
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*
*. No blocking
*
*. Lower half
        DO  J = 1, NDIM
          DO  I = J, NDIM
            MAT(I,J) =MAT(I,J) + FACTOR * MAT(J,I)
          END DO
        END DO
*. Upper half
        IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
          FAC2 = 1.0D0 - FACTOR**2
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
            END DO
          END DO
        ELSE IF(FACTOR .EQ. 1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = MAT(I,J )
            END DO
          END DO
        ELSE IF(FACTOR .EQ. -1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) =-MAT(I,J )
            END DO
          END DO
        END IF
      ELSE IF(IWAY .EQ. 2 ) THEN
*. Simple blocking of matrix
        LBLK = 40
        NBLK = NDIM/LBLK
        IF(NBLK*LBLK.LT.NDIM) NBLK = NBLK + 1
        IOFF = 1-LBLK
C?      write(6,*) 'NBLK ',nblk
        DO IBLK = 1, NBLK
          IF(IBLK.EQ.-1) write(6,*) 'IBLK = ',IBLK
          IOFF = IOFF + LBLK
          IEND = MIN(IOFF+LBLK-1,NDIM)
          JOFF = 1 - LBLK
          DO JBLK = 1, IBLK
            JOFF = JOFF + LBLK
            JEND = MIN(JOFF+LBLK-1,NDIM)
*. Lower half
            DO  I = IOFF,IEND
              IF(IBLK.EQ.JBLK) JEND = I   
              DO J = JOFF,JEND
                MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
              END DO
            END DO
*. Upper half
            IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
              FAC2 = 1.0D0 - FACTOR**2
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND 
                  MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
                 END DO
               END DO
            ELSE IF(FACTOR .EQ. 1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I -1
                DO J = JOFF, JEND 
                  MAT(J,I) = MAT(I,J )
                END DO
              END DO
            ELSE IF(FACTOR .EQ. -1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND 
                  MAT(J,I) = - MAT(I,J )
                END DO
              END DO
            END IF
*. ENd of loop over blocks
          END DO
        END DO
*. End of IWAY branching
      END IF
      RETURN
      END
      SUBROUTINE TRPMT3(XIN,NROW,NCOL,XOUT)
*
* XOUT(I,J) = XIN(J,I)
*
*. With a few considerations for large scale cases with cache minimization
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XIN(NROW,NCOL),XOUT(NCOL,NROW)
*
      INCLUDE 'rou_stat.inc'
*
      NCALL_TRPMT = NCALL_TRPMT + 1
      XOP_TRPMT = XOP_TRPMT + NCOL*NROW
*
*. To get rid of annoying and incorrect compiler warnings
      IROFF = 0
      ICOFF = 0
*
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*. Straightforward, no blocking
        DO IROW =1, NROW
          DO ICOL = 1, NCOL
            XOUT(ICOL,IROW) = XIN(IROW,ICOL)
          END DO   
        END DO
      ELSE IF(IWAY.EQ.2) THEN
*. Simple blocking of matrix
        LRBLK = 40 
        LCBLK = 40 
        NRBLK = NROW/LRBLK
        NCBLK = NCOL/LCBLK
        IF(LRBLK*NRBLK.NE.NROW) NRBLK = NRBLK + 1
        IF(LCBLK*NCBLK.NE.NCOL) NCBLK = NCBLK + 1
*
        DO IRBLK = 1,NRBLK
          IF(IRBLK.EQ.1) THEN
            IROFF = 1
          ELSE
            IROFF = IROFF + LRBLK
          END IF
          IREND = MIN(NROW,IROFF+LRBLK-1)
          DO ICBLK = 1, NCBLK
            IF(ICBLK.EQ.1) THEN
              ICOFF = 1
            ELSE
              ICOFF = ICOFF + LCBLK
            END IF
            ICEND = MIN(NCOL,ICOFF+LCBLK-1)
*
            DO IROW = IROFF,IREND
              DO ICOL = ICOFF,ICEND
                XOUT(ICOL,IROW) = XIN(IROW,ICOL)
              END DO
            END DO
*
          END DO
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE MXRESCPH(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                  NSMST,NSTFSMSPGP,
     &                  MXPNSMST,
     &                  NSMOB,MXPTOB,NTPOB,NTSOB,NTESTG,MXPKA,
     &                  NEL1234,
     &                  MXCJ,MXCIJA,
     &                  MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,
     &                  NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,
     &                  MX_NSPII)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
* Find largest dimension of matrix S(P,Ia,I,K) for a single K-string
*
* Particle hole version : hole electrons added, particle elec removed
*
* Largest block of single excitations MXSXBL

*. Input 
* IAB :allowed combination of alpha and beta supergroups 
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups
*
* Version of Jan 98 : IPHGAS added
* June 2003, playing a bit around with multisymmetry blocking

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXPNSMST,*)
      DIMENSION NTSOB(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
      DIMENSION IPHGAS(*)
      INTEGER NHLFSPGP(*)
* 
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from MXRESCPH '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' MXRESCPH : MXPKA ', MXPKA
        WRITE(6,*) ' NTEST = ', NTEST
      END IF
*
* matrix C(j,Ka,Ib)
*
      MXKAB = 0
      MXCJ = 0
      MXCJ_ALLSYM = 0
      MXADKBLK = 0
      MXADKBLK_AS = 0
      MX_NSPII = 0
      DO IAORC= 1,2
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          
          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100) 
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            ITOTB = 0
            DO 210 ISM = 1, NSMST
              MXB =MAX(MXB,NSTFSMSPGP(ISM,IBTPABS))
              ITOTB = ITOTB + NSTFSMSPGP(ISM,IBTPABS)
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB,ITOTB = ', MXB,ITOTB
            DO 300 IOBTP = 1, NTPOB
*. No K strings obtained from creation in particle space 
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) GOTO 300
*. type of K string obtained 
              CALL NEWTYP(IATPABS,IAORC,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
*. addi constraint to avoid calc with long columns and few rows
*. Works only in connection with active advice routine !
              IF(KATP.GT.0) THEN
                IF(IAORC.EQ.1.AND.IADVICE.EQ.1.AND.
     &          NHLFSPGP(IBTPABS)+NHLFSPGP(KATP).LT.MNHL.AND.
     &          NHLFSPGP(IATPABS).GT.(NHLFSPGP(IBTPABS)+1)) THEN
C                 WRITE(6,*) ' N-1 hole space eliminated '
C                 WRITE(6,*) ' IOBTP,IBTPABS,KATP',
C    &            IOBTP,IBTPABS,KATP
C                 KATP = 0
*. Eliminated Jan 2008, to avoid problem...
                END IF
              END IF
*
              IF(KATP.GT.0) THEN 
C                    DIM_SPII(IASPGRP,IBSPGRP,IOBTP,IAB,IAC,NSPII)
C               CALL DIM_SPII(IATPABS,IBTPABS,IOBTP,1,IAORC,NSPII)
C               MX_NSPII = MAX(MX_NSPII,NSPII)
                MX_NSPII = 0
              END IF
                
              IF(KATP.GT.0) THEN
                MXKA = 0         
                DO 310 KSM = 1, NSMST   
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE 
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) 
     &          MXKA= MXPKA
                MXSOB = 0
                NSOB_AS = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                  NSOB_AS = NSOB_AS + NTSOB(IOBTP,ISMOB) 
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                MXADKBLK_AS = MAX(MXADKBLK_AS,NSOB_AS*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                LCJBLK_ALLSYM = NSOB_AS*MXKA*ITOTB
*. June27, 2003 : Trying to get a multisymmetry version working, 
*. that only has one free symmetry index to 
                LCJBLK_ALLSYM = NSOB_AS*MXKA*MXB
*
                IF(LCJBLK.GT.MXCJ) THEN
                  MXCJ = LCJBLK
                  IATP_MX = IATP
                  IBTP_MX = IBTP
                  KATP_MX = KATP
                  IOBTP_MX = IOBTP
                END IF
                MXCJ_ALLSYM = MAX(MXCJ_ALLSYM,LCJBLK_ALLSYM)
                MXKAB = MAX(MXKAB,MXKA)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
      END DO
*     ^ End of anni/crea map
*
* matrix C(j,Ia,Kb)
*
      DO IAORC = 1, 2
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          
          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100) 
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            ITOTA = 0
            DO ISM = 1, NSMST
              MXA =MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
              ITOTA = ITOTA + NSTFSMSPGP(ISM,IATPABS)
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) GOTO 2812
              CALL NEWTYP(IBTPABS,IAORC,IOBTP,1,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN 
                IF(IAORC.EQ.1.AND.IADVICE.EQ.1.AND.
     &          NHLFSPGP(IATPABS)+NHLFSPGP(KBTP).LT.MNHL.AND.
     &          NHLFSPGP(IBTPABS).GT.NHLFSPGP(IATPABS)+1) THEN
C                 WRITE(6,*) ' N-1 hole space eliminated '
C                 WRITE(6,*) ' IOBTP,IATPABS,KBTP',
C    &            IOBTP,IATPABS,KBTP
                  KBTP = 0
                END IF
              END IF
              IF(KBTP.GT.0) THEN 
C               CALL DIM_SPII(IATPABS,IBTPABS,IOBTP,2,IAORC,NSPII)
C               MX_NSPII = MAX(MX_NSPII,NSPII)
                MX_NSPII = 0
              END IF
*
              IF(KBTP.GT.0) THEN
                MXKB = 0         
                DO KSM = 1, NSMST   
                  MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                END DO   
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA) 
     &          MXKB= MXPKA
                MXSOB = 0
                NSOB_AS = 0
                DO ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                  NSOB_AS = NSOB_AS + NTSOB(IOBTP,ISMOB)
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
                MXADKBLK_AS = MAX(MXADKBLK_AS,NSOB_AS*MXKBO)
                LCJBLK = MXSOB*MXKB*MXA
                LCJBLK_ALLSYM = NSOB_AS*MXKB*ITOTA
*. June27, 2003 : Trying to get a multisymmetry version working, 
*. that only has one free symmetry index to 
                LCJBLK_ALLSYM = NSOB_AS*MXKB*MXA
                MXCJ = MAX(MXCJ,LCJBLK)
                MXCJ_ALLSYM = MAX(MXCJ_ALLSYM,LCJBLK_ALLSYM)
                MXKAB = MAX(MXKAB,MXKB)
*
              END IF
 2812       CONTINUE
            END DO   
          END IF
        END DO   
      END DO   
      END DO
*     ^ End of loop over creation/annihilation
      IF(NTEST.GT.100) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
        WRITE(6,*) ' MXCJ_ALLSYM = ', MXCJ_ALLSYM
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      IF(NTEST.GE.100) WRITE(6,*) ' Dim for C(ij,Ka,Ib) to be constr.'
      MXCIJA = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1 
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100) THEN
             WRITE(6,*) ' TEST: IOCTPA, IOCTPB ', IOCTPA, IOCTPB
             WRITE(6,*) ' IATP, IBTP = ', IATP, IBTP
             WRITE(6,*) ' IATPABS, IBTPABS = ', IATPABS, IBTPABS
            END IF
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO IAORC = 1, 2
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,IAORC,IOBTP,1,K1ATP)
*. No N+1 mappings for particle spaces
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) K1ATP = 0
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IAORC,IOBTP K1ATP ',IAORC,IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JAORC = 1, 2
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1ATP,JAORC,JOBTP,1,KATP)
                  IF(JAORC.EQ.2.AND.IPHGAS(JOBTP).EQ.1) KATP = 0
                  IF(NTEST.GE.100)
     &            WRITE(6,*) ' JAORC,JOBTP KATP ',IAORC,IOBTP,KATP
                  IF(KATP.GT.0) THEN
                    MXKA = 0         
                    DO KSM = 1, NSMST   
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) 
     &              WRITE(6,*) ' KATP, MXKA = ', KATP ,MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) 
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    IF(NTEST.GE.100) 
     &              WRITE(6,*) ' MXISOB, MXJSOB, MXKA, MXIB = ',
     &                           MXISOB, MXJSOB, MXKA, MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
                END DO
*               ^ End of loop over JOBTP, JAORC
              END IF
            END DO
            END DO
*           ^ End of loop over IOBTP, IAORC
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' Dim for C(ij,Ka,Ib) constructed.'
*
*
* matrix C(ij,Ia,kb)
* both Ka and Ib blocked
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1 
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO IAORC = 1, 2
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP(IBTPABS,IAORC,IOBTP,1,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1)K1BTP = 0
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JAORC = 1, 2
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1BTP,JAORC,JOBTP,1,KBTP)
                  IF(JAORC.EQ.2.AND.IPHGAS(JOBTP).EQ.1) KBTP = 0
                  IF(KBTP.GT.0) THEN
                    MXKB = 0         
                    DO KSM = 1, NSMST   
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA) 
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
                END DO
*               ^ End of loop over JOBTP,JAORC
              END IF
            END DO
            END DO
*           ^ End of loop over IOBTP,IAORC
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.10) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,kb)
* both Ka and Kb blocked
*
*. Modified : Only used if atmost two elecs in i and j
*.            No batching 
*.            Used for hardwired few electron code
      MXCIJAB = 0
      MXKACTEL = 1
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1 
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
C           NEL1234(JOBTP,IATPABS)
              IF(NEL1234(IOBTP,KATP).GT.MXKACTEL) KATP = 0
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
*. No partitioning
C               IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from IBTP
                  CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
                  IF(NEL1234(JOBTP,KBTP).GT.MXKACTEL) KBTP = 0
                  IF(KBTP.GT.0) THEN
                    MXKB = 0         
                    DO KSM = 1, NSMST   
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
*. No partitioning
C                   IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA) MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*. Huge and not used pt so 
      MXCIJAB = 0
*
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is 
* obtained as single annihilations
      MXSXBL = 0
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??                                             
          CALL NEWTYP(IATPABS,1,JOBTP,1,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron 
            DO IOBTP = 1, NTPOB
*  Allowed ? 
              CALL NEWTYP(KATP,2,IOBTP,1,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??                                             
          CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron 
            DO IOBTP = 1, NTPOB
*  Allowed ? 
              CALL NEWTYP(KBTP,2,IOBTP,1,JBTP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.10) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL 
        WRITE(6,*) ' MXRESC_PH : MXKAB = ', MXKAB
        WRITE(6,*) ' Info on largest C(Ka,j,Jb) block'
        WRITE(6,*) ' IATP_MX, IBTP_MX, KATP_MX, IOBTP_MX ',
     &               IATP_MX, IBTP_MX, KATP_MX, IOBTP_MX 
        WRITE(6,*) ' MX_NSPII = ', MX_NSPII
      END IF
*
      RETURN
      END
      SUBROUTINE MXRESC(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                  NSMST,NSTFSMSPGP,
     &                  MXPNSMST,
     &                  NSMOB,MXPTOB,NTPOB,NTSOB,NTESTG,MXPKA,
     &                  NEL1234,
     &                  MXCJ,MXCIJA,
     &                  MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
*
* Largest block of single excitations MXSXBL

*. Input 
* IAB :allowed combination of alpha and beta supergroups 
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXPNSMST,*)
      DIMENSION NTSOB(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
* 
      NTESTL = 1000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) WRITE(6,*) ' MXRESC : MXPKA ', MXPKA
*
* matrix C(j,Ka,Ib)
*
*. Note : Only done for alpha-strings, problems when transposing
*         constructing C(J,Ia,Kb)
      MXCJ = 0
      MXADKBLK = 0
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          
          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100) 
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            DO 210 ISM = 1, NSMST
              MXB =MAX(MXB,NSTFSMSPGP(ISM,IBTPABS))
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB = ', MXB
            DO 300 IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0         
                DO 310 KSM = 1, NSMST   
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE 
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) 
     &          MXKA= MXPKA
                MXSOB = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
*
* matrix C(j,Ia,Kb)
*
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          
          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100) 
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            DO ISM = 1, NSMST
              MXA =MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IBTPABS,1,IOBTP,1,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN
                MXKB = 0         
                DO KSM = 1, NSMST   
                  MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                END DO   
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) 
     &          MXKB= MXPKA
                MXSOB = 0
                DO ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
                LCJBLK = MXSOB*MXKB*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
            END DO   
          END IF
        END DO   
      END DO   
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1 
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,K1ATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1ATP,1,JOBTP,1,KATP)
                  IF(KATP.GT.0) THEN
                    MXKA = 0         
                    DO KSM = 1, NSMST   
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) 
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
*
*
* matrix C(ij,Ia,kb)
* both Ka and Ib blocked
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1 
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP(IBTPABS,1,IOBTP,1,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1BTP,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0         
                    DO KSM = 1, NSMST   
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA) 
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,kb)
* both Ka and Kb blocked
*
      MXCIJAB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1 
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from IBTP
                  CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0         
                    DO KSM = 1, NSMST   
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA) 
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is 
* obtained as single annihilations
      MXSXBL = 0
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??                                             
          CALL NEWTYP(IATPABS,1,JOBTP,1,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron 
            DO IOBTP = 1, NTPOB
*  Allowed ? 
              CALL NEWTYP(KATP,2,IOBTP,1,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??                                             
          CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron 
            DO IOBTP = 1, NTPOB
*  Allowed ? 
              CALL NEWTYP(KBTP,2,IOBTP,1,JBTP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL 
      END IF
*
      RETURN
      END
      SUBROUTINE NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
*
* an input  supergroup is given.                    
* apply an string of elementary operators to this supergroup and
* obtain supergroup mumber of new string
*
* Jeppe Olsen, October 1993
* GAS-version : July 95 
*     
* ------
* Input
* ------
*     
* INSPGP  : input super group ( given as absolute supergroup number)
* IACOP(I) = 1 : operator I is an annihilation operator
*          = 2 : operator I is a  creation   operator
* ITPOP(I) : orbitals space of operator I
* NOP : Number of operators
*
* Output
* ------
* OUTSPGP  : supergroup of resulting string
*
*
c      IMPLICIT REAL*8(A-H,O-Z) 
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER ITPOP(*),IACOP(*)
*. Number of active spaces  (NGAS )
      INCLUDE 'cgas.inc'
      INCLUDE 'strbas.inc'
*. Local scratch
      DIMENSION IEL(MXPNGAS)
*. output
      INTEGER OUTSPGP
*
      INEW_OR_OLD = 1
      IF(INEW_OR_OLD.EQ.1) THEN
        CALL NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &       NGAS,int_mb(KSPGPAN),int_mb(KSPGPCR),OUTSPGP)
C     NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
      ELSE IF(INEW_OR_OLD.EQ.2) THEN
*. Number of electrons in AS1,AS2, .... for input supergroup
        CALL  GTSPGP(IEL,INSPGP,2)
*
        IDELTA = 0
        DO IOP = 1, NOP
*. Change in number of orbitals
          IF(IACOP(IOP).EQ.1) THEN
            IDELTA = IDELTA - 1
            IDELTAI = -1
          ELSE 
            IDELTA = IDELTA + 1
            IDELTAI = +1
          END IF
          IEL(ITPOP(IOP)) = IEL(ITPOP(IOP)) + IDELTAI
        END DO
*. output supergroup 
        CALL  GTSPGP(IEL,OUTSPGP,1)
      END IF
*     /\ End  of NEW_OR_OLD switch
*
      NTEST = 0 
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NEWTYP ,  OUTSPGP ', OUTSPGP
      END IF
*
      RETURN
      END
      SUBROUTINE NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &           NGAS,ISPGPAN,ISPGPCR,OUTSPGP)
*
* Strinf of operators X supergroup => new supergroup
*
* Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ISPGPAN(NGAS,*),ISPGPCR(NGAS,*)
*. Specific input
      DIMENSION IACOP(NOP),ITPOP(NOP)
*. Output
      INTEGER OUTSPGP 
*
      OUTSPGP = INSPGP
      DO IOP = 1, NOP
        IF(IACOP(IOP).EQ.1) THEN
          OUTSPGP = ISPGPAN(ITPOP(IOP),OUTSPGP)
        ELSE
          OUTSPGP = ISPGPCR(ITPOP(IOP),OUTSPGP)
        END IF
*
        IF(OUTSPGP.EQ.0) THEN
C         WRITE(6,*) ' NEWTYPS, cul de sac'
C         WRITE(6,*) ' undefined supergroup type '
C         WRITE(6,*) ' String of operator : IAC,ITPOP'
C         CALL IWRTMA(IACOP,1,NOP,1,NOP)
C         CALL IWRTMA(ITPOP,1,NOP,1,NOP)
          GOTO 1001
        END IF
      END DO
*
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NEWTYPS : IAC and ITPOP'
        CALL IWRTMA(IACOP,1,NOP,1,NOP)
        CALL IWRTMA(ITPOP,1,NOP,1,NOP)
        WRITE(6,*) ' Input and Output group ',INSPGP,OUTSPGP
      END IF
*
      RETURN
      END
C      CALL GTSTTP(KCLS,KEL1,KEL3,KTYPE)
      SUBROUTINE GTSPGP(IEL,ISPGP,IWAY)
*
*
* Relation between number of electrons in AS1, AS2 ... and 
* supergoup number   
*
* IWAY = 1 : 
* Get ISPGP : Supergroup of strings 
*             with IEL(*)  electrons in each AS
* IWAY = 2 :
* GET IEL(*)  : Number of electrons in each AS for supergroup ISPGP 
*
*
* Jeppe Olsen, Another lonely night in Lund
*               GAS version July 1995
*
c      IMPLICIT REAL*8 (A-H,O-Z)
*. Generel input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc' !dma: wrkspc not used.
      INCLUDE 'strbas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
C     COMMON/GASSTR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),
C    &              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),
C    &              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS),
C    &              NELFTP(MXPSTT),NSPGPFTP(MXPSTT),IBSPGPFTP(MXPSTT),
C    &              ISPGPFTP(MXPNGAS,MXPSTT),NELFSPGP(MXPNGAS,MXPSTT),
C    &              NGRP,NSTTP,MXNSTR,NTSPGP

C     COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
C    &            NGSOB(MXPOBS,MXPNGAS),
C    &            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,
C    &            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC
*. input(IWAY = 2 ), output (IWAY = 1 )
      INTEGER IEL(*)
*
      IF(IWAY.EQ.1) THEN
*. Occupation => Number 
        ISPGP = -1
        DO JSPGP = 1, NTSPGP
          IF(ISPGP.EQ.-1) THEN
            IEQUAL = 1
            DO IGAS = 1, NGAS  
              IF(NELFSPGP(IGAS,JSPGP).NE.IEL(IGAS))  IEQUAL= 0 
            END DO
            IF(IEQUAL.EQ.1) ISPGP = JSPGP
          END IF
        END DO
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Number => Occupation
        DO IGAS = 1, NGAS
         IEL(IGAS) = NELFSPGP(IGAS,ISPGP)
        END DO
      END IF
*
      NTEST  = 000
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from GTSPGP '
        WRITE(6,*) 
     &   ' IWAY ISPGP IEL ', IWAY,ISPGP,(IEL(IGAS),IGAS = 1, NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE RSBB2B(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                  JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,
     &                  IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,
     &                  ADSXA,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  ICJKAIB,CJRES,SIRES,S2,SCLFAC,NTEST,
     &                  NSEL2E,ISEL2E)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ICJKAIB =1 =>  construct C(Ka,Jb,j) and S(Ka,Ib,i) as intermediate 
*                 matrices in order to reduce overhead
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + ICJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks 
*         
*. October 96 : New routines for accessing annihilation information
*.              Cleaned and shaved, only IROUTE = 3 option active
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
*
      INTEGER ISEL2E(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION S2(*)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
*
*
      CALL LUCIAQENTER('RS2B ')
C?    NTEST = 0
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Warning : RSBB2B, SCLFAC .ne 1 '
C?    END IF
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' RSBB2B speaking '
        WRITE(6,*) ' =============== '
      END IF
      IROUTE = 3
* IROUTE = 1 : Normal (i.e. old) route,
* IROUTE = 2 : New route with j first
* IROUTE = 3 : C(Ka,j,Jb)
*
*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC)
      CALL SXTYP_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC)           
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
*
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all K strings
*. a+j!ka> = +/-/0 * !Ja>
          CALL ADSTN_GAS(JSM,JTYP,JATP,JASM,IAGRP,
     &                   I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFAC)
*. a+i!ka> = +/-/0 * !Ia>
          ONE = 1.0D0
          CALL ADSTN_GAS(ISM,ITYP,IATP,IASM,IAGRP,
     &                   I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE)
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
C             COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
              CALL COMPRS2LST(I1,XI1S,NJ,I3,XI3S,NI,NKASTR,NKAEFF)
          ELSE 
              NKAEFF = NKASTR
          END IF
            
*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK   
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
C         IF(NKABTC.GT.1) THEN
C           WRITE(6,*) ' RSBB2B : NKABTC .GT. 1 '
C           WRITE(6,*) ' I am not prepared for this '
C           STOP 'RSBB2B :  NKABTC .GT. 1'
C         END IF
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, NJ
              CALL GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*
            ZERO = 0.0D0
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC*NI)
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*. Should this group of excitations be included 
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
                KOFF = IOBPTS(KTYP,KSM)
                LOFF = IOBPTS(LTYP,LSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0  
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
                CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (j i k l )
*
                IXCHNG = 0
                ICOUL = 1
                ONE = 1.0D0
                CALL LUCIAGETINT(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,
     &                     LTYP,LSM,IXCHNG,0,0,ICOUL,ONE,ONE)
*
* S(Ka,j,Ib) = sum(k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                CALL SKICKJ(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,NI,NJ,NK,NL,NKBSTR,
     &               I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE )
                FACS = 1.0D0
*
 1930         CONTINUE
 2000       CONTINUE
*
*. Scatter out from s(Ka,Ib,i)
*
            DO II = 1, NI
              CALL ADD_SKAIIB(SB,NI,NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*. End of loop over partitioning of alpha strings
 1940   CONTINUE
 2001 CONTINUE
*
 9999 CONTINUE
*
*
      CALL LUCIAQEXIT('RS2B ')
      RETURN
      END
      SUBROUTINE SKICKJ(SKII,CKJJ,NKA,NIB,NJB,NKB,XIJKL,
     &                  NI,NJ,NK,NL,MAXK,
     &                  KBIB,XKBIB,KBJB,XKBJB,IKORD,
     &                  FACS,IROUTE)
*
*
* Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i) 
*          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
*
*
*
* Jeppe Olsen, Spring of 94
*
* : Note : Route 1 has retired, March 97
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKJJ(NKA,*) 
      DIMENSION XIJKL(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*. Input and output
      DIMENSION SKII(NKA,*)
*. Scratch
      INCLUDE 'mxpdim.inc'
      DIMENSION XIJILS(MXPTSOB)
*. To get rid of annoying and incorrect compiler warnings
      JKINTOF = 0
      IKINTOF = 0
*
      IF(NI.GT.MXPTSOB.OR.NJ.GT.MXPTSOB.OR.NK.GT.MXPTSOB
     &   .OR.NL.GT.MXPTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : > MXPTSOB '
         WRITE(6,*) ' N, MXPTSOB ',MAX(NI,NJ,NK,NL),MXPTSOB
         STOP ' Redim MXPTSOB '
      END IF
*
      NTEST =  000
*
      CALL LUCIAQENTER('SKICK')
      IF(IROUTE.EQ.3) THEN
* S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
        DO KB = 1, NKB
*. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
*
         IF(KK.NE.0.AND.LL.NE.0) THEN
           DO K = 1, NK
             IB = KBIB(KB,K)
             IF(IB.NE.0) THEN
               SGNK = XKBIB(KB,K)
               DO L = 1, NL
                 JB = KBJB(KB,L)
                 IF(NTEST.GE.100) 
     &           WRITE(6,'(a,5I4)') ' KB,K,L,IB,JB',KB,K,L,IB,JB
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
*. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
                   IMAX = NI
*
                   IF(IKORD.NE.0) THEN
*. Restrict so (ij) .le. (kl)
                     IMAX  = K
                     JKINTOF = INTOF + (K-1)*NJ
C                    CALL COPVEC(XIJKL(JKINTOF),XIJILS,NJ)
                     DO J = L,NL
                       XIJILS(J) = XIJKL(JKINTOF-1+J)  
                     END DO
                     XIJKL(JKINTOF-1+L) = 0.5D0*XIJKL(JKINTOF-1+L)
                     DO J = L+1, NL
                      XIJKL(JKINTOF-1+J) = 0.0D0
                     END DO
                   END IF
C                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF,1),  CKJJ(ICOFF,1),
     &                         XIJKL(INTOF),NKA,IMAX,NKA,NJ,
     &                         NJ,IMAX,FACS,FACTOR ,0)
                   IF(IKORD.NE.0) THEN
                      DO J = L,NL
                        XIJKL(JKINTOF-1+J) =  XIJILS(J) 
                      END DO
C                    CALL COPVEC(XIJILS,XIJKL(JKINTOF),NJ)
                   END IF
*
                 END IF
               END DO
*
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
      ELSE IF(IROUTE.EQ.2) THEN
* S(I,Ka,Ib) = S(I,Ka,Ib) + sum(j) (ij!kl) C(j,Ka,Jb)
        DO KB = 1, NKB
*. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
*
         IF(KK.NE.0.AND.LL.NE.0) THEN
           DO K = 1, NK
             IB = KBIB(KB,K)
             IF(IB.NE.0) THEN
               SGNK = XKBIB(KB,K)
               DO L = 1, NL
                 JB = KBJB(KB,L)
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
*. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
*
                   JMAX = NJ
                   IF(IKORD.NE.0) THEN
*. Restrict so (ji) .le. (kl)
                     JMAX  = K
                     IKINTOF = INTOF + (K-1)*NI
                     CALL COPVEC(XIJKL(IKINTOF),XIJILS,NI)
                     XIJKL(IKINTOF-1+L) = 0.5D0*XIJKL(IKINTOF-1+L)
                     DO I = L+1, NL
                      XIJKL(IKINTOF-1+I) = 0.0D0
                     END DO
                   END IF
*
C                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF,1), XIJKL(INTOF),
     &                         CKJJ(ICOFF,1),NI,NKA,NI,NJ,
     &                         NJ,NKA,FACS,FACTOR,0)
*
                 IF(IKORD.NE.0) THEN
                   CALL COPVEC(XIJILS,XIJKL(IKINTOF),NI)
                 END IF
*
                 END IF
               END DO
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
*

      ELSE IF (IROUTE.EQ.1) THEN
         WRITE(6,*) ' Sorry route 1 has retired, March 1997'
         STOP'SKICKJ:Invalid route=1'
C     DO 1000 KB = 1, NKB
*. Number of nonvanishing a+lb !Kb>
C       LL = 0
C       DO L = 1, NL
C         IF(KBJB(KB,L).NE.0) LL = LL + 1
C       END DO
*
C       IKEFF = 0
C       DO 900 K = 1, NK
C         IB = KBIB(KB,K)
C         IF(IB.EQ.0) GOTO 900
C         SGNK = XKBIB(KB,K)
*
C         IF(IKORD.EQ.0) THEN
C            LI = NI
C            IMIN = 1
C         ELSE
C            LI = NI-K+1
C            IMIN = K
C         END IF
*
C         DO 700 I = IMIN, NI
C           IKEFF = IKEFF + 1
C           IOFF = (IKEFF-1)*NJ*LL
*. Offset for S(1,IB,i)
C           IBOFF(IKEFF)  = (I-1)*NIB+IB
C           LEFF = 0
C           DO 800 L = 1, NL  
C             JB = KBJB(KB,L)
C             IF(JB.EQ.0) GOTO 800
C             LEFF = LEFF + 1
C             SGNL = XKBJB(KB,L)
C             IF(IKORD.EQ.1.AND.I.EQ.K)THEN
C                FACTOR = 0.5D0*SGNK*SGNL
C             ELSE
C                FACTOR =       SGNK*SGNL
C             END IF
C             JL0 = (LEFF-1)*NJ
C             JLIK0 = (K-1)*NJ*NL*NI
C    &              + (I-1)*NJ*NL
C    &              + (L-1)*NJ
C             DO 600 J = 1, NJ
C               JL = JL0 + J
*. Offsets for C(1,JB,j)
C               JBOFF(JL) = (J-1)*NJB + JB
*. integral * signs in SCR(jl,ik)
*. Integrals are stored as (j l i k )
C               SCR((IKEFF-1)*NJ*LL+JL) = FACTOR*XIJKL(JLIK)
C               SCR(IOFF+JL) = FACTOR*XIJKL(JLIK0+J)
C 600         CONTINUE
C 800       CONTINUE
C 700     CONTINUE
C 900   CONTINUE
*
C       CALL GSAXPY(SKII,CKJJ,SCR,IKEFF,NJ*LL,NKA,IBOFF,JBOFF)
C1000 CONTINUE
      END IF
*. End of IROUTE branchning 
*
      CALL LUCIAQEXIT('SKICK')
      RETURN
      END
      SUBROUTINE GSAXPY(AB,A,B,NABCOL,NACOL,NROW,IABCOL,IACOL)
*
* AB(I,IABCOL(J)) = AB(I,IABCOL(J)) + A(I,IACOL(K))*B(K,J)
*
*
* Jeppe Olsen, Spring of 94 Daughter of MSAXPY*
*
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*. Input
      DIMENSION A(NROW,NACOL),B(NACOL,NABCOL)
      DIMENSION IACOL(NACOL),IABCOL(NABCOL)
*. Output
      DIMENSION AB(NROW,NABCOL)
*
CT    CALL LUCIAQENTER('GSAXP')
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*. Straightforward sequence of SAXPY's
         DO 1100 J = 1, NABCOL
           DO 1200 K = 1, NACOL
             JACT = IABCOL(J)
             KACT = IACOL(K)
             FACTOR = B(K,J)
             DO 1000 I = 1, NROW
               AB(I,JACT) = AB(I,JACT) + FACTOR*A(I,KACT)
 1000        CONTINUE
 1200      CONTINUE
 1100    CONTINUE
*
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Unrolling over columns of A
       NROL = 5
       NRES = MOD(NACOL,NROL)
*. overhead
       IF(NRES.EQ.1) THEN
         DO 2201 J = 1, NABCOL
          JACT = IABCOL(J)
          K1ACT = IACOL(1)
          B1J = B(1,J)
          DO 2101 I = 1, NROW
            AB(I,JACT) = AB(I,JACT) + A(I,K1ACT)*B1J
 2101     CONTINUE
 2201    CONTINUE
       ELSE IF (NRES .EQ. 2 ) THEN
         DO 2202 J = 1, NABCOL
           K1ACT = IACOL(1)
           K2ACT = IACOL(2)
           B1J   = B(1,J)
           B2J   = B(2,J)
           JACT =  IABCOL(J)
           DO 2102 I = 1, NROW
            AB(I,JACT) = AB(I,JACT)
     &    + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
 2102      CONTINUE
 2202    CONTINUE
       ELSE IF (NRES .EQ. 3 ) THEN
         DO 2203 J = 1, NABCOL
           K1ACT = IACOL(1)
           K2ACT = IACOL(2)
           K3ACT = IACOL(3)
           JACT =  IABCOL(J)
           B1J   = B(1,J)
           B2J   = B(2,J)
           B3J   = B(3,J)
           DO 2103 I = 1, NROW
            AB(I,JACT) = AB(I,JACT)
     &    + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
     &    + A(I,K3ACT)*B3J
 2103      CONTINUE
 2203    CONTINUE
       ELSE IF (NRES .EQ. 4 ) THEN
         DO 2204 J = 1, NABCOL
           K1ACT = IACOL(1)
           K2ACT = IACOL(2)
           K3ACT = IACOL(3)
           K4ACT = IACOL(4)
           JACT =  IABCOL(J)
           B1J   = B(1,J)
           B2J   = B(2,J)
           B3J   = B(3,J)
           B4J   = B(4,J)
           DO 2104 I = 1, NROW
            AB(I,JACT) = AB(I,JACT)
     &    + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
     &    + A(I,K3ACT)*B3J + A(I,K4ACT)*B4J
 2104      CONTINUE
 2204    CONTINUE
        END IF
*. ( End of Overhead )
        DO 2305 K = NRES+1,NACOL,NROL
          DO 2205 J = 1, NABCOL
            K1ACT = IACOL(K)
            K2ACT = IACOL(K+1)
            K3ACT = IACOL(K+2)
            K4ACT = IACOL(K+3)
            K5ACT = IACOL(K+4)
            JACT =  IABCOL(J)
            B1J   = B(K,J)
            B2J   = B(K+1,J)
            B3J   = B(K+2,J)
            B4J   = B(K+3,J)
            B5J   = B(K+4,J)
            DO 2105 I = 1, NROW
             AB(I,JACT) = AB(I,JACT)
     &     + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
     &     + A(I,K3ACT)*B3J + A(I,K4ACT)*B4J
     &     + A(I,K5ACT)*B5J
 2105       CONTINUE
 2205     CONTINUE
 2305   CONTINUE
      END IF
*( End of IWAY branching )
*
CT    CALL LUCIAQEXIT('GSAXP')
      RETURN
      END
      

*
  
      SUBROUTINE ZIJPNT(IJPNT,IJSM,IGEJ,NSMOB,IOBSM,LOBSM,ISTOB,NTOTOB)
*
* Construct pointer IJPNT, for accessing h(i,j) 
* in symmetry packed matrix.
* h(i,j) is assumed in symmetry order , so all orbitals of the same 
* symmetry are stored consecutively
*
* In pointer IJPNT IJPNT(i,j) indices i and j are in  type order  
*
*. Input
* IJSM : Symmetry of h(i,j)
* IGEJ : matrix is packed so h(i,j) is included only for i. ge. j 
*        elements IJPNT (i,j) and IJPNT (j,i) are then
*        pointing to the same elements in h(i,j)
* NSMOB : Number of symmetries of orbitals
* LOBSM : Number of orbitals per symmetry
* IOBSM : offset for orbitals of given symmetry
* ISTOB : Symmetry => Type reorder for orbitals
* NTOTOB : Total number of orbitals.
*
* Output
* IJPNT : Pointer
*
* Jeppe Olsen, April 1 1994
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER LOBSM(NSMOB),IOBSM(*),ISTOB(NTOTOB)
*. Output
      INTEGER IJPNT(NTOTOB,NTOTOB)
*
      IJEFF = 0
C?    WRITE(6,*) ' ZIJPNT : NSMOB IJSM',NSMOB,IJSM
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
C?      write(6,*) ' ISM JSM ', ISM,JSM
        IF(JSM.NE.0) THEN
          DO JOB = IOBSM(JSM),IOBSM(JSM) + LOBSM(JSM)-1
            DO IOB = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
              IOBP = ISTOB(IOB)
              JOBP = ISTOB(JOB)
C?            write(6,*) ' JOB IOB JOBP IOPB',
C?   &        JOB,IOB,JOBP,IOBP
              IF(IGEJ.EQ.0) THEN
                IJEFF = IJEFF + 1
                IJPNT(IOBP,JOBP) = IJEFF
              ELSE IF(IGEJ.NE.0 .AND. IOB.GE.JOB) THEN
                IJEFF = IJEFF + 1
                IJPNT(IOBP,JOBP) = IJEFF
                IJPNT(JOBP,IOBP) = IJEFF
              END IF
            END DO
          END DO
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
        WRITE(6,*) 
     &  ' ZIJPNT delivering pointer for IJSM = ', IJSM
        WRITE(6,*)
        WRITE(6,*) ' Pointer matrix as delivered '
        CALL IWRTMA(IJPNT,NTOTOB,NTOTOB,NTOTOB,NTOTOB)
      END IF
*
      RETURN
      END
      SUBROUTINE IJKK(XIJKK,IGEJ,
     &           IJPNT,IJSM,NSMOB,LOBSM,IOBSM,ISTOB,NTOTOB)
*
* Construct matrix XIJKK(K,IJ) = (IJ!KK) - (IK!KJ)
* for pair of orbitals ij  with total symmetry IJSM
* stored in symmetry ordered form
* If IgeJ .ne. 0 , I. ge. j ordering is used
*
* Jeppe Olsen, April 1, 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION LOBSM(NSMOB),IOBSM(NSMOB)
      DIMENSION ISTOB(*)
*.Output
      DIMENSION XIJKK(NTOTOB,*)
*
      IJ = 0   
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        IF(JSM.NE.0) THEN
          DO J = IOBSM(JSM),IOBSM(JSM) + LOBSM(JSM)-1
            DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
              IP = ISTOB(I)
              JP = ISTOB(J)
              IF(IGEJ.EQ.0.OR.(IGEJ.NE.0.AND.I.GE.J)) THEN
                IJ= IJ + 1
                DO K = 1, NTOTOB
                  XIJKK(K,IJ) =
     &            GTIJKL(IP,JP,K,K)-GTIJKL(IP,K,K,JP)
                END DO
              END IF
            END DO
          END DO
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' XIJKK fresk from GTIJKK'
        IF(IGEJ.EQ.0) THEN
          NCOL = NTOTOB ** 2
        ELSE 
          NCOL = NTOTOB*(NTOTOB+1)/2
        END IF
        CALL WRTMAT(XIJKK,NTOTOB,NCOL,NTOTOB,NCOL)
      END IF
* 
      RETURN
      END
      SUBROUTINE MSAXPYN6(AX,A,X,NDIM,NVEC,INDEX)
*
* AX(I) = AX(I) + SUM(L=1,NVEC) A(L)*X(I,INDEX(L))
*
* Initial RISC/6000 version
*
* Jeppe Olsen
*
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AX(*),X(NDIM,*)
      DIMENSION A(*) ,INDEX(*)
*
* level of unrolling of column loop
      NROL = 5
      NRES = MOD(NVEC,NROL)
*
*. part that cannot be column unrolled
*
      IF(NRES.EQ.0.OR.NVEC.EQ.0) THEN
*. Nothing to do !
      ELSE IF (NRES.EQ.1) THEN
        I1 = INDEX(1)
        A1 = A(1)
        DO 20 I = 1, NDIM
          AX(I) = AX(I) + A1*X(I,I1)
   20   CONTINUE
      ELSE IF (NRES.EQ.2) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        DO 30 I = 1, NDIM
         AX(I) = AX(I) + A1*X(I,I1)+A2*X(I,I2)
   30   CONTINUE
      ELSE IF (NRES.EQ.3) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        I3 = INDEX(3)
        A3 = A(3)
        DO 40 I = 1, NDIM
         AX(I) = 
     &   AX(I) + A1*X(I,I1)+A2*X(I,I2)+A3*X(I,I3)
   40   CONTINUE
      ELSE IF (NRES.EQ.4) THEN
        I1 = INDEX(1)
        A1 = A(1)
        I2 = INDEX(2)
        A2 = A(2)
        I3 = INDEX(3)
        A3 = A(3)
        I4 = INDEX(4)
        A4 = A(4)
        DO 50 I = 1, NDIM
         AX(I) = 
     &   AX(I) + A1*X(I,I1)+A2*X(I,I2)+A3*X(I,I3)+A4*X(I,I4)
   50   CONTINUE
      ELSE IF (NRES.GE.5) THEN
        WRITE(6,*) ' NRES to large in MSAXPY so STOP '
        STOP ' MSAXPY '
      END IF
*. unrolled column loops
      DO 285 L = 1+NRES,NVEC,NROL
        AL= A(L)
        IL= INDEX(L)
        AL1= A(L+1)
        IL1 = INDEX(L+1)
        AL2= A(L+2)
        IL2 = INDEX(L+2)
        AL3= A(L+3)
        IL3 = INDEX(L+3)
        AL4= A(L+4)
        IL4 = INDEX(L+4)
        DO 205 I = 1, NDIM
          AX(I) = AX(I) + AL*X(I,IL)+AL1*X(I,IL1)
     &          + AL2*X(I,IL2)+AL3*X(I,IL3)+AL4*X(I,IL4)
  205   CONTINUE
  285 CONTINUE
*
      RETURN
      END
  
      SUBROUTINE SCAMT2(MATIN,MATUT,NROWIN,NROWUT,NCOL, 
     &                  LRSCA,IRSCA,SRSCA,
     &                  LCSCA,ICSCA,SCSCA)
*
*
* Scatter-add  rows/columns of a matrix
* 
* If LRSCA .ne.0 row    scattering is performed,
* If LCSCA .ne.0 column scattering is performed,
*      
*
* Note : NCOL is the number of columns in input matrix
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MATIN,MATUT
*.Input
      DIMENSION IRSCA(*),SRSCA(*),MATIN(NROWIN,*)
      DIMENSION ICSCA(*),SCSCA(*)
*.Output
      DIMENSION MATUT(NROWUT,*)
*
      IWAY = 3
      IF(LRSCA.NE.0) THEN
        IF(LCSCA.EQ.0) THEN
*. Row scattering, no column scattering
          IF(IWAY.EQ.1) THEN
          DO 200 ICOL = 1, NCOL
            DO 100 IROW = 1, NROWIN
              IF(IRSCA(IROW).NE.0)  THEN 
                MATUT(IRSCA(IROW),ICOL) = 
     &          MATUT(IRSCA(IROW),ICOL) + MATIN(IROW,ICOL)*SRSCA(IROW)
              END IF
  100       CONTINUE
  200     CONTINUE
          ELSE IF (IWAY .EQ. 2 ) THEN
            DO 300 IROW = 1, NROWIN
              IIROW = IRSCA(IROW)
              IF(IIROW.NE.0) THEN
                FACTOR = SRSCA(IROW)
                DO 400 ICOL = 1, NCOL
                  MATUT(IIROW,ICOL) = MATUT(IIROW,ICOL) 
     &          + FACTOR*MATIN(IROW,ICOL)
  400           CONTINUE
              END IF
  300       CONTINUE
          ELSE IF (IWAY .EQ. 3 ) THEN
*. row index as inner most loop, unrolled over columns
            NROL = 5
            NRES = MOD(NCOL,NROL)
            IF(NRES.EQ.1) THEN
              DO 1101 IROW = 1, NROWIN
                IF(IRSCA(IROW).NE.0)  THEN 
                  MATUT(IRSCA(IROW),1) = 
     &            MATUT(IRSCA(IROW),1) + MATIN(IROW,1)*SRSCA(IROW)
                END IF
 1101         CONTINUE
            ELSE IF (NRES.EQ.2) THEN
              DO 1102 IROW = 1, NROWIN
                IR = IRSCA(IROW)
                IF(IR.NE.0) THEN
                  SR = SRSCA(IROW)
                  MATUT(IR,1) =  MATUT(IR,1) + MATIN(IROW,1)*SR
                  MATUT(IR,2) =  MATUT(IR,2) + MATIN(IROW,2)*SR
                END IF
 1102         CONTINUE
            ELSE IF (NRES.EQ.3) THEN
              DO 1103 IROW = 1, NROWIN
                IR = IRSCA(IROW)
                IF(IR.NE.0) THEN
                  SR = SRSCA(IROW)
                  MATUT(IR,1) =  MATUT(IR,1) + MATIN(IROW,1)*SR
                  MATUT(IR,2) =  MATUT(IR,2) + MATIN(IROW,2)*SR
                  MATUT(IR,3) =  MATUT(IR,3) + MATIN(IROW,3)*SR
                END IF
 1103         CONTINUE
            ELSE IF (NRES.EQ.4) THEN
              DO 1104 IROW = 1, NROWIN
                IR = IRSCA(IROW)
                IF(IR.NE.0) THEN
                  SR = SRSCA(IROW)
                  MATUT(IR,1) =  MATUT(IR,1) + MATIN(IROW,1)*SR
                  MATUT(IR,2) =  MATUT(IR,2) + MATIN(IROW,2)*SR
                  MATUT(IR,3) =  MATUT(IR,3) + MATIN(IROW,3)*SR
                  MATUT(IR,4) =  MATUT(IR,4) + MATIN(IROW,4)*SR
                END IF
 1104         CONTINUE
            END IF
*. Main body
            DO 1200 ICOL = NRES+1, NCOL-NROL+1, NROL
            DO 1100 IROW = 1, NROWIN
              IR = IRSCA(IROW)
              IF(IR.NE.0) THEN
                SR = SRSCA(IROW)
                MATUT(IR,ICOL  ) = 
     &          MATUT(IR,ICOL  ) + MATIN(IROW,ICOL  )*SR
                MATUT(IR,ICOL+1) = 
     &          MATUT(IR,ICOL+1) + MATIN(IROW,ICOL+1)*SR
                MATUT(IR,ICOL+2) = 
     &          MATUT(IR,ICOL+2) + MATIN(IROW,ICOL+2)*SR
                MATUT(IR,ICOL+3) = 
     &          MATUT(IR,ICOL+3) + MATIN(IROW,ICOL+3)*SR
                MATUT(IR,ICOL+4) = 
     &          MATUT(IR,ICOL+4) + MATIN(IROW,ICOL+4)*SR
              END IF
 1100       CONTINUE
 1200       CONTINUE
          END IF
*. (end if IWAY loop )
*. Row scattering, Column scattering
        ELSE IF(LCSCA.NE.0) THEN
          DO 201 ICOL = 1, NCOL
            IICOL = ICSCA(ICOL)
            IF(IICOL.NE.0) THEN
              FACTOR = SCSCA(ICOL)
              DO 101 IROW  = 1, NROWIN
               IF(IRSCA(IROW).NE.0) THEN
                 MATUT(IRSCA(IROW),IICOL) =
     &           MATUT(IRSCA(IROW),IICOL)   
     &       +   FACTOR*SRSCA(IROW)*MATIN(IROW,ICOL)
               END IF
  101         CONTINUE
            END IF
  201     CONTINUE
        END IF
      ELSE IF (LRSCA.EQ.0) THEN
*
* No row scattering
*
        IF(LCSCA.EQ.0) THEN
*. No row scattering, no column scattering
          DO 202 ICOL = 1, NCOL
            DO 102 IROW = 1, NROWIN
              MATUT(IROW,ICOL) = 
     &         MATUT(IROW,ICOL) + MATIN(IROW,ICOL)
  102       CONTINUE
  202     CONTINUE
*.  No row scattering, column scattering
        ELSE IF(LCSCA.NE.0) THEN
          DO 203 ICOL = 1, NCOL
            IICOL = ICSCA(ICOL)
            IF(IICOL.NE.0) THEN
              FACTOR = SCSCA(ICOL)
              DO 103 IROW  = 1, NROWIN
               MATUT(IROW,IICOL) 
     &       = MATUT(IROW,IICOL) + FACTOR*MATIN(IROW,ICOL)
  103         CONTINUE
            END IF
  203     CONTINUE
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE GATMT2(MATIN,MATUT,NROWIN,NROWUT,NCOL,
     &                  LRGAT,IRGAT,SRGAT,
     &                  LCGAT,ICGAT,SCGAT)
*
* Gather rows of a matrix :  
* 
* If LRGAT .ne.0 row    gathering is performed,
* If LRCOL .ne.0 column gathering is performed,
*      
*
* If LCGAT = 0 :
*   MATUT(I,J) = SRGAT(I)*MATIN(IRGAT(I),J) if IRGAT(I).NE.0
*              = 0                          if IRGAT(I) = 0
*
* IF LCGAT .ne. 0 :
* MATUT(I,J) = SRGAT(I)*SRGAT(J)*MATIN(IRGAT(I),ICGAT(J)) IF IRGAT(I)*ICGAT(I)
*                                                           .ne. 0
*              = 0                                        IF IRGAT(I)*ICGAT(I)=0
*
*
* Note : NCOL is the number of columns in output matrix
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MATIN,MATUT
*.Input
      DIMENSION IRGAT(*),SRGAT(*),MATIN(NROWIN,*)
      DIMENSION ICGAT(*),SCGAT(*)
*.Output
      DIMENSION MATUT(NROWUT,*)
*
      IWAY = 3
      IF(LRGAT.NE.0) THEN
        IF(LCGAT.EQ.0) THEN
*. Row gathering, no column gathering 
          IF(IWAY.EQ.1) THEN
            DO 200 ICOL = 1, NCOL
              DO 100 IROW = 1, NROWUT
                IF(IRGAT(IROW).EQ.0)  THEN 
                  MATUT(IROW,ICOL) = 0.0D0
                ELSE
                  MATUT(IROW,ICOL) 
     &          = MATIN(IRGAT(IROW),ICOL)*SRGAT(IROW)
                END IF
  100         CONTINUE
  200       CONTINUE
          ELSE IF (IWAY .EQ. 2 ) THEN
            DO 300 IROW = 1, NROWUT
              IIROW = IRGAT(IROW)
              IF(IIROW.EQ.0) THEN
                DO 500 ICOL = 1, NCOL
                  MATUT(IROW,ICOL) = 0.0D0
  500           CONTINUE
              ELSE 
                FACTOR = SRGAT(IROW)
                DO 400 ICOL = 1, NCOL
                  MATUT(IROW,ICOL) = MATIN(IIROW,ICOL)*FACTOR
  400           CONTINUE
              END IF
  300       CONTINUE
          ELSE IF (IWAY .EQ. 3 ) THEN
*. row index as inner most loop, unrolled over columns
            NROL = 5
            NRES = MOD(NCOL,NROL)
            IF(NRES.EQ.1) THEN
              DO 1101 IROW = 1, NROWUT
                IF(IRGAT(IROW).EQ.0)  THEN 
                  MATUT(IROW,1) = 0.0D0
                ELSE
                  MATUT(IROW,1) = MATIN(IRGAT(IROW),1)*SRGAT(IROW)
                END IF
 1101         CONTINUE
            ELSE IF (NRES.EQ.2) THEN
              DO 1102 IROW = 1, NROWUT
                IR = IRGAT(IROW)
                IF(IR.EQ.0) THEN
                  MATUT(IROW,1) = 0.0D0          
                  MATUT(IROW,2) = 0.0D0          
                ELSE
                  SR = SRGAT(IROW)
                  MATUT(IROW,1) = MATIN(IR,1)*SR
                  MATUT(IROW,2) = MATIN(IR,2)*SR
                END IF
 1102         CONTINUE
            ELSE IF (NRES.EQ.3) THEN
              DO 1103 IROW = 1, NROWUT
                IR = IRGAT(IROW)
                IF(IR.EQ.0) THEN
                  MATUT(IROW,1) = 0.0D0          
                  MATUT(IROW,2) = 0.0D0          
                  MATUT(IROW,3) = 0.0D0          
                ELSE
                  SR = SRGAT(IROW)
                  MATUT(IROW,1) = MATIN(IR,1)*SR
                  MATUT(IROW,2) = MATIN(IR,2)*SR
                  MATUT(IROW,3) = MATIN(IR,3)*SR
                END IF
 1103         CONTINUE
            ELSE IF (NRES.EQ.4) THEN
              DO 1104 IROW = 1, NROWUT
                IR = IRGAT(IROW)
                IF(IR.EQ.0) THEN
                  MATUT(IROW,1) = 0.0D0
                  MATUT(IROW,2) = 0.0D0
                  MATUT(IROW,3) = 0.0D0
                  MATUT(IROW,4) = 0.0D0
                ELSE
                  SR = SRGAT(IROW)
                  MATUT(IROW,1) = MATIN(IR,1)*SR
                  MATUT(IROW,2) = MATIN(IR,2)*SR
                  MATUT(IROW,3) = MATIN(IR,3)*SR
                  MATUT(IROW,4) = MATIN(IR,4)*SR
                END IF
 1104         CONTINUE
            END IF
*. Main body
            DO 1200 ICOL = NRES+1, NCOL-NROL+1, NROL
            DO 1100 IROW = 1, NROWUT
              IR = IRGAT(IROW)
              IF(IR.EQ.0) THEN
                MATUT(IROW,ICOL  ) = 0.0D0                
                MATUT(IROW,ICOL+1) = 0.0D0                
                MATUT(IROW,ICOL+2) = 0.0D0                
                MATUT(IROW,ICOL+3) = 0.0D0                
                MATUT(IROW,ICOL+4) = 0.0D0                
              ELSE
                SR = SRGAT(IROW)
                MATUT(IROW,ICOL  ) = MATIN(IR,ICOL  )*SR
                MATUT(IROW,ICOL+1) = MATIN(IR,ICOL+1)*SR
                MATUT(IROW,ICOL+2) = MATIN(IR,ICOL+2)*SR
                MATUT(IROW,ICOL+3) = MATIN(IR,ICOL+3)*SR
                MATUT(IROW,ICOL+4) = MATIN(IR,ICOL+4)*SR
              END IF
 1100       CONTINUE
 1200       CONTINUE
          END IF
*. (end if IWAY loop )
*. Row gathering, Column gathering
        ELSE IF(LCGAT.NE.0) THEN
          DO 201 ICOL = 1, NCOL
            IICOL = ICGAT(ICOL)
            IF(IICOL.EQ.0) THEN
              DO  90 IROW = 1, NROWUT
               MATUT(IROW,ICOL) = 0.0D0
   90         CONTINUE
            ELSE
              FACTOR = SCGAT(ICOL)
              DO 101 IROW  = 1, NROWUT
               IF(IRGAT(IROW).EQ.0) THEN
                  MATUT(IROW,ICOL) = 0.0D0
               ELSE 
                 MATUT(IROW,ICOL) 
     &       =   FACTOR*SRGAT(IROW)*MATIN(IRGAT(IROW),IICOL)
               END IF
  101         CONTINUE
            END IF
  201     CONTINUE
        END IF
      ELSE IF (LRGAT.EQ.0) THEN
*
* No row gathering
*
        IF(LCGAT.EQ.0) THEN
*. No row gathering, no column gathering 
          DO 202 ICOL = 1, NCOL
            DO 102 IROW = 1, NROWUT
              MATUT(IROW,ICOL) = MATIN(IROW,ICOL)
  102       CONTINUE
  202     CONTINUE
*.  No row gathering, column gathering
        ELSE IF(LCGAT.NE.0) THEN
          DO 203 ICOL = 1, NCOL
            IICOL = ICGAT(ICOL)
            IF(IICOL.EQ.0) THEN
              DO  93 IROW = 1, NROWUT
               MATUT(IROW,ICOL) = 0.0D0
   93         CONTINUE
            ELSE
              FACTOR = SCGAT(ICOL)
              DO 103 IROW  = 1, NROWUT
               MATUT(IROW,ICOL) 
     &       = FACTOR*MATIN(IROW,IICOL)
  103         CONTINUE
            END IF
  203     CONTINUE
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE RSBB1E(ISCSM,ISCTP,ISBCTP,ICCSM,ICCTP,ICBCTP,
     &                  IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,
     &                  MXPNGASX,NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBSX,MOC,
     &                  NSCOL,MXSXST,IEXSTR,FACSTR,IH2TRM,SCLFAC,
     &                  IUSE_PH,IPHGAS,NTESTG)
*
* One electron excitations on column strings
*. If IH2TRM .ne. 0 then the diagonal and one-electron
*  excitations arising from the two body operator is also 
* included
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ISBCTP : Base for sigma column types 
* ICCSM,ICCTP : Symmetry and type of C     columns
* ICBCTP : Base for C     column types 
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active sets 
* ISEL : Occupation in each active set for sigma block 
* ICEL : Occupation in each active set for C     block 
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* MOC  : Use MOC method ( instead of N-1 resolution method )
* 
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for one electron integrals
*
* Jeppe Olsen, Winter of 1991
*              IUSE_PH added winter of 97
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. MAX dimensions 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER IPHGAS(NGAS)
*.Specific Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),H(*)
      DIMENSION I2(*),XI2S(*)
      DIMENSION IEXSTR(MXSXST,*),FACSTR(MXSXST,*)
*.Local arrays ( assume MPNGAS = 16 ) !!! 
      DIMENSION ITP(16),JTP(16)
      DIMENSION ISGRP(16),ICGRP(16)
*. For transposing integral block
      DIMENSION HSCR(MXPTSOB*MXPTSOB)
*
      DIMENSION IJ_REO(2),IJ_DIM(2),IJ_SM(2),IJ_TP(2),IJ_AC(2)
      DIMENSION ISCR(2)
      CALL LUCIAQENTER('RS1   ')
* Type of single excitations that connects the two column strings
C     MOC = 1
      NTESTL = 0000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500)THEN
        WRITE(6,*) 
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Information from RSBB1E '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' RSBB1E : MOC,IH2TRM,IUSE_PH ', MOC,IH2TRM,IUSE_PH
        WRITE(6,*) ' ISEL : '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' ICEL : '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*. Obtain groups
C     GET_SPGP_INF(ISPGP,ITP,IGRP)
      CALL GET_SPGP_INF(ICCTP,IGRP,ICGRP)
      CALL GET_SPGP_INF(ISCTP,IGRP,ISGRP)
*
      IFRST = 1
      JFRST = 1
*. Types of single excitations that connect ISEL and ICEL
      CALL SXTYP2_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL,IPHGAS)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.2000)
     &  write(6,*) ' ITYP JTYP ', ITYP,JTYP
*. Is this combination of types allowed
        IJ_ACT = I_SX_ACT(ITYP,JTYP) 
        IF(IJ_ACT.EQ.0) GOTO 900
*. Hvilken vej skal vi valge, 
        NOP = 2
        IJ_AC(1) = 2
        IJ_AC(2) = 1
        IJ_TP(1) = ITYP
        IJ_TP(2) = JTYP
        IF(IUSE_PH.EQ.1) THEN
          CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TP,IJ_AC,IJ_REO,SIGNIJ)
        ELSE
          IJ_REO(1) = 1
          IJ_REO(2) = 2
          SIGNIJ = 1.0D0
        END IF
*
        ISCR(1) = IJ_AC(1)
        ISCR(2) = IJ_AC(2)
        IJ_AC(1) = ISCR(IJ_REO(1))
        IJ_AC(2) = ISCR(IJ_REO(2))
*
        ISCR(1) = ITYP
        ISCR(2) = JTYP
        IJ_TP(1) = ISCR(IJ_REO(1))
        IJ_TP(2) = ISCR(IJ_REO(2))
*
        DO 800 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
*. New intermediate strings will be accessed so
          KFRST = 1
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.2000)
     &    write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
*. Reorder 
* 
          ISCR(1) = ISM
          ISCR(2) = JSM
          IJ_SM(1) = ISCR(IJ_REO(1))
          IJ_SM(2) = ISCR(IJ_REO(2))
*
          ISCR(1) = NIORB
          ISCR(2) = NJORB
          IJ_DIM(1) = ISCR(IJ_REO(1))
          IJ_DIM(2) = ISCR(IJ_REO(2))
*
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*.Fetch integrals : For CI-transformations using RSBB1E 
*.most of the blocks vanishes 
*.Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
           IF(IJ_REO(1).EQ.1) THEN
*. obtain integrals h(j,i)
             CALL GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
             CALL TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
C?           WRITE(6,*) ' RSBB1E: One-electron integrals h(i,j)(T) '
C?           CALL WRTMAT(H,IJ_DIM(2),IJ_DIM(1),IJ_DIM(2),IJ_DIM(1))
           ELSE
*. Obtain integrals h(i,j)
             CALL GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
           END IF
C          XNORM = INPROD(H,H,IJ_DIM(1)*IJ_DIM(2))
C          IF(XNORM.EQ.0) GOTO 800
          IF(MOC.EQ.0) THEN
*
*
* ======================================================================
*.                   Use N-1 resolution method
* ======================================================================
*
*
*. Obtain annihilation/creation maps for all K strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          SCLFACS = SIGNIJ*SCLFAC
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
          CALL ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC(1))
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          ONE = 1.0D0
          CALL ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,
     &         I2,XI2S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC(1))
*. Compress list to common nonvanishing elements
          IF(NTEST.GE.2000) WRITE(6,*)
     &    ' NKASTR efter ADAST_GAS ', NKASTR
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I2,XI2S,IJ_DIM(1),
     &             NKASTR,NKAEFF)
          ELSE 
              NKAEFF = NKASTR
          END IF
*. Loop over partitionings of the row strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = MIN(KTOP + MAXK,NKAEFF)
              IF(KTOP.EQ.NKAEFF) THEN
                KEND = 1
              ELSE
                KEND = 0
              END IF
              LKABTC = KTOP - KBOT +1
*. This is the place to start over partitioning of I strings
              DO 701 IIPART = 1, NIPART
                IBOT = (IIPART-1)*MAXI+1
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP - IBOT + 1
* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
                DO JJORB = 1,IJ_DIM(2)
                  ICGOFF = 1 + (JJORB-1)*LKABTC*NIBTC
                  CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                 LKABTC,I1(KBOT+(JJORB-1)*NKASTR),
     &                 XI1S(KBOT+(JJORB-1)*NKASTR) )
                END DO
*.Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
C               CALL GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
C               CALL TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
*. Problems when HOLE switches blocks around ?
C               CALL GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' RSBB1E H BLOCK '
                  CALL WRTMAT(H,IJ_DIM(2),IJ_DIM(1),IJ_DIM(2),IJ_DIM(1))
                END IF
*.Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
                NIK = NIBTC*LKABTC
                IF(NTEST.GE.2000) WRITE(6,*) ' NIBTC, LKABTC = ',
     &                                         NIBTC, LKABTC
                FACTORC = 0.0D0
                FACTORAB = 1.0D0
                IF(NTEST.GE.2000) THEN 
                  WRITE(6,*) ' CSCR array,NIK X NJORB array '
                  CALL WRTMAT(CSCR,NIK,IJ_DIM(2),NIK,IJ_DIM(2))
                END IF
                CALL MATML7(SSCR,CSCR,H,
     &               NIK,IJ_DIM(1),NIK,IJ_DIM(2),IJ_DIM(2),IJ_DIM(1),
     &               FACTORC,FACTORAB,0)
                IF(NTEST.GE.2000) THEN 
                  WRITE(6,*) ' SSCR array,NIK X NIORB array '
                  CALL WRTMAT(SSCR,NIK,IJ_DIM(1),NIK,IJ_DIM(1))
                END IF
*.S(I,a+ K) =  S(I, a+ K) + sgn*Sscr(I,K,i)
                DO IIORB = 1,IJ_DIM(1)
                  ISBOFF = 1+(IIORB-1)*LKABTC*NIBTC
                  CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                 LKABTC,I2(KBOT+(IIORB-1)*NKASTR),
     &                 XI2S(KBOT+(IIORB-1)*NKASTR))
                END DO
  701       CONTINUE
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
*. End of loop over I partitioninigs
          END IF
*.(End of algorithm switch)
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
      CALL LUCIAQEXIT('RS1  ')
C!    WRITE(6,*) ' Enforced stop in RSBB1E '
C!    STOP' Enforced stop in RSBB1E '
      RETURN
      END
      SUBROUTINE H1STRO(H,ISTRSM,ISTRTP,ISTRGP,
     &                  IOBSM,IOBTP,JOBSM,JOBTP,MXSXST,
     *                  IEXSTR,FACSTR,NEX,IH2TRM,NTESTG)
*                                                                       
* Obtain one-electron string matrix elements sum(ij) h(ij) <J!e(ij)!I>
* with orbital i and j restricted to given orbital TS subsets
*
* from strings of symmetry ISTRSM and class ISTRTP
*
*
* =====                                                             
* Input                                                             
* =====                                                             
*

*   ISTRSM,ISTTP,ISTRGP :  symmetry, type and group  of input string 
*   IOBSM,IOBTP : symmetry and type of orbital i
*   JOBSM,JOBTP : symmetry and type of orbital j
*   MXSXST      : Max number of single excitations for given string
*   H           : one-electron integrals for subsets of orbitals 
*   IH2TRM      : ne 0 => include diagonal and one-electron 
*                 operators from two body operator 
*                                                                       
* ======                                                            
* Output                                                            
* ======                                                            
*
*     NEX(ISTR)        : Number of excitations
*     IEXSTR(IEX,ISTR) : Number of excited string ( relative to offset )
*     FACSTR(IEX,ISTR)     : Phase factor of excitation.
*                                                                       
*     Jeppe Olsen, March 1994, LUCIA version 
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. a few include blocks
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,
C    &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR,
C    &              KSBIA,KSBIB,KVEC3,KPNIJ,KIJKK
      
      DIMENSION H(*)
*. A bit of local scratch
      DIMENSION IACAR(2),ITPAR(2)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ***************** '
        WRITE(6,*) ' Output from H1STR '
        WRITE(6,*) ' ***************** '
      END IF
*
*. Number of orbitals and offsets
      NI   = NOBPTS(IOBTP,IOBSM)
      IOFF = IBTSOB(IOBTP,IOBSM)
*
      NJ   = NOBPTS(JOBTP,JOBSM)
      JOFF = IBTSOB(JOBTP,JOBSM)
      IF(NTEST.GE.1000) THEN
        write(6,*) ' iobtp,iobsm,jobtp,jobsm',
     &               iobtp,iobsm,jobtp,jobsm
        write(6,*) ' NI NJ IOFF JOFF ', NI,NJ,IOFF,JOFF
      END IF
*
*. Offset for input string
      ISTROF =int_mb(KISTSO(ISTRGP)+(ISTRSM-1)*NOCTYP(ISTRGP)+ISTRTP-1)
*. Number of input strings  
      NIST = int_mb(KNSTSO(ISTRGP)
     &        + (ISTRSM-1)*NOCTYP(ISTRGP)+ISTRTP-1)
*. Group of output strings
      IACAR(1) = 2
      IACAR(2) = 1
*
      ITPAR(1) = IOBTP
      ITPAR(2) = JOBTP
      CALL NEWTYP(ISTRTP,IACAR,ITPAR,2,KSTRTP)
C     NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
*. Symmetry of output strings
      CALL SYMCOM(3,0,IOBSM,JOBSM,IJSXSM)
      CALL SYMCOM(3,0,ISTRSM,IJSXSM,KSTRSM)
C     WRITE(6,*) ' ISTRSM,IJSXSM KSTRSM',ISTRSM,IJSXSM,KSTRSM
C     SYMCOM(ITASK,IOBJ,I1,I2,I12)
*. Offset for output strings
      KSTROF = int_mb(KISTSO(ISTRGP)+(KSTRSM-1)*NOCTYP(ISTRGP)+KSTRTP-1)
C     WRITE(6,*) ' off set for output strings ', KSTROF
*
*. Type of mappings
*
*. Nel => Nel - 1 electrons
      IF(ISTAC(ISTRGP,1).NE.0.AND.ISTAC(ISTRGP,2).NE.0) THEN
*. full list
        IM1FUL = 1
        LM1    = NACOB
      ELSE
*. Truncated list
        IM1FUL = 0
        LM1 = NELEC(ISTRGP)
      END IF
*. Nel -1 => Nel mapping, must contain creation so
      IP1FUL = 1
      LP1 = NACOB
      NIEL = NELEC(ISTRGP)
      CALL  H1STRSO(H,ISTROF,NIST,KSTROF,             
     &            NI,IOFF,NJ,JOFF,
     &            int_mb(KSTSTM(ISTRGP,1)),int_mb(KSTSTM(ISTRGP,2)),
     &            LM1,IM1FUL,
     &            int_mb(KSTSTM(ISTRGP+1,1)),int_mb(KSTSTM(ISTRGP+1,2)),
     &            LP1,IP1FUL,
     &            MXSXST,IEXSTR,FACSTR,NEX,IH2TRM,WORK(KPNIJ),
     &            WORK(KIJKK),int_mb(KOCSTR(ISTRGP)),NIEL,
     &            NTOOB,NTEST)

*
      RETURN
      END
      SUBROUTINE H1STRSO(H,ISTROF,NIST,KSTROF,
     &                   NI,IOFF,NJ,JOFF,
     &                   IAMAPO,IAMAPS,LAMAP,IAMPFL,
     &                   ICMAPO,ICMAPS,LCMAP,ICMPFL,
     &                   MXSXST,IEXSTR,FACSTR,NEX,
     &                   IH2TRM,IPIJKK,XIJKK,IOCSTR,
     &                   NEL,NTOOB,NTEST)
*
* Slave routine mastered by H1STR 
* ( See my master for further information about my role in life )
*
* ==================
*. Additional input ( compared to SXSTR)
* ==================
*
* ISTROF : Absolute number of first string to be excited from
* NIST   : Number of strings to be excited from
* KSTROF : Offset of strings in resulting type-symmetry block
* N*,*OFF,*=I,J : Number and offset for each orbital set
*
* IAMAPO : Annihilation mapping, orbital part
* IAMAPS : Annihilation mapping, string part
* LAMAP  : Row dimension of Annihilation map
* IAMPFL : Annihilation map complete ?
*
* ICMAPO : Creation     mapping, orbital part
* ICMAPS : Creation     mapping, string part
* LCMAP  : Row dimension of Creation     map
*
* IH2TRM : ne. 0 : include zero and one-electron excitations from
*                  twobody operator
* IPIJKK : Pointer to symmetry adapted integral h(ij)
* XIJKK  : List of integrals  (ij!kk) - (ik!kj)
* NTOOB  : Number of orbitals
* IOCSTR : Occupation of input strings
* NEL    : Number of electrons in input string
*
* Jeppe Olsen, March 1994
      IMPLICIT REAL*8(A-H,O-Z)
*
*. Input
*
      INTEGER IAMAPO(LAMAP,*), IAMAPS(LAMAP,*)
      INTEGER ICMAPO(LCMAP,*), ICMAPS(LCMAP,*)
      DIMENSION H(NI,NJ)
      DIMENSION IPIJKK(NTOOB,NTOOB),XIJKK(NTOOB,*)
      DIMENSION IOCSTR(NEL,*)
*. Output
      INTEGER IEXSTR(MXSXST,*)
      INTEGER NEX(*)
      DIMENSION FACSTR(MXSXST,*)
*
C?    WRITE(6,*) ' SXFSTS : NTEST = ', NTEST
C?    WRITE(6,*) ' LCMAP LAMAP ', LCMAP,LAMAP
C?    WRITE(6,*) ' IAMPFL, ICMPFL ',IAMPFL, ICMPFL
C?    WRITE(6,*) ' MXSXST ', MXSXST 
*. To get rid of annoying and incorrect compiler warnings
      JISTR = 0
      SJ = 0.0D0
      IJISTR = 0
      SIJ = 0.0D0
*
      DO 1100 ISTR = ISTROF,ISTROF+NIST-1
C?      write(6,*) ' ISTR = ', ISTR 
        LEX = 0
        DO 1002 JORB = JOFF,JOFF+NJ-1
C?        write(6,*) ' JORB = ', JORB
*
* =====================================
* 1 :  Remove orbital JORB from ISTR
* =====================================
*
          JOCC = 0
          IF(IAMPFL.EQ.1) THEN
*. Read from full map
            IF(IAMAPO(JORB,ISTR).EQ.-JORB) THEN
              JOCC = 1
              IF(IAMAPS(JORB,ISTR).GT.0) THEN
                JISTR = IAMAPS(JORB,ISTR)
                SJ = 1.0D0
               ELSE 
                JISTR = -IAMAPS(JORB,ISTR)
                SJ = -1.0D0
              END IF
            END IF
          ELSE
*. Read from compact map
            DO JELEC = 1, LAMAP
              IF(IAMAPO(JELEC,ISTR).EQ.-JORB) THEN
                JOCC = 1
                IF(IAMAPS(JELEC,ISTR).GT.0) THEN
                   JISTR = IAMAPS(JELEC,ISTR)
                   SJ = 1.0D0
                ELSE 
                   JISTR = -IAMAPS(JELEC,ISTR)
                   SJ = -1.0D0
                END IF
              END IF
            END DO
          END IF
C?        WRITE(6,*) ' JOCC = ', JOCC
          IF(JOCC.EQ.0) GOTO 1002
*
* ==================================
*. 2 : Add orbital I to string JISTR
* ==================================
*
          DO 1001 IORB = IOFF,IOFF+NI-1
C?          WRITE(6,*) ' IORB = ', IORB
            IJACT = 0
            IF(ICMAPO(IORB,JISTR).EQ.+IORB) THEN
              IJACT  = 1
              IF(ICMAPS(IORB,JISTR).GT.0) THEN
                IJISTR = ICMAPS(IORB,JISTR)
                SIJ = SJ
              ELSE 
                IJISTR = -ICMAPS(IORB,JISTR)
                SIJ = -SJ 
              END IF
            END IF
C?          WRITE(6,*) ' IJACT = ', IJACT 
            IF(IJACT.EQ.0) GOTO 1001
*
*. A new excitation has been born, enlist it !
*
            LEX = LEX + 1
            IEXSTR(LEX,ISTR-ISTROF+1) = IJISTR-KSTROF+1
            FACSTR(LEX,ISTR-ISTROF+1) = 
     &      SIJ*H(IORB-IOFF+1,JORB-JOFF+1)
*. If IH2TRM .ne. 0 add
*
* sum (k)  a+i aj a+k a k /1+delta(i,j) (ij!kk)-(ik!kj)
           IF(IORB.EQ.JORB) THEN 
             FACIJ = 0.5D0*SIJ
           ELSE 
             FACIJ = 1.0D0*SIJ
           END IF
           IJEFF = IPIJKK(IORB,JORB)
*
C          WRITE(6,*) ' TESTING in H1STRS '
C          WRITE(6,*) ' IORB JORB FACIJ IJEFF ',
C    &                  IORB,JORB,FACIJ,IJEFF
           DO KEL = 1, NEL
             KORB = IOCSTR(KEL,ISTR)
             FACSTR(LEX,ISTR-ISTROF+1) = 
     &       FACSTR(LEX,ISTR-ISTROF+1) +
     &       FACIJ*XIJKK(KORB,IJEFF)
*
C            WRITE(6,*) ' TESTING in H1STRS '
C            WRITE(6,*) ' ISTR KEL KORB',ISTR,KEL,KORB
C            WRITE(6,*) ' XIJKK(K,IJEFF) and explicit '
C            XIJKK2 = GTIJKL(IORB,JORB,KORB,KORB)-
C    &                GTIJKL(IORB,KORB,KORB,JORB)
C            WRITE(6,*) XIJKK2, XIJKK(KORB,IJEFF) 
*
           END DO
            
             
 1001     CONTINUE
 1002   CONTINUE
        NEX(ISTR-ISTROF+1) = LEX
 1100 CONTINUE
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' **************** '
         WRITE(6,*) ' H1STRS reporting '
         WRITE(6,*) ' **************** '
         WRITE(6,*)
         DO ISTR = ISTROF,ISTROF+NIST-1
           WRITE(6,*) ' excitations from string ',ISTR
           WRITE(6,*)
           WRITE(6,*) '    exc.string      factor '
           WRITE(6,*) ' =============================='
           DO LEX = 1, NEX(ISTR-ISTROF+1)
             WRITE(6,'(3X,I8,F13.8)')
     &       IEXSTR(LEX,ISTR-ISTROF+1),FACSTR(LEX,ISTR-ISTROF+1)
           END DO
           WRITE(6,*)
         END DO
      END IF
*
      RETURN
      END
      SUBROUTINE SXSTRSO(ISTROF,NIST,KSTROF,
     &                   NI,IOFF,NJ,JOFF,
     &                   IAMAPO,IAMAPS,LAMAP,IAMPFL,
     &                   ICMAPO,ICMAPS,LCMAP,ICMPFL,
     &                   MXSXST,ISXSTR,JSXSTR,IEXSTR,FACSTR,NEX,
     &                   IOFFDG,NTEST)
*
* Obtain single excitations from string ISTROF-ISTROF+NIST-1, 
* Slave routine mastered by SXSTR 
*
* ==================
*. Additional input ( compared to SXSTR)
* ==================
*
* ISTROF : Absolute number of first string to be excited from
* NIST   : Number of strings to be excited from
* KSTROF : Offset of strings in resulting type-symmetry block
* N*,*OFF,*=I,J : Number and offset for each orbital set
*
* IAMAPO : Annihilation mapping, orbital part
* IAMAPS : Annihilation mapping, string part
* LAMAP  : Row dimension of Annihilation map
* IAMPFL : Annihilation map complete ?
*
* ICMAPO : Creation     mapping, orbital part
* ICMAPS : Creation     mapping, string part
* LCMAP  : Row dimension of Creation     map
*
* Jeppe Olsen, March 1994
      IMPLICIT REAL*8(A-H,O-Z)
*
*. Input
*
      INTEGER IAMAPO(LAMAP,*), IAMAPS(LAMAP,*)
      INTEGER ICMAPO(LCMAP,*), ICMAPS(LCMAP,*)
*. Output
      INTEGER ISXSTR(MXSXST,*),JSXSTR(MXSXST,*)
      INTEGER IEXSTR(MXSXST,*)
      DIMENSION FACSTR(MXSXST,*)
      INTEGER NEX(*)
*. To get rid of annoying and incorrect compiler warnings
      JISTR = 0
      SJ = 0.0D0
      IJISTR = 0
      SIJ = 0.0D0
*
C     WRITE(6,*) ' SXFSTS : NTEST = ', NTEST
C     WRITE(6,*) ' LCMAP LAMAP ', LCMAP,LAMAP
C     WRITE(6,*) ' IAMPFL, ICMPFL ',IAMPFL, ICMPFL
C     WRITE(6,*) ' MXSXST ', MXSXST 
      DO 1100 ISTR = ISTROF,ISTROF+NIST-1
C?      write(6,*) ' ISTR = ', ISTR 
        LEX = 0
        DO 1002 JORB = JOFF,JOFF+NJ-1
C?        write(6,*) ' JORB = ', JORB
*
* =====================================
* 1 :  Remove orbital JORB from ISTR
* =====================================
*
          JOCC = 0
          IF(IAMPFL.EQ.1) THEN
*. Read from full map
            IF(IAMAPO(JORB,ISTR).EQ.-JORB) THEN
              JOCC = 1
              IF(IAMAPS(JORB,ISTR).GT.0) THEN
                JISTR = IAMAPS(JORB,ISTR)
                SJ = 1.0D0
               ELSE 
                JISTR = -IAMAPS(JORB,ISTR)
                SJ = -1.0D0
              END IF
            END IF
          ELSE
*. Read from compact map
            DO JELEC = 1, LAMAP
              IF(IAMAPO(JELEC,ISTR).EQ.-JORB) THEN
                JOCC = 1
                IF(IAMAPS(JELEC,ISTR).GT.0) THEN
                   JISTR = IAMAPS(JELEC,ISTR)
                   SJ = 1.0D0
                ELSE 
                   JISTR = -IAMAPS(JELEC,ISTR)
                   SJ = -1.0D0
                END IF
              END IF
            END DO
          END IF
C?        WRITE(6,*) ' JOCC = ', JOCC
          IF(JOCC.EQ.0) GOTO 1002
*
* ==================================
*. 2 : Add orbital I to string JISTR
* ==================================
*
          DO 1001 IORB = IOFF,IOFF+NI-1
C?          WRITE(6,*) ' IORB = ', IORB
            IF(IOFFDG.NE.0 .AND. IORB.EQ.JORB) GOTO 1001
            IJACT = 0
            IF(ICMAPO(IORB,JISTR).EQ.+IORB) THEN
              IJACT  = 1
              IF(ICMAPS(IORB,JISTR).GT.0) THEN
                IJISTR = ICMAPS(IORB,JISTR)
                SIJ = SJ
              ELSE 
                IJISTR = -ICMAPS(IORB,JISTR)
                SIJ = -SJ 
              END IF
            END IF
C?          WRITE(6,*) ' IJACT = ', IJACT 
            IF(IJACT.EQ.0) GOTO 1001
*
*. A new excitation has been born, enlist it !
*
            LEX = LEX + 1
            ISXSTR(LEX,ISTR-ISTROF+1) = IORB
            JSXSTR(LEX,ISTR-ISTROF+1) = JORB
            IEXSTR(LEX,ISTR-ISTROF+1) = IJISTR-KSTROF+1
            FACSTR(LEX,ISTR-ISTROF+1) = SIJ
C?          WRITE(6,*) 'FACSTR = ',  FACSTR(LEX,ISTR-ISTROF+1)
 1001     CONTINUE
 1002   CONTINUE
        NEX(ISTR-ISTROF+1) = LEX
 1100 CONTINUE
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' **************** '
         WRITE(6,*) ' SXSTRS reporting '
         WRITE(6,*) ' **************** '
         WRITE(6,*)
         DO ISTR = ISTROF,ISTROF+NIST-1
           WRITE(6,*) ' excitations from string ',ISTR
           WRITE(6,*)
           WRITE(6,*) ' iorb jorb exc.string phase '
           WRITE(6,*) ' =========================='
           DO LEX = 1, NEX(ISTR-ISTROF+1)
             WRITE(6,'(2I4,I8,F8.3)')
     &       ISXSTR(LEX,ISTR-ISTROF+1),JSXSTR(LEX,ISTR-ISTROF+1),
     &       IEXSTR(LEX,ISTR-ISTROF+1),FACSTR(LEX,ISTR-ISTROF+1)
           END DO
           WRITE(6,*)
         END DO
      END IF
*
      RETURN
      END
      SUBROUTINE DXFSTR(IKIN1,IKIN2,IKST,IKSG,NIK,
     &                  KJIN1,KJIN2,KJST,KJSG,NKJ,
     &                  NI,NJ,NK,NL,
     &                  IOFF,JOFF,KOFF,LOFF,
     &                  ISTR,MXSXST,XIJKL,
     *                  IEXSTR,FACSTR,NEX,NTESTG)
*. Two list of single excitations occurs 
*
*  IK : String I => String K (IKIN1,IKIN2,IKST,IKSG)
*  KJ : String K => String J (KJIN1,KJIN2,KJST,KJSG)
*
* Obtain effect of two-electron excitation part of two-electron operator 
*                                                                       
*     
* g = 1/2sum(i,j,k,l) a+ia+ka l a j (ij!kl)
*   = sum(i>k,j>l)    a+ia+ka l a j ((ij!kl) -(il!jk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   + sum(i>k,k>l)    a+ia k a+k a l ((ik!kl) -(il!kk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   - sum(i>k,k>l)    a+ia+k a k a l ((ik!kl) -(il!kk))
*   = sum(i>k,j>l,j.ne.k)    a+ia j a+k a l ((ij!kl) -(il!jk))
*   - sum(i>k,k>j)    a+ia j a+k a k ((ik!kj) -(ij!kk))
*
* COLD :
*   = sum(i>k,j.ne.k, j.gt.l)
*     a+ia j a+k a l ((ij!kl) -(il!jk))
*   + sum(i>k, j>l)
*     a+ia j a+j a l ((ij!jl) -(il!jj))
* 
* as 
*   G !ISTR> = sum(i>k,j.ne.k,if(k.ne.l ) j.gt.l)
*     a+ia ji!KSTR><KSTR! a+k a l!ISTR> ((ij!kl) -(il!jk))
*
* Version only employing single annihilation tables
*
* =====                                                             
* Input                                                             
* =====                                                             
*
* IKIN1 : Index one of <K STR! SX ! ISTR > i.e k
* IKIN2 : Index two of <K STR! SX ! ISTR > i.e l
* IKST  : Corresponding K string
* IKSG :  Phase shift of excitation
* NIK  : Number of excitations per string I=> K
*
* KJIN1 : Index one of <J STR! SX ! KSTR > i.e i
* KJIN2 : Index two of <J STR! SX ! KSTR > i.e j
* KJST  : Corresponding J string
* KJSG :  Phase shift of excitation
* NKJ  : Number of excitations per string K => J
*
* NI,NJ,NK,NL : Number of orbitals in each group
* IOFF,JOFF,KOFF,LOFF : Offset for each group of orbitals.
*
*                                                                       
* XIJKL : Integrals stored as (j,l,i,k)
* ======                                                            
* Output                                                            
* ======                                                            
*      
*     NEX    : Number of excitations
*                                                                       
*     IEXSTR : IEXSTR(K) is string number of excited string K           
*              relative to start of TS block
*     FACSTR : FACSTR(K) is the factor that string INDSTR(K) is multiplied 
*              with
*
*     Jeppe Olsen, March 1994, LUCIA version of DXFSTR
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IKIN1(MXSXST,*),IKIN2(MXSXST,*)
      INTEGER IKST(MXSXST,*),NIK(*)
      REAL*8  IKSG(MXSXST,*)
*
      INTEGER KJIN1(MXSXST,*),KJIN2(MXSXST,*)
      INTEGER KJST(MXSXST,*),NKJ(*)
      REAL*8  KJSG(MXSXST,*)
*
      DIMENSION XIJKL(NJ,NL,NI,NK)
*
*. Output
      INTEGER IEXSTR(*)
      DIMENSION FACSTR(*)
*
      CALL LUCIAQENTER('DXFST')
*
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ****************************************'
        WRITE(6,*) ' DXFSTR at your command : ISTR = ',ISTR
        WRITE(6,*) ' ****************************************'
      END IF
*   
      NEX = 0
C?    WRITE(6,*) ' Number of excitations from string ', NIK(ISTR)
      DO 1000 KLEX = 1, NIK(ISTR)
        K = IKIN1(KLEX,ISTR)
        L = IKIN2(KLEX,ISTR)
        KSTR = IKST(KLEX,ISTR)
        SIGN1 = IKSG(KLEX,ISTR)
C?      WRITE(6,*) ' KLEX K L ', KLEX,K,L
C?      WRITE(6,*) ' KSTR SIGN1 ',KSTR,SIGN1
        DO 900 IJEX = 1, NKJ(KSTR)
          I = KJIN1(IJEX,KSTR)
          J = KJIN2(IJEX,KSTR)
C?        WRITE(6,*) ' I J ', I,J
          ICASE = 0
          IF(I.GT.K .AND. J.GT.L .AND.
     &       I.NE.L.AND.I.NE.J.AND.K.NE.L.AND.K.NE.J) THEN
            ICASE = 1
          END IF
          IF(ICASE.EQ.0) GOTO 900
          JSTR = KJST(IJEX,KSTR)
C?        write(6,*) ' KSTR = ', KSTR
C?        write(6,'(A,5I3)') 
C?   &   ' Allowed excitation i,j,k,l  JSTR ',I,J,K,L,JSTR
          IF(KJSG(IJEX,KSTR).EQ.1.0D0) THEN
            SIGN =   SIGN1
          ELSE
            SIGN = - SIGN1
          END IF
C?        write(6,*) ' jstr,kstr sign ', jstr,kstr,sign
*
          NEX = NEX + 1
          IEXSTR(NEX) = JSTR
          IF(ICASE.EQ.1) THEN
C    &      SIGN*XIJKL(J-JOFF+1,L-LOFF+1,I-IOFF+1,K-KOFF+1)
            FACSTR(NEX) = 
     &      SIGN*(GTIJKL(I,J,K,L)-GTIJKL(I,L,K,J))
C    &      SIGN*XIJKL(J-JOFF+1,L-LOFF+1,I-IOFF+1,K-KOFF+1)
          END IF
C         IF(FORM1.NE.FACSTR(NEX) ) THEN
C           write(6,*) ' problemo in DXFSTR'
C           write(6,*) ' i j k l Form1 FACSTR '
C           write(6,*) i,j,k,l,FORM1,FACSTR(NEX)
C           write(6,*) ' ijkl and ilkj '
C           write(6,*) GTIJKL(I,J,K,L),GTIJKL(I,L,K,J)
C         END IF
  900   CONTINUE
 1000 CONTINUE
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 'Number of excitations obtained ', NEX
        WRITE(6,*)
        WRITE(6,*) ' iex      exc.string         factor '
        WRITE(6,*) ' =================================='
        DO IEX = 1, NEX
          WRITE(6,'(I4,4X,I8,6X,F13.5)')
     &    IEX,IEXSTR(IEX),FACSTR(IEX)
        END DO
      END IF
*
      CALL LUCIAQEXIT ('DXFST')
      RETURN
      END 
      SUBROUTINE MATML7_ORIG(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C =  FACTORC*C + FACTORAB* A * B             FOR ITRNSP = 0
C
C     C =  FACTORC*C + FACTORAB* A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C =  FACTORC*C + FACTORAB* A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN,
C
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
*
c     CALL LUCIAQENTER('MATML')
*
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' A, B and C MATRIX FROM MATML7 ' 
      WRITE(6,*)
      CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL)
      CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL ' 
      WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      WRITE(6,*) ' ITRNSP : ', ITRNSP
      END IF
*
      IF(NAROW*NACOL*NBROW*NBCOL*NCROW*NCCOL .EQ. 0 ) THEN
        IZERO = 1
      ELSE
        IZERO = 0
      END IF

      IESSL = 0
      ICONVEX = 0
      IF ((ICONVEX.EQ.1 .OR.IESSL.EQ.1).AND. IZERO.EQ.0 ) THEN
*. DGEMM from CONVEX/ESSL  lib
        LDA = MAX(1,NAROW)
        LDB = MAX(1,NBROW)
* 
        LDC = MAX(1,NCROW)
        IF(ITRNSP.EQ.0) THEN
        CALL DGEMM('N','N',NAROW,NBCOL,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF (ITRNSP.EQ.1) THEN
        CALL DGEMM('T','N',NACOL,NBCOL,NAROW,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF(ITRNSP.EQ.2) THEN
        CALL DGEMM('N','T',NAROW,NBROW,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        END IF
      ELSE
* Use Jeppes version ( it should be working )
        IF( ITRNSP .EQ. 0 ) THEN     
* ======
* C=A*B
* ======
          DO  J = 1, NCCOL
            DO  I = 1, NCROW
              T = 0.0D0
              DO  K = 1, NBROW
                T = T  + A(I,K)*B(K,J)
              END DO  
              C(I,J) = FACTORC*C(I,J) + FACTORAB*T
            END DO  
          END DO  
C         CALL SGEMMX@ (NCROW, NCCOL, NBROW, FACTORAB, 
C    &                 A(1,1), 1, NAROW, B(1,1), 
C    1       1   , NBROW, FACTORC , C(1,1), 1, NCROW)
        END IF
        IF ( ITRNSP .EQ. 1 ) THEN      
*
* =========
* C=A(T)*B
* =========
*
          CALL SCALVE(C,FACTORC,NCROW*NCCOL)
C         CALL SETVEC(C,0.0D0,NCROW*NCCOL)
          DO J = 1, NCCOL
            DO K = 1, NBROW
              DO I = 1, NCROW
                C(I,J)= C(I,J) + FACTORAB*A(K,I)*B(K,J)
              END DO   
            END DO   
          END DO   
C         CALL SGEMMX@ (NCROW, NCCOL, NBROW, FACTORAB, A(1,1),
C    &          1*NAROW, 1, B(1,1)
C    1   , 1    , NBROW, FACTORC , C(1,1), 1, NCROW)
        END IF   
C
        IF ( ITRNSP .EQ. 2 ) THEN     
C ===========
C. C = A*B(T)
C ===========
          DO J = 1,NCCOL
            DO I = 1, NCROW
              T = 0.0D0
              DO K = 1,NBCOL
                T  = T  + A(I,K)*B(J,K)
              END DO   
            C(I,J) = FACTORC*C(I,J) + FACTORAB*T
            END DO   
          END DO   
C         CALL SGEMMX@ (NCROW, NCCOL, NBCOL, FACTORAB, A(1,1),
C    &   1,             NAROW, B(1,1), 1
C    1   *              NBROW, 1, FACTORC , C(1,1), 1, NCROW)
        END IF   
      END IF
C
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' C MATRIX FROM MATML7 ' 
      WRITE(6,*)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
C
c     CALL LUCIAQEXIT('MATML')
      RETURN
      END
      SUBROUTINE PR_MATML_STAT
*
* Print statistics from matrix multiplier
* and other statistics !
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT
      COMMON/COPVECST/XNCALL_COPVEC,XNMOVE_COPVEC
*
      CHARACTER*6 MARKC, IIDENT,MAX_BLKC
      PARAMETER(MAXLVL = 1024 )
      PARAMETER(MAXMRK = 1024)
      COMMON/CMEMO/NWORD,KFREES,KFREEL,NS,NL,NM,IBASE(MAXLVL),
     &             LENGTH(MAXLVL),IIDENT(MAXLVL),IMARK(MAXMRK),
     &             MARKL(MAXMRK),MARKS(MAXMRK),MARKC(MAXMRK),
     &             MAX_MEM,MAX_BLK,MAX_BLKC
*
      WRITE(6,*)
      WRITE(6,*) ' Information about COPVEC calls: '
      WRITE(6,*) ' ==============================='
      WRITE(6,*)
      WRITE(6,'(A,E8.3)') ' Number of calls = ', XNCALL_COPVEC
      WRITE(6,'(A,E8.3)') ' Number of R*8 words copied = ', 
     &                    XNMOVE_COPVEC
*
      WRITE(6,*)
      WRITE(6,*) ' Information about MATML7 calls: '
      WRITE(6,*) ' ================================'
      WRITE(6,*)
      WRITE(6,*) 'Number of calls= ', XNCALL
      WRITE(6,'(A,E8.3)') ' Number of flops executed = ', XNFLOP
      WRITE(6,'(A,E8.3)') ' Average row length of C = ',
     &           XLCROWCOL/XLCCOL
      WRITE(6,'(A,E8.3)') ' Average column  length of C = ',
     &           XLCROWCOL/XLCROW
      WRITE(6,'(A,E8.3)') 
     &           ' Average number of operations per element of  C = ',
     &           XNFLOP/XLCROWCOL
      WRITE(6,'(A,E8.3)') 
     &           ' Average number of operations per per call = ',
     &           XNFLOP/XNCALL
      WRITE(6,'(A,E8.3)') ' Number of seconds spent in MATML7 = ', TMULT
*
      WRITE(6,'(A,E8.3)') ' Average MFLOPS = ',XNFLOP/TMULT/1000000.0D0
*
c now output by MEMMAN, TASK 'STATI'
c      WRITE(6,*) ' Large size of active memory     : ', MAX_MEM
c      WRITE(6,*) ' Large allocated block of memory : ', MAX_BLK
c      WRITE(6,'(A,A)') ' Identifier of largest block : ', MAX_BLKC
      RETURN
      END


      SUBROUTINE MATML7_COPY(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C =  FACTORC*C + FACTORAB* A * B             FOR ITRNSP = 0
C
C     C =  FACTORC*C + FACTORAB* A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C =  FACTORC*C + FACTORAB* A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN,
C
*. Notice : If the summation index has dimension zero nothing
*           is performed
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT
*
c     CALL LUCIAQENTER('MATML')
*
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' A, B and initial C MATRIX FROM MATML7 ' 
      WRITE(6,*)
      CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL)
      WRITE(6,*)
      CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL)
      WRITE(6,*)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL ' 
      WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      WRITE(6,*) ' FACTORC, FACTORAB, ITRNSP : ', 
     &             FACTORC, FACTORAB, ITRNSP
      END IF
*. Statistics
      XNCALL = XNCALL + 1
      XLCROW = XLCROW + NCROW
      XLCCOL = XLCCOL + NCCOL
      XLCROWCOL = XLCROWCOL + DBLE(NCROW)*DBLE(NCCOL)
      T_INI = 0.01*LTIMEX()
*
      IF(ITRNSP.EQ.1) THEN
        XNFLOP = XNFLOP + 2*DBLE(NCROW)*DBLE(NCCOL)*DBLE(NAROW)
      ELSE
        XNFLOP = XNFLOP + 2*DBLE(NCROW)*DBLE(NCCOL)*DBLE(NACOL)
      END IF
*
      IF(NAROW*NACOL*NBROW*NBCOL*NCROW*NCCOL .EQ. 0 .OR.
     &   FACTORC.EQ.0.0D0) THEN
        IZERO = 1
      ELSE
        IZERO = 0
      END IF
*
      IF(IZERO.EQ.1.AND.NCROW*NCCOL.NE.0) THEN 
        IF(FACTORC.NE.0.0D0) THEN
         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
        ELSE 
         ZERO = 0.0D0
         CALL SETVEC(C,ZERO,NCROW*NCCOL)
        END IF
C?        WRITE(6,*) ' Scaled C '
C?        CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
*
      IESSL = 0
      ICONVEX = 0
      IF ((ICONVEX.EQ.1 .OR.IESSL.EQ.1).AND. IZERO.EQ.0 ) THEN
*. DGEMM from CONVEX/ESSL  lib
        LDA = MAX(1,NAROW)
        LDB = MAX(1,NBROW)
* 
        LDC = MAX(1,NCROW)
        IF(ITRNSP.EQ.0) THEN
        CALL DGEMM('N','N',NAROW,NBCOL,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF (ITRNSP.EQ.1) THEN
        CALL DGEMM('T','N',NACOL,NBCOL,NAROW,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF(ITRNSP.EQ.2) THEN
        CALL DGEMM('N','T',NAROW,NBROW,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        END IF
      ELSE
* Use Jeppes version ( it should be working )
        IF( ITRNSP .EQ. 0 ) THEN     
* ======
* C=A*B
* ======
*
C         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
          DO J =1, NCCOL
*. Initialize with FACTORC*C(I,J) + FACTORAB*A(I,1)*B(1,J)
            IF(NBROW.GE.1) THEN
              B1J = FACTORAB*B(1,J)
              DO I = 1, NCROW
                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1)
              END DO
            END IF
*. and the major part
            DO K =2, NBROW
              BKJ = FACTORAB*B(K,J)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BKJ*A(I,K)
              END DO
            END DO
          END DO
*
        END IF
        IF ( ITRNSP .EQ. 1 ) THEN      
*
* =========
* C=A(T)*B
* =========
*
          DO J = 1, NCCOL
            DO I = 1, NCROW
              T = 0.0D0         
              DO K = 1, NBROW
                T = T  + A(K,I)*B(K,J)
              END DO   
              C(I,J) = FACTORC*C(I,J) + FACTORAB*T
            END DO   
          END DO   
        END IF   
C
        IF ( ITRNSP .EQ. 2 ) THEN     
C ===========
C. C = A*B(T)
C ===========
          DO J = 1,NCCOL
*. Initialization
            IF(NBCOL.GE.1) THEN
              BJ1 = FACTORAB*B(J,1)
              DO I = 1, NCROW
                C(I,J) = FACTORC*C(I,J) + BJ1*A(I,1)
              END DO   
            END IF
*. And the rest
            DO K = 2,NBCOL
              BJK = FACTORAB*B(J,K)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BJK*A(I,K)
              END DO   
            END DO   
          END DO   
        END IF   
      END IF
C
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' C MATRIX FROM MATML7 ' 
      WRITE(6,*)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
C
      T_END = 0.01*LTIMEX()
      TMULT = TMULT + T_END - T_INI 
C
C     CALL LUCIAQEXIT('MATML')
      RETURN
      END
      SUBROUTINE MICDV6(MV7,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,
     &                  MAXIT,NVAR,
     &                  LU3,LU4,LU5,LU6,LU7,LU8,LUDIA,NROOT,MAXVEC,
     &                  NINVEC,
     &                  APROJ,AVEC,WORK,IPRTXX,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,IROOTHOMING,NCNV_RT,RTCNV,IPRECOND,
     &                  CONVER,RNRM_CNV,ISBSPPR_ACT)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version 
*
* Version with 2 or 3 vectors in subspace per root
*
* From MICDV4, Jeppe Olsen, April 1997
*                           Oct 97: root homing added
*                           May 99: New preconditioner switch 
*                           Feb.13: Extensive recoding for improved stability and readability..
*
* Input :
* =======
*        MV7: Sigma routine 
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENLVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENLVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       LOGICAL RTCNV(NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
       DIMENSION RNRM_CNV(NROOT)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*MAXVEC                               
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
*      IROOTHOMING : Do roothoming, i.e. select the
*      eigenvectors in iteration n+1 as the approximations
*      with largest overlap with the previous space
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER
       REAL*8 INPRDD, INPROD
*. Notice XJEP is also used for ROOTHOMING, should be allocated
* outside (for roothoming :dim = 3*MAXVEC )
       DIMENSION XJEP(10000)
       INTEGER   IXJEP(10000)
*. Dimension of Work2: 3*MAXVEC**2
       DIMENSION WORK2(10000)
*. For working with routines accessing several files
       DIMENSION LU_FILE1(5),NVC_FILE1(5)
       DIMENSION LU_FILE2(5),NVC_FILE2(5)
*
       EXTERNAL MV7
*
       IPICO = 0
       I_CHECK_OVLAP = 1
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
       IPRT = 10
       IPRT = MAX(IPRTXX,IPRT)
       IF(IPRT.GE.0) WRITE(6,*) ' MICDV6 in action, NROOT, MAXVEC =',
     &                            NROOT,MAXVEC
       WRITE(6,*) ' TESTY, ISBSPPR_ACT = ', ISBSPPR_ACT
*
       IOLSTM = 0
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop on MICDV6 '
         STOP 20
       END IF
*
C?     WRITE(6,*) ' MICDV6 : NROOT, IPRECOND ', NROOT, IPRECOND
C?     WRITE(6,*) ' LU1 LU2 LU3 LU4 LU5 LU6',
C?   &              LU1,LU2,LU3,LU4,LU5,LU6 
       IF(IROOTHOMING.EQ.1) THEN
         WRITE(6,*) ' Root homing performed '
       END IF
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = -1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       CALL GFLUSH(6)
       ITER = 1
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       DO IVEC = 1,NINVEC
         CALL REWINO(LU5)
         CALL REWINO(LU6)
         CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         CALL MV7(VEC1,VEC2,LU5,LU6,0,0)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINO(LU6)
         CALL COPVCD(LU6,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINO(LU2)
         DO JVEC = 1, IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
         END DO
*        ^ End of loop over JVEC
       END DO
*      ^ End of loop over IVEC
*
       IF(IPRT.GE.10000) THEN
         CALL REWINO(LU1)
         CALL REWINO(LU2)
         DO IVEC = 1, NINVEC
           WRITE(6,*)
           WRITE(6,*) ' ==============================='
           WRITE(6,*) ' Initial C and S vector ', IVEC
           WRITE(6,*) ' ==============================='
           WRITE(6,*)
           CALL WRTVCD(VEC1,LU1,0,LBLK)
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
     
       IF( IPRT .GE.5 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGENL(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
       END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,*) ' Eigenvalues of initial iteration (with shift)'
         WRITE(6,'(5F21.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(6,*) ' Initial set of eigen values (no shift) '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
*. Transform vectors on LU1 so they become the actual
*. eigenvector approximations
       CALL REWINO(LU3)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINO(LU3)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,LBLK,LU4,LU5)
*
       IF(IPRT.GE.10000) THEN
         WRITE(6,*) ' Initial set of eigenvectors '
         CALL REWINO(LU1)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
*
         WRITE(6,*) ' Initial set of sigma vectors '
         CALL REWINO(LU2)
         DO IROOT = 1, NROOT
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
*. And the corresponding Hamiltonian matrix, no problems
*. with numerical stabilities, so just use eigenvalues
       ZERO = 0.0D0
       CALL SETVEC(APROJ,ZERO,NROOT*(NROOT+1)/2)
       DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
       END DO
*
       NVEC = NROOT 
CJAN20 IF (MAXIT .EQ. 1 ) GOTO  901
C      IF (MAXIT .EQ. 1 ) GOTO  1001
*
* ======================
*. Loop over iterations
* ======================
*
*
 1000 CONTINUE
       ITER = ITER + 1
       IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*) ' Info from iteration .... ', ITER
       END IF
       CALL GFLUSH(6)
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       CALL REWINO(LU1)
       CALL REWINO(LU2)
*
       DO 100 IROOT = 1, NROOT
         IF(IPRT.GE.100) WRITE(6,*) ' Loop 100, IROOT = ', IROOT
*
*
*. Save current eigenvector IROOT on LU7
         CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
         CALL REWINO(LU7)
         CALL COPVCD(LU1,LU7,VEC1,0,LBLK)
*. Calculate (HX - EX ) and store on LU5
*. Current eigenvector is  on LU7, corresponding sigma vector
*. on LU2
         EIGAPR = EIG(ITER-1,IROOT)
         ONE = 1.0D0
*
         CALL REWINO(LU7)
         CALL REWINO(LU5)
         FACTOR = - EIGAPR
         CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU2,LU7,LU5,0,LBLK)
*
         IF ( IPRT  .GE. 10000 ) THEN
           WRITE(6,*) '  ( HX - EX ) for IROOT =', IROOT
           CALL WRTVCD(VEC1,LU5,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         RNORM = SQRT( ABS(INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)) )
         RNRM(ITER-1,IROOT) = RNORM
C?       WRITE(6,*) ' MICDV6: ITER = ', ITER
*
         IF(IPRT.GE.1000) WRITE(6,*) ' RNORM = ', RNORM
         IF(RNORM.LT. TEST .OR. 
     &      (ITER.GT.2.AND.
     &      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.E_CONV)) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            IF(IROOT.LE.NCNV_RT) CONVER = .FALSE.
         END IF
*
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU6
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.10) THEN
             WRITE(6,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINO(LUDIA)
           CALL REWINO(LU5)
           CALL REWINO(LU6)
           CALL H0M1TD(LU6,LUDIA,LU5,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND,ISBSPPR_ACT)
           IF ( IPRT  .GE. 1000) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU6,1,LBLK)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU5
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             CALL REWINO(LUDIA)
*
             CALL H0M1TD(LU5,LUDIA,LU7,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND,ISBSPPR_ACT)
*. Get X back on LU7
             CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
             CALL REWINO(LU7)
             CALL COPVCD(LU1,LU7,VEC1,0,LBLK)

*
* Gamma = X(T) * (H0 - E) ** -1 * X
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             GAMMA = INPRDD(VEC1,VEC2,LU5,LU7,0,LBLK)
             IF(IPRT.GE.10) WRITE(6,*) ' Gamma = ', Gamma
* is X an eigen vector for (H0 - 1 ) - 1
             CALL REWINO(LU5)
             CALL REWINO(LU7)
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU7,LU5,0,LBLK))
              IF(IPRT.GE.1000) write(6,*) ' VNORM = ', VNORM
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                FACTOR = -(DELTA/GAMMA)
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU6,LU5,LU7,1,LBLK)
                CALL COPVCD(LU7,LU6,VEC1,1,LBLK)
                XNORM2 = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
                XNORM_INI = SQRT(ABS(XNORM2))
*
                IF(IPRT.GE.10000) THEN
                  WRITE(6,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU6,1,LBLK)
                END IF
*
              END IF
            END IF
*
*. 1.3 Orthogonalize to all previous vectors
*
*
* MGS- more I/O intensive and stable
*
C              ADD_ORTN_2SUBSPC(LU1,LU2,N1,N2,LUNEW,IMNEW,
C    &                          VEC1,VEC2,LUSC1,LUSC2,THRES)
            THRES = 1.0D-3 
            CALL ADD_ORTN_2SUBSPC(LU1,LU3,NROOT,NVEC+IADD-1-NROOT,
     &           LU6,IMNEW,VEC1,VEC2,LU7,LU5,THRES)
            I_SUPER_ORT = 0
            IF(I_SUPER_ORT.EQ.1) THEN
*. A second orthonormalization 
              THRES = 1.0D-3
              CALL ADD_ORTN_2SUBSPC(LU1,LU3,NROOT,NVEC+IADD-1-NROOT,
     &             LU6,IMNEW,VEC1,VEC2,LU7,LU5,THRES)
            END IF
*
            IF(IMNEW.EQ.1) THEN
*. Vector is new, copy to LU3- which was left at EOF by add
              CALL REWINO(LU6)
              CALL COPVCD(LU6,LU3,VEC1,0,LBLK)
            ELSE
*. Nothing new 
              IADD = IADD - 1
              IF(IPRT.GE.10) THEN
               WRITE(6,*) ' Correction vector removed for root',IROOT
              END IF
            END IF
         END IF ! Root is not converged
  100  CONTINUE
*. IADD new vectors have been obtained, NROOT are on LU1, IADD + NVEC-NROOT on LU3
       IF(IPRT.GE.3) 
     &  WRITE(6,*) ' Number of orthogonal corrections added ',  IADD
       NADD1 = IADD
       NVEC3 = NVEC-NROOT+IADD
*. Notice: NVEC is not updated yet, it is still the number of vectors from prev. it
*
       IF(I_CHECK_OVLAP.EQ.1) THEN
*. Overlap LU1 LU2
         WRITE(6,*) ' Overlap roots'
C             GET_OVLAP_FOR_VECTORS(LU1,N1,LU2,N2,LUSCR,OVLAP,VEC1,VEC2)
         CALL GET_OVLAP_FOR_VECTORS(LU1,NROOT,LU1,NROOT,LU7,
     &        XJEP,VEC1,VEC2)
         CALL CHECK_UNIT_MAT(XJEP,NROOT,XMAX_DIA,XMAX_OFD,1)
C             CHECK_UNIT_MAT(UNI,NDIM,XMAX_DIFF_DIAG,XMAX_DIFF_OFFD,ISYM)
 
         WRITE(6,*) ' Overlap roots and corrections '
         CALL GET_OVLAP_FOR_VECTORS(LU1,NROOT,LU3,NVEC3,LU7,
     &        XJEP,VEC1,VEC2)
         XMAX = XMNMX(XJEP,NROOT*NVEC3,3)
         WRITE(6,*) ' Largest element of <roots!Corrections> ', XMAX
         WRITE(6,*) ' Overlap corrections '
         CALL GET_OVLAP_FOR_VECTORS(LU3,NVEC3,LU3,NVEC3,LU7,
     &        XJEP,VEC1,VEC2)
         CALL CHECK_UNIT_MAT(XJEP,NVEC3,XMAX_DIA,XMAX_OFD,1)
       END IF
       IF(CONVER) THEN
*. Well, we converged     
C         ITER = ITER-1
          GOTO  1001
       END IF
       IF( ITER.GT. MAXIT) THEN
          ITER = MAXIT
          GOTO 1001
       END IF
*
* ====================================================
*  2 : Optimal combination of new and old directions
* ====================================================
*
*  2.1: Multiply new directions with matrix
*
       IF(IPRT.GE.10000) THEN
         WRITE(6,*) ' Vectors on LU3'
         WRITE(6,*) ' Number of vectors = ', NVEC3
         CALL REWINO(LU3)
         DO IVEC = 1, NVEC3
           CALL WRTVCD(VEC1,LU3,0,LBLK)
         END DO 
       END IF
*
       CALL SKPVCD(LU3,NVEC-NROOT,VEC1,1,LBLK)
       CALL SKPVCD(LU4,NVEC-NROOT,VEC1,1,LBLK)
       DO IVEC = 1, NADD1
         CALL REWINE(LU5,LBLK)
         CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         CALL MV7(VEC1,VEC2,LU5,LU6,0,0)
         CALL REWINE(LU6,LBLK)
         CALL COPVCD(LU6,LU4,VEC1,0,LBLK)
*.2.2 Augment projected matrix
         CALL REWINE( LU1,LBLK)
         DO JVEC = 1, NROOT
           CALL REWINE(LU6,LBLK)
           IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
         END DO
*
         CALL REWINE(LU3,LBLK)
         DO JVEC = NROOT+1, NVEC+IVEC
          CALL REWINE(LU6,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
         END DO
       END DO
*      /\ End do over new trial vectors
*. 2.3 Diagonalize projected matrix
       NVEC = NVEC + IADD
       IF(IPRT.GE.10) THEN
         WRITE(6,*) ' New Subspace matrix '
         CALL PRSYM(APROJ,NVEC)
       END IF
*
       CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
       CALL EIGENL(WORK(KAPROJ),AVEC,NVEC,0,1)
*. Test : transform projected matrix
C TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
C      CALL TRAN_SYM_BLOC_MAT(APROJ,AVEC,1,NVEC,XJEP(1000),XJEP(1))
C      WRITE(6,*) ' Explicitly transformed matrix'
C      CALL PRSYM(XJEP(1000),NVEC)

       IF(IPICO.NE.0) THEN
         E0VAR = WORK(KAPROJ)
         C0VAR = AVEC(1)
         C1VAR = AVEC(2)
         C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
         AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
         AVEC(2) = -(C1NRM/SQRT(1.0D0+C1NRM**2))
         E0PERT = AVEC(1)**2*APROJ(1)
     &          + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &          + AVEC(2)**2*APROJ(3)
         WORK(KAPROJ) = E0PERT
         IF(IPRT.GE.10) THEN
          WRITE(6,*) ' Var and Pert solution, energy and coefficients'
          WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
          WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
         END IF
       END IF
*
       IF(IROOTHOMING.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first 
*
*. Norm of wavefunction in previous space
        DO IVEC = 1, NVEC
          XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                 AVEC(1+(IVEC-1)*NROOT),NROOT)
        END DO
        IF(IPRT.GT.10) WRITE(6,*) 
     & ' Norm of projections to previous vector space '
        CALL WRTMAT(XJEP,1,NVEC,1,NVEC)
*. My sorter arranges in increasing order, multiply with minus 1
*  so the eigenvectors with largest overlap comes out first
        ONEM = -1.0D0
        CALL SCALVE(XJEP,ONEM,NVEC)
        CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
        IF(NSORT.LT.NVEC) THEN
          WRITE(6,*) ' Warning : Some elements lost in sorting '
          WRITE(6,*) ' NVEC,NSORT = ', NSORT,NVEC
        END IF
        IF(IPRT.GE.3) THEN
          WRITE(6,*) ' New roots choosen as vectors '
          CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
        END IF
*. Reorder
        DO INEW = 1, NVEC
          IOLD = IXJEP(INEW)
          CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
        END DO
        CALL COPVEC(XJEP,AVEC,NROOT*NVEC)
        DO INEW = 1, NVEC
          IOLD = IXJEP(INEW)
          XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
        END DO
        DO INEW = 1, NVEC
          WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
        END DO
*
        IF(IPRT.GE.3) THEN
          WRITE(6,*) ' Reordered WORK and AVEC arrays '
          CALL PRSYM(WORK,NVEC)
          CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
        END IF
*
       END IF
*      ^ End of root homing procedure
*
       DO IROOT = 1, NROOT
         EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
       END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F21.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(6,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(6,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
       IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,F22.13)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,NVEC,NROOT)
       END IF
*
**  Reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset      
*
*
* case 1 : Only NROOT vectors can be stored
*          save current eigenvector approximations
* Case 2 : Atleast 2*NROOT eigenvectors can be saved
*          Current eigenvactor approximations+
*          vectors allowing generation of previous approxs.
*
*
C       IF(IPRT.GE.1000) THEN 
        IF(IPRT.GE.000) THEN 
          write(6,*) ' I am going to reset '
          write(6,*) ' nroot, nvec ', nroot,nvec
        END IF
C
        IF(IPRT.GE.10000) THEN
          WRITE(6,*) ' Initial vectors on LU1'
          CALL REWINO(LU1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
          WRITE(6,*) ' Initial vectors on LU3'
          CALL REWINO(LU3)
          DO IVEC = 1, NVEC3
             CALL WRTVCD(VEC1,LU3,0,LBLK)
          END DO
        END IF
*. Orthogonalization of March 2013
*
*. First orthonormalization, based on subspace matrices
*
C       GET_CNEWCOLD_BAS(CN,CNO,NVEC,NROOT,SCR,NVECUT,THRES)
        THRES = 1.0D-7
        CALL GET_CNEWCOLD_BAS(AVEC,WORK,NVEC,NROOT,WORK2,NVECUT,
     &       THRES)
        NVEC3 = NVECUT - NROOT
        NADD = NVEC3
*. Obtain roots on LU5
        LU_FILE1(1) = LU1
        LU_FILE1(2) = LU3
        NVC_FILE1(1) = NROOT
        NVC_FILE1(2) = NVEC-NROOT
        CALL TRAVCMF(2,LU_FILE1,NVC_FILE1,NVEC,NROOT,LU5,WORK,
     &       LU6,LU7,VEC1,VEC2)
*. Obtain also correction vectors on LU8
        IF(MAXVEC.EQ.3*NROOT) THEN
          KC = NROOT*NVEC+1
          CALL TRAVCMF(2,LU_FILE1,NVC_FILE1,NVEC,NVEC3,LU8,WORK(KC),
     &         LU6,LU7,VEC1,VEC2)
          CALL COPNVCD(LU8,LU3,NVEC3,VEC1,1,LBLK)
        END IF
        CALL COPNVCD(LU5,LU1,NROOT,VEC1,1,LBLK)
*. The sigma vectors for new roots
        LU_FILE1(1) = LU2
        LU_FILE1(2) = LU4
        NVC_FILE1(1) = NROOT
        NVC_FILE1(2) = NVEC-NROOT
        CALL TRAVCMF(2,LU_FILE1,NVC_FILE1,NVEC,NROOT,LU5,WORK,
     &       LU6,LU7,VEC1,VEC2)
*. The sigma-vectors for corrections
        IF(MAXVEC.EQ.3*NROOT) THEN
          KC = NROOT*NVEC+1
          CALL TRAVCMF(2,LU_FILE1,NVC_FILE1,NVEC,NVEC3,LU8,WORK(KC),
     &         LU6,LU7,VEC1,VEC2)
          CALL COPNVCD(LU8,LU4,NVEC3,VEC1,1,LBLK)
        END IF
        CALL COPNVCD(LU5,LU2,NROOT,VEC1,1,LBLK)
        NNVEC = NROOT + NADD
*
        NVEC = NNVEC 
        WRITE(6,*) ' End of new reset '
*
        IF(IPRT.GE.1000) THEN
         WRITE(6,*) ' Updated eigenvector approximations '
         WRITE(6,*) ' ================================== '
         WRITE(6,*)
         CALL REWINO(LU1)
         DO JVEC = 1, NROOT
            CALL WRTVCD(VEC1,LU1,0,LBLK)
            WRITE(6,*)
         END DO
*
         WRITE(6,*) ' Updated Sigma for eigenvector approximations '
         WRITE(6,*) ' ============================================ '
         CALL REWINO(LU2)
         DO JVEC = 1, NADD
            CALL WRTVCD(VEC2,LU2,0,LBLK)
            WRITE(6,*)
         END DO
*
         WRITE(6,*) ' Updated Sigma for Additional vectors '
         WRITE(6,*) ' ==================================== '
         CALL REWINO(LU4)
         DO JVEC = 1, NADD
            CALL WRTVCD(VEC2,LU4,0,LBLK)
            WRITE(6,*)
         END DO
        END IF  ! Print
*
* We now have:
* NROOT current approximations to the eigenvectors and their sigmavectors on LU1 and LU3, respectively
* NVEC-NROOT additional vectors and their sigmavectors on LU2 and LU4, respectively
*
*. New subspace
*
*. Calculate subspace Hamiltonian from actual vectors 
*. on disc
        I_NEW_OR_OLD = 1
*
        IF(I_NEW_OR_OLD.EQ.2) THEN
         IF(IPRT.GE.1000) write(6,*) ' Subspace hamiltonian' 
         CALL REWINO(LU1)
         CALL REWINO(LU3)
         DO IVEC = 1, NVEC
*
          CALL REWINO(LU5)
          IF(IVEC.LE.NROOT) THEN
            CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
          ELSE
            CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
          END IF
*
          CALL REWINO(LU2)
          DO JVEC = 1, MIN(IVEC,NROOT)
            CALL REWINO(LU5)
            IJ = IVEC*(IVEC-1)/2+JVEC
            APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
          END DO
*
          CALL REWINO(LU4)
          DO JVEC = NROOT+1,IVEC
            CALL REWINO(LU5)
            IJ = IVEC*(IVEC-1)/2+JVEC
            APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
          END DO
         END DO ! Loop over IVEC
         IF(IPRT.GE.10) THEN
           write(6,*) ' Reset hamiltonian'
           call prsym(aproj,nvec)
         END IF
*. Test : Orthogonality of new vectors
         CALL REWINO(LU1)
         CALL REWINO(LU3)
         DO IVEC = 1, NVEC
*
          CALL REWINO(LU5)
          IF(IVEC.LE.NROOT) THEN
            CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
          ELSE
            CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
          END IF
*
          CALL REWINO(LU1)
          DO JVEC = 1, MIN(IVEC,NROOT)
            CALL REWINO(LU5)
            IJ = IVEC*(IVEC-1)/2+JVEC
            XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU1,0,LBLK)
          END DO
          CALL REWINO(LU3)
          DO JVEC = NROOT+1,IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU3,0,LBLK)
          END DO
         END DO ! Loop over IVEC
         IF(IPRT.GE.10) THEN
          write(6,*) ' Overlap matrix    '
          call prsym(xjep,nvec)
         END IF
        ELSE ! Switch between old and new subspace matrices
*
*. Alternative calculation of subspace Hamiltonian '
*
         NVC_FILE1(1) = NROOT
         LU_FILE1(1) = LU1
         NVC_FILE1(2) = NADD
         LU_FILE1(2) = LU3
*
         NVC_FILE2(1) = NROOT
         LU_FILE2(1) = LU2
         NVC_FILE2(2) = NADD
         LU_FILE2(2) = LU4
*
         CALL SUBSPC_MF(2,LU_FILE1,NVC_FILE1,
     &                  2,LU_FILE2,NVC_FILE2,
     &                  VEC1,VEC2,APROJ,1,LU6)
         IF(IPRT.GE.10) THEN
           write(6,*) ' Reset hamiltonian (new route)'
           call prsym(aproj,nvec)
         END IF
*. Alternative calculation of subspace Overlap '
         IF(I_CHECK_OVLAP.EQ.1) THEN
           NVC_FILE1(1) = NROOT
           LU_FILE1(1) = LU1
           NVC_FILE1(2) = NADD
           LU_FILE1(2) = LU3
           CALL SUBSPC_MF(2,LU_FILE1,NVC_FILE1,
     &                    2,LU_FILE1,NVC_FILE1,
     &                    VEC1,VEC2,XJEP,1,LU6)
C          SUBSPC_MF(NFL_L,LUFL_L,NVCFL_L,NFL_R,LUFL_R,NVCFL_R,
C    &                       VEC1,VEC2,SUBMAT,ISYM,LUSCR)
           IF(IPRT.GE.10) THEN
             write(6,*) ' Overlap matrix (alt. route) '
            call prsym(xjep,nvec)
           END IF
           CALL CHECK_UNIT_MAT(XJEP,NVEC,XMAX_DIA,XMAX_OFD,1)
         END IF ! overlap should be calculated
        END IF ! Switch between reset types
*. End of resetting business
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000 ! End of loop over iterations
 1001 CONTINUE
*
* ( End of loop over iterations )
*
C?    WRITE(6,*) ' ITER, MAXIT = ', ITER, MAXIT
      IF(CONVER.AND.ITER.LE.MAXIT) THEN
        ITER_LAST = ITER-1
      ELSE
        ITER_LAST = MAXIT
      END IF
C?    WRITE(6,*) ' ITER, MAXIT, ITER_LAST, CONVER = ',
C?   &             ITER, MAXIT, ITER_LAST, CONVER
      DO 1601 IROOT = 1, NROOT
         FINEIG(IROOT) = EIG(ITER_LAST,IROOT)+EIGSHF
         RNRM_CNV(IROOT) = RNRM(ITER_LAST,IROOT)
C?    WRITE(6,*) ' ITER, MAXIT, ITER_LAST, RNRM_CNV(1) = ', 
C?   &             ITER, MAXIT, ITER_LAST, RNRM_CNV(1)
 1601 CONTINUE
C?    WRITE(6,*) ' MICDV6(still), FINEIG(1) = ', FINEIG(1)
 
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      DO IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(ITER,IROOT) + EIGSHF
      END DO
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
C?    WRITE(6,*) ' MICDV6(c), FINEIG(1) = ', FINEIG(1)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.1000) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
CM      DO 1601 IROOT = 1, NROOT
CM         FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
CM         RNRM_CNV(IROOT) = RNRM(ITER,IROOT)
CM1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE TRAVCD(VEC1,VEC2,X,NVECIN,NVECOUT,LUIN,LUOUT,
     &                  ICOPY,LBLK,LUSCR1,LUSCR2)
*
* NVECIN vectors reside on LU1, Transform these vectors,
*  using LUSCR1 and  LUSCR2  as 
* scratch files,  with matrix X to produce output file LUOUT.
* 
* Since LUIN is accessed several times it is always 
* REWINDED. LUOUT is written to from current start.
*
* I ICOPY .ne. 0 the transformed vectors are written back to LUIN
*
* Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input 
      DIMENSION X(NVECIN,NVECOUT)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*
C             MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
      DO IVECOUT = 1, NVECOUT
         CALL MVCSMD(LUIN,X(1,IVECOUT), LUSCR1,LUSCR2,VEC1,VEC2,
     &               NVECIN,1,LBLK)
         CALL REWINO(LUSCR1)
         CALL COPVCD(LUSCR1,LUOUT,VEC1,0,LBLK)
      END DO
*
      IF(ICOPY.EQ.1) THEN
        CALL REWINO(LUIN)
        CALL REWINO(LUOUT)
        DO IVECOUT = 1, NVECOUT
          CALL COPVCD(LUOUT,LUIN,VEC1,0,LBLK)
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE MVCSMD2(LUIN,FAC,FACLUOUT,LUOUT,LUSCR,
     &           VEC1,VEC2,NVEC,IREW,LBLK)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data 
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1)
*
      IF(IREW .NE. 0 ) CALL REWINE(LUIN,LBLK)
*
      LLUOUT = LUOUT
      LLUSCR = LUSCR
      DO IVEC = 1, NVEC
        CALL REWINE(LLUSCR,LBLK)
        CALL REWINE(LLUOUT,LBLK)
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,
     &                  LUOUT,LUSCR,0,LBLK)
          ELSE
            CALL SCLVCD(LUIN,LUSCR,FAC(1),VEC1,0,LBLK)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &                0,LBLK)
        END IF
*
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
      END DO
*
      IF(LLUOUT.NE.LUOUT) THEN
C            COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
        CALL COPVCD(LLUOUT,LUOUT,VEC1,1,LBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE SEL_ROOT(SUBSPCVC,SUBSPCMT,ISEL_MET,NVEC,NROOT,LUC,
     &                    VEC1)
*
* A subspace of vectors is defined by coefficientc SUBPSVEC, SUBPSCMAT and
* Vectors on LUC.
*
*. Determine those vectors that should be the new roots and saved
*  as the first subspace expansion in SUBSP*.
*
* The method used for defining the roots are given by ISEL_MET
*
* ISEL_MET = 0 => No selection, just continue
* ISEL_MET = 1 => Root homing
* ISEL_MET = 2 => Super-symmetry
*
*. Jeppe Olsen, february 2013
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Input and output
      DIMENSION SUBSPCVC(NVEC*NVEC),SUBSPCMT(NVEC*(NVEC+1)/2)
*. Scratch
      DIMENSION VEC1(*)
*. Local scratch
      PARAMETER (MXLNVEC = 100)
      DIMENSION SUBSPCM2(MXLNVEC*MXLNVEC)
      DIMENSION SUBSPCV2(2*MXLNVEC)
      DIMENSION ISUBSPCV2(MXLNVEC)
*
      NTEST = 10
*
      IF(ISEL_MET.LT.0.OR.ISEL_MET.GT.2) THEN
          WRITE(6,*) ' Unknown ISEL_MET parameter, = ', ISEL_MET
          STOP ' SEL_ROOT: Unknown ISEL_MET parameter '
      END IF
*
      IF(NVEC.GT.MXLNVEC) THEN
        WRITE(6,*) ' Problem in SEL_ROOT: NVEC > MXLNVEC '
        WRITE(6,*) ' NVEC, MXLNVEC = ', NVEC, MXLNVEC
        WRITE(6,*) ' Increase local parameter MXLNVEC and recompile'
        STOP ' Problem in SEL_ROOT: NVEC > MXLNVEC '
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from SEL_ROOT: '
      END IF
      IF(NTEST.GE.100) THEN
        IF(ISEL_MET.EQ.0) THEN
          WRITE(6,*) ' No root selection '
        ELSE IF (ISEL_MET.EQ.1) THEN
          WRITE(6,*) ' Root homing '
        ELSE IF (ISEL_MET.EQ.2) THEN
          WRITE(6,*) ' Super-symmetry '   
        END IF
      END IF
*
      IF(ISEL_MET.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first 
*
*. Norm of wavefunction in previous space
       DO IVEC = 1, NVEC
         SUBSPCV2(IVEC) = INPROD(SUBSPCVC(1+(IVEC-1)*NVEC),
     &                    SUBSPCVC(1+(IVEC-1)*NVEC),NROOT)
       END DO
       IF(NTEST .GT.20) THEN 
         WRITE(6,*)  ' Norm of projections to previous vector space '
         CALL WRTMAT(SUBSPCV2,1,NVEC,1,NVEC)
       END IF
*. My sorter arranges in increasing order, multiply with minus 1
*  so the eigenvectors with largest overlap comes out first
       ONEM = -1.0D0
       CALL SCALVE(SUBSPCV2,ONEM,NVEC)
       CALL SORLOW(SUBSPCV2,SUBSPCV2(1+NVEC),ISUBSPCV2,
     &      NVEC,NVEC,NSORT,0)
       IF(NSORT.LT.NVEC) THEN
         WRITE(6,*) ' Warning : Some elements lost in sorting '
         WRITE(6,*) ' NVEC,NSORT = ', NSORT,NVEC
       END IF
* Resort?
       IRESORT = 0
       DO I = 1, NROOT
         IF(ISUBSPCV2(I).NE.I) IRESORT = 1
       END DO
       IF(NTEST.GE.10) THEN
        IF(IRESORT.EQ.0) THEN
         WRITE(6,*) ' No resort '
        ELSE
         WRITE(6,*) ' Resorted roots '
         WRITE(6,*) ' New roots choosen as vectors '
         CALL IWRTMA(ISUBSPCV2,1,NROOT,1,NROOT)
        END IF
       END IF
      END IF !  End of switch between different selection methods
*
* We have now obtained, by some method, the new order of the vectors
* reorder the coeffiecients and energies
*
      DO INEW = 1, NVEC
        IOLD = ISUBSPCV2(INEW)
        CALL COPVEC(SUBSPCVC(1+(IOLD-1)*NVEC),
     &       SUBSPCM2(1+(INEW-1)*NVEC),NVEC)
      END DO
      CALL COPVEC(SUBSPCM2,SUBSPCVC,NVEC**2)
      DO INEW = 1, NVEC
        IOLD = ISUBSPCV2(INEW)
        SUBSPCV2(INEW) = SUBSPCMT(IOLD*(IOLD+1)/2)
      END DO
      DO INEW = 1, NVEC
        SUBSPCMT(INEW*(INEW+1)/2) = SUBSPCV2(INEW)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reordered SUBSPCMT and SUBSPCVC arrays '
        CALL PRSYM(SUBSPCMT,NVEC)
        CALL WRTMAT(SUBSPCVC,NVEC,NVEC,NVEC,NVEC)
      END IF
*
      RETURN
      END
      SUBROUTINE ADD_ORTN_2SUBSPC(LU1,LU2,N1,N2,LUNEW,IMNEW,
     &           VEC1,VEC2,LUSC1,LUSC2,THRES)
*
*. The files of orthonormal vectors, LU1, LU2  are given
*  Orthogonalize vector on LUNEW to these and save on LUNEW
*
* On output IMNEW = 1 => New vector is on LUNEW 
*.                = 0 => No new vector
*
* Using Modified Gram Schmidt and plenty of I/O
*
*. Jeppe Olsen; Feb. 26, 2013
*
      INCLUDE 'implicit.inc'
      REAL*8 INPRDD
*. block of vector
      DIMENSION VEC1(*), VEC2(*)
*. scratch for subspace 
*
      ONE = 1.0D0
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from ADD_ORTN_2SUBSPC '
        WRITE(6,*) ' ========================= '
        WRITE(6,*)
        WRITE(6,*) ' N1, N2 = ', N1, N2
        WRITE(6,*) ' LU1, LU2, LUNEW,LUSC1,LUSC2= ',
     &               LU1, LU2, LUNEW,LUSC1,LUSC2
      END IF
*
*. Initial norm of New vector 
*
      CALL REWINO(LUNEW)
      LBLK = -1
      XNORM2 = INPRDD(VEC1,VEC2,LUNEW,LUNEW,0,LBLK)
      XNORM_INI = SQRT(ABS(XNORM2))
*
*. Orthogonalize to files on LU1
*
      LUIN = LUNEW
      LUOUT = LUSC1
*.
      CALL REWINO(LU1)
      DO I1 = 1, N1
        IF(NTEST.GE.1000) WRITE(6,*) ' I1 = ', I1
        CALL REWINO(LUSC2)
        CALL COPVCD(LU1,LUSC2,VEC1,0,LBLK)
        OVLAP =  INPRDD(VEC1,VEC2,LUIN,LUSC2,1,LBLK) 
        FACTOR = -OVLAP
        CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LUIN,LUSC2,LUOUT,1,LBLK)
C            VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        LUX = LUIN
        LUIN = LUOUT
        LUOUT = LUX
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) '  Finished with ortho to LU1 '
*
*. Orthogonalize to files on LU2
*
      CALL REWINO(LU2)
      DO I2 = 1, N2
        CALL REWINO(LUSC2)
        CALL COPVCD(LU2,LUSC2,VEC1,0,LBLK)
        OVLAP =  INPRDD(VEC1,VEC2,LUIN,LUSC2,1,LBLK) 
        FACTOR = -OVLAP
        IF(NTEST.GE.1000) WRITE(6,*) ' I2, OVLAP = ', I2, OVLAP
        CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LUIN,LUSC2,LUOUT,1,LBLK)
C            VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        LUX = LUIN
        LUIN = LUOUT
        LUOUT = LUX
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) '  Finished with ortho to LU2 '
*. Test: ortho to LU2 a second time... 
      I_DO_TEST = 0
      IF(I_DO_TEST.EQ.1) THEN
      CALL REWINO(LU2)
      DO I2 = 1, N2
        CALL REWINO(LUSC2)
        CALL COPVCD(LU2,LUSC2,VEC1,0,LBLK)
        OVLAP =  INPRDD(VEC1,VEC2,LUIN,LUSC2,1,LBLK) 
        FACTOR = -OVLAP
        IF(NTEST.GE.1000) WRITE(6,*) ' I2(2), OVLAP = ', I2, OVLAP
        CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LUIN,LUSC2,LUOUT,1,LBLK)
C            VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        LUX = LUIN
        LUIN = LUOUT
        LUOUT = LUX
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) '  Finished with ortho(2) to LU2 '
      END IF ! I_DO_TEST = 1


*
* Norm of new vector and check if it is below threshold
*
      XNORM2 = INPRDD(VEC1,VEC2,LUIN,LUIN,1,LBLK)
      XNORM_FINAL = SQRT(ABS(XNORM2))
      RATIO = XNORM_FINAL/XNORM_INI
      IF(NTEST.GE.10) 
     &  WRITE(6,*) ' Ratio between initial and final norm ', RATIO
      IF(RATIO.GT.THRES) THEN
*. New vector is more than numerical noise
        IMNEW = 1
        FACTOR = 1.0D0/XNORM_FINAL
C SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
        CALL SCLVCD(LUIN,LUOUT,FACTOR,VEC1,1,LBLK)
        IF(LUOUT.NE.LUNEW) THEN
          CALL COPVCD(LUOUT,LUNEW,VEC1,1,LBLK)
        END IF
      ELSE
        IMNEW = 0
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' IMNEW = ', IMNEW
      END IF
   
      RETURN
      END
      SUBROUTINE GET_OVLAP_FOR_VECTORS(LU1,N1,LU2,N2,LUSCR,OVLAP,
     &                                 VEC1,VEC2)
*
* Obtain overlap between vectors on files LU1, LU2
* If LU1 = LU2, then only the lower half matrix is constructed
*
*. Jeppe Olsen; Feb. 26, 2013; tracking some numerical instabilities of eigenvalue problem
*
      INCLUDE 'implicit.inc'
      REAL*8 INPRDD
*. Output
      DIMENSION OVLAP(*)
*. Scratch
      DIMENSION VEC1(*), VEC2(*)
*
      LBLK = -1 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_OVLAP_FOR_VECTORS '
        WRITE(6,*) ' ================================'
        WRITE(6,*) 
        WRITE(6,*) ' LU1, LU2 = ', LU1, LU2 
        WRITE(6,*) ' N1 and N2 = ', N1, N2
      END IF
*
      IF(NTEST.GE.1000) THEN 
        WRITE(6,*) ' Vectors on LU2 '
        CALL REWINO(LU2)
        DO IVEC = 1, N2
          CALL WRTVCD(VEC1,LU2,0,LBLK)
        END DO
      END IF
*
      IF(LU1.NE.LU2) THEN
*
* Different files 
*
        CALL REWINO(LU1)
        DO IVEC1 = 1, N1
          CALL REWINO(LUSCR)
          CALL COPVCD(LU1,LUSCR,VEC1,0,LBLK)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Vector copied to LUSCR '
            CALL WRTVCD(VEC1,LUSCR,1,LBLK)
          END IF
          CALL REWINO(LU2)
          DO IVEC2 = 1, N2
            CALL REWINO(LUSCR)
            S = INPRDD(VEC1,VEC2,LUSCR,LU2,0,LBLK)
            OVLAP((IVEC2-1)*N1+IVEC1) = S
          END DO
        END DO
      ELSE
*
* Identical files
*
        CALL REWINO(LU1)
        DO IVEC1 = 1, N1
          CALL REWINO(LUSCR)
          CALL COPVCD(LU1,LUSCR,VEC1,0,LBLK)
          OVLAP(IVEC1*(IVEC1-1)/2+IVEC1) = S11
          CALL REWINO(LU1)
          DO IVEC2 = 1, IVEC1
            CALL REWINO(LUSCR)
            S12 = INPRDD(VEC1,VEC2,LUSCR,LU1,0,LBLK)
            OVLAP(IVEC1*(IVEC1-1)/2+IVEC2) = S12
          END DO
        END DO
      END IF ! files are different
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap matrix '
        IF(LU1.NE.LU2) THEN
          CALL WRTMAT(OVLAP,N1,N2,N1,N2)
        ELSE
          CALL PRSYM(OVLAP,N1)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE CNEW_AS_COLD_PLUS_CORR(CNEW,NROOT,NVEC,DELTA)
*
* A expansion CNEW of a set of vectors in a basis is given. 
*. Make the rather trivial rewrite to the original vectors plus a set of corrections
*. No reordering this time around.
*
* CNEW = 1 + DELTA 
*
* Jeppe Olsen, improving stability of multi-root diagonalizers, Feb. 2013
* 
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION CNEW(NVEC,NROOT)
*. Output
      DIMENSION DELTA(NVEC,NROOT)
*
      CALL COPVEC(CNEW,DELTA,NVEC*NROOT)
      DO IROOT = 1, NROOT
        DELTA(IROOT,IROOT) =  DELTA(IROOT,IROOT) - 1.0D0
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DELTA = CNEW - 1 '
        CALL WRTMAT(DELTA,NVEC,NROOT,NVEC,NROOT)
      END IF
*
      RETURN
      END 
      SUBROUTINE MGSCD(LUIN,NVECIN,LUOUT,NVECUT,X,VEC1,VEC2,
     &                LUSCR1,LUSCR2,LUSCR3,THRES,XC)
*
* Modified Gram Shmidt for NVECIN vectors on LUIN to produce NVECUT 
*.orthonormal vectors on LUOUT
*
*. LUIN = LUOUT is allowed
*
*. Matrix giving Final vectors in terms of initial vectors is given in X 
*  (complete NVECIN*NVECIN matrix)
*
* LUSCR1, LUSCR2 will only contain a single vector, whereas LUSCR3 will 
* contain NVEC vectors
*
*. Jeppe Olsen, Febr. 2013; 
*
      INCLUDE 'implicit.inc'
      REAL*8 INPRDD
*. Output
      DIMENSION X(NVECIN,NVECIN)
*.  Scratch
      DIMENSION VEC1(*), VEC2(*), XC(NVECIN)
*
      NTEST = 100
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Info from MGSCD (Modified Gram-Schmidt) '
       WRITE(6,*) ' ======================================= '
       WRITE(6,*)
       WRITE(6,*) ' Number of input vectors = ', NVECIN
       WRITE(6,*)
       WRITE(6,*) ' LUIN, LUOUT = ', LUIN, LUOUT
      END IF
*
      LBLK = -1
      ONE = 1.0D0
      ZERO = 0.0D0
*
*. Initial X-matrix
*
      CALL SETVEC(X,ZERO,NVECIN*NVECIN)
C          SETDIA(MATRIX,VALUE,NDIM,IPACK)
      CALL SETDIA(X,ONE,NVECIN,0)
*. 
      CALL REWINO(LUOUT)
*
*. Initial norms of vectors
      CALL REWINO(LUIN)
      DO IVEC = 1, NVECIN
        S2 = ABS(INPRDD(VEC1,VEC2,LUIN,LUIN,0,LBLK))
        XC(IVEC) = SQRT(S2)
      END DO
*
      LLUSCR1 = LUSCR3
      LLUSCR2 = LUOUT
*. Copy all vectors to LLUSCR1
      CALL REWINO(LUIN)
      CALL REWINO(LLUSCR1)
      DO IVEC = 1, NVECIN
        CALL COPVCD(LUIN,LLUSCR1,VEC1,0,LBLK)
      END DO 
*
      NVECUT = 0
      DO IVEC = 1, NVECIN
*. Vector NVECUT + 1 at LLUSCR1 is the new (not normalized) orthogonal vector
*. Copy the first NVECUT vectors to LLUSCR2
        CALL REWINO(LLUSCR1)
        CALL REWINO(LLUSCR2)
        DO JVEC = 1, NVECUT
          CALL COPVCD(LLUSCR1,LLUSCR2,VEC1,0,LBLK)
        END DO
*. And normalize vector NVECUT + 1
        CALL REWINO(LUSCR1)
        CALL COPVCD(LLUSCR1,LUSCR1,VEC1,0,LBLK)
        S2 = ABS(INPRDD(VEC1,VEC2,LUSCR1,LUSCR1,1,LBLK))
        S = SQRT(S2)
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Norm of vector, before and after orth ',
     &   XC(IVEC),S
        END IF
        IF(S/XC(IVEC).GE.THRES) THEN
*. The vector is independent, add
          NVECUT = NVECUT + 1
          IF(NVECUT.NE.IVEC) 
     &    CALL COPVEC(X(1,IVEC),X(1,NVECUT),NVECIN)
          FACTOR = 1.0D0/S
          CALL SCLVCD(LUSCR1,LUSCR2,FACTOR,VEC1,1,LBLK)
C              SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
          CALL COPVCD(LUSCR2,LUSCR1,VEC1,1,LBLK)
          CALL REWINO(LUSCR1)
          CALL COPVCD(LUSCR1,LLUSCR2,VEC1,0,LBLK)
          CALL SCALVE(X(1,NVECUT),FACTOR,NVECIN)
*. Orthogonalize the remaining vectors on LLUSCR1 to this and save on LLUSCR2
          DO JVEC = IVEC + 1, NVECIN
            CALL REWINO(LUSCR2)
            CALL COPVCD(LLUSCR1,LUSCR2,VEC1,0,LBLK)
            OVLAP = INPRDD(VEC1,VEC2,LUSCR1,LUSCR2,1,LBLK)
            FACTOR = -OVLAP
            CALL REWINO(LUSCR1)
            CALL REWINO(LUSCR2)
            CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LUSCR2,LUSCR1,LLUSCR2,0,
     &                  LBLK)
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
C  VECSUM(C,A,B,FACA,FACB,NDIM)
            CALL VECSUM(X(1,JVEC),X(1,JVEC),X(1,NVECUT),
     &                  ONE,FACTOR,NVECIN)
          END DO
        ELSE
          WRITE(6,*) ' Vector eliminated '
*. Just copy remaining vectors from LLUSCR1 to LLUSCR2
          DO JVEC = IVEC+1, NVECIN
            CALL COPVCD(LLUSCR1,LLUSCR2,VEC1,0,LBLK)
          END DO
        END IF
*. Prepare files for next round
        LX = LLUSCR1
        LLUSCR1 = LLUSCR2
        LLUSCR2 = LX
      END DO
*
      IF(LLUSCR1.NE.LUOUT) THEN
        CALL REWINO(LLUSCR1)
        CALL REWINO(LUOUT)
        DO IVEC = 1, NVECUT
            CALL COPVCD(LLUSCR1,LUOUT,VEC1,0,LBLK)
        END DO
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Number of vectors obtained ', NVECUT
        WRITE(6,*) 
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformation matrix '
        CALL WRTMAT(X,NVECIN,NVECUT,NVECIN,NVECUT)
      END IF

*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' And the orthonormal vectors '
        WRITE(6,*) ' --------------------------- '
        WRITE(6,*)
        CALL REWINO(LUOUT)
        DO JVEC = 1, NVECUT
          CALL WRTVCD(VEC1,LUOUT,0,LBLK)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE ORTVECTORS_TO_ORTNVECTORS(LU1,NVEC1,LU2,NVEC2,
     &           THRES,VEC1,VEC2,LU3,NVEC3,LUSC1,LUSC2,LUSC3,X,XC,IORIG)
*
* A set of NVEC1 vectors are given at LU1.
* Orthonormalize these vectors to the  NVEC2 orthonormal vectors on LU2 and 
* save the NVEC3 independent vectors on LU3.
* X is the transformation matrix and IORIG(IVEC) is the original vector on LU1
* connected with independent vector IVEC
*
* LUSC1, LUSC2 are single vector files, LUSC3 is multiple vector file
*
* LUSC1, LUSC2 may be equal to LUOUT
*
*. Jeppe Olsen; Feb. 2013; Improving and stabilizing diagonalers
*
      INCLUDE 'implicit.inc'
      REAL*8 INPRDD
*. Output
      DIMENSION X(NVEC2,NVEC1),XC(NVEC1), IORIG(NVEC1)
*. Scratch
      DIMENSION VEC1(*), VEC2(*)
*
      LBLK = -1
      ONE = 1.0D0
*
      NTEST = 100
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from ORTVECTORS_TO_ORTNVECTORS'
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
        WRITE(6,*) ' Input files, LU1, LU2 = ', LU1, LU2
        WRITE(6,*) ' Output file, LU3 = ', LU3
        WRITE(6,*) ' Scratch files, LUSC1, LUSC2, LUSC3 = ',
     &                              LUSC1, LUSC2, LUSC3
        WRITE(6,*) ' Number of vectors on LU1, LU2 = ', NVEC1, NVEC2
        WRITE(6,*) ' Threshold ', THRES
      END IF
*
      CALL REWINO(LU3)
      CALL REWINO(LU1)
      CALL REWINO(LU2)
*
*. Overlap between vectors on LU1 and LU2
*
      DO IVEC1 = 1, NVEC1
C COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
        CALL REWINO(LUSC1)
        CALL COPVCD(LU1,LUSC1,VEC1,0,LBLK)
        S2 = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
        XC(IVEC1) = SQRT(ABS(S2))
        CALL REWINO(LU2)
        DO IVEC2 = 1, NVEC2
          CALL REWINO(LUSC1)
          OVLAP = INPRDD(VEC1,VEC2,LUSC1,LU2,0,LBLK)
          X(IVEC2,IVEC1) = -OVLAP
        END DO
      END DO
      WRITE(6,*) ' Overlap between LU1 and LU2 determined '
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthogonalization matrix 1 '
        CALL WRTMAT(X,NVEC2,NVEC1,NVEC2,NVEC1)
        WRITE(6,*) ' Initial norm of vectors on LU1 '
        CALL WRTMAT(XC,NVEC1,1,NVEC1,1)
      END IF
*
*. Check for singularities using X
*
      NVEC3 = 0
      DO IVEC3 = 1, NVEC1
       XNORM2 = XC(IVEC3)**2
       DO IVEC2 = 1, NVEC2
         XNORM2 = XNORM2 - X(IVEC2,IVEC3)**2
       END DO
       XNORM = SQRT(ABS(XNORM2))
       IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Norm of vector, before and after orth.',
     &   XC(IVEC3),XNORM
       END IF
       IF(XNORM/XC(IVEC3).GE.THRES) THEN
*. Independent vector
         NVEC3 = NVEC3 + 1
         IORIG(NVEC3) = IVEC3
         IF(NVEC3.NE.IVEC3) THEN
           CALL COPVEC(X(1,IVEC3),X(1,NVEC3),NVEC2)
         END IF
       ELSE
         WRITE(6,*) ' Vector eliminated '
       END IF
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthogonalization matrix 2 '
        CALL WRTMAT(X,NVEC2,NVEC3,NVEC2,NVEC3)
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' The IORIG array '
        CALL IWRTMA(IORIG,1,NVEC3,1,NVEC3)
      END IF
*
* Transform the vectors to the basis orthogonal
*
C     TRAVCD(VEC1,VEC2,X,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
      CALL REWINO(LUSC3)
      CALL TRAVCD(VEC1,VEC2,X,NVEC2,NVEC3,LU2,LUSC3,0,
     &            LBLK,LUSC1,LUSC2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Vectors after TRAVCD '
        CALL REWINO(LUSC3)
        DO IVEC3 = 1, NVEC3
          CALL WRTVCD(VEC1,LUSC3,0,LBLK)
        END DO
      END IF
      CALL REWINO(LU1)
      CALL REWINO(LUSC3)
      CALL REWINO(LU3)
      IVEC3 = 1
      DO IVEC1 = 1, NVEC1
        IF(IORIG(IVEC3).EQ.IVEC1) THEN
          CALL VECSMD(VEC1,VEC2,ONE,ONE,LU1,LUSC3,LU3,0,LBLK)
          IVEC3 = IVEC3 + 1
        ELSE
C              SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
          CALL SKPVCD(LU1,1,VEC1,0,LBLK)
        END IF
      END DO
CE    DO IVEC3 = 1, NVEC3
CE      CALL VECSMD(VEC1,VEC2,ONE,ONE,LU1,LUSC3,LU3,0,LBLK)
CE    END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The X(IVEC2,IVEC2) transformation matrix '
        WRITE(6,*) ' ========================================='
        WRITE(6,*)
        CALL WRTMAT(X,NVEC2,NVEC3,NVEC2,NVEC3)
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The orthogonalized vectors '
        CALL REWINO(LU3)
        DO IVEC3 = 1, NVEC3
          CALL WRTVCD(VEC1,LU3,0,LBLK)
        END DO
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Leaving ORT....  '
      END IF
*
      RETURN
      END
      SUBROUTINE TRAVCMF(NFILE,LU_FILE,NVC_FILE,NIN,NOUT,LUOUT,X,
     &           LUSC1, LUSC2,VEC1,VEC2)
*
* A set of files is defined by NFILE with NVCFILE(IFILE) vectors on FILE IFILE
* Perform the transformation indicated by transformation matrix X.
* NIN is total number of vectors in the input files.
*
*. Jeppe Olsen, March 1, Cleaning Diagonalization codes
*
      INCLUDE 'implicit.inc'
*. FIles
      INTEGER LU_FILE(NFILE),NVC_FILE(NFILE)
*. Transformation matrix
      DIMENSION X(NIN,NOUT)
*. Scratch
      DIMENSION VEC1(*), VEC2(*)
*
      LLUSC1 = LUSC1
      LLUSC2 = LUSC2
*
      ONE = 1.0D0
      LBLK = -1
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from TRAVCMF '
        WRITE(6,*) ' =================='
      END IF
*
      CALL REWINO(LUOUT)
      DO IVECOUT = 1, NOUT
        IVECIN = 0
        DO IIFILE = 1, NFILE
          IFILE = LU_FILE(IIFILE)
          CALL REWINO(IFILE)
          DO IIVECIN = 1, NVC_FILE(IIFILE)
            IVECIN = IVECIN + 1
            FACTOR = X(IVECIN,IVECOUT)
            IF(IVECIN.EQ.1) THEN
C SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
              CALL REWINO(LLUSC2)
              CALL SCLVCD(IFILE,LLUSC2,FACTOR,VEC1,0,LBLK)
            ELSE
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
              CALL REWINO(LLUSC1)
              CALL REWINO(LLUSC2)
              CALL VECSMD(VEC1,VEC2,FACTOR,ONE,IFILE,LLUSC1,LLUSC2,
     &             0,LBLK)
            END IF
*
            LUX = LLUSC1
            LLUSC1 = LLUSC2
            LLUSC2 = LUX
          END DO
        END DO
*. And transfer to final destination
        CALL REWINO(LLUSC1)
        CALL COPVCD(LLUSC1,LUOUT,VEC1,0,LBLK)
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The transformed files on LUOUT '
        WRITE(6,*) ' ==============================='
        WRITE(6,*) 
        CALL REWINO(LUOUT)
        DO IVEC = 1, NOUT
          CALL WRTVCD(VEC1,LUOUT,0,LBLK)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE COPNVCD(LUIN,LUOUT,NVEC,SEGMNT,IREW,LBLK)
*
* Copy NVEC vectors from LUIN to LUOUT
*
*. Jeppe Olsen, Feb. 2013
*
      INCLUDE 'implicit.inc'
*. Scratch
      DIMENSION SEGMNT(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Entering COPNVCD '
      END IF
*
      IF(IREW.EQ.1) THEN
        CALL REWINO(LUIN)
        CALL REWINO(LUOUT)
      END IF 
*
      DO IVEC = 1, NVEC
        CALL COPVCD(LUIN,LUOUT,SEGMNT,0,LBLK)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Leaving COPNVCD '
      END IF
*
      RETURN
      END
      SUBROUTINE SUBSPC_MF(NFL_L,LUFL_L,NVCFL_L,NFL_R,LUFL_R,NVCFL_R,
     &                     VEC1,VEC2,SUBMAT,ISYM,LUSCR)
*
* Obtain subspace matrix for vectors residing on multiple files
*
*. Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
      REAL*8 INPRDD
*. Input
      DIMENSION LUFL_L(NFL_L),NVCFL_L(NFL_L)
      DIMENSION LUFL_R(NFL_R),NVCFL_R(NFL_R)
*. Output
      DIMENSION SUBMAT(*)
* Scratch
      DIMENSION VEC1(*),VEC2(*)
*
      LBLK = -1
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from SUBSPC_MF '
        WRITE(6,*) ' ===================='
        WRITE(6,*)
        WRITE(6,*) ' LUSCR = ', LUSCR
      END IF
*
      NTVC_L = 0
      DO IFL_L = 1, NFL_L
        NTVC_L = NTVC_L + NVCFL_L(IFL_L)
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NTVC_L = ',NTVC_L
*
      IVC_L = 0
      DO IFL_L = 1, NFL_L
       LU_L = LUFL_L(IFL_L)
C?     WRITE(6,*) ' IFL_L, LU_L = ', IFL_L, LU_L
       CALL REWINO(LU_L)
       IF(ISYM.EQ.0) THEN
         MXFL_R = NFL_R
       ELSE
         MXFL_R = IFL_L
       END IF
       DO IIVC_L = 1, NVCFL_L(IFL_L)
        IVC_L = IVC_L + 1
C?      WRITE(6,*) ' IVC_L = ', IVC_L
C?      WRITE(6,*) ' LUSCR = ', LUSCR
        CALL REWINO(LUSCR)
        CALL COPVCD(LU_L,LUSCR,VEC1,0,LBLK)
C            COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK
        IVC_R = 0
        DO IFL_R = 1, MXFL_R
          LU_R = LUFL_R(IFL_R)
C?        WRITE(6,*) ' IFL_R, LU_R = ', IFL_R, LU_R
          CALL REWINO(LU_R)
          IF(ISYM.EQ.1.AND.IFL_R.EQ.IFL_L) THEN
           MXVC_R = IIVC_L
          ELSE
           MXVC_R = NVCFL_R(IFL_R)
          END IF
          DO IIVC_R = 1, MXVC_R
            IVC_R = IVC_R + 1
C?          WRITE(6,*) ' IIVC_R, IVC_R = ', IIVC_R, IVC_R
            IF(ISYM.EQ.0) THEN
              IJ = (IVC_R-1)*NTVC_L + IVC_L
            ELSE
              IJ = IVC_L*(IVC_L-1)/2 + IVC_R
            END IF
C?          WRITE(6,*) ' IVC_L, IVC_R, IJ = ', IVC_L, IVC_R, IJ 
            CALL REWINO(LUSCR)
            SIJ = INPRDD(VEC1,VEC2,LU_R,LUSCR,0,LBLK)
C                 INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
            SUBMAT(IJ) = SIJ
          END DO
        END DO! End of right files
       END DO
      END DO! End of left files
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Subspace matrix '
       IF(ISYM.EQ.0) THEN
         CALL WRTMAT(SUBMAT,IVC_L,IVC_R,IVC_L,IVC_R)
       ELSE
         CALL PRSYM(SUBMAT,IVC_L)
       END IF
      END IF
*

      END
      SUBROUTINE MGSCD_MF(NFL,LUFL_IN,NVCFL_IN,LUFL_OUT,NVCFL_OUT,X,
     &           NTVEC_IN,VEC1,VEC2,LUSCR1,LUSCR2,LUSCR3,THRES,XC)
*
* Modified Gram Shmidt for vectors on multiple files 
* Files on LUFL_IN(IFL) is after orthonormalization stored on LUFL_OUT(IFL)
*
*. Intermediate files LUSCR1, LUSCR2 will be full files.
*
*. Matrix giving Final vectors in terms of initial vectors is given in X 
*  (complete NTVEC_IN,* matrix)
*
*
*. Jeppe Olsen, Febr. 2013; 
*
      INCLUDE 'implicit.inc'
      REAL*8 INPRDD
*. Input
       INTEGER LUFL_IN(NFL),NVCFL_IN(NFL),NVCFL_OUT(NFL)
*. Output
      DIMENSION X(NTVEC_IN,*)
*.  Scratch
      DIMENSION VEC1(*), VEC2(*), XC(NTVEC_IN)
*
      NTEST = 100
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Info from MGSCD_MF (Modified Gram-Schmidt) '
       WRITE(6,*) ' ========================================== '
       WRITE(6,*)
      END IF
*
      LBLK = -1
      ONE = 1.0D0
      ZERO = 0.0D0
*
*
*. Initial X-matrix
*
      CALL SETVEC(X,ZERO,NTVEC_IN*NTVEC_IN)
*
*. Initial norms of vectors
      IVEC = 0
      DO IFL = 1, NFL
        LUIN = LUFL_IN(IFL)
        NVECIN  = NVCFL_IN(IFL)
        CALL REWINO(LUIN)
        DO IIVEC = 1, NVECIN
          IVEC = IVEC + 1
          S2 = ABS(INPRDD(VEC1,VEC2,LUIN,LUIN,0,LBLK))
          XC(IVEC) = SQRT(S2)
        END DO
      END DO
*
      LLUSCR1 = LUSCR3
      LLUSCR2 = LUOUT
*. Copy all vectors to LLUSCR1
      CALL REWINO(LUIN)
      CALL REWINO(LLUSCR1)
      DO IVEC = 1, NVECIN
        CALL COPVCD(LUIN,LLUSCR1,VEC1,0,LBLK)
      END DO 
*
      NVECUT = 0
      DO IVEC = 1, NVECIN
*. Vector NVECUT + 1 at LLUSCR1 is the new (not normalized) orthogonal vector
*. Copy the first NVECUT vectors to LLUSCR2
        CALL REWINO(LLUSCR1)
        CALL REWINO(LLUSCR2)
        DO JVEC = 1, NVECUT
          CALL COPVCD(LLUSCR1,LLUSCR2,VEC1,0,LBLK)
        END DO
*. And normalize vector NVECUT + 1
        CALL REWINO(LUSCR1)
        CALL COPVCD(LLUSCR1,LUSCR1,VEC1,0,LBLK)
        S2 = ABS(INPRDD(VEC1,VEC2,LUSCR1,LUSCR1,1,LBLK))
        S = SQRT(S2)
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Norm of vector, before and after orth ',
     &   XC(IVEC),S
        END IF
        IF(S/XC(IVEC).GE.THRES) THEN
*. The vector is independent, add
          NVECUT = NVECUT + 1
          IF(NVECUT.NE.IVEC) 
     &    CALL COPVEC(X(1,IVEC),X(1,NVECUT),NVECIN)
          FACTOR = 1.0D0/S
          CALL SCLVCD(LUSCR1,LUSCR2,FACTOR,VEC1,1,LBLK)
C              SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
          CALL COPVCD(LUSCR2,LUSCR1,VEC1,1,LBLK)
          CALL REWINO(LUSCR1)
          CALL COPVCD(LUSCR1,LLUSCR2,VEC1,0,LBLK)
          CALL SCALVE(X(1,NVECUT),FACTOR,NVECIN)
*. Orthogonalize the remaining vectors on LLUSCR1 to this and save on LLUSCR2
          DO JVEC = IVEC + 1, NVECIN
            CALL REWINO(LUSCR2)
            CALL COPVCD(LLUSCR1,LUSCR2,VEC1,0,LBLK)
            OVLAP = INPRDD(VEC1,VEC2,LUSCR1,LUSCR2,1,LBLK)
            FACTOR = -OVLAP
            CALL REWINO(LUSCR1)
            CALL REWINO(LUSCR2)
            CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LUSCR2,LUSCR1,LLUSCR2,0,
     &                  LBLK)
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
C  VECSUM(C,A,B,FACA,FACB,NDIM)
            CALL VECSUM(X(1,JVEC),X(1,JVEC),X(1,NVECUT),
     &                  ONE,FACTOR,NVECIN)
          END DO
        ELSE
          WRITE(6,*) ' Vector eliminated '
*. Just copy remaining vectors from LLUSCR1 to LLUSCR2
          DO JVEC = IVEC+1, NVECIN
            CALL COPVCD(LLUSCR1,LLUSCR2,VEC1,0,LBLK)
          END DO
        END IF
*. Prepare files for next round
        LX = LLUSCR1
        LLUSCR1 = LLUSCR2
        LLUSCR2 = LX
      END DO
*
      IF(LLUSCR1.NE.LUOUT) THEN
        CALL REWINO(LLUSCR1)
        CALL REWINO(LUOUT)
        DO IVEC = 1, NVECUT
            CALL COPVCD(LLUSCR1,LUOUT,VEC1,0,LBLK)
        END DO
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Number of vectors obtained ', NVECUT
        WRITE(6,*) 
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformation matrix '
        CALL WRTMAT(X,NVECIN,NVECUT,NVECIN,NVECUT)
      END IF

*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' And the orthonormal vectors '
        WRITE(6,*) ' --------------------------- '
        WRITE(6,*)
        CALL REWINO(LUOUT)
        DO JVEC = 1, NVECUT
          CALL WRTVCD(VEC1,LUOUT,0,LBLK)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE MERGE_VECFILES(NFL_IN,LUFL_IN, NVCFL_IN,LUOUT,
     &           IREW,LBLK,VEC1)
*
* A set of NFL_IN files LUFL_IN with NVCFL_IN vectors per file is given
* merge these vectors to a single file
*
*. Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LUFL_IN(NFL_IN),NVCFL_IN(NFL_IN)
*. Scratch
      DIMENSION VEC1(*)
*
      IF(IREW.EQ.1) THEN
        CALL REWINO(LUOUT)
      END IF
*
      DO IFL = 1, NFL_IN
        LUIN = LUFL_IN(IFL)
        NVEC = NVCFL_IN(IFL)
        CALL REWINO(LUIN)
        DO IVEC = 1, NVEC
C              COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
          CALL COPVCD(LUIN,LUOUT,VEC1,0,LBLK)
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE SPLIT_VECFILES(LUIN,NFL_OUT,LUFL_OUT, NVCFL_OUT,
     &           IREW,LBLK,VEC1)
*
* Split a single vector file to NFL_OUT files defined by LUFL_OUT, NVCFL_OUT
*
*. Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LUFL_OUT(NFL_OUT),NVCFL_OUT(NFL_OUT)
*. Scratch
      DIMENSION VEC1(*)
*
      IF(IREW.EQ.1) THEN
        CALL REWINO(LUIN)
      END IF
*
      DO IFL = 1, NFL_OUT
        LUOUT = LUFL_OUT(IFL)
        NVEC = NVCFL_OUT(IFL)
        CALL REWINO(LUOUT)
        DO IVEC = 1, NVEC
C              COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
          CALL COPVCD(LUIN,LUOUT,VEC1,0,LBLK)
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GET_CNEWCOLD_BAS(CN,CNO,NVEC,NROOT,SCR,NVECUT,
     &           THRES)
*
* A subspace expansion is given by CNEW
* Obtain orthogonal expansion of basis for New + old vectors
*
*. Jeppe Olsen, March 2013- reducing orthogonalization errors
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION CN(NVEC,NROOT)
*. Output
      DIMENSION CNO(NVEC,*)
*. Scratch: Min size: MAX(4*2*NROOT*NVEC,4*(2*NROOT)**2)
      DIMENSION SCR(3*NVEC**2)
*. Partitioning of SCR
      KS = 1
      KFREE = KS + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KX = KFREE
      KFREE = KX + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KSCR = KFREE
      KFREE = KSCR + (2*NROOT)*MAX(2*NROOT,NVEC)
*
      KVEC = KFREE
      KFREE = KVEC + 2*NROOT
*
      ZERO = 0.0D0
      ONE = 1.0D0
*
      NTEST = 1000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GET_CNEWCOLD_BAS '
        WRITE(6,*) ' ==========================='
      END IF
*. CNO: First Nroot vectors are the new vectors, the last NROOT original first roots
      ZERO = 0.0D0
      CALL SETVEC(CNO,ZERO,2*NROOT*NVEC)
      CALL COPVEC(CN,CNO,NROOT*NVEC)
*. And the original first roots
      DO IROOT = 1, NROOT
        CNO(IROOT,NROOT+IROOT) = 1.0D0
      END DO
      NDO = 2*NROOT
*
* Iterate over orthogonalizations
*
      MAXIT = 2
      DO ITER = 1, MAXIT
        IF(NTEST.GE.10) WRITE(6,*) ' Info from orth. iter. ', ITER
        IF(ITER.GT.1) THEN
          NDO = NVECUT
        END IF
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The CNO matrix '
          CALL WRTMAT(CNO,NVEC,NDO,NVEC,NDO)
        END IF
*. Overlap matrix S = CNO^T CNO
C         MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                    NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
        CALL MATML7(SCR(KS),CNO,CNO,NDO,NDO,NVEC,NDO,NVEC,NDO,ZERO,ONE,
     &              1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Overlap of NDO basis '
          CALL WRTMAT(SCR(KS),NDO,NDO,NDO,NDO)
        END IF
*. Perform Gram-Schmidt orthogonalization
C       MGS4(X,S,NDIM,SCR1,THRES,NVECUT)
        CALL MGS4(SCR(KX),SCR(KS),NDO,SCR(KVEC),THRES,NVECUT)
C   TRNMAD(A,X,SCR,NDIMI,NDIMO)
*.   Check that delivered X orthonormalizes
        CALL TRNMAD(SCR(KS),SCR(KX),SCR(KSCR),NDO,NVECUT)
        IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' New overlap matrix '
           CALL WRTMAT(SCR(KS),NVECUT,NVECUT,NVECUT)
        END IF
        CALL CHECK_UNIT_MAT(SCR(KS),NVECUT,XMAX_DIA,XMAX_OFD,0)
*. Transformation from original to orthonormal basis
        CALL MATML7(SCR(KS),CNO,SCR(KX),NVEC,NVECUT,
     &              NVEC,NDO,NDO,NVECUT,ZERO,ONE,0)
        CALL COPVEC(SCR(KS),CNO,NVEC*NVECUT)
      END DO
*
      RETURN
      END
    


      




c $Id$
