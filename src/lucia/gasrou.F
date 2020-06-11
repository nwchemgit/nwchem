*
* Some GAS routines
*
* Nomenclature
*
*. Group of strings : set of strings with a given number of orbitals
*                    in a given GASspace
*
*. Supergroup of strings  : product of NGAS groups, i.e. a string with a
*                    given numb er of electrons in each GAS space
*
*. Type of string : Type is defined by the total number of electrons
*                   in the string. A type will therefore in general
*                   consists of several supergroups
*
      SUBROUTINE SPGP_AC(INSPGRP,NINSPGRP,IOUTSPGRP,NOUTSPGRP,
     &           NGAS,MXPNGAS,IAC,ISPGRP_AC,IBASEIN,IBASEOUT)
*
* Annihilation/Creation mapping of strings
*
* Jeppe Olsen, April 1, 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input : Number of electrons in each gasspace
      INTEGER INSPGRP(MXPNGAS,*),IOUTSPGRP(MXPNGAS,*)
*. Output
      INTEGER ISPGRP_AC(NGAS,*)
*. Check first that supergroups + IAC information is consistent
      NELIN = 0
      NELOUT = 0
      DO IGAS = 1, NGAS
        NELIN  = NELIN + INSPGRP(IGAS,IBASEIN)
        NELOUT = NELOUT + IOUTSPGRP(IGAS,IBASEOUT)
      END DO
C?    WRITE(6,*) ' JEPTEST: NELIN, NELOUT = ', NELIN, NELOUT
      IF(.NOT.((IAC.EQ.1.AND.NELIN.EQ.NELOUT+1).OR.
     &         (IAC.EQ.2.AND.NELIN.EQ.NELOUT-1))) THEN
        WRITE(6,*) ' Inconsistent data provided to SPGRP_AC'
        WRITE(6,*) ' NELIN NELOUT IAC=',NELIN,NELOUT,IAC
        WRITE(6,*) ' IBASEIN, IBASEOUT = ', IBASEIN, IBASEOUT 
        WRITE(6,*) ' INSPGRP: ', (INSPGRP(I,IBASEIN),I= 1,NGAS)
        WRITE(6,*) ' IOUTSPGRP: ', (IOUTSPGRP(I,IBASEOUT),I= 1,NGAS)
        STOP' Inconsistent data provided to SPGRP_AC'
      END IF
*
      DO ISPGRP = IBASEIN, IBASEIN+NINSPGRP-1
        DO IGAS = 1, NGAS
          IF(IAC.EQ.1) THEN
            INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) - 1
          ELSE IF (IAC.EQ.2) THEN
             INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) + 1
          END IF
          IF(INSPGRP(IGAS,ISPGRP).EQ.-1) THEN
*. Trivial zero
           ISPGRP_AC(IGAS,ISPGRP) = 0
          ELSE 
*. Nonzero or nontrivial zero:
*. Find corresponding supergroup 
           I_NEW_OR_OLD = 1
           IF(I_NEW_OR_OLD.EQ.2) THEN
*. Explicit check of all arrays until match
            ITO = 0    
            DO JSPGRP = IBASEOUT, IBASEOUT+NOUTSPGRP-1
              IAMOKAY = 1
              DO JGAS = 1, NGAS
                IF( INSPGRP(JGAS,ISPGRP).NE.IOUTSPGRP(JGAS,JSPGRP))
     &          IAMOKAY=0
              END DO
              IF(IAMOKAY.EQ.1) THEN
                ITO = JSPGRP              
                GOTO 3006
              END IF
            END DO
 3006       CONTINUE
            ISPGRP_AC(IGAS,ISPGRP) = ITO    
           ELSE
C     FIND_INTARR_IN_ORD_INTARR_NLIST(INTARR,INTARR_LIST,
C                NELMNT,NELMNT_MAX,NARR,IREO,IMET,INUM)
*. Using bisection
            IMET = 2
            CALL FIND_INTARR_IN_ORD_INTARR_NLIST(
     &           INSPGRP(1,ISPGRP),IOUTSPGRP(1,IBASEOUT),NGAS,MXPNGAS,
     &           NOUTSPGRP,0,IMET,INUM)
            IF(INUM.EQ.0) THEN
              ISPGRP_AC(IGAS,ISPGRP) = 0
            ELSE 
             ISPGRP_AC(IGAS,ISPGRP) = INUM + IBASEOUT - 1
            END IF !inum = 0
           END IF !switch between old and new method
          END IF !switch betwen trivial zero or not
*. And clean up
          IF(IAC.EQ.1) THEN
            INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) + 1
          ELSE IF (IAC.EQ.2) THEN
             INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) - 1
          END IF
        END DO
      END DO
*
      NTEST = 000
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input supergroups '
        CALL IWRTMA(INSPGRP(1,IBASEIN),NGAS,NINSPGRP,
     &  MXPNGAS,NINSPGRP)
        WRITE(6,*) ' Output supergroups '
        CALL IWRTMA(IOUTSPGRP(1,IBASEOUT),NGAS,NOUTSPGRP,
     &  MXPNGAS,NOUTSPGRP)
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from SPGP_AC '
       WRITE(6,*) ' ===================='
       WRITE(6,*)
       WRITE(6,'(A, 2I5)') 
     & ' First and last IN supergroup: ',IBASEIN, IBASEIN+NINSPGRP-1
       WRITE(6,*) ' IAC = ', IAC
       WRITE(6,*) ' Mapping : '
       CALL IWRTMA(ISPGRP_AC(1,IBASEIN),NGAS,NINSPGRP,NGAS,NINSPGRP)
      END IF
*
      RETURN
      END
      SUBROUTINE ADD_STR_GROUP(NSTADD,IOFADD,ISTADD,NSTB,NSTA,
     &                         ISTRING,IELOF,NELADD,NELTOT)
*
* Part of assembling strings in individual types to
* super group of strings
*
*. Copying strings belonging to a given type to supergroup of strings
*
* Jeppe Olsen, for once improving performance of LUCIA
*
*.Input
* =====
* NSTADD : Number of strings to be added
* IOFADD : First string to be added
* ISTADD : Strings to be added
* NSTB   : Number of strings belonging to lower gasspaces
* NSTA   : Number of strings belonging to higher gasspaces
* ISTRING: Supergroup of strings under construction
* IELOF  : Place of first electron to be added
* NELADD : Number of electrons to be added
* NELTOT : Total number of electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input 
C     DIMENSION ISTADD(NELADD,*)
      DIMENSION ISTADD(*)
*. Input and output
C     DIMENSION ISTRING(NELTOT,*)
      DIMENSION ISTRING(*)
* 
C?    WRITE(6,*) '  Inside ADD ... '
      IF(NSTA.GT.1) THEN
        DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
          IOFFY = (IOFADD-2+IISTR)*NELADD
          IOFF1 = (IISTR-1)*NSTA + 1
          IADD2 = NSTADD*NSTA
          IOFF2 = IOFF1 - IADD2
          DO ISTB = 1, NSTB
*. Address of A(1,IISTR,ISTB)
C           IOFF2 = IOFF1 + (ISTB-1)*NSTADD*NSTA
            IOFF2 = IOFF2 + IADD2
            IOFFX = IELOF-1+(IOFF2-2)*NELTOT
            DO ISTA = 1, NSTA
              IOFFX = IOFFX + NELTOT
              DO JEL = 1, NELADD
                ISTRING(JEL+IOFFX) 
     &        = ISTADD(JEL+IOFFY)                          
C               ISTRING(IELOF-1+JEL,IOFF2-1+ISTA) 
C    &        = ISTADD(JEL,IOFADD-1+IISTR)
              END DO
            END DO
          END DO
        END DO 
      ELSE IF (NSTA .EQ. 1 ) THEN
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
        DO ISTB = 1, NSTB
          IOFF0 = (ISTB-1)*NSTADD
          IOFFY  = (IOFADD-2)*NELADD
          IOFFX = IELOF-1+(IOFF0-1)*NELTOT
          DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,ISTB)
C           IOFF2 = IISTR  + IOFF0             
C           IOFFX  = IELOF-1+(IOFF2-1)*NELTOT
            IOFFX  = IOFFX + NELTOT
            IOFFY  = IOFFY + NELADD            
            DO JEL = 1, NELADD
              ISTRING(JEL+IOFFX) 
     &      = ISTADD( JEL+IOFFY)
C             ISTRING(IELOF-1+JEL+(IOFF2-1)*NELTOT) 
C    &      = ISTADD(JEL+(IOFADD-1+IISTR-1)*NELADD)
C             ISTRING(IELOF-1+JEL,IOFF2) 
C    &      = ISTADD(JEL,IOFADD-1+IISTR)
            END DO
C?          WRITE(6,*) ' New string from ADD '
C?          CALL IWRTMA(ISTRING(IOFFX+1),1,NELADD,1,NELADD)
C?          WRITE(6,*) ' IOFFX, IOFFY, NELADD',
C?   &                   IOFFX, IOFFY, NELADD
          END DO
        END DO 
      END IF
*
      RETURN
      END
      SUBROUTINE ZNELFSPGP(NTESTG)
*
* Generate for each supergroup the number of electrons in each active 
* orbital space and store in NELFSPGP
*
* Jeppe Olsen, July 1995
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. input  
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
c..dongxia there's no WORK and MEMMAN involved.
      INCLUDE 'strbas.inc'
      INCLUDE 'cgas.inc'
*. Input and Output ( NELFSPGP(MXPNGAS,MXPSTT) )
      INCLUDE 'gasstr.inc'
#include "global.fh"
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      DO ITP = 1, NSTTP
        NSPGP = NSPGPFTP(ITP)
        IBSPGP = IBSPGPFTP(ITP) 
        DO ISPGP = IBSPGP,IBSPGP + NSPGP - 1
          DO IGAS = 1, NGAS
            NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        DO ITP = 1, NSTTP
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE ZSPGPIB(NSTSO,ISTSO,NSPGP,NSMST)
*
* Offset for supergroups of strings with given sym.
*. Each supergroup is given start address 1
*
* Jeppe Olsen, Still summer of 95
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER NSTSO(NSMST,NSPGP)
*. Output
      INTEGER ISTSO(NSMST,NSPGP)
*
      DO ISPGP = 1, NSPGP
        ISTSO(1,ISPGP) = 1
        DO ISMST = 2, NSMST
          ISTSO(ISMST,ISPGP) = ISTSO(ISMST-1,ISPGP) + NSTSO(ISMST,ISPGP) 
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ZSPGPIB '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        CALL IWRTMA(ISTSO,NSMST,NSPGP,NSMST,NSPGP)
      END IF
*
      RETURN
      END 
      SUBROUTINE ICOPMT(MATI,NRI,NCI,MATO,NRO,NCO)
*
* Copy integer matrix MATI to MATO 
*
*. Input
      INTEGER  MATI(NRI,NCI)
*. Output
      INTEGER MATO(NRO,NCO)
*
      NREFF = MIN(NRI,NRO)
      NCEFF = MIN(NCI,NCO)
*
      DO IC = 1, NCEFF
        DO IR = 1, NREFF
          MATO(IR,IC) = MATI(IR,IC)
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE NSTPTP_GAS_NEW(NGAS,ISPGRP,NSTSGP,NSMST,
     &                      NSTSSPGP,IGRP,MXNSTR,
     &                      NSMCLS,NSMCLSE,NSMCLSE1,NSTR_AS)
*
* Find number of strings per symmetry for the supergroup defined
* by the groups of ISPGRP. The obtained number of strings per sym 
* is stored in NSTSSPGP(*,IGRP)
*
* Jeppe Olsen, Winter 2011 - old version too slow for many gaspaces
*                            (new version simpler and quicker)
* Last Modification; Jeppe Olsen; Apr. 29, 2013; NSTR_AS added
* 
*. Also delivered:
*
* NSMCLS : MAX Number of symmetry classes for given supergroup,
*          i.e. number of combinations of symmetries of groups
*          containing strings
* NSMCLSE : Number of symmetry classes for given supergroup 
*          obtained by restricting allowed symmetries in 
*          a given group by a max and min.
* NSMCLSE1 : As NSMCLSE, but the symmetry of the last active 
*            orbital space where there is more than one symmetry 
*            is left out
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ISPGRP(NGAS),NSTSGP(NSMST,*)
*. Input and Output (column IGRP updated)
      DIMENSION NSTSSPGP(NSMST,IGRP)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'multd2h.inc'
#include "global.fh"
*. Scratch 
      INTEGER MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER MSM1(MXPNSMST),MSM2(MXPNSMST)
      INTEGER ISM1(MXPNSMST),ISM2(MXPNSMST)
*
      NTEST = 0
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' ==========================='
        WRITE(6,*) ' NSTPTP_GAS_NEW is speaking '
        WRITE(6,*) ' ==========================='
*
        WRITE(6,*) ' Supergroup in action '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
      END IF
*
*. The NSMCLS* parameters
*
*. Max and min for each GASpace
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
      END DO
*. Last space with more than one symmetry
      NGASL = 1
      DO IGAS = 1, NGAS
        IF(MXSM(IGAS).NE.MNSM(IGAS)) NGASL = IGAS
      END DO
*. NSMCLSE
      NSMCLSE = 1
      DO IGAS = 1, NGAS
        NSMCLSE = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE
      END DO
*. NSMCLSE1
      NSMCLSE1 = 1
      DO IGAS = 1, NGASL-1
        NSMCLSE1 = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE1
      END DO
*
      IZERO = 0
* 
      IF(NGAS.EQ.0) THEN
*. special treatment for (presumable zero electrons in zero orbitals
        CALL ISETVC(ISM2,IZERO,NSMST)
        ISM2(1) = 1
        CALL ISETVC(MSM2,IZERO,NSMST)
        MSM2(1) = 1
      END IF
*
      DO IGAS = 1, NGAS
*. In ISM1, the number of strings per symmetry for the first 
*  IGAS-1 spaces are given, obtain in ISM2 the number of strings per sym
*  for the first IGAS spaces
*. Also: in MSM1, MSM2, counts the number of nontrivial combinations per
*  sym
        IF(IGAS.EQ.1) THEN
*. ISM1: The number of strings per symmetry for zero electrons
         CALL ISETVC(ISM1,IZERO,NSMST)
         ISM1(1) = 1
         CALL ISETVC(MSM1,IZERO,NSMST)
         MSM1(1) = 1
        ELSE
*. copy from the ISM2 obtained for preceeding IGAS
         CALL ICOPVE(ISM2,ISM1,NSMST)
         CALL ICOPVE(MSM2,MSM1,NSMST)
        END IF
        CALL ISETVC(ISM2,IZERO,NSMST)
        CALL ISETVC(MSM2,IZERO,NSMST)
        DO ISM_IGASM1 = 1, NSMST
         DO ISM_IGAS = 1, NSMST
           ISM = MULTD2H(ISM_IGASM1,ISM_IGAS)
           ISM2(ISM) = ISM2(ISM) +
     &     ISM1(ISM_IGASM1)*NSTSGP(ISM_IGAS,ISPGRP(IGAS))
           IF(ISM1(ISM_IGASM1)*NSTSGP(ISM_IGAS,ISPGRP(IGAS)).NE.0)
     &     MSM2(ISM) = MSM2(ISM) + MSM1(ISM_IGASM1)
         END DO  
        END DO
      END DO !loop over IGAS
      CALL ICOPVE(ISM2,NSTSSPGP(1,IGRP),NSMST)
*
      MXNSTR = 0 
      NSMCLS = 0
      DO ISTRSM = 1, NSMST
        MXNSTR = MAX(MXNSTR,NSTSSPGP(ISTRSM,IGRP))
        NSMCLS = MAX(NSMCLS,MSM2(ISTRSM))
      END DO
      NSTR_AS = IELSUM(NSTSSPGP(1,IGRP),NSMST)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
     &  ' Number of strings per symmetry for supergroup',IGRP
        CALL IWRTMA10(NSTSSPGP(1,IGRP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Largest number of strings of given sym ',MXNSTR
*
        WRITE(6,'(A,3(2X,I8))') ' NSMCLS,NSMCLSE,NSMCLSE1=',
     &                       NSMCLS,NSMCLSE,NSMCLSE1
        WRITE(6,'(A,I9)') ' Number of strings, all symmetries=', NSTR_AS
      END IF
*
      RETURN
      END
      SUBROUTINE NSTPTP_GAS(NGAS,ISPGRP,NSTSGP,NSMST,
     &                      NSTSSPGP,IGRP,MXNSTR,
     &                      NSMCLS,NSMCLSE,NSMCLSE1)
*
* From number of strings per group and sym to number of strings
* per supergroup and sym for given super group.
*
* Jeppe Olsen , Fall of 94
* 
* Last Revision : Aug 97, NSMCLS,NSMCLSE added
*
* NSMCLS : MAX Number of symmetry classes for given supergroup,
*          i.e. number of combinations of symmetries of groups
*          containing strings
* NSMCLSE : Number of symmetry classes for given supergroup 
*          obtained by restricting allowed symmetries in 
*          a given group by a max and min.
* NSMCLSE1 : As NSMCLSE, but the symmetry of the last active 
*            orbital space where there is more than one symmetry 
*            is left out
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ISPGRP(NGAS),NSTSGP(NSMST,*)
*. Output
      DIMENSION NSTSSPGP(NSMST,IGRP)
*. Scratch 
      INCLUDE 'mxpdim.inc'
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER MSMCLS(MXPNSMST)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NSTPTP_GAS is speaking '
        WRITE(6,*) ' ======================='
*
        WRITE(6,*) ' Supergroup in action '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
      END IF
*
      IZERO = 0
      CALL ISETVC(NSTSSPGP(1,IGRP),IZERO,NSMST)
      CALL ISETVC(MSMCLS,IZERO,NSMST)
      
*. Max and Min for allowed symmetries
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
      END DO
*. Last space with more than one symmetry
      NGASL = 1
      DO IGAS = 1, NGAS
        IF(MXSM(IGAS).NE.MNSM(IGAS)) NGASL = IGAS
      END DO
*. NSMCLSE
      NSMCLSE = 1
      DO IGAS = 1, NGAS
        NSMCLSE = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE
      END DO
*. NSMCLSE1
      NSMCLSE1 = 1
      DO IGAS = 1, NGASL-1
        NSMCLSE1 = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE1
      END DO
*. First symmetry combination
      DO IGAS = 1, NGAS
         ISM(IGAS) = MNSM(IGAS)
      END DO
*. Loop over symmetries in each gas space
      IFIRST = 1
 1000 CONTINUE
      IF(IFIRST.EQ.1) THEN
        CALL ISETVC(ISM,1,NGAS)
        IFIRST = 0
        NONEW = 0
      ELSE
C       NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in this supergroup
        IF(NGAS.GT.0) THEN
          ISMSPGP = ISM(1)
          NST = NSTSGP(ISM(1),ISPGRP(1))
        ELSE
          NST = 1
          ISMSPGP = 1
        END IF
        DO JGRP = 2, NGAS
          CALL SYMCOM(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
          ISMSPGP = ISMSPGPO
          NST = NST * NSTSGP(ISM(JGRP),ISPGRP(JGRP))
        END DO
C?      WRITE(6,*) ' Symmetry of groups '
C?      CALL IWRTMA(ISM,1,NGAS,1,NGAS)
C?      WRITE(6,*) ' Symmetry of super group and number of strings '
C?      WRITE(6,*) ' ISMSPGP , NST = ', ISMSPGP , NST
        NSTSSPGP(ISMSPGP,IGRP) =   NSTSSPGP(ISMSPGP,IGRP) + NST
        IF(NST.NE.0) MSMCLS(ISMSPGP) = MSMCLS(ISMSPGP) + 1
        GOTO 1000
      END IF         
*
      MXNSTR = 0 
      NSMCLS = 0
      DO ISTRSM = 1, NSMST
        MXNSTR = MAX(MXNSTR,NSTSSPGP(ISTRSM,IGRP))
        NSMCLS = MAX(NSMCLS,MSMCLS(ISTRSM))
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
     &  ' Number of strings per symmetry for supergroup',IGRP
        CALL IWRTMA10(NSTSSPGP(1,IGRP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Largest number of strings of given sym ',MXNSTR
*
        WRITE(6,'(A,3(2X,I8))') ' NSMCLS,NSMCLSE,NSMCLSE1=',
     &                       NSMCLS,NSMCLSE,NSMCLSE1
      END IF
*
      RETURN
      END
      SUBROUTINE GASSPC2(I_IADP,IOC)
*
*
* Divide orbital spaces ( I_IADX ) into 
*
*  Inactive spaces : Orbitals that are doubly occupied in all CI spaces
*  Active orbitals : Orbitals that have variable occ in atleast some spaces.
*  Secondary spaces: Orbitals that are unoccupied in all spaces
*
* Division based upon occupation in CI space IOC
*
* Jeppe Olsen, Summer of 98 ( not much of an summer !)
*
* Slight modification of GASSPC
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER IOC(MXPNGAS,2)
*. Output
      INTEGER I_IADP(NGAS)
*
      NEL_REF = NELEC(1) + NELEC(2)
*
      DO IGAS = 1, NGAS
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS) 
       END IF
*
       IF(IOC(IGAS,1) .EQ. NEL_MAX  .AND.
     &    IOC(IGAS,2) .EQ. NEL_MAX       ) THEN 
*. Inactive  space
          I_IADP(IGAS) = 1
       ELSE IF(IGAS.GT.1.AND.IOC(IGAS-1,1) .EQ. NEL_REF ) THEN
*. Secondary space
          I_IADP(IGAS) = 3
       ELSE 
*. Active space
          I_IADP(IGAS) = 2
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Division of orbitals '
        WRITE(6,*) ' ======================= ' 
        WRITE(6,*)
        WRITE(6,*) ' Inactive = 1, Active = 2, Delete = 3 '
        WRITE(6,*)
        CALL IWRTMA(I_IADP,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE GASSPC
*
*
* Divide orbital spaces into 
*
*  Inactive spaces : Orbitals that are doubly occupied in all CI spaces
*  Active orbitals : Orbitals that have variable occ in atleast some spaces.
*  Secondary spaces: Orbitals that are unoccupied in all spaces
*
* I_IAD : Division based upon occupation in Compound CI spaces IGSOCC
* I_IADX: Division based upon occupation in First CI space 
*
* Jeppe Olsen, Summer of 98 ( not much of an summer !)
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
*
      NEL_REF = NELEC(1) + NELEC(2)
*
* For compound space 
*
      NEL_MAX = -2303
      DO IGAS = 1, NGAS
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS) 
       END IF
*
       IF(IGSOCC(IGAS,1) .EQ. NEL_MAX  .AND.
     &    IGSOCC(IGAS,2) .EQ. NEL_MAX       ) THEN 
*. Inactive  space
          I_IAD(IGAS) = 1
       ELSE IF(IGAS.GT.1.AND.IGSOCC(IGAS-1,1) .EQ. NEL_REF ) THEN
*. Delete space
          I_IAD(IGAS) = 3
       ELSE 
*. Active space
          I_IAD(IGAS) = 2
       END IF
*
      END DO
*
* For First CI space 
*
      DO IGAS = 1, NGAS
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS) 
       END IF
*
       IF(IGSOCCX(IGAS,1,1) .EQ. NEL_MAX  .AND.
     &    IGSOCCX(IGAS,2,1) .EQ. NEL_MAX       ) THEN 
*. Inactive  space
          I_IADX(IGAS) = 1
       ELSE IF(IGAS.GT.1.AND.IGSOCCX(IGAS-1,1,1) .EQ. NEL_REF ) THEN
*. Delete space
          I_IADX(IGAS) = 3
       ELSE 
*. Active space
          I_IADX(IGAS) = 2
       END IF
*
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) 
     &  ' Division of orbitals according to compound CI space'
        WRITE(6,*) 
     &  ' ================================================== ' 
        WRITE(6,*)
        WRITE(6,*) ' Inactive = 1, Active = 2, Delete = 3 '
        WRITE(6,*)
        CALL IWRTMA(I_IAD,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Division of orbitals according to first CI space'
        WRITE(6,*) 
     &  ' ================================================== ' 
        WRITE(6,*)
        WRITE(6,*) ' Inactive = 1, Active = 2, Delete = 3 '
        WRITE(6,*)
        CALL IWRTMA(I_IADX,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE STRTYP_GAS(IPRNT)
*
* Find groups of strings in each GA space
*
* Output : /GASSTR/
*
* Jeppe Olsen, Oct 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'spinfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Local scratch
      DIMENSION IOCTYP(MXPSTT),IREOSPGP(MXPSTT),ISCR(MXPSTT)
*
      CALL LUCIAQENTER('STRTY')
*
      NTESTL = 000
      NTEST = MAX(IPRNT,NTESTL)
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Info from routine STRTYP_GAS '
        WRITE(6,*) ' ============================='
      END IF
*. As input NCISPC GAS spaces IGSOCCX are given.
* Obtain space that cantains all these as special cases
*
C?    WRITE(6,*) ' NCISPC ', NCISPC 
      DO IGAS = 1, NGAS
       MINI = IGSOCCX(IGAS,1,1)
       MAXI = IGSOCCX(IGAS,2,1)
C?     WRITE(6,*) ' MINI and MAXI for ISPC = 1 ',MINI,MAXI
       DO ICISPC = 2, NCISPC
        MINI = MIN(MINI,IGSOCCX(IGAS,1,ICISPC))
        MAXI = MAX(MAXI,IGSOCCX(IGAS,2,ICISPC))
C?     WRITE(6,*) ' MINI and MAXI for ISPC =  ',ICISPC,MINI,MAXI
       END DO
       IGSOCC(IGAS,1) = MINI
       IGSOCC(IGAS,2) = MAXI
      END DO
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Compound GAS space : '
        WRITE(6,*) ' ====================='
        WRITE(6,'(A)')
        WRITE(6,'(A)') '         Min. occ    Max. occ '
        WRITE(6,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(6,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCC(IGAS,1),IGSOCC(IGAS,2)
        END DO
      END IF
*. Allowed number of electrons in ensemble gas for compound statw
      IF(NELVAL_IN_ENSGS(1).NE.-1) CALL Z_ENSOCC_CMP
*. Max number of unpaired electrons and occupation classes in 
*. compound space
      CALL MAX_NOPEN_NOCCLS_FOR_CISPAC(IGSOCC(1,1),IGSOCC(1,2),
     &     NOPEN_MAX,NOCCLS_MAX)
*. Number of AB supergroups in expansion of the compound space
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,'(A,I4)') 
     &  ' Max number of unpaired electrons ', NOPEN_MAX
        WRITE(6,'(A,I8)')
     &  ' Max number of occupation classes ', NOCCLS_MAX
      END IF
      MAXOP = NOPEN_MAX
      NMXOCCLS = NOCCLS_MAX
*. MINOP - from MS2 and multiplicity if relevant
      MINOP = ABS(MS2)
      IF(NOCSF.EQ.0) MINOP =  MULTS - 1
      if (ga_nodeid().eq.0.and.NTEST.ge.5) WRITE(6,*) ' MINOP == ', MINOP
*
*. Find min and max number of elecs in each subspace
*
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          MNGSOC(IGAS) = IGSOCC(IGAS,1)
          MXGSOC(IGAS) = IGSOCC(IGAS,2)
        ELSE
          MXGSOC(IGAS) = IGSOCC(IGAS,2)-IGSOCC(IGAS-1,1)
          MNGSOC(IGAS) = MAX(0,IGSOCC(IGAS,1)-IGSOCC(IGAS-1,2))
        END IF
      END DO
*
*. Number of occupation subclasses
*
      NOCSBCLST = 0
      DO IGAS = 1, NGAS
        NOCSBCLST = NOCSBCLST + MXGSOC(IGAS)-MNGSOC(IGAS)+1
      END DO
      IF (NTEST.GT.5.and.ga_nodeid().eq.0)
     &WRITE(6,*) ' Number of occupation subclasses ', NOCSBCLST
*
*. Particle and hole spaces  : 
*  Hole spaces are always more than half occupied
*
      IPHGASL = 0
      NPHGAS = 0
      DO IGAS = 1, NGAS
        IF(IUSE_PH.EQ.1) THEN 
*. P/H separation for compound space 
*. Below is correct and my standard def
          IF(MNGSOC(IGAS).GE.NOBPT(IGAS)) THEN
            IPHGAS(IGAS) = 2
            IPHGASL = IGAS
            NPHGAS = NPHGAS + 1
          ELSE
             IPHGAS(IGAS) = 1
          END IF
*. P/H separation for initial  space in IPHGAS1 and strict P/H
*. separation in IPHGASS - in the strict separation, only
*. doubly occupied spaces are hole spaces
          IF(IGAS.EQ.1) THEN
            MIN_OC1 = IGSOCCX(IGAS,1,1)
          ELSE
            MIN_OC1 = IGSOCCX(IGAS,1,1)-IGSOCCX(IGAS-1,2,1)
          END IF
*
          IF(MIN_OC1.GT.NOBPT(IGAS)) THEN
            IPHGAS1(IGAS) = 2
          ELSE
            IPHGAS1(IGAS) = 1
          END IF
*
          IF(MIN_OC1.EQ.2*NOBPT(IGAS)) THEN
            IPHGASS(IGAS) = 2
          ELSE
            IPHGASS(IGAS) = 1
          END IF
*
        ELSE IF(IUSE_PH.EQ.0) THEN
*. The strict division
          IF(IGAS.EQ.1) THEN
            MIN_OC1 = IGSOCCX(IGAS,1,1)
          ELSE
            MIN_OC1 = IGSOCCX(IGAS,1,1)-IGSOCCX(IGAS-1,2,1)
          END IF
*
          IF(MIN_OC1.EQ.2*NOBPT(IGAS)) THEN
            IPHGASS(IGAS) = 2
          ELSE
            IPHGASS(IGAS) = 1
          END IF
*
          IPHGAS(IGAS) = 1
          IPHGAS1(IGAS) = 1
        END IF
      END DO
*. Determine also HPV spaces based on occupation of reference space (1)
      CALL CC_AC_SPACES(1,IREFTYP)
*. Large number of particle and hole orbitals of given type
      MXTSOB_P = 0
      MXTSOB_H = 0
      DO IGAS = 1, NGAS
        IF(IPHGAS1(IGAS).EQ.1) THEN
          MXTSOB_P = MAX(MXTSOB_P,NOBPT(IGAS))
        ELSE
          MXTSOB_H = MAX(MXTSOB_H,NOBPT(IGAS))
        END IF
      END DO
C?    WRITE(6,*) ' MXTSOB_H, MXTSOB_P = ', MXTSOB_H, MXTSOB_P
*
C?    IF(IUSE_PH.EQ.1) THEN
C?    IPHGAS(1) = 2
C?    DO I = 1, 100
C?      WRITE(6,*) ' First space enforced to hole space'
C?    END DO
C?    END IF
*
*. In the following I assume that the hole spaces are the first NPGAS spaces
* ( only used when calculating min number of electrons, so it can be modified
*   easily )
*. Min number of electrons in hole spaces
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) 
     &   WRITE(6,*) ' IPHGASL, NPHGAS ',  IPHGASL, NPHGAS
      IF(IPHGASL.NE.NPHGAS) THEN
        WRITE(6,*) ' The hole spaces are not the first orbital spaces'
        WRITE(6,*) ' The hole spaces are not the first orbital spaces'
        WRITE(6,*) ' The hole spaces are not the first orbital spaces'
        WRITE(6,*) ' The hole spaces are not the first orbital spaces'
C       STOP       ' The hole spaces are not the first orbital spaces'
      END IF
      MNHL = IGSOCC(IPHGASL,1)
C     MNHL  = 0
C     DO IGAS = 1, NGAS 
C       IF(IPHGAS(IGAS).EQ.2) THEN
C         MNHL = MNHL + MNGSOC(IGAS)
C       END IF
C     END DO
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) WRITE(6,*) ' MNHL' , MNHL 
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Min and Max occupation in each GAS space: '
        WRITE(6,'(A)') ' ========================================= '
        WRITE(6,*)
        DO IGAS = 1,  NGAS
          WRITE(6,'(A,I2,4X,2I3)') 
     &    '  GAS',IGAS,MNGSOC(IGAS),MXGSOC(IGAS)
        END DO
*
        WRITE(6,*)' Particle(1) or hole(2) spaces (for compound space)'
        CALL IWRTMA(IPHGAS,1,NGAS,1,NGAS)
        WRITE(6,*)' Particle(1) or hole(2) spaces (for initial space)'
        CALL IWRTMA(IPHGAS1,1,NGAS,1,NGAS)
        WRITE(6,*)' Strict Particle(1) or hole(2) spaces '
        CALL IWRTMA(IPHGASS,1,NGAS,1,NGAS)
       END IF
*
* Split into alpha and beta parts
*
*. Number of alpha. and beta electrons
*
      NAEL = (MS2 + NACTEL ) / 2
      NBEL = (NACTEL - MS2 ) / 2
*
      IF(NAEL + NBEL .NE. NACTEL ) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I4)')   MS2,NACTEL,NAEL,NBEL
        WRITE(6,*)
     &  ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT '
        STOP ' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY '
      END IF
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I6)')   MS2,NACTEL,NAEL,NBEL
      END IF
*
      IF(NAEL + NBEL .NE. NACTEL ) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I4)')   MS2,NACTEL,NAEL,NBEL
        WRITE(6,*)
     &  ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT '
          STOP ' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY '
      END IF
*
* --------------------------------------------------
* Construct groups : occupations for a each GASpace
* --------------------------------------------------
*
*. Increase for MRCC to include N+3, N+4, N-3, N-4
cVOG  I_MAYBE_DO_MRCC = 1
      I_MAYBE_DO_MRCC = 0
*. Number of electrons to be subtracted or added
      MAXSUB = 2
      MAXADD = 2
*. electrons are only added for systems that atleast have halffilled 
*. shells or are declared valence by IHPVGAS
      IGRP = 0
      MXAL = NAEL
      MNAL = NAEL
      MXBL = NBEL
      MNBL = NBEL
      NORBL = NTOOB
      DO IGAS = 1, NGAS
*. occupation constraints 1 
       MXA1 = MIN(MXGSOC(IGAS),NOBPT(IGAS),MXAL)
       MXB1 = MIN(MXGSOC(IGAS),NOBPT(IGAS),MXBL)
       MNA1 = MAX(0,MNGSOC(IGAS)-MXA1)    
       MNB1 = MAX(0,MNGSOC(IGAS)-MXB1)    
*. Additional checks can be made here
       MXA = MXA1
       MXB = MXB1
       MNA = MNA1
       MNB = MNB1
*
       MXAL = MXAL - MNA
       MNAL = MAX(0,MNAL-MXA)
       MXBL = MXBL - MNB
       MNBL = MAX(0,MNBL-MXB)
*
       IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Occupation numbers for IGAS = ', IGAS
        WRITE(6,*) ' MXAL MNAL MXBL MNBL ',MXAL,MNAL,MXBL,MNBL
        WRITE(6,*) ' MXA MNA MXB MNB ',MXA,MNA,MXB,MNB
       END IF
*
       MNAB = MIN(MNA,MNB)
       MXAB = MAX(MXA,MXB)
*. Additional holes only allowed in particle spaces or valence space
       IF(IPHGAS(IGAS).EQ.1.OR.IHPVGAS(IGAS).EQ.3) THEN 
         MNAB = MAX(0,MNAB-MAXSUB)
       ELSE IF(IPHGAS(IGAS).EQ.2) THEN
         MNAB = MNAB
       END IF
C. For coupled cluster- could be refined ...
*. I do not think this is needed ..( Aug. 03)
*. Is needed as excitation strings must also be constructed,
* but something must be done soon JEPPE !!!!!
       IF(I_MAYBE_DO_MRCC.EQ.1) THEN
*. Is actually also needed for standard CC
       MNAB = 0
C?     WRITE(6,*) ' MNAB = 0 statement eliminated !! '
C?     WRITE(6,*) ' MNAB = 0 statement eliminated !! '
C?     WRITE(6,*) ' MNAB = 0 statement eliminated !! '
C?     WRITE(6,*) ' (Problems for CC) '
*. Additional electrons allowed in hole spaces
       IF(IPHGAS(IGAS).EQ.2.OR.IHPVGAS(IGAS).EQ.3)
     & MXAB = MIN(MXAB + MAXADD,NOBPT(IGAS))
       END IF
*
       IF(NTEST.GE.100.and.ga_nodeid().eq.0) 
     &    WRITE(6,*) ' MNAB,MXAB',MNAB,MXAB
       NGPSTR(IGAS) = MXAB-MNAB+1
       IBGPSTR(IGAS) = IGRP + 1
       MNELFGP(IGAS) = MNAB
       MXELFGP(IGAS) = MXAB
*
       IADD = 0
       DO JGRP = IGRP+1,IGRP+NGPSTR(IGAS)
         IF(JGRP.GT.MXPSTT) THEN
           WRITE(6,*) ' Too many string groups '
           WRITE(6,*) ' Current limit ', MXPSTT 
           WRITE(6,*) ' STOP : GASSTR, Too many string groups'
                        STOP ' GASSTR, Too many string groups'
          END IF
*
         IADD = IADD + 1
         IEL = MNAB-1+IADD
         NELFGP(JGRP) = IEL           
         IGSFGP(JGRP) = IGAS
         NSTFGP(JGRP) = IBION(NOBPT(IGAS),IEL)
       END DO
       IGRP = IGRP + NGPSTR(IGAS)
      END DO
      NGRP = IGRP
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Information about Groups of strings '
        WRITE(6,'(A)') ' =================================== '
        WRITE(6,*)
        WRITE(6,*) '     GAS  MNEL  MXEL IBGRP   NGRP'
        WRITE(6,*) '    ============================='
        DO IGAS = 1, NGAS
          WRITE(6,'(2X,5(2X,I4))') IGAS,MNELFGP(IGAS),
     &          MXELFGP(IGAS),IBGPSTR(IGAS),NGPSTR(IGAS)
        END DO
        WRITE(6,'(A,I3)')
     &  ' Total number of groups generated ', NGRP
*
        WRITE(6,'(A)') ' Information about each string group '
        WRITE(6,'(A)') ' ===================================='
        WRITE(6,*)
        IITYPE = 0
        WRITE(6,'(A)') '   GROUP  GAS   NEL      NSTR '
        WRITE(6,'(A)') ' ============================='
        DO IGRP = 1, NGRP
          IITYPE = IITYPE + 1
          WRITE(6,'(3(2X,I4),2X,I8)')
     &    IITYPE,IGSFGP(IGRP),NELFGP(IGRP),NSTFGP(IGRP)
        END DO
      END IF
*
*. Creation-annihilation connections between groups
*
      DO IGRP = 1, NGRP
        ISTAC(IGRP,1) = 0
        ISTAC(IGRP,2) = 0
        DO JGRP = 1, NGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)-1) ISTAC(IGRP,2) = JGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)+1) ISTAC(IGRP,1) = JGRP
        END DO
      END DO
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================================='
        WRITE(6,*) ' Annihilation / Creation connections of groups'
        WRITE(6,*) ' ============================================='
        WRITE(6,*)
        CALL IWRTMA(ISTAC,NGRP,2,MXPSTT,2)
      END IF
*.
*
* Construct number of type ( combinations of groups ) with nael and nbel strings
*
*
* Type 1 : NAEL electrons
*      2 : NBEL ELECTRONS
*      3 : NAEL -1 ELECTRONS
*      4 : NBEL -1 ELECTRONS
*      5 : NAEL -2 ELECTRONS
*      6 : NBEL -2 ELECTRONS
*      
      NSTTYP = 6
      NSTTP = 6
*
      IF(IUSE_PH.EQ.1) THEN
*. allow N+1,N+2 resolution string
C?      NSTTYP = 10
C?      NSTTP  = 10
      END IF
*. alpha 
      NELEC(1) = NAEL
      NELFTP(1) = NAEL
*. beta
      NELEC(2) = NBEL
      NELFTP(2) = NBEL
      NSTTYP = 2
*. alpha -1 
      NSTTYP = NSTTYP + 1
      IF(NAEL.GE.1) THEN 
        NELEC(NSTTYP) = NAEL-1
        NELFTP(NSTTYP) = NAEL-1
      ELSE 
        NELEC(NSTTYP) = 0 
        NELFTP(NSTTYP) = 0
      END IF
*. beta  -1 
      NSTTYP = NSTTYP + 1
      IF(NBEL.GE.1) THEN
        NELEC(NSTTYP) = NBEL-1
        NELFTP(NSTTYP) = NBEL-1
      ELSE 
        NELEC(NSTTYP) = 0 
        NELFTP(NSTTYP) = 0
      END IF
*. alpha -2 
      NSTTYP = NSTTYP + 1
      IF(NAEL.GE.2) THEN
        NELEC(NSTTYP) = NAEL-2
        NELFTP(NSTTYP) = NAEL-2
      ELSE 
        NELEC(NSTTYP) = 0 
        NELFTP(NSTTYP) = 0
      END IF
*. beta  -2  
      NSTTYP = NSTTYP + 1
      IF(NBEL.GE.2) THEN
        NELEC(NSTTYP) = NBEL-2
        NELFTP(NSTTYP) = NBEL-2
      ELSE 
        NELEC(NSTTYP) = 0 
        NELFTP(NSTTYP) = 0
      END IF
*. alpha-3
      IF( I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = MAX(0,NAEL - 3)
       NELFTP(NSTTYP) = MAX(0,NAEL - 3)
      END IF
*. beta-3
      IF( I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = MAX(0,NBEL - 3)
       NELFTP(NSTTYP) = MAX(0,NBEL - 3)
      END IF
*. alpha-4
      IF( I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = MAX(0,NAEL - 4)
       NELFTP(NSTTYP) = MAX(0,NAEL - 4)
      END IF
*. beta-4
      IF( I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = MAX(0,NBEL - 4)
       NELFTP(NSTTYP) = MAX(0,NBEL - 4)
      END IF
*
      IF((IUSE_PH.EQ.1).OR.(I_MAYBE_DO_MRCC.EQ.1)) THEN
*. Alpha + 1
        NSTTYP = NSTTYP + 1
        NELEC(NSTTYP) = NAEL+1
        NELFTP(NSTTYP) = NAEL+1
*. beta  + 1
        NSTTYP = NSTTYP + 1
        NELEC(NSTTYP) = NBEL+1
        NELFTP(NSTTYP) = NBEL+1
*. Alpha + 2
        NSTTYP = NSTTYP + 1
        NELEC(NSTTYP) = NAEL+2
        NELFTP(NSTTYP) = NAEL+2
*. beta  + 2
        NSTTYP = NSTTYP + 1
        NELEC(NSTTYP) = NBEL+2
        NELFTP(NSTTYP) = NBEL+2
      END IF
* alpha + 3
      IF(I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = NAEL + 3
       NELFTP(NSTTYP) = NELEC(NSTTYP)
      END IF
* beta + 3
      IF(I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = NBEL + 3
       NELFTP(NSTTYP) = NELEC(NSTTYP)
      END IF
* alpha + 4
      IF(I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = NAEL + 4
       NELFTP(NSTTYP) = NELEC(NSTTYP)
      END IF
* beta + 4
      IF(I_MAYBE_DO_MRCC.EQ.1) THEN
       NSTTYP = NSTTYP + 1
       NELEC(NSTTYP) = NBEL + 4
       NELFTP(NSTTYP) = NELEC(NSTTYP)
      END IF
*
      NSTTP = NSTTYP
*. Can easily be extended to relativistic case !!
      DO ITP = 1, NSTTYP
        NOCTYP(ITP) = 0
        NSPGPFTP(ITP) = 0
      END DO 
*
*  ------------------------------------------
*. Generate now the supergroups for each type 
*  ------------------------------------------
*
*
* Loop over types, i.e.  given number of electrons
*
      IOFF = 1
      NABEL = NAEL + NBEL 
      NSPGP_TOT = 0
      DO 2000 ITYP = 1, NSTTYP
*. Number of electrons in reference space ( alpha or beta )
        IF(MOD(ITYP,2).EQ.1) THEN
*. alpha type
          NELEC_REF = NELEC(1)     
          IAB = 1
        ELSE 
*. Beta type 
          NELEC_REF = NELEC(2)
          IAB = 2
        END IF
*. If we are studying beta type, and number of alpha and beta 
* electrons are identical, just refer to alpha
        IF(NAEL.EQ.NBEL.AND.MOD(ITYP,2).EQ.0) THEN
          IBSPGPFTP(ITYP) =  IBSPGPFTP(ITYP-1)
          NOCTYP(ITYP) =   NOCTYP(ITYP-1)
          NSPGPFTP(ITYP) =  NSPGPFTP(ITYP-1)
        ELSE
*
*
*. Number of electrons removed compared to reference 
        IDEL = NELEC(ITYP) - NELEC_REF
C?      WRITE(6,*) '  GASSPC : ITYP IDEL ', ITYP,IDEL 
*. Initial type of strings, relative to offset for given group
        DO IGAS = 1, NGAS
          IOCTYP(IGAS) = 1
        END DO
        NSPGP = 0
        IBSPGPFTP(ITYP) = IOFF 
        IF(NELEC(ITYP).LT.0) THEN
          NOCTYP(ITYP) = 0     
          NSPGPFTP(ITYP) =  0                
          GOTO 2000
        END IF
*. Number of electrons in present type
*. Loop over  SUPER GROUPS with current nomenclature! 
*. Temp max for loop 
        MXLOOP = 10000
        IONE = 1
        XNLOOP = 0.0D0
 1000   CONTINUE
*. Number of electrons in present supergroup
          NEL = 0
          DO IGAS = 1, NGAS
            NEL = NEL + NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
          END DO
*
          IF(NEL.GT.NELEC(ITYP)) THEN
*. If the number of electrons is to large find next number that 
* can be correct.
* The following uses that within a given GAS space 
* the number of elecs increases as the type number increases
*
*. First integer  that can be reduced
            IRED = 0
            DO IGAS = 1, NGAS
              IF(IOCTYP(IGAS).NE.1) THEN
                IRED = IGAS
                GOTO 888
              END IF
            END DO
  888       CONTINUE
            IF(IRED.EQ.NGAS) THEN
              NONEW = 1
            ELSE IF(IRED.LT.NGAS) THEN
              IOCTYP(IRED) = 1
*. Increase remanining part
              CALL NXTNUM2(
     &        IOCTYP(IRED+1),NGAS-IRED,IONE,NGPSTR(IRED+1),NONEW)
            END IF
            GOTO 2803 
          END IF

          IF(NEL.EQ.NELEC(ITYP)) THEN
*. test 1 has been passed, check additional occupation constraints
*
           I_AM_OKAY = 1 
*. Number of extra holes in hole spaces
CE         IF(IUSE_PH.EQ.1) THEN
CE           IDELP = 0
CE           IDELM = 0
CE           DO IGAS = 1, NGAS
CE             IF(IPHGAS(IGAS).EQ.2) THEN
CE               NELH =  NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
CE               IF(NELH.LT.MNGSOC(IGAS)) 
CE    &          IDELM = IDELM +MNGSOC(IGAS)-NELH
CE               IF(NELH.GT.MXGSOC(IGAS)) 
CE    &          IDELP = IDELP + NELH-MXGSOC(IGAS)
CE             END IF
CE           END DO
CE           IF(IDELM.GT.0.OR.IDELP.GT.MAX(0,IDEL)) THEN
CE             I_AM_OKAY = 0
CE             WRITE(6,*) ' P/H rejected supergroup '
CE             CALL IWRTMA(IOCTYP,1,NGAS,1,NGAS)
CE             WRITE(6,*) ' IDELM, IDELP ', IDELM, IDELP         
CE           END IF
CE         END IF
*
*. Check from above 
*
           DO IGAS = NGAS, 1, -1 
*. Number of electrons when all electrons of AS IGAS have been added
             IF(IGAS.EQ.NGAS ) THEN
               IEL = MAX(NABEL,NABEL+IDEL)
             ELSE
               IEL = IEL-NELFGP(IOCTYP(IGAS+1)+IBGPSTR(IGAS+1)-1)
               IF(IEL.LT.MAX(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL))
     &         I_AM_OKAY = 0
             END IF
           END DO
*
* Check from below
*
           IEL = 0
           IOELMX = 0
           DO IGAS = 1, NGAS
             IEL = IEL + NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
             IOELMX = IOELMX+NOBPT(IGAS)
             IF(IEL+IOELMX.LT.MIN(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL))
     &       I_AM_OKAY = 0
           END DO
*
* electrons should never be removed from hole spaces, so
*
* Added aug 03 : 
           IF(IDEL.LT.0) THEN
             DO IGAS = 1, NGAS 
               IF(IPHGAS(IGAS).EQ.2) THEN
*. Min number of electrons in reference supergroups 
                 JREFMIN = NELFGP(ISPGPFTP(IGAS,IBSPGPFTP(IAB)))
                 DO JJSPGP = 2, NSPGPFTP(IAB)
                   JREFMIN = MIN(JREFMIN,
     &             NELFGP(ISPGPFTP(IGAS,IBSPGPFTP(IAB)-1+JJSPGP)) ) 
                 END DO
                 IF(NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1).LT.
     &              JREFMIN) I_AM_OKAY = 0
               END IF
             END DO
           END IF
*
           IF(I_AM_OKAY.EQ.1) THEN
*. passed !!!
             NSPGP = NSPGP + 1
*. Copy supergroup to ISPGPFTP with absolute group numbers
             DO IGAS = 1, NGAS
               ISPGPFTP(IGAS,IOFF-1+NSPGP) 
     &       = IOCTYP(IGAS)+IBGPSTR(IGAS)-1
             END DO
           END IF
*
          END IF
*. Next type of strings
          IONE = 1
          CALL NXTNUM2(IOCTYP,NGAS,IONE,NGPSTR,NONEW)
          XNLOOP = XNLOOP + 1.0D0
 2803   CONTINUE
        IF(NONEW.EQ.0) GOTO 1000
C?      WRITE(6,*) ' ITYP, XNLOOP  = ', ITYP, XNLOOP
*. End of loop over possible supergroups, save information about current type
        IOFF = IOFF + NSPGP
        NOCTYP(ITYP) = NSPGP
        NSPGPFTP(ITYP) =  NSPGP            
        NSPGP_TOT =  NSPGP_TOT +  NSPGP
      END IF
 2000 CONTINUE
      NTSPGP = NSPGP_TOT
C        IMNMX(IVEC,NDIM,MINMAX)
*. Largest number of supergroups of any type
      NSPGPFTP_MAX = IMNMX(NSPGPFTP,NSTTYP,2)
      if (NTEST.GE.10.AND.ga_nodeid().eq.0)
     &   WRITE(6,*) ' Largest number of supergroups, any type ',
     &             NSPGPFTP_MAX
*
      IF(NSPGP_TOT .GT. MXPSTT ) THEN
        WRITE(6,*) ' Too many super groups = ', NSPGP_TOT               
        WRITE(6,*) ' Increase MXPSTT to this value'
        WRITE(6,*) ' See you later '              
        WRITE(6,*)
        WRITE(6,*) ' STOP Increase MXPSTT '
        STOP' Increase MXPSTT'
      END IF
*
*. Reorder supergroups according to dimension
*
* Reordering eliminated by Jeppe, Dec 2011 to allow for search in 
* ordered list
*. START OF CHANGE, Jeppe, December 2011
      I_DO_REO = 0
      IF(I_DO_REO.EQ.0) THEN
        IF (NTEST.GT.2.and.ga_nodeid().eq.0)
     &  WRITE(6,*) ' No reordering of supergroups'
      ELSE
       IF (NTEST.GT.2.and.ga_nodeid().eq.0)
     & WRITE(6,*) ' Reordering of supergroups'
       DO ITYP= 1, NSTTP
         IBTYP = IBSPGPFTP(ITYP)            
         NSPGP = NSPGPFTP(ITYP)
*.Dimension of supergroups
         DO ISPGP = 1, NSPGP
           IDIM = 1
           DO JGAS = 1, NGAS
             IDIM = IDIM * NSTFGP(ISPGPFTP(JGAS,ISPGP+IBTYP-1))
           END DO
           IOCTYP(ISPGP) = IDIM
         END DO
C            ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
         CALL ORDINT(IOCTYP,ISCR,NSPGP,IREOSPGP,NTEST)
*. And reorder the definition of supergroups
         DO ISPGP = 1, NSPGP
          CALL ICOPVE(ISPGPFTP(1,ISPGP+IBTYP-1),NELFSPGP(1,ISPGP),NGAS)
         END DO
         DO ISPGP_N = 1, NSPGP
           ISPGP_O = IREOSPGP(ISPGP_N)
           CALL ICOPVE(NELFSPGP(1,ISPGP_O),ISPGPFTP(1,ISPGP_N+IBTYP-1),
     &                 NGAS) 
         END DO
*
       END DO ! Loop over types
      END IF! reordering
*  END OF CHANGE, JEPPE, DEC. 2011
*
*
      IF(NTEST.GE.2.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' Total number of super groups ', NTSPGP
      END IF
      
      if (ga_nodeid().eq.0) WRITE(6,*)
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Information about types of strings'
        WRITE(6,*) ' =================================='
        WRITE(6,*)
        DO ITYP = 1, NSTTYP
          WRITE(6,*)
          WRITE(6,*) '      Type : ', ITYP
          WRITE(6,*) '      ==============='
          WRITE(6,*) '      Number of electrons  ',NELFTP(ITYP)
          WRITE(6,*) '      Number of super groups ', NSPGPFTP(ITYP)
          IF(NTEST.GE.10) THEN
          WRITE(6,*) '      Supergroups '
          DO ISPGP = 1, NSPGPFTP(ITYP)
            IOFF = IBSPGPFTP(ITYP)
            CALL IWRTMA(ISPGPFTP(1,IOFF-1+ISPGP),1,NGAS,1,NGAS)
          END DO
          END IF
        END DO
      END IF
*
*. Division of Orbitals and spinorbitals into holes, particle and valence
* (IHPVGAS and IHPVGAS_AB)
*. (Done down here as NELEC arrays must be defined )
      CALL CC_AC_SPACES(1,IREFTYP)
C?    WRITE(6,*) ' IREFTYP after call to CC_AC = ', IREFTYP
*
* Number of AB supergroups in compound CI space
*
      IATP = 1
      IBTP = 2
      NASPGP = NSPGPFTP(IATP)
      NBSPGP = NSPGPFTP(IBTP)
      IB_A = IBSPGPFTP(IATP)
      IB_B = IBSPGPFTP(IBTP)
C?    WRITE(6,*) ' IB_A, IB_B = ', IB_A, IB_B
      CALL ABSPGP_TO_CISPACE(ISPGPFTP(1,IB_A),NASPGP,
     &     ISPGPFTP(1,IB_B),NBSPGP,NGAS,IGSOCC,NELFGP,
     &     N_ABSPGP_MAX,IDUM,1)
*
C     ABSPGP_TO_CISPACE(
C    &           IASPGP,NASPGP,IBSPGP,NBSPGP,
C    &           NGAS,IGSOCC_MNMX,NELFGP,
C    &           N_ABSPGP,I_ABSPGP,
C    &           IONLY_NABSPGP)
      CALL LUCIAQEXIT('STRTY')
*
      RETURN
      END 
CNW   SUBROUTINE STRINF_GAS(STIN,IPRNT)
      SUBROUTINE STRINF_GAS(IPRNT)
*
* Obtain string information for GAS expansion
*
* Last modification; Jeppe Olsen; April 29, 2013; MXNSTR_AS added
*
* =====
*.Input
* =====
*
* /LUCINP/,/ORBINP/,/CSM/, /CGAS/, /GASSTR/
*
* =====
*.Output
* =====
*
* /STRINP/,/STINF/,/STRBAS/ and string information in STIN
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. Input
*     (and /LUCINP/ not occuring here )
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'irat.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*
CNW   DIMENSION STIN(*)
*. A bit of scratch
C     DIMENSION IOCTYP(MXPNGAS)
*
      CALL LUCIAQENTER('STRIN')
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
         WRITE(6,*) ' Output from STRINF_GAS '
         WRITE(6,*) ' ======================='
      END IF
*
**.2 : Number of classes per string type and mappings between
**.    string types (/STINF/)
*
      CALL ZSTINF_GAS(IPRNT)
*
**.3 : Static memory for string information
*
       CALL MEMSTR_GAS
*
** 4 : Info about group of strings
*
CBERT FIRST FREE ADDRESS NNOT USED, ALLOCATE BLOCK
*.First free address
CNW   CALL MEMMAN(KFREEL,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
c     IDUMMY = 2*NOCOB+(NACTEL+1)*(NOCOB+1)+NSMST
      IEL = 0
      DO I = 1, NGRP
         IEL = max(IEL,NELFGP(I))
      ENDDO
      IDUMMY = max((2+IEL+1)*(NOCOB+1),NSMST+IEL)
      CALL MEMMAN(KFREEL,IDUMMY,'ADDL  ',1,'KFREEL')
      DO IGRP = 1, NGRP  
*. A gas group can be considered as a RAS group with 0 electrons in
*  RAS1, RAS3 !
        IGAS = IGSFGP(IGRP)
        IF(IGAS.EQ.1) THEN
          NORB1 = NINOB
        ELSE
          NORB1 = NINOB + IELSUM(NOBPT,IGAS-1)
        END IF
        NORB2 = NOBPT(IGAS)
        NORB3 = NOCOB-NORB1-NORB2
        MNRS1X = 0
        MXRS1X = 0
        MNRS3X = 0
        MXRS3X = 0
        IEL = NELFGP(IGRP)
        IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) ' IGRP, NORB1, NORB2,NORB3, IEL = ',
     &                 IGRP, NORB1, NORB2,NORB3, IEL
        END IF
        IOCTYPX = 1
*. Reverse lexical adresing schemes for each group of string
        CALL WEIGHT(int_mb(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &              MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              int_mb(KFREEL),IPRNT )
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS(IEL,NORB1,NORB2,NORB3,
     &              MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              int_mb(KFREEL),NACOB,
     &              int_mb(KNSTSGP(1)),
     &              int_mb(KISTSGP(1)),
     &              IOCTYPX,NSMST,IGRP,IPRNT)
*. Construct the strings ordered by symmetry
        CALL GENSTR_GAS(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              int_mb(KISTSGP(1)),IGRP,
     &              IOCTYPX,NSMST,int_mb(KZ(IGRP)),int_mb(KFREEL),
     &              int_mb(KSTREO(IGRP)),int_mb(KOCSTR(IGRP)),
     &              int_mb(KFREEL+IOCTYPX*NSMST),IGRP,IPRNT)
*
       CALL ICOPVE2(int_mb(KNSTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              NSTFSMGP(1,IGRP))
       CALL ICOPVE2(int_mb(KISTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              ISTFSMGP(1,IGRP))
      END DO
*
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
      WRITE(6,*) ' Number of strings per group and symmetry '
      CALL IWRTMA10(int_mb(KNSTSGP(1)),NSMST,NGRP,NSMST,NGRP)
      WRITE(6,*) ' Number of strings per group and symmetry(2) '
      CALL IWRTMA10(NSTFSMGP,NSMST,NGRP,MXPNSMST,NGRP)
      END IF
*
*. Min and max of sym for each group
*
      DO IGP = 1, NGRP
       MX = 1
       DO ISM = 1, NSMST
         IF(NSTFSMGP(ISM,IGP).GT.0) MX = ISM
       END DO
*
       MN = NSMST
       DO ISM = NSMST,1,-1
         IF(NSTFSMGP(ISM,IGP).GT.0) MN = ISM
       END DO
*
       MINMAX_SM_GP(1,IGP) = MN
       MINMAX_SM_GP(2,IGP) = MX
*
      END DO
      IF(NTEST.GT.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' MINMAX array for sym of groups '
        WRITE(6,*) ' =============================='
        CALL IWRTMA(MINMAX_SM_GP,1,NSMST,1,NSMST)
      END IF
*
*       
* 4.5 : Creation/Annihilation mappings between different
*       types of strings 
*
      DO IGRP = 1, NGRP
*
        IGAS = IGSFGP(IGRP)
        NGSOBP = NOBPT(IGAS)
*. First orbital in GAS spacce
        IGSOB = NINOB + IELSUM(NOBPT,IGAS-1)+1
        IEL = NELFGP(IGRP)
        NSTINI = NSTFGP(IGRP)
*
*. Type of mapping : Only creation                  (LAC = 1)
*                    Only annihilation              (LAC = 2)
*                    Both annihilation and creation (LAC = 3)
* If only annihilation is present the string mapping arrays
* will only be over electronns
        IF(     ISTAC(IGRP,1).NE.0.AND.ISTAC(IGRP,2).NE.0) THEN
          LAC = 3
          IEC = 1
          LROW = NGSOBP
        ELSE IF(ISTAC(IGRP,1).NE.0.AND.ISTAC(IGRP,2).EQ.0) THEN
          LAC = 1
          IEC = 2
          LROW = IEL
        ELSE IF(ISTAC(IGRP,1).EQ.0.AND.ISTAC(IGRP,2).NE.0) THEN
          LAC = 2
          IEC = 0 
          LROW = NGSOBP
        ELSE IF(ISTAC(IGRP,1).EQ.0.AND.ISTAC(IGRP,2).EQ.0) THEN
          LAC = 0
          IEC = 0 
          LROW = 0
        END IF
*. Zero
        IF(LAC.NE.0) THEN
          IZERO = 0
          CALL ISETVC(int_mb(KSTSTM(IGRP,1)),IZERO,LROW*NSTINI)
          CALL ISETVC(int_mb(KSTSTM(IGRP,2)),IZERO,LROW*NSTINI)
        END IF
*
        IF(ISTAC(IGRP,2).NE.0) THEN
          JGRP = ISTAC(IGRP,2)
          IF (NTEST.GT.5) WRITE(6,*) ' Creation group = ', JGRP
          CALL CRESTR_GAS(int_mb(KOCSTR(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,int_mb(KZ(JGRP)),
     &         int_mb(KSTREO(JGRP)),0,IDUM,IDUM,
     &         int_mb(KSTSTM(IGRP,1)),int_mb(KSTSTM(IGRP,2)),NOCOB,
     &         IPRNT)
        END IF
        IF(ISTAC(IGRP,1).NE.0) THEN
          JGRP = ISTAC(IGRP,1)
          CALL ANNSTR_GAS(int_mb(KOCSTR(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,int_mb(KZ(JGRP)),
     &         int_mb(KSTREO(JGRP)),0,IDUM,IDUM,
     &         int_mb(KSTSTM(IGRP,1)),int_mb(KSTSTM(IGRP,2)),NOCOB,IEC,
     /         LROW,IPRNT)
        END IF
      END DO   
*
      CALL LUCIAQENTER('STIN2')
*
*. Now to supergroups , i.e. strings of with given number of elecs in 
*  each GAspace
*
      CALL ISETVC(NSTFSMSPGP,0,MXPNSMST*NTSPGP)
      MXNSTR = -1 
      MXNSTR_AS = -1
      DO ITP = 1, NSTTYP
*. Loop over supergroups of given type . i.e. strings
*  with given occupation in each GAS space
        DO IGRP = 1, NSPGPFTP(ITP)
          IGRPABS = IGRP-1 + IBSPGPFTP(ITP)
          CALL NSTPTP_GAS_NEW(NGAS,ISPGPFTP(1,IGRPABS),
     &                    int_mb(KNSTSGP(1)),NSMST,
     &                    int_mb(KNSTSO(ITP)),IGRP,MXNSTRFSG,
     &                    NSMCLS,NSMCLSE,NSMCLSE1,NSTR_AS)
*
          MXSMCLS   = MAX(MXSMCLS,NSMCLS)
          MXSMCLSE  = MAX(MXSMCLSE,NSMCLSE)
          MXSMCLSE1 = MAX(MXSMCLSE1,NSMCLSE1)
          MXNSTR_AS = MAX(MXNSTR_AS, NSTR_AS)
*
          MXNSTR = MAX(MXNSTR,MXNSTRFSG)
        END DO
        IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) ' MXNSTR, MXNSTR_AS = ',
     &                 MXNSTR, MXNSTR_AS  
        END IF
*        
        CALL ICOPMT(int_mb(KNSTSO(ITP)),NSMST,NSPGPFTP(ITP),
     &              NSTFSMSPGP(1,IBSPGPFTP(ITP)),MXPNSMST,NSPGPFTP(ITP))
*. Corresponding offset array : Each supergroup is generated individually
*. so each supergroup starts with offset 1 !
        CALL ZSPGPIB(int_mb(KNSTSO(ITP)),int_mb(KISTSO(ITP)),
     &                NSPGPFTP(ITP),NSMST)
*
        IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) 
     &    ' Info on strings for type = ', ITP
          IF(NTEST.GE.10) THEN
          WRITE(6,*) 
     &    ' Number of strings per sym (row) and supergroup(column)' 
          CALL IWRTMA(int_mb(KNSTSO(ITP)),NSMST,NSPGPFTP(ITP),
     &                NSMST,NSPGPFTP(ITP))
          END IF
          WRITE(6,'(A,3(2X,I8))') ' NSMCLS,NSMCLSE,NSMCLSE1=',
     &                         NSMCLS,NSMCLSE,NSMCLSE1
          WRITE(6,*)
        END IF
*
      END DO
*. Number of electron in each AS for each supergroup
      CALL ZNELFSPGP(IPRNT)
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) 
     &  ' Number of electrons in the GAS spaces for each supergroup '
        WRITE(6,*) 
     &  ' ========================================================= '
        WRITE(6,*)
        DO I = 1, NTSPGP
        WRITE(6,'(A,I4,A, 16I3)')
     &  ' Supergroup: ', I, ' Occupation: ', (NELFSPGP(J,I),J=1,NGAS)
        END DO
      END IF
*. Number of PH-particles per supergroups
      DO ISPGP = 1, NTSPGP
        NPHEL = 0
        IBTYP = 1
        ISPGP_A = ISPGP + IBTYP -1
        DO IGAS = 1, NGAS
          IF(IPHGAS(IGAS).EQ.1) THEN
             NPHEL = NPHEL + NELFSPGP(IGAS,ISPGP)
          ELSE
             NPHEL = NPHEL + NOBPT(IGAS) - NELFSPGP(IGAS,ISPGP)
          END IF
        END DO
        NPHELFSPGP(ISPGP) = NPHEL
      END DO
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
       WRITE(6,*)
       WRITE(6,*) ' Number of particle-hole operators per supergroup '
       WRITE(6,*)
       WRITE(6,*) '  Supergroup      Number of ph operators '
       WRITE(6,*) ' ----------------------------------------'
       DO ISPGP = 1, NTSPGP
         WRITE(6,'(6X,I3,12X,I3)') ISPGP, NPHELFSPGP(ISPGP)
       END DO
      END IF
*
* Number of holes per supergroup
      DO IISPGP = 1, NTSPGP
        NHOLE = 0
        DO IGAS = 1, NGAS
          IF(IPHGAS(IGAS).EQ.2) NHOLE = NHOLE + NELFSPGP(IGAS,IISPGP)
        END DO
        NHLFSPGP(IISPGP) = NHOLE
      END DO
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
      WRITE(6,*) ' Number of electrons in hole spaces per supergroup '
      CALL IWRTMA(NHLFSPGP,1,NTSPGP,1,NTSPGP)
      END IF
*. Largest number of strings belonging to given supergroup
*. Largest Occupation block for given supergroup and sym
      MAX_STR_OC_BLK = -1
      MAX_STR_PHOC_BLK = -1
      MAX_STR_SPGP = 0  
      DO ISPGP = 1, NTSPGP
        NSTR = IELSUM(NSTFSMSPGP(1,ISPGP),NSMST)
        MAX_STR_SPGP = MAX(MAX_STR_SPGP,NSTR)
        NEL = IELSUM(NELFSPGP(1,ISPGP),NGAS)
        NPHEL = NPHELFSPGP(ISPGP)
        DO ISTSM = 1, NSMST
          MAX_STR_OC_BLK 
     &  = MAX(MAX_STR_OC_BLK,NEL*NSTFSMSPGP(ISTSM,ISPGP))
C?   &  = MAX(MAX_STR_OC_BLK,(NEL+4)*NSTFSMSPGP(ISTSM,ISPGP))
          MAX_STR_PHOC_BLK = 
     &    MAX(MAX_STR_PHOC_BLK,NPHEL*NSTFSMSPGP(ISTSM,ISPGP))
        END DO
      END DO 
*     
      IF(NTEST.GE.2.and.ga_nodeid().eq.0) THEN
      WRITE(6,*)
     & ' Largest number of strings of given supergroup   ',
     &   MAX_STR_SPGP
      WRITE(6,*) 
     & ' Largest block of string occupations             ',
     &   MAX_STR_OC_BLK
      WRITE(6,*) 
     & ' Largest block  of string ph occupations         ',
     &    MAX_STR_PHOC_BLK
*
      WRITE(6,*)
     & ' Largest number of strings of given supergroup and sym', MXNSTR
      END IF
C?    WRITE(6,'(A,3I6)') ' MXSMCLS,MXSMCLSE,MXSMCLSE1 = ',
C?   &                     MXSMCLS,MXSMCLSE,MXSMCLSE1
*
*
* Possible occupation classes
*
COLD  CALL OCCLS(2,NMXOCCLS,WORK(KIOCLS),NACTEL,NGAS,
COLD &           IGSOCC(1,1),IGSOCC(1,2),
COLD &           0,0,NOBPT)
*
* Maps creation/annihilation of given gas orb from given supergroup
* gives new supergroup. 
*
C?    WRITE(6,*) ' IBSPGPFTP: '
C?    CALL IWRTMA(IBSPGPFTP,1,NSTTYP,1,NSTTYP) 
      IZERO = 0
      CALL ISETVC(int_mb(KSPGPCR),IZERO,NGAS*NTSPGP)
      CALL ISETVC(int_mb(KSPGPAN),IZERO,NGAS*NTSPGP)
*
      CALL LUCIAQENTER('SPGPA')
      DO ISTTYP = 1,NSTTYP
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) 
     &WRITE(6,*)  ' Info for SPGP_AC for ISTTYP = ', ISTTYP
      IF(NSPGPFTP(ISTTYP).GT.0) THEN
        IIEL = NELFTP(ISTTYP)
*
*. Creation map from this type
*
*. Type of string with one elec more
        ISTTYPC = 0
        DO JSTTYP = 1, NSTTYP
          IF(MOD(ISTTYP,2).EQ.MOD(JSTTYP,2).AND.
     &       NELFTP(JSTTYP).EQ.IIEL+1           ) ISTTYPC = JSTTYP
        END DO
        IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) ' ISTTYP and ISTTYPC ',ISTTYP,ISTTYPC
          WRITE(6,'(A, 2I4)') ' NEL(ISTTYP),NEL(ISTTYPC) = ',
     &    NELFTP(ISTTYP), NELFTP(ISTTYPC)
        END IF
*
        IF(ISTTYPC.GE.1.AND.NSPGPFTP(ISTTYPC).GT.0) THEN
           CALL SPGP_AC(NELFSPGP(1,1), NSPGPFTP(ISTTYP),
     &                  NELFSPGP(1,1),NSPGPFTP(ISTTYPC),
     &                  NGAS,MXPNGAS,2,int_mb(KSPGPCR),
     &                  IBSPGPFTP(ISTTYP),IBSPGPFTP(ISTTYPC))
        END IF
*
*. Annihilation maps
*
        ISTTYPA = 0       
        DO JSTTYP = 1, NSTTYP
          IF(MOD(ISTTYP,2).EQ.MOD(JSTTYP,2).AND.
     &       NELFTP(JSTTYP).EQ.IIEL-1           ) ISTTYPA = JSTTYP
        END DO
        IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) 'ISTTYP, ISTTYPA', ISTTYP,ISTTYPA
          WRITE(6,'(A, 2I4)') ' NEL(ISTTYP),NEL(ISTTYPA) = ',
     &    NELFTP(ISTTYP), NELFTP(ISTTYPA)
        END IF
        IF(ISTTYPA.GE.1 .AND.NSPGPFTP(ISTTYPA).GT.0) THEN
           CALL SPGP_AC(NELFSPGP(1,1), NSPGPFTP(ISTTYP),
     &                  NELFSPGP(1,1),NSPGPFTP(ISTTYPA),
     &                  NGAS,MXPNGAS,1,int_mb(KSPGPAN),
     &                  IBSPGPFTP(ISTTYP),IBSPGPFTP(ISTTYPA))
        END IF
      END IF ! type has a nonvanishing number of strings
      END DO !loop over types
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        NSPGP_TOT = IELSUM(NSPGPFTP,NSTTYP)
        WRITE(6,*) ' NSPGP_TOT = ', NSPGP_TOT
        WRITE(6,*) ' Annihilation mappings of supergroups: '
        CALL IWRTMA(int_mb(KSPGPAN),NSPGP_TOT,1,NSPGP_TOT,1)
        WRITE(6,*) ' Creation mappings of supergroups: '
        CALL IWRTMA(int_mb(KSPGPCR),NSPGP_TOT,1,NSPGP_TOT,1)
      END IF
*
      CALL LUCIAQEXIT('SPGPA')
*.AB supergroups in compound space
      IATP = 1
      IBTP = 2
      NASPGP = NSPGPFTP(IATP)
      NBSPGP = NSPGPFTP(IBTP)
      IB_A = IBSPGPFTP(IATP)
      IB_B = IBSPGPFTP(IBTP)
      CALL ABSPGP_TO_CISPACE(ISPGPFTP(1,IB_A),NASPGP,
     &     ISPGPFTP(1,IB_B),NBSPGP,NGAS,IGSOCC,NELFGP,
     &     N_ABSPGP_MAX,int_mb(KIABSPGP_FOR_CMPSPC),0)
*. Obtain mapping from occupation classes to supergroups
      IATP = 1
      IBTP = 2
      NASPGP = NSPGPFTP(IATP)
      NBSPGP = NSPGPFTP(IBTP)
      IB_A = IBSPGPFTP(IATP)
      IB_B = IBSPGPFTP(IBTP)
*
*. Initiallize two variables with MEMMAN related to the parallelization of 
*. ABSPGP_TO_OCCLS
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ABSPGP')
      CALL MEMMAN(KI_ABSPGP_A,N_ABSPGP_MAX,'ADDL  ',1,'ABSPGP_A')
      CALL MEMMAN(KI_ABSPGP_B,N_ABSPGP_MAX,'ADDL  ',1,'ABSPGP_B')
      CALL MEMMAN(KI_ABSPGP_VOG,NOCCLS_MAX,'ADDL  ',1,'ABSPGP_V')
      CALL MEMMAN(KI_ABSPGP,ga_nnodes(),'ADDL  ',1,'NABSPGP')
*
!      if (ga_nodeid().eq.0) write(*,*) 'VOG: N_ABSPGP_MAX', N_ABSPGP_MAX
!      if (ga_nodeid().eq.0) write(*,*) 'VOG: KI_ABSPGP_VOG', NOCCLS_MAX
*. First: offsets
      CALL ABSPGP_TO_OCCLS(
     &     ISPGPFTP(1,IB_A),NASPGP,ISPGPFTP(1,IB_B),NBSPGP,
     &     NGAS,IGSOCC,NELFGP,
     &     NOCCLS_MAX,int_mb(KIOCCLS),N_ABSPGP_TOT2,
     &     int_mb(KIABSPGP_FOR_OCCLS),int_mb(KNABSPGP_FOR_OCCLS),
     &     int_mb(KIBABSPGP_FOR_OCCLS),1,ga_nnodes(),N_ABSPGP_MAX,
     &     int_mb(KI_ABSPGP_A),int_mb(KI_ABSPGP_B),
     &     int_mb(KI_ABSPGP),int_mb(KI_ABSPGP_VOG))
*. And then the actual mappings
      CALL ABSPGP_TO_OCCLS(
     &     ISPGPFTP(1,IB_A),NASPGP,ISPGPFTP(1,IB_B),NBSPGP,
     &     NGAS,IGSOCC,NELFGP,
     &     NOCCLS_MAX,int_mb(KIOCCLS),N_ABSPGP_TOT2,
     &     int_mb(KIABSPGP_FOR_OCCLS),int_mb(KNABSPGP_FOR_OCCLS),
     &     int_mb(KIBABSPGP_FOR_OCCLS),0,ga_nnodes(),N_ABSPGP_MAX,
     &     int_mb(KI_ABSPGP_A),int_mb(KI_ABSPGP_B),
     &     int_mb(KI_ABSPGP),int_mb(KI_ABSPGP_VOG))
C?     WRITE(6,*) ' First two elements of KIABSPGP_FOR_OCCLS as int'
C?     CALL IWRTMA(WORK(KIABSPGP_FOR_OCCLS),1,2,1,2)
C ABSPGP_TO_OCCLS(
C      &           IASPGP,NASPGP,IBSPGP,NBSPGP,
C      &           NGAS,IGSOCC_MNMX,NELFGP,
C      &           NOCCLS,IOCCLS,
C      &           N_ABSPGP_TOT,I_ABSPGP_FOR_OCCLS,N_ABSPGP_FOR_OCCLS,
C      &           IB_ABSPGP_FOR_OCCLS,IONLY_NABSPGP)
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'ABSPGP')
*
*. Dimensions of PL strings
*
C     GET_DIM_PLSTRINGS(MXNSTRP,MXNSTRP_AS,MXNSTRL,
C    &           MXNSTRL_AS,MXNSTRSTRP_AS,MXNSTRSTRL_AS)
      CALL GET_DIM_PLSTRINGS(MXNSTRP,MXNSTRP_AS,MXNSTRL,
     &           MXNSTRL_AS,MXNSTRSTRP_AS,MXNSTRSTRL_AS)
      IF (NTEST.GT.10.and.ga_nodeid().eq.0) THEN
      WRITE(6,*) ' Info on PL expansion of strings and CI-space:'
      WRITE(6,'(A,2I6)') ' MXNSTRP,MXNSTRP_AS = ', MXNSTRP,MXNSTRP_AS
      WRITE(6,'(A,2I6)') ' MXNSTRL,MXNSTRL_AS = ', MXNSTRL,MXNSTRL_AS
      WRITE(6,'(A,2I6)') 
     &' MXNSTRSTRP_AS,MXNSTRSTRL_AS = ',MXNSTRSTRP_AS,MXNSTRSTRL_AS 
*


*
      WRITE(6,*) ' Memory Check at end of STRINF_GAS  '
      CALL MEMCHK
      ENDIF
*
      CALL LUCIAQEXIT('STIN2')
      CALL LUCIAQEXIT('STRIN')
      RETURN
      END

      SUBROUTINE CRESTR_GAS(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,
     &                  Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,NOCOB,
     &                  IPRNT)
*
* A type of strings containing NEL electrons are given
* set up all possible ways of adding an electron to this type of strings
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* IORBOF : Number of first orbital
* Z      : Lexical ordering matrix for output strings containing
*          NEL + 1 electrons
* NEWORD : Reordering array for N+1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL+1 strings
*
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the string have a negative sign
*          then the phase equals - 1
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(NORB,NSTINI),TTO(NORB,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTEST0 =  0
      NTEST = MAX(IPRNT,NTEST0)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' CRESTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
C?    WRITE(6,*) ' Reorder array NEWORD '
C?    CALL IWRTMA(NEWORD,1,NSTINO,1,NSTINO)
      LUOUT = 6
*
      DO 1000 ISTRIN = 1,NSTINI
      IF(NTEST.GE.100) 
     &write(6,*) ' Input string ',istrin,(string(i,istrin),i=1,nel)
        DO 100 IORB = IORBOF, IORBOF-1+NORB
        IF(NTEST.GE.100) write(6,*) ' orbital ',iorb
           IPLACE = 0
           IF(NEL.EQ.0) THEN
             IPLACE = 1
             GOTO 11
           ELSE IF ( NEL .NE. 0 ) THEN
            DO 10 IEL = 1, NEL
              IF(IEL.EQ.1.AND.STRING(1,ISTRIN).GT.IORB) THEN
                IPLACE = 1
                GOTO 11
              ELSE IF( (IEL.EQ.NEL.AND.IORB.GT.STRING(IEL,ISTRIN)) .OR.
     &                 (IEL.LT.NEL.AND.IORB.GT.STRING(IEL,ISTRIN).AND.
     &                  IORB.LT.STRING(IEL+1,ISTRIN)) ) THEN
                IPLACE = IEL+1
                GOTO 11
              ELSE IF(STRING(IEL,ISTRIN).EQ.IORB) THEN
                IPLACE = 0
                GOTO 11
              END IF
   10       CONTINUE
           END IF
   11     CONTINUE
*
C?        write(6,*) ' iplace = ', iplace   
          IF(IPLACE.NE.0) THEN
*. Generate next string
            DO 30 I = 1, IPLACE-1
   30       STRIN2(I) = STRING(I,ISTRIN)
            STRIN2(IPLACE) = IORB
            DO 40 I = IPLACE,NEL
   40       STRIN2(I+1) = STRING(I,ISTRIN)
C?          write(6,*) ' updated string (STRIN2) '
C?          call iwrtma(STRIN2,1,NEL+1,1,NEL+1)
            JSTRIN = ISTRNM(STRIN2,NOCOB,NEL+1,Z,NEWORD,1)
C?          write(6,*) ' corresponding number ', JSTRIN
*
            TTO(IORB-IORBOF+1,ISTRIN) = JSTRIN
            IIISGN = (-1)**(IPLACE-1)
            IF(LSGSTR.NE.0)
     &      IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
            IF(IIISGN .EQ. -1 )
     &      TTO(IORB-IORBOF+1,ISTRIN) = - TTO(IORB-IORBOF+1,ISTRIN)
            TI(IORB-IORBOF+1,ISTRIN ) = IORB
          END IF
  100   CONTINUE
*
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(LUOUT,*) ' Output from CRESTR : '
        WRITE(LUOUT,*) '==================='
*
        WRITE(6,*)
        WRITE(LUOUT,*) ' Strings with an electron added  '
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' New strings.. ',
     &     (TTO(I,ISTRIN),I = 1,NORB)
        END DO
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' orbitals added or removed ' ,
     &     (TI(I,ISTRIN),I = 1,NORB)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE ZSTINF_GAS(IPRNT)
*
* Set up common block /STINF/ from information in /STINP/
*
*=========
* Input
*=========
* Information in /CGAS/ and /GASSTR/  
*
*======================
* Output ( in /STINF/ )
*======================
* ISTAC (MXPSTT,2) : string type obtained by creating (ISTAC(ITYP,2))
*                    or annihilating (ISTAC(ITYP,1)) an electron
*                    from a string of type  ITYP . A zero indicates
*                    that this mapping is not included
*                    Only strings belonging to the same 
*                    Orbital group are mapped
*                    mapped
*. Input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
C     COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
C    &            NGSOB(MXPOBS,MXPNGAS),
C    &            NGSOBT(MXPNGAS),IGSOCC(2,MXPNGAS),IGSINA,IGSDEL
C     COMMON/GASSTR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),
C    &              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),
C    &              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS)
*. Output
      INCLUDE 'stinf.inc'
*./STINF/ 
C     COMMON/STINF/ISTAC(MXPSTT,2),NOCTYP(MXPSTT),NSTFTP(MXPSTT),
C    &             INUMAP(MXPSTT),INDMAP(MXPSTT)
*. Only the first element, i.e. ISTAC  is defined

*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
* ******************************************************************
* Mappings between strings with the same type ISTTP index , +/- 1 el
* ******************************************************************
      CALL ISETVC(ISTAC,0,2*MXPSTT)
      DO  IGAS = 1, NGAS     
*. groups for a given gas spaces goes with increasing number of orbitals,
*  so the first space does not have any creation mapping
*  and the last space does not have any annihilation mapping
*
        MGRP = NGPSTR(IGAS)
        DO IGRP = 1, MGRP
          IIGRP = IGRP + IBGPSTR(IGAS) -1 
          IF(IGRP.NE.1) THEN
*. Annihilation map is present : IIGRP => IIGRP - 1
            ISTAC(IIGRP,1) = IIGRP -1
          END IF
          IF(IGRP.NE.MGRP) THEN
*. Creation map is present : IIGRP => IIGRP + 1
             ISTAC(IIGRP,2) = IIGRP + 1
          END IF
        END DO
      END DO
*
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Group-group mapping array ISTAC '
        WRITE(6,*) ' =============================== '
        CALL IWRTMA(ISTAC,NGRP  ,2,MXPSTT,2)
      END IF
*
      RETURN
      END
      SUBROUTINE MEMSTR_GAS
*
*
* Construct pointers for saving information about strings and
* their mappings
*
* GAS version 
*
*========
* Input :
*========
* Number and groups of strings defined by /GASSTR/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in common block /STRBAS/
*
* Jeppe Olsen , Winter of 1994
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'stinf.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
*. Start of string information
      CALL MEMMAN(KSTINF,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      NTEST = 00
      IF(NTEST.GE.10)
     &WRITE(6,*) ' First word with string information',KSTINF
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)

        IF(NTEST.GE.10) THEN
        WRITE(6,*) ' IGRP, NSTRIN, NELFGP = ',
     &               IGRP, NSTRIN, NELFGP(IGRP)
        END IF

        CALL MEMMAN(KOCSTR(IGRP),LSTRIN,'ADDS  ',1,'OCSTR ')
        CALL MEMMAN(KSTREO(IGRP),NSTRIN,'ADDS  ',1,'STREO ')
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      CALL MEMMAN(KNSTSGP(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
      CALL MEMMAN(KISTSGP(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
*
*. Number of strings per symmetry and offset for strings of given sym
*. for types 
*
      DO  ITP  = 1, NSTTP  
        CALL MEMMAN(KNSTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'NSTSO ')
        CALL MEMMAN(KISTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'ISTSO ')
      END DO
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP  
        CALL MEMMAN(KZ(IGRP),NOCOB*NELFGP(IGRP),'ADDS  ',1,'Zmat  ')
      END DO
*
*. Mappings between different groups
*
      DO  IGRP = 1, NGRP   
        IEL = NELFGP(IGRP)
        IGAS = IGSFGP(IGRP)
        IORB = NOBPT(IGAS)
        ISTRIN = NSTFGP(IGRP)
*. IF creation is involve : Use full orbital notation
*  If only annihilation is involved, compact form will be used
        IF(ISTAC(IGRP,2).NE.0) THEN
          LENGTH = IORB*ISTRIN 
          CALL MEMMAN(KSTSTM(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
          CALL MEMMAN(KSTSTM(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
        ELSE IF(ISTAC(IGRP,1).NE.0) THEN
*. Only annihilation map so
          LENGTH = IEL*ISTRIN
          CALL MEMMAN(KSTSTM(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
          CALL MEMMAN(KSTSTM(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
        ELSE
*. Neither annihilation nor creation (?!)
          KSTSTM(IGRP,1) = -1
          KSTSTM(IGRP,2) = -1
        END IF
      END DO
*
*. Symmetry of conjugated orbitals and orbital excitations
*
*     KCOBSM,KNIFSJ,KIFSJ,KIFSJO
      CALL MEMMAN(KCOBSM,NACOB,'ADDS  ',1,'Cobsm ')
      CALL MEMMAN(KNIFSJ,NACOB*NSMSX,'ADDS  ',1,'Nifsj ')
      CALL MEMMAN(KIFSJ,NACOB**2 ,'ADDS  ',1,'Ifsj  ')
      CALL MEMMAN(KIFSJO,NACOB*NSMSX,'ADDS  ',1,'Ifsjo ')
*
*. Symmetry of excitation connecting  strings of given symmetry
*
      CALL MEMMAN(KSTSTX,NSMST*NSMST,'ADDS  ',1,'Ststx ')
*
*. Occupation classes
*
      CALL MEMMAN(KIOCLS,NMXOCCLS*NGAS,'ADDS  ',1,'IOCLS ')
*. Annihilation/Creation map of supergroup types
      CALL MEMMAN(KSPGPAN,NTSPGP*NGAS,'ADDS  ',1,'SPGPAN')
      CALL MEMMAN(KSPGPCR,NTSPGP*NGAS,'ADDS  ',1,'SPGPCR')
*
*. Last word of string information
*
      CALL MEMMAN(KSTINE,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      IF(NTEST.GE.10)
     &WRITE(6,*) ' Last word with string information',KSTINE-1
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' R*8 words used for storing Strinformation ',
     &               KSTINE - KSTINF
      END IF
*
      RETURN
      END
      SUBROUTINE NSTRSO_GAS(NEL,NORB1,NORB2,NORB3,
     &                  NELMN1,NELMX1,NELMN3,NELMX3,
     &                  IOC,NORB,NSTASO,ISTASO,
     &                  NOCTYP,NSMST,IOTYP,IPRNT)
*
* Number of strings per symmetry for group IOTYP
*
* Gas version, no check of type : set to 1
*
* Jeppe Olsen Winter of 1994
*
      IMPLICIT REAL*8           ( A-H,O-Z)
C     DIMENSION IOC(*),NSTASO(NOCTYP,NSMST)
      DIMENSION IOC(*),NSTASO(NSMST,*),ISTASO(NSMST,*)
*
      CALL ISETVC(NSTASO(1,IOTYP),0,NSMST)
      NTEST0 = 0
      NTEST = MAX(IPRNT,NTEST0)
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
* Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
* Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( NTEST .GE.500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
         IFRST2 = 1
         IFRST3 = 1
* Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
           IF( NTEST .GE.500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
           IFRST3 = 1
* Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
             IF( NTEST .GE. 500) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry of string
             ISYM = ISYMST(IOC,NEL)
C                   ISYMST(STRING,NEL)
*. occupation type of string
COLD         ITYP = IOCTP2(IOC,NEL,IOTYP)
C                   IOCTP2(STRING,NEL)
*
             NSTASO(ISYM,IOTYP) = NSTASO(ISYM,IOTYP)+ 1
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
*. The corresponding offset
*
      DO ISM = 1, NSMST
        IF(ISM.EQ.1) THEN
          ISTASO(ISM,IOTYP) = 1
        ELSE
          ISTASO(ISM,IOTYP) = ISTASO(ISM-1,IOTYP)+NSTASO(ISM-1,IOTYP)
        END IF
      END DO
 
      IF(NTEST.GE.5)
     &WRITE(6,*) ' Number of strings generated   ', NSTRIN
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Number of strings per sym for group = ', IOTYP 
        WRITE(6,*) '================================================'
        CALL IWRTMA(NSTASO(1,IOTYP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Offset for given symmetry for group = ', IOTYP 
        WRITE(6,*) '================================================'
        CALL IWRTMA(ISTASO(1,IOTYP),1,NSMST,1,NSMST)
      END IF
C
      RETURN
      END
      SUBROUTINE NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL,MAXVAL(I). 
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen Oct 1994
*
*. Input
      DIMENSION MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
       NONEW = 1
       GOTO 1001
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
C08-08-01   INUM(JPLACE) = 1               
            INUM(JPLACE) = MINVAL         
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
      SUBROUTINE GENSTR_GAS(NEL,NELMN1,NELMX1,NELMN3,NELMX3,
     &                  ISTASO,IGRP,NOCTYP,NSMST,Z,LSTASO,
     &                  IREORD,STRING,IOC,IOTYP,IPRNT)
*
* Generate strings consisting of  NEL electrons fullfilling
*   1 : Between NELMN1 AND NELMX1 electrons in the first NORB1 orbitals
*   2 : Between NELMN3 AND NELMX3 electrons in the last  NORB3 orbitals
*
* In the present version the strings are directly ordered into
* symmetry and occupation type .
*
* Jeppe Olsen Winter of 1990
*
* Special GAS version, Winter of 94 All strings of group IGRP
*
* ========
* Output :
* ========
* STRING(IEL,ISTRIN) : Occupation of strings.
* IREORD             : Reordering array going from lexical
*                      order to symmetry and occupation type order.
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Input
      DIMENSION ISTASO(NSMST,*)
      INTEGER Z(NACOB,NEL)
*.Orbinp
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*
*.Output
      INTEGER STRING(NEL,*),IREORD(*)
*.Scratch arrays
      DIMENSION IOC(*),LSTASO(NOCTYP,NSMST)
*
      CALL ISETVC(LSTASO,0,NOCTYP*NSMST)
      NTEST0 = 00
      NTEST = MAX(NTEST0,IPRNT)
      IF( NTEST .GE. 10) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' GENSTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
        WRITE(6,*) ' NELMN1, NELMX1,  = ', NELMN1, NELMX1
        WRITE(6,*) ' NELMN3, NELMX3,  = ', NELMN3, NELMX3
        WRITE(6,*) ' NEL, NORB1, NORB2, NORB3 = ', 
     &               NEL, NORB1, NORB2, NORB3
        WRITE(6,*)
      END IF
*
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
* Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
* Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( NTEST .GE. 500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
         IFRST2 = 1
         IFRST3 = 1
* Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
           IF( NTEST .GE. 500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
           IFRST3 = 1
* Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
             IF( NTEST .GE. 500 ) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry
*                   ISYMST(STRING,NEL)
             ISYM = ISYMST(IOC,NEL)
*. Occupation type
C            ITYP = IOCTP2(IOC,NEL,IOTYP)
             ITYP = 1
*
             IF(ITYP.NE.0) THEN
               LSTASO(ITYP,ISYM) = LSTASO(ITYP,ISYM)+ 1
C                      ISTRNM(IOCC,NACTOB,NEL,Z,NEWORD,IREORD)
               LEXCI = ISTRNM(IOC,NOCOB,NEL,Z,IREORD,0)
               LACTU = ISTASO(ISYM,IGRP)-1+LSTASO(ITYP,ISYM)
               IREORD(LEXCI) = LACTU
               IF(NTEST.GT.10) WRITE(6,*) ' LEXCI,LACTU',
     &         LEXCI,LACTU
               CALL ICOPVE(IOC,STRING(1,LACTU),NEL)
             END IF
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
      IF(NTEST.GE.1 ) THEN
        WRITE(6,*) ' Number of strings generated   ', NSTRIN
      END IF
      IF(NTEST.GE.10)  THEN
        IF(NTEST.GE.100) THEN
          NPR = NSTRIN
        ELSE
          NPR = MIN(NSTRIN,50)
        END IF
        WRITE(6,*) ' Strings generated '
        WRITE(6,*) ' =================='
        ISTRIN = 0
        DO 100 ISYM = 1, NSMST
        DO 100 ITYP = 1,NOCTYP
          LSTRIN = MIN(LSTASO(ITYP,ISYM),NPR-ISTRIN)
          IF(LSTRIN.GT.0) THEN
            WRITE(6,*) ' Strings of type and symmetry ',ITYP,ISYM
            DO 90 KSTRIN = 1,LSTRIN
              ISTRIN = ISTRIN + 1
              WRITE(6,'(2X,I5,A,8X,(10I5))')
     &        ISTRIN,':',(STRING(IEL,ISTRIN),IEL = 1,NEL)
   90       CONTINUE
          END IF
  100   CONTINUE
*
        WRITE(6,*) ' Array giving actual place from lexical place'
        WRITE(6,*) ' ============================================'
        CALL IWRTMA(IREORD,1,NPR,1,NPR)
      END IF
 
      RETURN
      END
      SUBROUTINE LCISPC(IPRNT)
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1994/1995 ( woops !)
*               MXSOOB_AS added,MXSB removed May 1999
*
* GAS VERSION
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
* ===================
*.Input common blocks
* ===================
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
      INCLUDE 'cicisp.inc'
      CALL LUCIAQENTER('LCISP')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LCISPC')
*
*
*. Number of spaces
      NICISP = NCMBSPC
C?    write(6,*) ' LCISPC : NICISP ', NICISP
*. Type of alpha- and beta strings
      IATP = 1
      IBTP = 2
*
      NOCTPA =  NOCTYP(IATP)
      NOCTPB =  NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*.Local memory
      CALL MEMMAN(KLBLTP,NSMST,'ADDL  ',2,'KLBLTP') !done
      IF(IDC.EQ.3 .OR. IDC .EQ. 4 )
     &CALL MEMMAN(KLCVST,NSMST,'ADDL  ',2,'KLCVST') !done
      CALL MEMMAN(KLIOIO,NOCTPA*NOCTPB,   'ADDL  ',2,'KLIOIO') !d
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
      IF(IDC.EQ.3.OR.IDC.EQ.4)
     &CALL SIGVST(dbl_mb(KLCVST),NSMST)
 
*. Array defining symmetry combinations of internal strings
*. Number of internal dets for each symmetry
        CALL SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
*. MXSB is not calculated anymore, set to 0
      MXSB = 0
*
      MXSOOB = 0
      MXSOOB_AS = 0
      XMXSOOB = 0.0D0
      XMXSOOB_AS = 0.0D0
      DO 100 ICI = 1, NICISP
*. allowed combination of types
      CALL IAIBCM(ICI,dbl_mb(KLIOIO))

      DO  50 ISYM = 1, NSMCI
          CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,dbl_mb(KLBLTP),
     $               dbl_mb(KLCVST))
          CALL NGASDT(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),NGAS,ISYM,
     &                NSMST,NOCTPA,NOCTPB,int_mb(KNSTSO(IATP)),
     &                int_mb(KNSTSO(IBTP)),
     &                ISPGPFTP(1,IBSPGPFTP(IATP)),
     &                ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NELFGP,
     &                NCOMB,XNCOMB,MXS,MXSOO,dbl_mb(KLBLTP),NTTSBL,
     &                LCOL,dbl_mb(KLIOIO),MXSOO_AS,XMXSOO,XMXSOO_AS)
      
          XISPSM(ISYM,ICI) = XNCOMB
          MXSOOB = MAX(MXSOOB,MXSOO)
          XMXSOOB = MAX(XMXSOOB,XMXSOO)
          MXSB = MAX(MXSB,MXS)
          MXSOOB_AS = MAX(MXSOO_AS,MXSOOB_AS)
          XMXSOOB_AS = MAX(XMXSOO_AS,XMXSOOB_AS)
          
          NBLKIC(ISYM,ICI) = NTTSBL
          LCOLIC(ISYM,ICI) = LCOL
   50 CONTINUE
  100 CONTINUE
*
      NTEST = 10
      NTEST = MAX(NTEST,IPRNT)
      if (NTEST.GE.10.and.ga_nodeid().eq.0) then
      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,*) 
     &' Dimension of CI-expansions per symmetry '
      WRITE(6,*)
     & ' ====================================== '
*
      DO 200 ICI = 1, NCMBSPC
          WRITE(6,*) ' CI space ', ICI
          WRITE(6,'(1H , 4E22.15)') (XISPSM(II,ICI),II=1,NSMCI)
C         CALL WRTMAT(XISPSM(1,ICI),1,NSMCI,1,NSMCI)
  200 CONTINUE
      WRITE(6,*)
      WRITE(6,*) ' Largest Symmetry-type-type block ',MXSOOB
      WRITE(6,*) ' Largest type-type block (all symmetries) ',MXSOOB_AS
      WRITE(6,*)
      WRITE(6,'(A,E22.13)')  
     & ' Largest Symmetry-type-type block(as real)', XMXSOOB
      WRITE(6,'(A,E22.13)')  
     & ' Largest type-type block (all symmetries)(as real) ',XMXSOOB_AS
      WRITE(6,*)
      endif
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) 
     &  ' Number of TTS subblocks per CI expansion '
        WRITE(6,*)
     &   ' ======================================== '
*
        DO  ICI = 1,  NCMBSPC
            WRITE(6,*) ' Internal CI space ', ICI
            CALL IWRTMA(NBLKIC(1,ICI),1,NSMCI,1,NSMCI)
        END DO
      END IF
*. Largest number of BLOCKS in a CI expansion
      MXNTTS = 0
      DO ICI = 1,NCMBSPC
       DO ISM =1, NSMCI
        MXNTTS = MAX(MXNTTS,NBLKIC(ISM,ICI))
       END DO
      END DO
*
      if (NTEST.GE.10.AND.ga_nodeid().eq.0) WRITE(6,*) 
     &  ' Largest number of blocks in CI expansion',MXNTTS
*
      IF(NTEST.GE.5.and.ga_nodeid().eq.0) THEN
      WRITE(6,*) 
     &' Number of columns per CI expansion '
      WRITE(6,*)
     & ' =================================== '
*
      DO  ICI = 1,  NCMBSPC
          WRITE(6,*) ' Internal CI space ', ICI
          CALL IWRTMA(LCOLIC(1,ICI),1,NSMCI,1,NSMCI)
      END DO
      END IF
*
      CALL LUCIAQEXIT('LCISP')
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'LCISPC')
*
      RETURN
      END
      SUBROUTINE NGASDT(IOCCMN,IOCCMX,NGAS,ITOTSM,
     &                  NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,
     &                  IAOCC,IBOCC,MXPNGAS,NELFGP,
     &                  NCOMB,XNCOMB,MXSB,MXSOOB,
     &                  IBLTP,NTTSBL,LCOL,IOCOC,MXSOOB_AS,
     &                  XMXSOOB,XMXSOOB_AS)
*
* Number of combimations with symmetry ITOTSM and
* occupation determined by IOCOC.
* IOCCMN, IOCCMX are included but are not active
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
* NTTSBL is number of TTS blocks in vector
* LCOL is the sum of the number of columns in each block
*
*
* Winter 94/95        
* May 1999 : Loops restructrured to sym,type,type (leftmost innerst)
*            MXSB not calculated
*
* Last Modification; Oct. 2012; Jeppe Olsen, Real versions of MXSOOB, MXSOOB_AS
* added
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Allowed combinations of alpha and beta types
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Occupation constraints
      DIMENSION IOCCMN(NGAS),IOCCMX(NGAS)
*. Occupation of alpha and beta strings
      DIMENSION IAOCC(MXPNGAS,*),IBOCC(MXPNGAS,*)
*. Number of strings per supergroup and symmetry
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*),NELFGP(*)
*. block types
      DIMENSION IBLTP(*)
*
      CALL LUCIAQENTER('NGASD')
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' NGASDT speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
        WRITE(6,*) ' ITOTSM ', ITOTSM
        WRITE(6,*) ' Upper and lower occupation constraints'
        CALL IWRTMA(IOCCMN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCCMX,1,NGAS,1,NGAS)
        WRITE(6,*) ' IOCOC matrix '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' Number of alpha and beta strings '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
      MXSB = 0
      MXSOOB = 0
      MXSOOB_AS = 0
      XMXSOOB = 0.0D0
      XMXSOOB_AS = 0.0D0
      NCOMB = 0
      XNCOMB = 0.0D0
      NTTSBL = 0
      LCOL = 0
*
      DO 200 IATP = 1, NOCTPA
        DO 100 IBTP = 1, NOCTPB 
*
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' Alpha super group and beta super group'
            CALL IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
          END IF
*
          IF(IOCOC(IATP,IBTP).EQ.1) THEN
*
            LTTS_AS = 0
            XLTTS_AS = 0.0D0
            DO 300 IASM = 1, NSMST
              IF(IBLTP(IASM).EQ.0) GOTO 300
              CALL SYMCOM(2,1,IASM,IBSM,ITOTSM)
              IF(IBSM.NE.0) THEN
                IF(IBLTP(IASM).EQ.2) THEN
                  ISYM = 1
                ELSE
                  ISYM = 0
                END IF
                IF(ISYM.EQ.1.AND.IBTP.GT.IATP) GOTO 300
                LASTR = NSSOA(IASM,IATP)
                LBSTR = NSSOB(IBSM,IBTP)
*. Size of unpacked block
                LTTSUP =  LASTR*LBSTR                         
                XLTTSUP =  DBLE(LASTR)*DBLE(LBSTR)
*. Size of packed block
                IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                  LTTSBL = LASTR*LBSTR          
                  XLTTSBL = DBLE(LASTR)*DBLE(LBSTR)
                  XNCOMB = XNCOMB + DBLE(LASTR)*DBLE(LBSTR)        
                ELSE
                  LTTSBL = LASTR*(LASTR+1)/2
                  XLTTSBL = DBLE(LASTR+1)*DBLE(LASTR)/2
                  XNCOMB = XNCOMB + DBLE(LASTR+1)*DBLE(LASTR)/2        
                END IF
                LTTS_AS = LTTS_AS + LTTSUP
                XLTTS_AS = XLTTS_AS + XLTTSUP
                NCOMB = NCOMB + LTTSBL
                MXSOOB = MAX(MXSOOB,LTTSUP)
                XMXSOOB =MAX(XMXSOOB,XLTTSUP)
                NTTSBL = NTTSBL + 1
                LCOL = LCOL + NSSOB(IBSM,IBTP)
              END IF
  300       CONTINUE
            MXSOOB_AS = MAX(MXSOOB_AS,LTTS_AS)
            XMXSOOB_AS = MAX(XMXSOOB_AS,XLTTS_AS)
          END IF
  100   CONTINUE
  200 CONTINUE
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' NGASDT : NCOMB XNCOMB,NTTSBL,MXSOOB', 
     &               NCOMB,XNCOMB,NTTSBL,MXSOOB
      END IF
*
      CALL LUCIAQEXIT('NGASD')
*
      RETURN
      END
      SUBROUTINE OCCLS_OLD(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in the 
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and 
* IOCCLS = allowed distributions of electrons
*
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IGSMIN(NGAS),IGSMAX(NGAS)
*. Output
      DIMENSION  IOCCLS(NGAS,*)
*. Local scratch 
      INCLUDE 'mxpdim.inc'
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' NGAS NEL ', NGAS,NEL
      END IF
*
      ISKIP = 0
      NOCCLS = 0
*. start with smallest allowed number 
      DO IGAS = 1, NGAS
        IOCA(IGAS) = IGSMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation 
          CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
        END IF
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
        IF(ISKIP.EQ.1) THEN
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IGSMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
        END IF
C?      WRITE(6,*) ' Another accumulated occupation: ' 
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space 
        NEGA=0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: ' 
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
        IFIRST = 0
*. Correct number of electrons
        IEL = IELSUM(IOC,NGAS)
        IF(IEL.EQ.NEL.AND.NEGA.EQ.0) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Another allowed class : ' 
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
            CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes '
           CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX,
     &                  I_DO_BASSPC,IBASSPC,NOBPT)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in the 
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and 
* IOCCLS = allowed distributions of electrons
*
* Added Oct 98 : IBASSPC
* The basespace of 
* a given class is the first space where this class occurs
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
*. Input
      DIMENSION IGSMIN(NGAS),IGSMAX(NGAS),NOBPT(NGAS)
*. Output
      DIMENSION IOCCLS(NGAS,*)
      DIMENSION IBASSPC(*)
*. Local scratch 
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' NGAS NEL ', NGAS,NEL
      END IF
*
      ISKIP = 1
      NOCCLS = 0
*. start with smallest allowed number 
      DO IGAS = 1, NGAS
        IOCA(IGAS) = IGSMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation 
          CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
        END IF
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
        IF(ISKIP.EQ.1) THEN
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IGSMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
        END IF
C?      WRITE(6,*) ' Another accumulated occupation: ' 
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space 
        NEGA=0
        IM_TO_STUFFED = 0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1 
            IF(IOC(IGAS).GT.2*NOBPT(IGAS)) IM_TO_STUFFED = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: ' 
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
        IFIRST = 0
*. Correct number of electrons
        IEL = IELSUM(IOC,NGAS)
        IF(IEL.EQ.NEL.AND.NEGA.EQ.0.AND.IM_TO_STUFFED.EQ.0) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Another allowed class : ' 
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
            CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
*
            IF(I_DO_BASSPC.EQ.1) THEN
              IBASSPC(NOCCLS) = IBASSPC_FOR_CLS(IOC)
            END IF
*
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes : '
           WRITE(6,*) ' ===================='
           WRITE(6,*)
           WRITE(6,*) ' Class    Occupation in GASpaces '
           WRITE(6,*) ' ================================'
           DO I = 1, NOCCLS
             WRITE(6,'(1H ,I5,3X,16I3)')
     &       I, (IOCCLS(IGAS,I),IGAS=1, NGAS)
           END DO
C          CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
*
      IF(I_DO_BASSPC.EQ.1) THEN
C       WRITE(6,*) ' Base CI spaces for the classes '
C       CALL IWRTMA(IBASSPC,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END 
      SUBROUTINE NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
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
      IF(NELMNT.EQ.0) THEN
        NONEW = 1
        GOTO 1001
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
      SUBROUTINE GASANA(C,NBLOCK,IBLOCK,IBLTP,LUC,ICISTR)
*
*
*
* Analyze CI vector 
*
* Jeppe Olsen, August 1995 
*
* Driven By IBLOCK, May 1997
*           String occupations added, Feb. 98
*
c      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      DIMENSION IBLOCK(8,NBLOCK),IBLTP(*)
      integer C
      CALL LUCIAQENTER('ANACI')
      CALL MEMMAN(KLOFF,DUMMY,'MARK  ',DUMMY,'GASANA')
*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
* Type of alpha and beta strings
      IATP = 1             
      IBTP = 2              
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
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' GASANA speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL '
        WRITE(6,*)   IATP,IBTP,NAEL,NBEL
        WRITE(6,*) ' NOCTPA NOCTPB ', NOCTPA,NOCTPB 
      END IF
*
**. Info on block structure of space
*
*
*. 
*. Analyze CI vector in terms of SD's and occupation classes
*. Number of terms to be printed 
C?    WRITE(6,*) ' GASANA : IPRNCIV = ', IPRNCIV
        IF(IPRNCIV.EQ.0) THEN
          THRES = 0.10 
          MAXTRM = 200
        ELSE
          THRES = 0.0D0
*. Well atmost 100000 coefs - to save disk ..
          MAXTRM = 100000
        END IF
*
      IF( ICISTR .NE. 1 ) CALL REWINO(LUC)
CM    IF(NOCSF.EQ.1) THEN
*. Analyze CI vector in terms of CSF's and configurations
CM    ELSE
        IUSLAB = 0
*. Number of occupation classes 
        IWAY = 1
        NEL = NAEL + NBEL
        CALL OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,
     &             IGSOCC(1,1),IGSOCC(1,2),
     &             0,0,NOBPT)
*. and then the occupation classes 
        CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL') !d
        IWAY = 2
        CALL OCCLS(IWAY,NOCCLS,int_mb(KLOCCLS),NEL,NGAS,
     &             IGSOCC(1,1),IGSOCC(1,2),
     &             0,0,NOBPT)
*
        CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR') !d
        CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR') !d
*
        LENGTH = NOCCLS*10
        CALL MEMMAN(KNCPMT,LENGTH,'ADDL  ',1,'KNCPMT') !done
        CALL MEMMAN(KWCPMT,LENGTH,'ADDL  ',2,'KWCPMT') !done
*
*. Occupation of strings of given sym and supergroup
        CALL GASANAS(C,LUC,int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &       NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &       NBLOCK,IBLOCK,
     &       THRES,MAXTRM,NAEL,NBEL,
     &       int_mb(KLASTR),int_mb(KLBSTR),
     &       IBLTP,NSMST,IUSLAB,
     &       IDUMMY,int_mb(KNCPMT),dbl_mb(KWCPMT),NELFSPGP,      
     &       NOCCLS,NGAS,int_mb(KLOCCLS),ICISTR,NTOOB,NINOB,IPRNCIV)
CM    END IF !Switch between CSF's and SD's
   
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'GASANA')
      CALL LUCIAQEXIT('ANACI')
*
C?    STOP ' Jeppe forced me to stop after INTANA '
      RETURN
      END
      SUBROUTINE GASANAS(C,LUC,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                 MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                 THRES,MAXTRM,NAEL,NBEL,
     &                 IASTR,IBSTR,IBLTP,NSMST,IUSLAB,
     &                 IOBLAB,NCPMT,WCPMT,NELFSPGP,NOCCLS,NGAS,
     &                 IOCCLS,ICISTR,NORB,NINOB,IPRNCIV)
*
* Analyze CI vector :
*
*      1) Print atmost MAXTRM  combinations with coefficients
*         larger than THRES
*         Currently the corresponding dets are not GIVEN !!
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given 
*         occupation class         
*
* Jeppe Olsen , Jan. 1989 ,   
*               Aug 1995 : GAS version
*               May 1997 : BLOCK driven
*                                                  

*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
      IMPLICIT REAL*8(A-H,O-Z)
CNW   DIMENSION C(*)
      integer C
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)  
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IBLTP(*)
      CHARACTER*6 IOBLAB(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
*. Output
      DIMENSION NCPMT(10,NOCCLS)                           
      DIMENSION WCPMT(10,NOCCLS)                          
*
      NTEST = 00
*
      IF(IUSLAB.NE.0) THEN 
       WRITE(6,*)
       WRITE(6,*) 
     & ' Labels for orbitals are of the type n l ml starting with n = 1' 
       WRITE(6,*) 
     & ' so the user should not be  alarmed by labels like 1 f+3 '  
       WRITE(6,*)
      END IF
     
C     WRITE(6,*) 'C(1) = ',C(1)
      MINPRT = 0
      ITRM = 0
      IDET = 0
      IIDET = 0
      ILOOP = 0
      NCIVAR = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
2001  CONTINUE
      IF( ICISTR .GE. 2 ) CALL REWINO(LUC)
      IIDET = 0
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = 1.0D0
        XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNCIV.EQ.1) THEN
*. Print in one shot
       XMAX = 3006.1956
       XMIN = 0.0D0
      END IF
      IDET = 0
      ITDET = 0
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,'(A,E10.4,A,E10.4)')
     &'  Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(6,'(A)')
     &'  ========================================================='
      WRITE(6,*)
*
C?    WRITE(6,*) ' GASANAS : NBLOCK = ',NBLOCK
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
        IF(NTEST.GE.100) THEN
        WRITE(6,'(A,4I4)') 
     &  ' IATP, IBTP, IASM, IBSM = ', IATP, IBTP, IASM, IBSM
        END IF
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN 
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(IDET,1,-1,LUC)
          NO_ZEROING = 1
          CALL FRMDSC2(C,IDET,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
C?        WRITE(6,*) ' Number of elements readin ',IDET
          IDET = 0
        END IF
        IF(IMZERO.NE.1) THEN
*. Obtain alpha strings of sym IASM and type IATP
          IDUM = 0
          CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                             NORB,0,IDUM,IDUM)
*. Obtain Beta  strings of sym IBSM and type IBTP
          IDUM = 0
          CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                             NORB,0,IDUM,IDUM)
*
          IF(IBLTP(IASM).EQ.2) THEN
            IRESTR = 1
          ELSE
            IRESTR = 0
          END IF
*
          NIA = NSSOA(IASM,IATP)
          NIB = NSSOB(IBSM,IBTP)
*
        IBBAS = 1
        IABAS = 1
*

        DO  IB = IBBAS,IBBAS+NIB-1
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB - IBBAS + IABAS
          ELSE
            MINIA = IABAS
          END IF
          DO  IA = MINIA,IABAS+NIA-1
*
!            IF(ILOOP .EQ. 1 ) NCIVAR = NCIVAR + 1
!            IDET = IDET + 1
!            ITDET = ITDET + 1
!            call ga_get(C,IDET,IDET,1,1,XVAL,1)
!CNW         IF( XMAX .GE. ABS(C(IDET)) .AND.
!CNW  &      ABS(C(IDET)).GT. XMIN ) THEN
!            IF( XMAX .GE. ABS(XVAL) .AND.
!     &      ABS(XVAL).GT. XMIN ) THEN
!              ITRM = ITRM + 1
!              IIDET = IIDET + 1
!              IF( ITRM .LE. MAXTRM ) THEN
!CNW             CNORM = CNORM + C(IDET) ** 2
!                CNORM = CNORM + XVAL ** 2
!                WRITE(6,'(A)')
!                WRITE(6,'(A)')
!     &          '                 =================== '
!                WRITE(6,*)
!*
!C                WRITE(6,'(A,I8,A,E14.8)')
!C    &          '  Coefficient of combination ',IDET,' is ',
!                WRITE(6,'(A,I8,A,E14.8)')
!     &          '  Coefficient of combination ',ITDET,' is ',
!     &          XVAL   
!CNW  &          C(IDET)
!                WRITE(6,'(A)')
!     &          '  Corresponding alpha - and beta string '
!                IF(IUSLAB.EQ.0) THEN
!                  IF(NINOB.EQ.0) THEN
!                    WRITE(6,'(4X,10I4)')
!     &              (IASTR(IEL,IA),IEL = 1, NAEL )
!                    WRITE(6,'(4X,10I4)')
!     &              (IBSTR(IEL,IB),IEL = 1, NBEL )
!                  ELSE
!                    WRITE(6,'(4X,A,10I4,(/,
!     &                        16X,10I4))') ' (Inactive) ',
!     &              (IASTR(IEL,IA),IEL = 1, NAEL )
!                    WRITE(6,'(4X,A,10I4,(/,
!     &                        16X,10I4))') ' (Inactive) ',
!     &              (IBSTR(IEL,IB),IEL = 1, NBEL )
!                  END IF
!                ELSE 
!                  WRITE(6,'(4X,10(1X,A6))')
!     &            (IOBLAB(IASTR(IEL,IA)),IEL = 1, NAEL )
!                  WRITE(6,'(4X,10(1X,A6))')
!     &            (IOBLAB(IBSTR(IEL,IB)),IEL = 1, NBEL )
!                END IF
!              END IF
!            END IF
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
        END IF
*       ^ End of if statement for nonvanishing blocks
        END DO
*       ^ End of loop over blocks
       IF(IIDET .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
       IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
*
       WRITE(6,'(A,E15.8)')
     & '  Norm of printed CI vector .. ', CNORM
*
*.Size of CI coefficients
*
*
      IDET = 0
      IF(ICISTR .GE. 2 ) CALL REWINO(LUC)
      CALL ISETVC(NCPMT,0    ,10*NOCCLS)
      CALL SETVEC(WCPMT,0.0D0,10*NOCCLS)
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
*. Occupation class corresponding to given occupation
        JOCCLS = 0
        DO JJOCCLS = 1, NOCCLS
          IM_THE_ONE = 1
          DO IGAS = 1, NGAS
            IF(NELFSPGP(IGAS,IATP-1+IOCTPA)+
     &         NELFSPGP(IGAS,IBTP-1+IOCTPB).NE.IOCCLS(IGAS,JJOCCLS))
     &         IM_THE_ONE = 0
          END DO
          IF(IM_THE_ONE .EQ. 1 ) JOCCLS = JJOCCLS
        END DO
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN 
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(IDET,1,-1,LUC)
          NO_ZEROING = 1
          CALL FRMDSC2(C,IDET,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
          IDET = 0
        END IF
        IF(IMZERO.EQ.0) THEN

        DO IB = IBBAS,IBBAS+NIB-1
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB - IBBAS + IABAS
          ELSE
            MINIA = IABAS
          END IF
          DO IA = MINIA,IABAS+NIA-1
            IDET = IDET + 1
            call ga_get(C,IDET,IDET,1,1,XVAL,1)
            DO IPOT = 1, 10
CNW           IF(10.0D0 ** (-IPOT+1).GE.ABS(C(IDET)).AND.
CNW  &           ABS(C(IDET)).GT. 10.0D0 ** ( - IPOT )) THEN
              IF(10.0D0 ** (-IPOT+1).GE.ABS(XVAL).AND.
     &           ABS(XVAL).GT. 10.0D0 ** ( - IPOT )) THEN
                 NCPMT(IPOT,JOCCLS)= NCPMT(IPOT,JOCCLS)+ 1  
CNW              WCPMT(IPOT,JOCCLS)= WCPMT(IPOT,JOCCLS)+ C(IDET) ** 2
                 WCPMT(IPOT,JOCCLS)= WCPMT(IPOT,JOCCLS)+ XVAL ** 2
              END IF
            END DO
*           ^ End of loop over powers of ten
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
        END IF
*       ^ End of test for novanishing blocks
      END DO
*     ^ End of loop over blocks
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') '   Magnitude of CI coefficients '
      WRITE(6,'(A)') '  =============================='
      WRITE(6,'(A)')
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, 10
        W = 0.0D0
        N = 0
        DO 290 JOCCLS = 1, NOCCLS 
            N = N + NCPMT(IPOT,JOCCLS)                    
            W = W + WCPMT(IPOT,JOCCLS)                    
  290   CONTINUE
        WACC = WACC + W
        NACC = NACC + N
        WRITE(6,'(A,I2,A,I2,3X,I9,X,E15.8,3X,E15.8)')
     &  '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
  300 CONTINUE
*
      WRITE(6,*) ' Number of coefficients less than  10-11',
     &           ' IS  ',NCIVAR - NACC
*
      IF(NOCCLS.NE.1) THEN                      
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Magnitude of CI coefficients for each excitation level '
      WRITE(6,'(A)') 
     & '  ========================================================='
      WRITE(6,'(A)')
      DO 400 JOCCLS = 1, NOCCLS  
          N = 0
          DO 380 IPOT = 1, 10
            N = N + NCPMT(IPOT,JOCCLS)                     
  380     CONTINUE
          IF(N .NE. 0 ) THEN
            WRITE(6,*)
            WRITE(6,'(A,15I3)')'       Occupation of active sets :',
     &      (IOCCLS(IGAS,JOCCLS),IGAS=1, NGAS)
            WRITE(6,'(A,I9)')  
     &      '         Number of coefficients larger than 10-11 ', N
            WRITE(6,*)
            WACC = 0.0D0
            DO 370 IPOT = 1, 10
              N =  NCPMT(IPOT,JOCCLS)                    
              W =  WCPMT(IPOT,JOCCLS)                    
              WACC = WACC + W
              WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
  370       CONTINUE
          END IF 
  400 CONTINUE
*
*. Total weight and number of dets per excitation level
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Total weight and number of SD''s (> 10 ** -11 )  : '          
      WRITE(6,'(A)') 
     & '  ================================================='
      WRITE(6,'(A)')
      WRITE(6,*) '        N      Weight      Acc. Weight   Occupation '
      WRITE(6,*) ' ==================================================='
      WACC = 0.0D0
      DO 500 JOCCLS = 1, NOCCLS
          N = 0
          W = 0.0D0
          DO 480 IPOT = 1, 10
            N = N + NCPMT(IPOT,JOCCLS)                   
            W = W + WCPMT(IPOT,JOCCLS)                   
  480     CONTINUE
          WACC = WACC + W
          IF(N .NE. 0 ) THEN
            WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,16(1X,I2))') 
     &      N,W,WACC,(IOCCLS(IGAS,JOCCLS),IGAS=1,NGAS)
          END IF
  500 CONTINUE
      END IF
*
      RETURN
      END
*
      SUBROUTINE GSDNBBO(RHO1,IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                  NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  NAEL,NBEL,
     &                  IJAGRP,IJBGRP,
     &                  SB,CB,IDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  MXPNGAS,NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,RHO1S)
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
* NGAS : Number of As'es
* IAOC : Occpation of each AS for alpha strings in L
* IBOC : Occpation of each AS for beta  strings in L
* JAOC : Occpation of each AS for alpha strings in R
* JBOC : Occpation of each AS for beta  strings in R
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* MXPNGAS : Largest number of As'es allowed by program
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
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
      DIMENSION C2(*)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      NTESTO= NTEST
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSDNBB :  R block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' GSDNBB :  L block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
*
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in L '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in L '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in R '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in R '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
      END IF
      IACTIVE = 0
*
      IF(NBEL.GE.1.AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
*  beta beta contribution
* =============================
*
        IACTIVE = 1
        CALL GSBBD1(RHO1,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &       NGAS,IBOC,JBOC,
     &       SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
C     GSBBD1(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
C    &                  NGAS,ISEL,ICEL,
C    &                  SB,CB,
C    &                  ADSXA,SXSTST,STSTSX,MXPNGAS
C    &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
C    &                  SSCR,CSCR,I1,XI1S,H,
C    &                  NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
      END IF
      IF(NAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* =============================
*  alpha alpha contribution
* =============================
*
        IACTIVE = 1
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
        CALL GSBBD1(RHO1,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &                   NGAS,IAOC,JAOC,
     &                   SB,CB,
     &                   ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                   NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                   NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMAT(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
      IF(IACTIVE.NE.0.AND. NTEST.GE.50) THEN
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' GSDNBB : updated RHO 1'
        WRITE(6,*) ' =========================='
        CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
      END IF
      NTESTO = NTEST
      RETURN
      END
      SUBROUTINE GSBBD1O(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen 
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces 
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1 : Updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), ITSOB(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*)
      DIMENSION I2(MAXK,*),XI2S(MAXK,*)
*
*.Local arrays
      DIMENSION ITP(16*16),JTP(16*16)
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GSBBD1 in action '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
      IFRST = 1
      JFRST = 1
*
* Type of single excitations that connects the two column strings
      CALL SXTYP_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(NTEST.GE.1000)    
     &WRITE(6,*) ' ISCSM,ICCSM IJSM ', ISCSM,ICCSM,IJSM
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
        DO 800 ISM = 1, NSMOB
*. new i and j so new intermediate strings
          KFRST = 1
*
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.1000) write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
          IBIORB = IOBPTS(ITYP,ISM)
          IBJORB = IOBPTS(JTYP,JSM)
          IF(NTEST.GE.2000)
     &    WRITE(6,*) ' NIORB NJORB ', NIORB,NJORB
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
COLD. Loop over partitionings of the row strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
C?        write(6,*) ' NROW MAXI NIPART ', NROW,MAXI,NIPART
C         DO 701 IPART = 1, NIPART
C           IBOT = (IPART-1)*MAXI+1
C           ITOP = MIN(IBOT+MAXI-1,NROW)
C           NIBTC = ITOP - IBOT + 1
*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = KTOP + MAXK
*. Single excitation information independent of I strings
*
*.set up I1(K) =  XI1S(K) a JORB !J STRING >
              J12 = 1
              K12 = 1
              IF(JFRST.EQ.1) KFRST=1
              SCLFAC = 1.0D0
              CALL ADST_GAS(IBJORB,NJORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &             I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,J12,K12,
     &             KKJACT,SCLFAC)
              JFRST = 0
              KFRST = 0
*.set up I2(K) =  XI1S(K) a JORB !J STRING >
              I12 = 2
              IF(IFRST.EQ.1) KFRST = 1
              SCLFAC = 1.0D0
              CALL ADST_GAS(IBIORB,NIORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &             I2,XI2S,MAXK,NKBTC,KEND,IFRST,KFRST,I12,K12,
     &             KKIACT,SCLFAC)
              IFRST = 0
              KFRST = 0
*. Appropriate place to start partitioning over I strings
*. Loop over partitionings of the row strings
          DO 701 IIPART = 1, NIPART
            IBOT = (IIPART-1)*MAXI+1
            ITOP = MIN(IBOT+MAXI-1,NROW)
            NIBTC = ITOP - IBOT + 1

* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
*.Gather  C Block
              DO JJORB = 1,NJORB
                ICGOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I1(1,JJORB),XI1S(1,JJORB))
              END DO    
*
* Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
              DO IIORB = 1,NIORB
*.Gather S Block
                ISGOFF = 1 + (IIORB-1)*NKBTC*NIBTC
                CALL MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I2(1,IIORB),XI2S(1,IIORB))
              END DO  
              NKI = NKBTC*NIBTC
*
              FACTORC = 0.0D0
              FACTORAB = 1.0D0
              CALL MATML7(RHO1S,SSCR,CSCR,NIORB,NJORB,NKI,NIORB,
     &                    NKI,NJORB,FACTORC,FACTORAB,1)
*. Scatter out to complete matrix 
              DO 610 JJORB = 1, NJORB
C               JORB = ITSOB(JBORB-1+JJORB)
                JORB = IBJORB-1+JJORB
                DO 605 IIORB = 1, NIORB
C                 IORB = ITSOB(IBORB-1+IIORB)
                  IORB = IBIORB-1+IIORB
                  RHO1((JORB-1)*NACOB+IORB) =
     &            RHO1((JORB-1)*NACOB+IORB) +
     &            RHO1S((JJORB-1)*NIORB+IIORB)
 605            CONTINUE
 610         CONTINUE

  701     CONTINUE
*. /\ end of this I partitioning  
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
C 701     CONTINUE
*. End of loop over I partitioninigs
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
C!    stop ' enforrced stop in RSBBD1 '
      RETURN
      END
      SUBROUTINE GSBBD1VO(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen 
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces 
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1 : Updated density block
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
* RHO1S : Space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), ITSOB(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),RHO1S(*)
*.Local arrays
      DIMENSION ITP(16*16),JTP(16*16)
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GSBBD1 in action '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
      IFRST = 1
      JFRST = 1
*
* Type of single excitations that connects the two column strings
      CALL SXTYP_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(NTEST.GE.1000)    
     &WRITE(6,*) ' ISCSM,ICCSM IJSM ', ISCSM,ICCSM,IJSM
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
        DO 800 ISM = 1, NSMOB
*. new i and j so new intermediate strings
          KFRST = 1
*
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.1000) write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
          IBIORB = IOBPTS(ITYP,ISM)
          IBJORB = IOBPTS(JTYP,JSM)
          IF(NTEST.GE.2000)
     &    WRITE(6,*) ' NIORB NJORB ', NIORB,NJORB
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
*. Loop over partitionings of the row strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
C?        write(6,*) ' NROW MAXI NIPART ', NROW,MAXI,NIPART
          DO 701 IIPART = 1, NIPART
            IBOT = (IIPART-1)*MAXI+1
            ITOP = MIN(IBOT+MAXI-1,NROW)
            NIBTC = ITOP - IBOT + 1
*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = KTOP + MAXK
*
* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
*
              DO 500 JJORB = 1,NJORB
                JORB = IBJORB -1 + JJORB    
*.set up I1(K) =  XI1S(K) a JORB !J STRING >
                J12 = 1
                K12 = 1
                IF(JFRST.EQ.1) KFRST=1
                SCLFAC = 1.0D0
                CALL ADST_GAS(JORB,1,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &               I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,J12,K12,
     &               KKJACT,SCLFAC)
                JFRST = 0
                KFRST = 0
*.Gather  C Block
                ICGOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I1,XI1S)
  500         CONTINUE
*
* Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
*
              DO 501 IIORB = 1,NIORB
                IORB = IBIORB-1+IIORB       
*.set up I2(K) =  XI1S(K) a JORB !J STRING >
                I12 = 2
                IF(IFRST.EQ.1) KFRST = 1
                SCLFAC = 1.0D0
                CALL ADST_GAS(IORB,1,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &               I1,XI1S,MAXK,NKBTC,KEND,IFRST,KFRST,I12,K12,
     &               KKIACT,SCLFAC)
                IFRST = 0
                KFRST = 0
*.Gather S Block
                ISGOFF = 1 + (IIORB-1)*NKBTC*NIBTC
                CALL MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I1,XI1S)
  501         CONTINUE
              NKI = NKBTC*NIBTC
              CALL MATML5(RHO1S,SSCR,CSCR,NIORB,NJORB,NKI,NIORB,
     &                    NKI,NJORB,1)
*. Scatter out to complete matrix 
              DO 610 JJORB = 1, NJORB
C               JORB = ITSOB(JBORB-1+JJORB)
                JORB = IBJORB-1+JJORB
                DO 605 IIORB = 1, NIORB
C                 IORB = ITSOB(IBORB-1+IIORB)
                  IORB = IBIORB-1+IIORB
                  RHO1((JORB-1)*NACOB+IORB) =
     &            RHO1((JORB-1)*NACOB+IORB) +
     &            RHO1S((JJORB-1)*NIORB+IIORB)
 605            CONTINUE
 610         CONTINUE
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
  701     CONTINUE
*. End of loop over I partitioninigs
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
C!    stop ' enforrced stop in RSBBD1 '
      RETURN
      END
      SUBROUTINE EXPCIV(ISM,ISPCIN,LUIN,ISPCUT,LUUT,LBLK,
     &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
*
* Expand CI vector in CI space ISPCIN to CI vector in ISPCUT
* Input vector is supposed to be on LUIN
* Output vector will be placed on unit LUUT
*. If ICOPY .ne. 0 the output vectors will be copied to LUIN
*
* Storage form is defined by ICISTR 
*
* Jeppe Olsen, February 1994
* GAS version August 1995
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      

*
      CALL LUCIAQENTER('EXPCV')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPCIV')
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.1) THEN
        WRITE(6,'(A,2(2X,I4))')
     &  ' CI-vector is expanded: in- and out-space = ', ISPCIN,ISPCUT
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from EXPCIV:'
        WRITE(6,*) ' LUIN, LUUT, LUSCR= ', LUIN,LUUT,LUSCR
        WRITE(6,*) ' ICOPY, ICISTR = ', ICOPY, ICISTR
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
*. Allowed combinations of strings types for input and output
*. spaces
*
      CALL IAIBCM(ISPCIN,int_mb(KCIOIO))
      CALL IAIBCM(ISPCUT,int_mb(KSIOIO))
*
* type of each symmetry block ( full, lower diagonal, absent )
*
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,int_mb(KCBLTP),IDUMMY)
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,int_mb(KSBLTP),IDUMMY)

*
*. Number of dets etc per TTS block
*
*
*. Allocate memory for blocks of input and output space in 
*  accordance with ICISTR
*
      NCOMBI = 0
      NCOMBU = 0
      IF(ICISTR.EQ.1) THEN
*. Real to integer should be okay
        NCOMBI = XISPSM(ISM,ISPCIN)
        NCOMBU = XISPSM(ISM,ISPCUT)
C?      WRITE(6,*) ' EXPCIV: NCOMBI, NCOMBU = ',  NCOMBI, NCOMBU
        LENGTHI = NCOMBI
        LENGTHU = NCOMBU
      ELSE 
        LENGTHI = MXSOOB
        LENGTHU = MXSOOB
      END IF
*
C?    WRITE(6,*) ' ICISTR,  MXSOOB = ',ICISTR,  MXSOOB
      CALL MEMMAN(KLBLI,LENGTHI,'ADDL  ',2,'KLBLI ') !done
      CALL MEMMAN(KLBLU,LENGTHU,'ADDL  ',2,'KLBLU ') !done
*
*. and now : Let another subroutine complete the taks
*
      CALL REWINO(LUIN)
      CALL REWINO(LUUT)
*
*. Print for testing initial vectors out
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial vectors in EXPCIV '
*
        DO IROOT = 1, NROOT
          WRITE(6,*) ' Root number ', IROOT 
          IF(ICISTR.EQ.1) THEN
            CALL FRMDSC(dbl_mb(KLBLI),NCOMBI,-1,LUIN,IMZERO,IAMPACK)
            CALL WRTMAT(dbl_mb(KLBLI),1,NCOMBI,1,NCOMBI)
          ELSE
            CALL WRTVCD(dbl_mb(KLBLI),LUIN,0,-1)
          END IF
        END DO
        CALL REWINO(LUIN)
      END IF
*
      DO IROOT = 1, NROOT
*. Input vector should be first vector on file so
        IF(IROOT.EQ.1) THEN
          LLUIN = LUIN
        ELSE
*. With the elegance of an elephant: Place IROOT as first vector on LUSCR
          CALL REWINO(LUIN)
          DO JROOT = 1, IROOT
            CALL REWINO(LUSCR)
            IF(ICISTR.EQ.1) THEN
              CALL FRMDSC(dbl_mb(KLBLI),NCOMBI,-1,LUIN,IMZERO,IAMPACK)
              CALL  TODSC(dbl_mb(KLBLI),NCOMBI,-1,LUSCR)
            ELSE
              CALL COPVCD(LUIN,LUSCR,dbl_mb(KLBLI),0,-1)
            END IF
          END DO
          CALL REWINO(LUSCR)
          LLUIN = LUSCR
        END IF! IROOT = 1
*. Expcivs may need the IAMPACK parameter ( in case it must write
*  a zero block before any blocks have been read in.
*  Use IDIAG to decide
       IF(IDIAG.EQ.1) THEN
         IAMPACK = 0
       ELSE
         IAMPACK = 1
       END IF
C?     WRITE(6,*) ' IAMPACK in EXPCIV ', IAMPACK
C?     WRITE(6,*) ' ISMOST before EXPCIVS for ISM = ', ISM
C?     CALL IWRTMA(ISMOST(1,ISM),1,NSMST,1,NSMST)

        CALL EXPCIVS(LLUIN,dbl_mb(KLBLI),NCOMBI,
     &       int_mb(KCIOIO),NOCTPA,NOCTPB,int_mb(KCBLTP),
     &       LUUT,dbl_mb(KLBLU),NCOMBU,
     &       int_mb(KSIOIO),
     &       int_mb(KSBLTP),
     &       ICISTR,IDC,NSMST,
     &       LBLK,IAMPACK,ISMOST(1,ISM),
     &       int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)))
C     EXPCIVS(LUIN,VECIN,NCOMBIN,IABIN,
C    &                   NOCTPA,NOCTPB,IBLTPIN,
C    &                   LUUT,VECUT,NCOMBUT,IABUT,
C    &                   IBLTPUT,
C    &                   ICISTR,IDC,NSMST,LBLK,IAMPACKED_IN,
C    &                   ISMOST,NSSOA,NSSOB)
*
      END DO
*
      IF(ICOPY.NE.0) THEN
*. Copy expanded vectors to LUIN
        CALL REWINO(LUIN)
        CALL REWINO(LUUT)
        DO IROOT = 1, NROOT
          IF(ICISTR.EQ.1) THEN
            CALL FRMDSC(dbl_mb(KLBLU),NCOMBU,-1,LUUT,IMZERO,IAMPACK)
            CALL  TODSC(dbl_mb(KLBLU),NCOMBU,-1,LUIN)
          ELSE
            CALL COPVCD(LUUT,LUIN,dbl_mb(KLBLU),0,-1)
          END IF
        END DO
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output  vectors in EXPCIV '
        WRITE(6,*) ' =========================='
*
        CALL REWINO(LUUT)
        CALL REWINO(LUIN)
        DO IROOT = 1, NROOT
        WRITE(6,*) ' Root number ', IROOT 
*
        IF(ICISTR.EQ.1) THEN
          IF(ICOPY.EQ.0) THEN
            CALL FRMDSC(dbl_mb(KLBLU),NCOMBU,-1,LUUT,IMZERO,IAMPACK)
            CALL WRTMAT(dbl_mb(KLBLU),1,NCOMBU,1,NCOMBUT)
          ELSE
            CALL FRMDSC(dbl_mb(KLBLU),NCOMBU,-1,LUIN,IMZERO,IAMPACK)
            CALL WRTMAT(dbl_mb(KLBLU),1,NCOMBU,1,NCOMBUT)
          END IF
        ELSE
          CALL WRTVCD(dbl_mb(KLBLU),LUUT,0,-1)
          END IF
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXPCIV')
      CALL LUCIAQEXIT('EXPCV')
*
      RETURN
      END 
      SUBROUTINE EXPCIVS(LUIN,VECIN,NCOMBIN,IABIN,
     &                   NOCTPA,NOCTPB,IBLTPIN,
     &                   LUUT,VECUT,NCOMBUT,IABUT,
     &                   IBLTPUT,
     &                   ICISTR,IDC,NSMST,LBLK,IAMPACKED_IN,
     &                   ISMOST,NSSOA,NSSOB)
*
* Obtain those part of vector in cispace UT , 
* that can be obtained from terms in cispace IN
*
* Input vector on LUIN, Output vector in LUUT
* Output vector is supposed on start of vector
*
* LUIN is assumed to be single vector file,
* so rewinding will place vector on start of vector
*
* Both files are assumed on start of vector
*
* Jeppe Olsen, February 1994
*              March 2012: Cleaned a bit
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER IABIN(NOCTPA,NOCTPB),IABUT(NOCTPA,NOCTPB)
      INTEGER IBLTPIN(NSMST),IBLTPUT(NSMST)
      DIMENSION VECIN(*)
*, Symmetry of other string, given total symmetry
      INTEGER ISMOST(NSMST)
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
*. Output
      DIMENSION VECUT(*)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' EXPCIVS in action'
        WRITE(6,*) ' ================='
        IF(ICISTR.EQ.1) THEN
          WRITE(6,*) '      Number of input  parameters  ',NCOMBIN
          WRITE(6,*) '      Number of output parameters  ',NCOMBUT
          WRITE(6,*) ' ICISTR = ', ICISTR
        END IF
      END IF
*
      IF(ICISTR.EQ.1) THEN
        CALL FRMDSC(VECIN,NCOMBIN,-1,LUIN,IMZERO,IAMPACK)
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ISMOST = '
        CALL IWRTMA(ISMOST,1,NSMST,1,NSMST)
        WRITE(6,*) ' Input vector '
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
        ELSE
          CALL WRTVCD(VECIN,LUIN,1,-1)
        END IF
      END IF
*
      ZERO = 0.0D0
*
*. Loop over TTS blocks of output vector
*
      IATPIN = 1
      IBTPIN = 1
      IASMIN = 0
*
      IATPUT = 1
      IBTPUT = 1
      IASMUT = 0
      IOFFUT = 1
      IOFFIN = 1
      CALL REWINO(LUIN)
* 
 1000 CONTINUE
*. Next output block 
        CALL NXTBLK(IATPUT,IBTPUT,IASMUT,NOCTPA,NOCTPB,NSMST,
     &              IBLTPUT,IDC,NONEWUT,IABUT,ISMOST,
     &              NSSOA,NSSOB,LBLOCK,LBLOCKP)
        IF(NONEWUT.EQ.0) THEN
         IF(NTEST.GE.1000) THEN
           WRITE(6,'(A,4(2X,I5))') ' Next output block, TTS = ',
     &     IATPUT, IBTPUT, IASMUT
         END IF
     
*. Corresponding input TTS block
         JATPIN = IATPUT
         JBTPIN = IBTPUT
         IF(IABIN(JATPIN,JBTPIN).EQ.0) THEN
           IZERO = 1
         ELSE 
           IZERO = 0
         END IF
         NELMNT = LBLOCKP
*
         IF(IZERO.NE.0) THEN
*. Zero-block
           IF(ICISTR.EQ.1) THEN
             IF(NTEST.GE.1000) WRITE(6,*) ' Input block is zero '
             CALL SETVEC(VECUT(IOFFUT),ZERO,NELMNT)
             IOFFUT = IOFFUT + NELMNT
           ELSE
             CALL ITODS(NELMNT,1,-1,LUUT)
             CALL ZERORC(-1,LUUT,IAMPACKED_IN)
           END IF
         ELSE ! IZERO switch
*. Obtain this block in input file
  999      CONTINUE
           CALL NXTBLK(IATPIN,IBTPIN,IASMIN,NOCTPA,NOCTPB,NSMST,
     &                 IBLTPIN,IDC,NONEWIN,IABIN,ISMOST,
     &                 NSSOA,NSSOB,LBLOCK,LBLOCKPIN)
           IF(NTEST.GE.1000) THEN
            WRITE(6,'(A,4(2X,I5))') ' Next input block, TTS = ',
     &      IATPIN, IBTPIN, IASMIN
*
           END IF
           IF(NONEWIN.NE.0) THEN
             CALL REWINO(LUIN)
             IATPIN = 1
             IBTPIN = 1
             IASMIN = 0
             IOFFIN = 1
             GOTO 999
           END IF
           IF(ICISTR.GT.1) THEN
             CALL IFRMDS(LENGTH,1,-1,LUIN)
             CALL FRMDSC(VECIN,LENGTH,-1,LUIN,IMZERO,IAMPACK)
           END IF
*. There was a input blocks, right one?
           IF(IATPIN.EQ.JATPIN.AND.IBTPIN.EQ.JBTPIN.AND.
     &        IASMIN.EQ.IASMUT) THEN
              IF(NTEST.GE.1000)
     &        WRITE(6,*) ' Match between input and output block '
*. Correct block, save it
             IF(ICISTR.EQ.1) THEN
               IF(NTEST.GE.1000) WRITE(6,*) ' IOFFIN, IOFFUT = ',
     &         IOFFIN, IOFFUT
               CALL COPVEC(VECIN(IOFFIN),VECUT(IOFFUT),NELMNT)
               IOFFIN = IOFFIN + NELMNT
               IOFFUT = IOFFUT + NELMNT
             ELSE
               CALL ITODS(LENGTH,1,-1,LUUT)
               IF(IMZERO.EQ.1) THEN
                 CALL ZERORC(-1,LUUT,IAMPACKED_IN)
               ELSE
                 IF(IAMPACK.EQ.0) THEN
                   CALL TODSC(VECIN,LENGTH,-1,LUUT)
                 ELSE
                   CALL TODSCP(VECIN,LENGTH,-1,LUUT)
                 END IF
               END IF! IMZERO switch
             END IF! ICISTR switch
           ELSE
*. This was not the block we were after so
             IF(ICISTR.EQ.1) IOFFIN = IOFFIN + LBLOCKPIN
             GOTO 999
           END IF! Input = Output
         END IF !IZERO switch
      GOTO 1000
        END IF !Nonewut = 0
*. End of loop over output blocks
      IF(ICISTR.EQ.1) THEN
        CALL TODSC(VECUT,NCOMBUT,-1,LUUT)
      ELSE 
        CALL ITODS(-1,1,-1,LUUT)
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' EXPCIVS Speaking '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' ============ '
        WRITE(6,*) ' Input Vector '
        WRITE(6,*) ' ============ '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
        ELSE
          CALL WRTVCD(VECIN,LUIN,1,LBLK)
        END IF
        WRITE(6,*)
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Output Vector '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECUT,1,NCOMBUT,1,NCOMBUT)
        ELSE
          CALL WRTVCD(VECUT,LUUT,1,LBLK)
        END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
     &                          NNSTSGP,IISTSGP)
*
* Obtain all superstrings containing  strings of given sym and type 
*
* ( Superstring :contains electrons belonging to all gasspaces  
*        string :contains electrons belonging to a given GAS space
* A super string is thus a product of NGAS strings )
*
* Jeppe Olsen, Summer of 95
*              Optimized version, october 1995
*
*. In this subroutine the ordering of strings belonging to a given type 
*  is defined !!
* Currently we are using the order 
* Loop over GAS 1 strings 
*  Loop over GAS 2 strings 
*   Loop over GAS 3 strings --
*
*     Loop over gas N strings 
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'csm.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input 
      INTEGER ITPFGS(*), ISMFGS(*)
      INTEGER NNSTSGP(MXPNSMST,*), IISTSGP(MXPNSMST,*)
*. Local scratch 
C     INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS), ISTRNM(MXPNGAS)
      INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS)
*. Output 
      INTEGER ISTROC(NEL,*)
*. Number of strings per GAS space
C?    write(6,*) ' entering problem child '
      DO IGAS = 1, NGAS
        NSTFGS(IGAS)  = NNSTSGP(ISMFGS(IGAS),IGAS)
        IBSTFGS(IGAS) = IISTSGP(ISMFGS(IGAS),IGAS)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GETSTR_GASSM_SPGP speaking '
        WRITE(6,*) '  =========================== '
        WRITE(6,*) ' ISMFGS,ITPFGS (input) '
        CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' NSTFGS, IBSTFGS ( intermediate results ) '
        CALL IWRTMA(NSTFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(IBSTFGS,1,NGAS,1,NGAS)
      END IF
*. Last gasspace with a nonvanishing number of electrons
      IGASL = 0
      DO IGAS = 1, NGAS
        IF( NELFGP(ITPFGS(IGAS)) .NE. 0 ) IGASL = IGAS
      END DO
C     WRITE(6,*) ' IGASL = ', IGASL
*
      NSTRTOT = 1
      DO IGAS = 1, NGAS
        NSTRTOT = NSTRTOT*NSTFGS(IGAS)
      END DO
C     WRITE(6,*) ' NSTRTOT = ', NSTRTOT
      IF(IGASL.EQ.0) GOTO 2810
*
      NELL = NELFGP(ITPFGS(IGASL))
      NELML = NEL - NELL
      NSTRGASL = NSTFGS(IGASL)
      IBGASL = IBSTFGS(IGASL)
*
      IF(NSTRTOT.EQ.0) GOTO 1001
*. Loop over GAS spaces
      DO IGAS = 1, IGASL
*. Number of electrons in GAS = 1, IGAS - 1
        IF(IGAS.EQ.1) THEN
          NELB = 0
        ELSE
          NELB = NELB +  NELFGP(ITPFGS(IGAS-1))
        END IF
*. Number of electron in IGAS 
        NELI = NELFGP(ITPFGS(IGAS))
C?      WRITE(6,*) ' NELI and NELB ', NELI,NELB
        IF(NELI.GT.0) THEN
         
*. The order of strings corresponds to a matrix A(I(after),Igas,I(before))
*. where I(after) loops over strings in IGAS+1 - IGASL and
*  I(before) loop over strings in 1 - IGAS -1
          NSTA = 1
          DO JGAS = IGAS+1, IGASL
            NSTA = NSTA * NSTFGS(JGAS)
          END DO
*
          NSTB =  1
          DO JGAS = 1, IGAS-1
            NSTB = NSTB * NSTFGS(JGAS)
          END DO
*
          NSTI = NSTFGS(IGAS)
         
C?        write(6,*) ' before call to add_str_group '
          IF(NTEST.GE.200) THEN
            WRITE(6,*) ' NSTI,NSTB,NSTA,NELB,NELI,NEL ',
     &                   NSTI,NSTB,NSTA,NELB,NELI,NEL
            WRITE(6,*) ' IBSTFGS(IGAS),KOC()',
     &                   IBSTFGS(IGAS),KOCSTR(ITPFGS(IGAS))
          END IF
*
          CALL ADD_STR_GROUP(NSTI,
     &          IBSTFGS(IGAS),
     &          int_mb(KOCSTR(ITPFGS(IGAS))),
     &          NSTB,NSTA,ISTROC,NELB+1,NELI,NEL)
*. Loop over strings in IGAS 
        END IF
      END DO
 1001 CONTINUE
 2810 CONTINUE
      NSTR = NSTRTOT  
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from  GETSTR_GASSM_SPGP ' 
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and type strings : '
        WRITE(6,*)
        WRITE(6,*) '   AS    Sym  Type '
        WRITE(6,*) ' =================='
        DO IGAS = 1, NGAS
          WRITE(6,'(3I6)') IGAS,ISMFGS(IGAS),ITPFGS(IGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Number of strings generated : ', NSTR
        WRITE(6,*) ' Strings generated '
        CALL PRTSTR(ISTROC,NEL,NSTR)
      END IF
*
      RETURN
      END 
      SUBROUTINE REORD_GENMAT(NDIM,LDIM,IFIRST,CIN,COUT)
*
* A matrix CIN is given as a multidimensional array
* CIN(I1,I2,I3.... INDIM)
* with the usual notation that the first index is the inner index.
*
* Reorganize so index IFIRST become the first index.
*
* Jeppe Olsen, July 95
* 
      RETURN
      END 
      SUBROUTINE SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
*
* Two supergroups are given. Find single excitations that connects 
* these supergroups 
*
* Jeppe Olsen, July 1995
*
* Dec 97 : IPHGAS added :  
*          Occupations of particle spaces (IPHGAS=2) are allowed to
*          have occupations less than zero in intermidiate steps
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ILTP(NGAS),IRTP(NGAS),IPHGAS(*)
*. Output
      DIMENSION ITP(*),JTP(*)
*. Differences between occupations :
      NCREA = 0
      NANNI = 0
      ICREA = -2303
      IANNI = -2303
      DO IAS = 1, NGAS
        IF(ILTP(IAS).GT.IRTP(IAS)) THEN
         NCREA = NCREA + ILTP(IAS) - IRTP(IAS)
         ICREA = IAS
        ELSE IF(IRTP(IAS).GT.ILTP(IAS)) THEN
         NANNI = NANNI + IRTP(IAS) - ILTP(IAS)
         IANNI = IAS 
        END IF
      END DO
*
      IF(NCREA.GT.1) THEN
*. Sorry : No single excitation connects
        NSXTYP = 0
      ELSE IF(NCREA .EQ. 1 ) THEN
*. Supergroups differ by one single excitation.
        NSXTYP = 1
        ITP(1) = ICREA
        JTP(1) = IANNI
      ELSE IF (NCREA.EQ.0 ) THEN
*. Supergroups are identical, connects with all
*  diagonal excitations.
        NSXTYP = 0
        DO IAS = 1, NGAS
          IF(IRTP(IAS).NE.0.OR.IPHGAS(IAS).EQ.2) THEN
            NSXTYP = NSXTYP + 1
            ITP(NSXTYP) = IAS 
            JTP(NSXTYP) = IAS
          END IF
        END DO
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SXTYP_GAS : '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input left  supergroup '
        CALL IWRTMA(ILTP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Input right supergroup '
        CALL IWRTMA(IRTP,1,NGAS,1,NGAS)
        WRITE(6,*)
     &  ' Number of connecting single excitations ', NSXTYP
        IF(NSXTYP.NE.0) THEN
          WRITE(6,*) ' Connecting single excitations '
          DO ISX = 1, NSXTYP
            WRITE(6,*) ITP(ISX),JTP(ISX)
          END DO
        END IF
      END IF
*
      RETURN
      END 


*
      SUBROUTINE SXTYP_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP)
*
* Two supergroups are given. Find single excitations that connects 
* these supergroups 
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ILTP(NGAS),IRTP(NGAS)
*. Output
      DIMENSION ITP(*),JTP(*)
*. Differences between occupations :
      NCREA = 0
      NANNI = 0
      ICREA = -2810
      IANNI = -2810
      DO IAS = 1, NGAS
        IF(ILTP(IAS).GT.IRTP(IAS)) THEN
         NCREA = NCREA + ILTP(IAS) - IRTP(IAS)
         ICREA = IAS
        ELSE IF(IRTP(IAS).GT.ILTP(IAS)) THEN
         NANNI = NANNI + IRTP(IAS) - ILTP(IAS)
         IANNI = IAS 
        END IF
      END DO
*
      IF(NCREA.GT.1) THEN
*. Sorry : No single excitation connects
        NSXTYP = 0
      ELSE IF(NCREA .EQ. 1 ) THEN
*. Supergroups differ by one sngle excitation.
        NSXTYP = 1
        ITP(1) = ICREA
        JTP(1) = IANNI
      ELSE IF (NCREA.EQ.0 ) THEN
*. Supergroups are identical, connects with all
*  diagonal excitations.
        NSXTYP = 0
        DO IAS = 1, NGAS
          IF(IRTP(IAS).NE.0) THEN
            NSXTYP = NSXTYP + 1
            ITP(NSXTYP) = IAS 
            JTP(NSXTYP) = IAS
          END IF
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SXTYP_GAS : '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input left  supergroup '
        CALL IWRTMA(ILTP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Input right supergroup '
        CALL IWRTMA(IRTP,1,NGAS,1,NGAS)
        WRITE(6,*)
     &  ' Number of connecting single excitations ', NSXTYP
        IF(NSXTYP.NE.0) THEN
          WRITE(6,*) ' Connecting single excitations '
          DO ISX = 1, NSXTYP
            WRITE(6,*) ITP(ISX),JTP(ISX)
          END DO
        END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE GETSTR_TOTSM_SPGP(ISTRTP,ISPGRP,ISPGRPSM,NEL,NSTR,ISTR,
     &                             NORBT,IDOREO,IZ,IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space 
*
*.If  IDOREO .NE. 0 THEN reordering array : lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete 
*                           orbital space, product of strings in
*                           each GAS space 
* =====
* Input 
* =====
*
* ISTRTP  : Type of of superstrings ( alpha => 1, beta => 2 )
* ISPGRP :  supergroup number, (relative to start of this type )
* ISPGRPSM : Total symmetry of superstrings 
* NEL : Number of electrons 
* IZ  : Reverse lexical ordering array for this supergroup
* 
*
* ======
* Output 
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0) 
*
*
* Jeppe Olsen, July 1995
* Last modification; May 3, 2013; Jeppe Olsen; Changed order of symmetry-
*                    distributions for AB-CONF reorder..
*
c      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INTEGER IZ(NORBT,NEL)
*. output
      INTEGER ISTR(NEL,*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
      NTEST = 000
      CALL LUCIAQENTER('GETST')
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A,3I3)')
     & ' Strings to be obtained : Type, supergroup, symmetry ',
     &   ISTRTP,ISPGRP,ISPGRPSM
        WRITE(6,*)
      END IF
*. The strings are generated as a double loop
*.    Loop over symmetrydistributions of strings with correct
*     total sym : symmetry of each gas
*.     Loop over strings of this symmetry distribution
*      End of loop over strings
*     End of loop over symmetrydistribution
* The symmetry distributions are ordered as
* 1: Loop over  gas-space 2- last active (occupied gas-space=
* 2: Determine the symmetry of gas-space 1, so total sym is correct 
*. Absolut number of this supergroup
      ISPGRPA = IBSPGPFTP(ISTRTP) - 1 + ISPGRP
* The approach with NGASL has been eliminated, but code is kept
*. Occupation per gasspace
*. Largest occupied space 
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPA)
        NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))          
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
      IF(NTEST.GE.200) WRITE(6,*) ' NGASL = ', NGASL
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        CALL ICOPVE2(int_mb(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(int_mb(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
* Largest and lowest active symmetries for each GAS space
      IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF 
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL 
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1 
         CALL NXTNUM3(ISMFGS(2),NGASL-1,MINVAL(2),MAXVAL(2),NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
        IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) ' next symmetry of spaces 2-NGASL '
          CALL IWRTMA(ISMFGS(2),NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of GASpaces 2- NGASL 
        ISTSMM1 = 1
        DO IGAS = 2, NGASL 
C         SYMCOM(ITASK,IOBJ,I1,I2,I12)
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. Required symmetry of GASpace 1
        CALL SYMCOM(2,1,ISTSMM1,ISMGS1,ISPGRPSM)
        ISMFGS(1) = ISMGS1
*. A test that ISFGS(1) is within bounds could be inserted here
*. Correct symmetry, so proceed
         DO IGAS = NGASL+1,NGAS
           ISMFGS(IGAS) = 1
         END DO
         IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
           WRITE(6,*) ' Next symmetry distribution '
           CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
         END IF
*. Obtain all strings of this symmetry 
CT       CALL LUCIAQENTER('GASSM')
         CALL GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTR(1,ISTRBS),NSTR,NEL,
     &                          NNSTSGP,IISTSGP)
CT       CALL LUCIAQEXIT('GASSM')
*. Reorder Info : Lexical => actual number 
         IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces 
          DO JSTR = ISTRBS, ISTRBS+NSTR-1
            LEX = 1
            DO IEL = 1, NEL 
              LEX = LEX + IZ(ISTR(IEL,JSTR),IEL)
            END DO
*
            IF(NTEST.GE.200.and.ga_nodeid().eq.0) THEN
              WRITE(6,*) ' string '
              CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
              WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
            END IF

*
            MAXLEX = MAX(MAXLEX,LEX)
            IREO(LEX) = JSTR
          END DO
         END IF ! Order array was required
         ISTRBS = ISTRBS + NSTR 
*. ready for next symmetry distribution 
        IF(NGAS-1.NE.0) GOTO 1000
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings: '
        WRITE(6,*)
        CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,MAXLEX,1,MAXLEX)
        END IF
      END IF
*
      CALL LUCIAQEXIT('GETST')
      RETURN
      END 
      SUBROUTINE PRTSTR(ISTR,NEL,NSTR)
*
* Print NSTR strings each containing NEL electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISTR(NEL,NSTR)
*
      DO JSTR = 1, NSTR
        WRITE(6,'(1H0,A,I6,A,4X,10(2X,I4),/,(1H ,19X,10(2X,I4)))' )
     &   ' String ',JSTR,' : ',(ISTR(IEL,JSTR),IEL=1,NEL)
      END DO
*
      RETURN
      END 
      SUBROUTINE GETSTR_GASSM_SPGP(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
     &                          NNSTSGP,IISTSGP)
*
* Obtain all superstrings containing  strings of given sym and type 
*
* ( Superstring :contains electrons belonging to all gasspaces  
*        string :contains electrons belonging to a given GAS space
* A super string is thus a product of NGAS strings )
*
* Jeppe Olsen, Summer of 95
*
*. In this subroutine the ordering of strings belonging to a given type 
*  is defined !!
* Currently we are using the order 
* Loop over GAS 1 strings 
*  Loop over GAS 2 strings 
*   Loop over GAS 3 strings --
*
*     Loop over gas N strings 
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'csm.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input 
      INTEGER ITPFGS(*), ISMFGS(*)
      DIMENSION NNSTSGP(MXPNSMST,*), IISTSGP(MXPNSMST,*)
*. Local scratch 
      INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS), ISTRNM(MXPNGAS)
*. Output 
      INTEGER ISTROC(NEL,*)
*. Number of strings per GAS space
      DO IGAS = 1, NGAS
        NSTFGS(IGAS)  = NNSTSGP(ISMFGS(IGAS),IGAS)
        IBSTFGS(IGAS) = IISTSGP(ISMFGS(IGAS),IGAS)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GETSTR_GASSM_SPGP speaking '
        WRITE(6,*) '  =========================== '
        WRITE(6,*) ' ISMFGS,ITPFGS (input) '
        CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' NSTFGS, IBSTFGS ( intermediate results ) '
        CALL IWRTMA(NSTFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(IBSTFGS,1,NGAS,1,NGAS)
      END IF
*
      NSTRGASN = NSTFGS(NGAS)
      IBGASN = IBSTFGS(NGAS)
      NELN = NELFGP(ITPFGS(NGAS))
      NELMN = NEL - NELN
*. Last gasspace with a nonvanishing number of electrons
      IGASL = 0
      DO IGAS = 1, NGAS
        IF( NELFGP(ITPFGS(IGAS)) .NE. 0 ) IGASL = IGAS
      END DO
*
      NELL = NELFGP(ITPFGS(IGASL))
      NELML = NEL - NELL
      NSTRGASL = NSTFGS(IGASL)
      IBGASL = IBSTFGS(IGASL)
*
*. Loop over strings : spaces 1 - IGASL - 1 are treated by
*  call to 'next number generator', while the strings in 
*  GAS  IGASL are added directly 

      IFIRST = 1
      ISTR = 1
*
      NSTRTOT = 1
      DO IGAS = 1, NGAS
        NSTRTOT = NSTRTOT*NSTFGS(IGAS)
      END DO
      IF(NSTRTOT.EQ.0) GOTO 1001
 1000 CONTINUE
*
        IF(IFIRST.EQ.1) THEN
          DO IGAS  = 1, IGASL-1
            ISTRNM(IGAS) = 1
            NONEW = 0 
          END DO
        ELSE
*. Next number 
          MINVAL = 1
          CALL NXTLEXNUM(ISTRNM,IGASL-1,MINVAL,NSTFGS,NONEW)
          IF(NONEW.NE.0) GOTO 1001
        END IF
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' String in IGASL-1 spaces (ISTRNM) '
          CALL IWRTMA(ISTRNM,1,IGASL-1,1,IGASL-1)
        END IF
        IFIRST = 0
*.IGASL-1 string
        IBEL = 1
        DO IGAS = 1, IGASL-1
          NELI = NELFGP(ITPFGS(IGAS))
C         ICOPVE2(IIN,IOFF,NDIM,IOUT)
          CALL ICOPVE2(int_mb(KOCSTR(ITPFGS(IGAS))),
     &                 (ISTRNM(IGAS)+IBSTFGS(IGAS)-2)*NELI+1,
     &                 NELI,ISTROC(IBEL,ISTR))
C         DO JEL = 1, NELI
C         IFRMR(WORK,IROFF,IELMNT)
C           ISTROC(IBEL-1+JEL,ISTR) = 
C    &      IFRMR(WORK(KOCSTR(ITPFGS(IGAS))),1,
C    &            (ISTRNM(IGAS)+IBSTFGS(IGAS)-2)*NELI+JEL)
C         END DO
          IBEL = IBEL + NELI
        END DO
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Next string in IGASL -1 spaces (ISTROC) '
          CALL IWRTMA(ISTROC(1,ISTR),1,NELML,1,NELML)
        END IF
*. copy  to remaining strings with the same GAS1-GASL-1 string
        DO JSTR = 1, NSTRGASL-1
          DO JEL = 1, NELML 
            ISTROC(JEL,ISTR+JSTR) = ISTROC(JEL,ISTR)
          END DO
        END DO
*. Add IGASL strings 
        CALL COPSTR(int_mb(KOCSTR(ITPFGS(IGASL))),IBGASL,NELL,NEL,1,
     &                 NELML+1,NELL,NSTRGASL,ISTROC(1,ISTR))
C     COPSTR(INSTR,IBSTR,NELI,NELO,IBELI,IBELO,NELAD,
C    &                  NSTR,IOUSTR)
        ISTR = ISTR + NSTRGASL
        IF(IGASL-1.NE.0) GOTO 1000
 1001 CONTINUE
*. END Loop over new strings
      NSTR = ISTR - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from  GETSTR_GASSM_SPGP ' 
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and type strings : '
        WRITE(6,*)
        WRITE(6,*) '   AS    Sym  Type '
        WRITE(6,*) ' =================='
        DO IGAS = 1, NGAS
          WRITE(6,'(3I6)') IGAS,ISMFGS(IGAS),ITPFGS(IGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Number of strings generated : ', NSTR
        WRITE(6,*) ' Strings generated '
        CALL PRTSTR(ISTROC,NEL,NSTR)
      END IF
*
      RETURN
      END 
      SUBROUTINE COPSTR(INSTR,IBSTR,NELI,NELO,IBELI,IBELO,NELAD,
     &                  NSTR,IOUSTR)
*
* IOUSTR(IEL+IBELO-1,ISTR) = INSTR(IEL+IBELI-1,IBSTR-1+ISTR),
*                             IEL = 1, NELAD, ISTR = 1, NSTR
*
* Jeppe Olsen, July 95
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER INSTR(NELI,*)
*. Input and Output 
      INTEGER IOUSTR(NELO,*)
*
      DO ISTR = 1, NSTR
        DO IEL = 1, NELAD
          IOUSTR(IBELO+IEL-1,ISTR) = INSTR(IBELI-1+IEL,IBSTR+ISTR-1)
        END DO
      END DO
*
      RETURN
      END 
      SUBROUTINE NXTLEXNUM(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL,MAXVAL(I).
*
* In this version the usual numbering is used, so the 
* rightmost digits are first changed !
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen July 1995
*
*. Input
      DIMENSION MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 000
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IPLACE = NELMNT + 1
 1000 CONTINUE
        IPLACE = IPLACE - 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.GT.1     ) THEN
          DO JPLACE = IPLACE, NELMNT 
            INUM(JPLACE) = MINVAL
          END DO
        ELSE IF ( IPLACE. EQ. 1) THEN
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
      SUBROUTINE GASDIAS(NAEL,IASTR,NBEL,IBSTR,
     &           NORB,DIAG,NSMST,H,XA,XB,SCR,RJ,RK,
     &           NSSOA,NSSOB,LUDIA,ECORE,
     &           PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,I12,
     &           IBLTP,NBLOCK,IBLKFO)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* Driven by IBLKFO, May 97
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
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
C     REAL * 8  INPROD
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
*
      INCLUDE 'cprnt.inc'
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*), buffer(100)
*. Output
CNW   DIMENSION DIAG(*)
      integer diag
*
      integer nxtask, next, itask, icount, nproc
      external nxtask
* 
      NTEST = 1000 !VOG 000
      NTEST = MAX(NTEST,IPRNT)
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
*
 
      IF( NTEST .GE. 20 .OR.IPRINTEGRAL.GE.100) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NORB,1,NORB)
        WRITE(6,*) ' Core energy ', ECORE
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
        END IF
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK               
          WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' IBLTP: '
        CALL IWRTMA(IBLTP,1,NSMST,1,NSMST)
*
        WRITE(6,*) ' I12 = ',I12
      END IF
*
*  Diagonal elements according to Handys formulae
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
      IF(LUDIA.GT.0) CALL REWINO(LUDIA)
*
      itask = 0
      nproc = ga_nnodes()
      next = nxtask(nproc, 1)
      DO IBLK = 1, NBLOCK
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
        IOFF =  1                 
        DO IA = 1, NSSOA(IASM,IATP)
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
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
        IBSTRT = 1                
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
          IBREL = IB - IBSTRT + 1
*
*. Terms depending only on IB
*
          HB = 0.0D0
          RJBB = 0.0D0
          CALL SETVEC(XB,0.0D0,NORB)
*
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
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1                 
          END IF
          IASTOP = NSSOA(IASM,IATP) 
*
          if (itask.eq.next) then
          icount = 0
          CALL SETVEC(buffer,0.0d0,100)
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            icount = icount + 1
            IF(NTEST.GE.1000) WRITE(6,*) ' IA IB,IDET = ',
     &      IA,IB,IDET
            ITDET = ITDET + 1
            buffer(icount) = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NAEL
              buffer(icount) = buffer(icount) +XB(IASTR(IEL,IA)) 
            END DO
CNW         DIAG(IDET) = X
CNW         IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
CNW  &         IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
            IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &         IB.EQ.IA) buffer(icount) = buffer(icount) + XADD
            if (icount.eq.100.OR.IA.eq.IASTOP) then
               call ga_put(DIAG,IDET-icount+1,IDET,1,1,buffer,1)
               CALL SETVEC(buffer,0.0d0,100)
               icount = 0
            endif
          END DO
*         ^ End of loop over alpha strings|
          next = nxtask(nproc,1)
          else
          IDET = IDET + IASTOP - IASTRT + 1
          endif
          itask = itask + 1
*         ^ End task loop
        END DO
*       ^ End of loop over betastrings
*. Yet a RAS block of the diagonal has been constructed
        IF(ICISTR.GE.2) THEN
          IF(NTEST.GE.100) THEN
            write(6,*) ' number of diagonal elements to disc ',IDET
CNW         CALL WRTMAT(DIAG,1,IDET,1,IDET)
            call ga_print(DIAG)
          END IF
          CALL ITODS(IDET,1,-1,LUDIA)
          CALL TODSC(DIAG,IDET,-1,LUDIA)
          IDET = 0
        END IF
      END DO
*        ^ End of loop over blocks
      next = nxtask(-nproc, 1)
      call ga_sync()
 
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) WRITE(6,*)
     &' Number of diagonal elements generated ',ITDET
*
      IF(NTEST .GE.100 .AND.ICISTR.LE.1 ) THEN
        WRITE(6,*) ' CIDIAGONAL ',diag
CNW     CALL WRTMAT(DIAG(1),1,IDET,1,IDET)
        call ga_print(diag)
      END IF
*
      IF ( ICISTR.GE.2 ) CALL ITODS(-1,1,-1,LUDIA)
*
      RETURN
      END
      SUBROUTINE WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,
     &           ISCR,NTEST)        
*
* construct vertex weights for given supergroup 
*
* Reverse lexical ordering is used 
*. Inactive orbitals added, June 2010
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Input
      INTEGER NELFTP(NORBTP),NORBFTP(NORBTP)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Ouput
      INTEGER Z(*)
*. Scratch length : 2 * NORB + (NEL+1)*(NORB+1)
       INTEGER ISCR(*)
*
       NORB = IELSUM(NORBFTP,NORBTP)+NINOB
       NEL  = IELSUM(NELFTP,NORBTP)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Subroutine WEIGHT_SPGP in action '
        WRITE(6,*) ' ================================='
        WRITE(6,*) 'NELFTP '
        CALL IWRTMA(NELFTP,1,NORBTP,1,NORBTP)
      END IF
*
      KLFREE = 1
      KLMAX = KLFREE
      KLFREE = KLFREE + NORB
*
      KLMIN = KLFREE
      KLFREE = KLFREE + NORB
*
      KW = KLFREE
      KLFREE = KW + (NEL+1)*(NORB+1)
*.Max and min arrays for strings
      CALL MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),NORBTP,NORBFTP,NELFTP,
     &                 NINOB,NTEST)
*. Arc weights
      CALL GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL,NTEST)
*
      RETURN
      END
      SUBROUTINE MXMNOC_SPGP(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,
     &                  NINOB,NTESTG)
*
* Construct accumulated MAX and MIN arrays for a GAS supergroup
*
      IMPLICIT REAL*8           ( A-H,O-Z)
      INCLUDE 'mxpdim.inc'
*. Output
      DIMENSION  MINEL(*),MAXEL(*)
*. Input
      INTEGER NORBFTP(*),NELFTP(*)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
C???  WRITE(6,*) ' NTEST, NTESTG = ', NTEST, NTESTG
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ===========' 
        WRITE(6,*) ' MXMNOC_SPGP'
        WRITE(6,*) ' ===========' 
        WRITE(6,*)
        WRITE(6,'(A,I5)') '  NINOB = ', NINOB
        WRITE(6,*) ' NORBFTP : '    
        CALL IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
      END IF
*
      IZERO = 0
      CALL ISETVC(MINEL,IZERO,NINOB)
      CALL ISETVC(MAXEL,IZERO,NINOB)
      DO IORBTP = 1, NORBTP
*. Max and min at start of this type and at end of this type
        IF(IORBTP.EQ.1) THEN
          IORB_START = NINOB + 1
          IORB_END = NORBFTP(1) + NINOB
          IF(NTEST.GE.1000) WRITE(6,'(A,2I4)') 
     &    ' IORB_START, IORB_END = ',IORB_START, IORB_END
          NEL_START = 0
          NEL_END   = NELFTP(1)
        ELSE
          IORB_START =  IORB_START + NORBFTP(IORBTP-1)
          IORB_END   =  IORB_START + NORBFTP(IORBTP)-1
          NEL_START = NEL_END 
          NEL_END   = NEL_START + NELFTP(IORBTP)
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,'(A,5I4)') 
     &    ' IORBTP,IORB_START,IORB_END,NEL_START,NEL_END ',
     &      IORBTP,IORB_START,IORB_END,NEL_START,NEL_END 
        END IF
*
        DO IORB = IORB_START, IORB_END
          MAXEL(IORB) = MIN(IORB-NINOB,NEL_END)
          MINEL(IORB) = NEL_START
          IF(NEL_END-MINEL(IORB).GT. IORB_END-IORB) 
     &    MINEL(IORB) = NEL_END - ( IORB_END - IORB ) 
        END DO
      END DO
*
      IF( NTEST .GE. 100 ) THEN
        NORB = IELSUM(NORBFTP,NORBTP) + NINOB
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
      SUBROUTINE ADADST_GAS(IOB,IOBSM,IOBTP,NIOB,
     &                      JOB,JOBSM,JOBTP,NJOB,
     &                      ISPGP,ISM,ITP,KMIN,KMAX,
     &                      I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12,
     &                      SCLFAC)
*
*
* Obtain mappings
* a+IORB a+ JORB !KSTR> = +/-!ISTR>
* In the form
* I1(KSTR) =  ISTR if a+IORB a+ JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,IGRP.
* (numbering relative to TS start)
*. Only excitations IOB. GE. JOB are included 
* The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
* and IOB (JOB) is the first orbital to be used, and the number of orbitals
* to be checked is NIOB ( NJOB).
*
* Only orbital pairs IOB .gt. JOB are included 
*
* The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , August of 95   
*
* ======
*. Input
* ======
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*./ORBINP/
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Local scratch
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/SSAVE/NELIS(4), NSTRKS(4)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' IOB,IOBSM,IOBTP ', IOB,IOBSM,IOBTP
        WRITE(6,*) ' JOB,JOBSM,JOBTP ', JOB,JOBSM,JOBTP
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) 'Problemo, ADADST '
C?      WRITE(6,*) ' SCLFAC = ',SCLFAC
C?    END IF

*
*. Internal affairs
*
      IF(I12.LE.4.AND.K12.LE.2) THEN
        KLLOC = KLOCSTR(K12)
        KLLZ = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*) ' ADST_GAS : Illegal value of I12 = ', I12
        STOP' ADST_GAS : Illegal value of I12  '
      END IF

*
*. Supergroup and symmetry of K strings
*
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      CALL NEWTYP(ISPGPABS,1,IOBTP,1,K1SPGPABS)
*. Added June 2012
      IF(K1SPGPABS.EQ.0) THEN
        NK = 0
        RETURN
      END IF
*. End added June 12
      CALL NEWTYP(K1SPGPABS,1,JOBTP,1,KSPGPABS)
*. Added June 2012
      IF(KSPGPABS.EQ.0) THEN
        NK = 0
        RETURN
      END IF
*. End added
C?    IF(K1SPGPABS.GT.100.OR.KSPGPABS.GT.100) THEN
C?      WRITE(6,*) ' Strange number of supergroup: '
C?      WRITE(6,*) ' ISPGPABS, K1SPGPABS, KSPGPABS',
C?   &               ISPGPABS, K1SPGPABS, KSPGPABS
C?      WRITE(6,*) ' ITP, ISPGP', ITP, ISPGP
C?      WRITE(6,*) ' IOBTP,  JOBTP = ', IOBTP, JOBTP
C?    END IF
      CALL SYMCOM(2,0,IOBSM,K1SM,ISM)
      CALL SYMCOM(2,0,JOBSM,KSM,K1SM)
      IF(NTEST.GE.100) WRITE(6,*)
     & ' K1SM,K1SPGPABS,KSM,KSPGPABS : ',
     &   K1SM,K1SPGPABS,KSM,KSPGPABS
* In ADADS1_GAS we need : Occupation of KSTRINGS
*                         lexical => Actual order for I strings
* Generate if required
*
      IF(IFRST.NE.0) THEN
*.. Generate information about I strings
*. Arc weights for ISPGP
        NTEST2 = NTEST
        CALL WEIGHT_SPGP(int_mb(KLLZ),NGAS,
     &                  NELFSPGP(1,ISPGPABS),
     &                  NOBPT,int_mb(KLZSCR),NTEST2)
        NELI = NELFTP(ITP)
        NELIS(I12) = NELI
*. Reorder array for I strings
        CALL GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,
     &                         int_mb(KLLOC),NOCOB,
     &                         1,int_mb(KLLZ),int_mb(KLLREO))
      END IF
      NELK = NELIS(I12) - 2
      IF(KFRST.NE.0) THEN
*. Generate occupation of K STRINGS
       CALL GETSTR_TOTSM_SPGP(1,KSPGPABS,KSM,NELK,NSTRK,
     &                        int_mb(KLLOC),NOCOB,
     &                        0,IDUM,IDUM)
       NSTRKS(K12) = NSTRK
      END IF
*
      NSTRK = NSTRKS(K12)
*
      IIOB = IOBPTS(IOBTP,IOBSM) + IOB - 1
      JJOB = IOBPTS(JOBTP,JOBSM) + JOB - 1
      CALL ADADS1_GAS(NK,I1,XI1S,LI1,IIOB,NIOB,JJOB,NJOB,
     &          int_mb(KLLOC),NELK,NSTRK,int_mb(KLLREO),int_mb(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,SCLFAC)
*
      RETURN
      END
      SUBROUTINE ADADS1_GAS(NK,I1,XI1S,LI1,IORB,NIORB,JORB,NJORB,
     &                KSTR,NKEL,NKSTR,IREO,IZ,
     &                NOCOB,KMAX,KMIN,IEND,SCLFAC)
*
* Obtain I1(KSTR) = +/- A+ IORB A+ JORB !KSTR>
* Only orbital pairs IOB .gt. JOB are included 
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : First I orbital to be added 
* NIORB : Number of orbitals to be added : IORB to IORB-1+NIORB
*        are used. They must all be in the same TS group
* JORB : First J orbital to be added 
* LORB : Number of orbitals to be added : JORB to JORB-1+NJORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 =>  a+IORB a+JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      NTEST = 000
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' ADADS1_GAS speaking '
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' IORB,NIORB ', IORB,NIORB       
       WRITE(6,*) ' JORB,NJORB ', JORB,NJORB       
      END IF
*
      IORBMIN = IORB
      IORBMAX = IORB + NIORB - 1
*
      JORBMIN = JORB
      JORBMAX = JORB + NJORB - 1
*
      NIJ = NIORB*NJORB
*
      KEND = MIN(NKSTR,KMAX)
      IF(KEND.LT.NKSTR) THEN
        IEND = 0
      ELSE
        IEND = 1
      END IF
      NK = KEND-KMIN+1
*
      DO KKSTR = KMIN,KEND 
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Occupation of string ', KKSTR
         CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
       END IF
*. Loop over electrons after which JORB can be added 
       DO JEL = 0, NKEL
*
         IF(JEL.EQ.0 ) THEN     
           JORB1 = JORBMIN - 1
         ELSE
           JORB1 = MAX(JORBMIN-1,KSTR(JEL,KKSTR))
         END IF
         IF(JEL.EQ.NKEL) THEN
           JORB2 = JORBMAX + 1
         ELSE
           JORB2 = MIN(JORBMAX+1,KSTR(JEL+1,KKSTR))
         END IF
         IF(NTEST.GE.1000)
     &    WRITE(6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
*
         IF(JEL.GT.0.AND.JORB1.GE.JORBMIN.AND.
     &                   JORB1.LE.JORBMAX) THEN
*. vanishing for any IORB
           IJOFF = (JORB1-JORBMIN)*NIORB 
           DO IIORB = 1, NIORB
             IJ = IJOFF + IIORB
             if(ij.gt.nij) then
               write(6,*) ' ij .gt. nij '
               write(6,*) ' JORB1 IIORB' , JORB1,IIORB 
               write(6,*) ' ijoff ', ijoff 
               stop 
             end if
             I1(KKSTR-KMIN+1,IJ) = 0   
             XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
           END DO
         END IF
*
         IF(JORB1.LT.JORBMAX.AND.JORB2.GT.JORBMIN) THEN
*. Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
           SIGNJ = (-1) ** JEL * SCLFAC
*. reverse lexical number of the first JEL ELECTRONS
           ILEX0 = 1
           DO JJEL = 1, JEL  
             ILEX0 = ILEX0 + IZ(KSTR(JJEL,KKSTR),JJEL)
           END DO
           DO JJORB = JORB1+1, JORB2-1
* And electron JEL + 1
             ILEX1 = ILEX0 + IZ(JJORB,JEL+1)
*. Add electron IORB
             DO IEL = JEL, NKEL
               IF(IEL.EQ.JEL) THEN
                 IORB1 = MAX(JJORB,IORBMIN-1)
               ELSE
                 IORB1 = MAX(IORBMIN-1,KSTR(IEL,KKSTR))
               END IF
               IF(IEL.EQ.NKEL) THEN
                 IORB2 = IORBMAX+1
               ELSE 
                 IORB2 = MIN(IORBMAX+1,KSTR(IEL+1,KKSTR))
               END IF
               IF(NTEST.GE.1000)
     &          WRITE(6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
               IF(IEL.GT.JEL.AND.IORB1.GE.IORBMIN.AND.
     &                           IORB1.LE.IORBMAX) THEN
                 IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
             if(ij.gt.nij) then
               write(6,*) ' ij .gt. nij '
               write(6,*) ' JJORB IORB1' , JJORB,IORB1 
               write(6,*) ' ijoff ', ijoff 
               stop 
             end if
                 I1(KKSTR-KMIN+1,IJ) = 0
                 XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
               END IF
               IF(IORB1.LT.IORBMAX.AND.IORB2.GT.IORBMIN) THEN
*. Orbitals IORB1+1 - IORB2 -1 can be added after ELECTRON IEL in KSTR
*. Reverse lexical number of the first IEL+1 electrons
                 ILEX2 = ILEX1
                 DO IIEL = JEL+1,IEL
                   ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+1)
                 END DO
*. add terms for the last electrons
                 DO IIEL = IEL+1,NKEL
                   ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+2)
                 END DO
                 IJOFF = (JJORB-JORBMIN)*NIORB 
                 SIGNIJ =  SIGNJ*(-1.0D0) ** (IEL+1)
                 DO IIORB = IORB1+1, IORB2-1
                   IJ = IJOFF + IIORB - IORBMIN + 1
                   IF(IJ.LE.0.OR.IJ.GT.NIJ) THEN
                     WRITE(6,*) ' PROBLEMO ADADS1 : IJ : ', IJ
                     WRITE(6,*) ' IJOFF IORBMIN ', IJOFF,IORBMIN
                     WRITE(6,*) ' IIORB JJORB ', IIORB,JJORB
                     stop
                     NTEST = 1000
                   END IF
                   ILEX = ILEX2 + IZ(IIORB,IEL+2)
                   IACT = IREO(ILEX)
                   IF(NTEST.GE.1000) 
     &             WRITE(6,*) ' IIORB JJORB ', IIORB,JJORB
                   IF(NTEST.GE.1000) 
     &             WRITE(6,*) ' IJ ILEX,IACT',IJ, ILEX,IACT
                   IF(NTEST.GE.1000) 
     &             WRITE(6,*) ' ILEX0 ILEX1 ILEX2 ILEX ',   
     &                          ILEX0,ILEX1,ILEX2,ILEX
                   I1(KKSTR-KMIN+1,IJ) = IACT
                   XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
                   IF(IJ.LT.0) THEN
                     STOP ' NEGATIVE IJ in ADADS1 '
                   END IF
                 END DO
               END IF
             END DO
           END DO
         END IF
       END DO
      END DO
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADADST1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          IJ = 0
          DO  JJORB = JORB,JORB+NJORB-1
            JJORBR = JJORB-JORB+1
            DO  IIORB = IORB, IORB + NIORB - 1
              IJ = IJ + 1
C?            WRITE(6,*) ' IJ = ', IJ
C?            IF(IIORB.GT.JJORB) THEN
                IIORBR = IIORB - IORB + 1
                WRITE(6,*)
     &          ' Info for orbitals (iorb,jorb) ', IIORB,JJORB
                WRITE(6,*) ' Excited strings and sign '
                CALL IWRTMA(I1(1,IJ),1,NK,1,NK)
                CALL WRTMAT(XI1S(1,IJ),1,NK,1,NK)
C?            END IF
            END DO
          END DO
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE ADST_GAS(IORB,NORB,ISPGP,ISM,ITP,KMIN,KMAX,
     &                    I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12,
     &                    KACT,SCLFAC)
*
*
* Obtain mappings
* a+JORB !KSTR> = +/-!ISTR> for orbitals IORB - IORB+NORB-1
* and KSTR in the range KMIN-KMAX
*. All orbitals are assumed to belong to the same TS class
*. and are type-symmetry ordered
* 
* The results are given in the form
* I1(KSTR,JORB) =  ISTR if A+JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,ITP .
*
* if some nonvanishinr excitations were found, KACT is set to 1,
* else it is zero

* (numbering relative to TS start)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* If all K strings have been searched IEND is set to 1
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version 
*
* ======
*. Input
* ======
*
*./BIGGY
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Local scratch
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. save internally
      COMMON/SSAVE/NELIS(4), NSTRKS(4)
*
* =======
*. Output
* =======
*
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      CALL LUCIAQENTER('ADST  ')
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' ADST_GAS in service '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Problemo, ADST_GAS'
C?      WRITE(6,*) 'SCLFAC = ',SCLFAC
C?    END IF
*. Internal affairs 
      IF(I12.LE.4.AND.K12.LE.2 ) THEN
        KLLOC = KLOCSTR(K12)
        KLLZ = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*)
     &   ' ADST_GAS : Illegal value of I12,K12 = ', I12,K12
        STOP' ADST_GAS : Illegal value of I12,K12  '
      END IF
*
*. Supergroup and symmetry of K strings
*
      IORBTP = ITPFTO(IORB)
      IORBSM = ISMFTO(IORB)
      IF(NTEST.GE.100) WRITE(6,*) 
     &' ADST : IORB IORBTP IORBSM : ', IORB,IORBTP,IORBSM
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      CALL NEWTYP(ISPGPABS,1,IORBTP,1,KSPGPABS)
      CALL SYMCOM(2,0,IORBSM,KSM,ISM)
      IF(NTEST.GE.100) WRITE(6,*) 
     & ' KSM and KSPGPABS : ', KSM,KSPGPABS
*. In order to construct mapping we need to generate 
*.      Reorder array : lexical order => actual order for I strings
*.      Occupation of K strings
*
        IF(IFRST.NE.0) THEN
*. Arc weights for ISPGP
          NTEST2 = NTEST
          CALL WEIGHT_SPGP(int_mb(KLLZ),NGAS,
     &                     NELFSPGP(1,ISPGPABS),
     &                     NOBPT,int_mb(KLZSCR),NTEST2)
*. Reorder array for I strings
         NELI = NELFTP(ITP)
         NELIS(I12) = NELI
C?       CALL LUCIAQENTER('GENSTR')
         CALL GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,
     &                          int_mb(KLLOC),NOCOB,
     &                          1,int_mb(KLLZ),int_mb(KLLREO))
C?       CALL LUCIAQEXIT('GENSTR')
       END IF
       NELK = NELIS(I12) - 1
       IF(KFRST.NE.0) THEN
*. Occupation of K strings
CT       CALL LUCIAQENTER('GENSTR')
         CALL GETSTR_TOTSM_SPGP(1,KSPGPABS,KSM,NELK,NSTRK,
     &                         int_mb(KLLOC),NOCOB,
     &                         0,IDUM,IDUM)
CT       CALL LUCIAQEXIT('GENSTR')
         NSTRKS(K12) = NSTRK
      END IF
*. a+i !kstr> = +/- !istr>
      NSTRK = NSTRKS(K12)
      CALL ADS1_GAS(NK,I1,XI1S,LI1,IORB,NORB,
     &          int_mb(KLLOC),NELK,NSTRK,int_mb(KLLREO),int_mb(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,KACT,SCLFAC)
*
C     SAVE NELK , NSTRK
      CALL LUCIAQEXIT('ADST  ')
      RETURN
      END
      SUBROUTINE ADS1_GAS(NK,I1,XI1S,LI1,IORB,LORB,
     &                KSTR,NKEL,NKSTR,IREO,IZ,
     &                NOCOB,KMAX,KMIN,IEND,KACT,SCLFAC)
*
* Obtain I1(KSTR) = +/- A+ IORB !KSTR>
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : Firat orbital to be added 
* LORB : Number of orbitals to be added : IORB to IORB-1+LORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 => a + JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      NTEST = 000
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ============== '
       WRITE(6,*) ' ADSTS_GAS speaking '
       WRITE(6,*) ' ============== '
       WRITE(6,*) ' IORB,LORB ', IORB,LORB       
      END IF
*
      KACT = 0
      IORBMIN = IORB
      IORBMAX = IORB + LORB - 1
*
      KEND = MIN(NKSTR,KMAX)
      IF(KEND.LT.NKSTR) THEN
        IEND = 0
      ELSE
        IEND = 1
      END IF
      NK = KEND-KMIN+1
*. Start by zeroing
      ZERO = 0.0D0
      IZERO = 0
      DO IIORB = 1, LORB
C       CALL SETVEC(XI1S(1,IIORB),ZERO,NK)
        CALL ISETVC(I1(1,IIORB),IZERO,NK)
      END DO
*
      DO KKSTR = KMIN,KEND 
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Occupation of string ', KKSTR
         CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
       END IF
*
       SGN = -1.0D0*SCLFAC
       DO KEL = 1, NKEL+1
         SGN = - SGN
         IF(KEL.LE.NKEL) THEN
           KORB = KSTR(KEL,KKSTR)
         ELSE
           KORB = IORB+LORB
         END IF
*          
         IF(KORB.GT.IORBMIN) THEN
*. Orbitals can be added before electron KEL
           IORB2 = MIN(KORB-1,IORBMAX)
           IF(KEL.EQ.1) THEN
             IORB1 = IORBMIN
           ELSE 
             IORB1 = MAX(IORBMIN,KSTR(KEL-1,KKSTR)+1)
           END IF
           IF(IORB1.GT.IORBMAX) GOTO 2912
           IF(IORB2.GE.IORB1) THEN
*
*. common piece of lexical ordering 
C            IF(KEL.EQ.1) THEN
C              IREVLEX = 1
C              DO JEL = 1, NKEL
C                 IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL+1)
C              END DO
C            ELSE
C                 IREVLEX = IREVLEX + IZ(KSTR(KEL-1,KKSTR),KEL-1)
C    &                              - IZ(KSTR(KEL-1,KKSTR),KEL)
C            END IF
C          IF(IORB2.GE.IORB1) THEN
*
             IREVLEX = 1
             DO JEL = 1, KEL-1
               IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL)
             END DO
             DO JEL = KEL, NKEL
               IREVLEX = IREVLEX + IZ(KSTR(JEL,KKSTR),JEL+1)
             END DO
C?           IF(NTEST.GE.1000) WRITE(6,*)
C?   &       'KEL IORB1,IORB2,IREVLEX ',KEL,IORB1,IORB2,IREVLEX
             DO IIORB = IORB1, IORB2
               I1(KKSTR-KMIN+1,IIORB-IORBMIN+1) = 
     &         IREO(IREVLEX  + IZ(IIORB,KEL))
               XI1S(KKSTR-KMIN+1,IIORB-IORBMIN+1) = SGN
*
               KACT = 1
             END DO
           END IF
         END IF
       END DO
 2912  CONTINUE
      END DO
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ASTR_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO 200 IIORB = IORB,IORB+LORB-1
            IIORBR = IIORB-IORB+1
            WRITE(6,*) ' Info for orbital ', IIORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(I1(1,IIORBR),1,NK,1,NK)
            CALL WRTMAT(XI1S(1,IIORBR),1,NK,1,NK)
  200     CONTINUE
        END IF
      END IF
*
C     WRITE(6,*) 'JEPPE forced me to stop in ADS1   '
C     STOP'JEPPE forced me to stop in ADS1   '
      RETURN
      END
      SUBROUTINE DXTYP2_GAS(NDXTP,ITP,JTP,KTP,LTP,
     &                     NOBTP,IL,IR,IPHGAS)
*
* Obtain types of I,J,K,l so
* <L!a+I a+K a L a J!R> is nonvanishing
* only combinations with type(I) .ge. type(K) and type(J).ge.type(L)
* are included
*
* Intermediate occupations less than zero allowed for particle spaces
* (IPHGAS=2)
*
*
      INTEGER IL(NOBTP),IR(NOBTP),IPHGAS(NOBTP)
      INTEGER ITP(*),JTP(*),KTP(*),LTP(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DXTYP_GAS in action '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' Occupation of left string '
        CALL IWRTMA(IL,1,NOBTP,1,NOBTP)
        WRITE(6,*) ' Occupation of right string '
        CALL IWRTMA(IR,1,NOBTP,1,NOBTP)
      END IF
*
*. Number of differing occupations
      NANNI = 0
      NCREA = 0
      NDIFT = 0
*
      ICREA1 = 0
      ICREA2 = 0
      IANNI1 = 0
      IANNI2 = 0
      DO IOBTP = 1, NOBTP
        NDIFT = NDIFT + ABS(IL(IOBTP)-IR(IOBTP))
        NDIF = IL(IOBTP)-IR(IOBTP)
        IF(NDIF.EQ.2) THEN
*. two electrons of type IOBTP must be created
          ICREA1 = IOBTP
          ICREA2 = IOBTP
          NCREA = NCREA + 2
        ELSE IF (NDIF .EQ. -2 ) THEN
*. Two electrons of type IOBTP must be annihilated
          IANNI1 = IOBTP
          IANNI2 = IOBTP
          NANNI = NANNI + 2
        ELSE IF (NDIF.EQ.1) THEN
*. one electron of type IOBTP must be created
          IF(NCREA.EQ.0) THEN
            ICREA1 = IOBTP
          ELSE
            ICREA2 = IOBTP
          END IF
          NCREA = NCREA + 1
        ELSE IF (NDIF.EQ.-1) THEN
* One electron of type IOBTP must be annihilated
          IF(NANNI.EQ.0) THEN
            IANNI1 = IOBTP
          ELSE
            IANNI2 = IOBTP
          END IF
          NANNI = NANNI + 1
        END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NCREA, NANNI ', NCREA, NANNI
        WRITE(6,*)  ' ICREA1, IANNI1 ', ICREA1,IANNI1
        WRITE(6,*)  ' ICREA2, IANNI2 ', ICREA2,IANNI2
      END IF
*
      NDXTP = 0
      IF(NDIFT.GT.4) THEN
        NDXTP = 0
      ELSE
      IF(NCREA.EQ.0.AND.NANNI.EQ.0) THEN
*. strings identical, include diagonal excitions  itp = jtp, ktp=ltp 
        DO IJTP = 1, NOBTP
          IF(IR(IJTP).GE.1.OR.IPHGAS(IJTP).EQ.2) THEN
            DO KLTP = 1, IJTP 
             IF((IJTP.NE.KLTP.AND.(IR(KLTP).GE.1.OR.IPHGAS(KLTP).EQ.2))
     &      .OR.(IJTP.EQ.KLTP.AND.(IR(KLTP).GE.2.OR.IPHGAS(KLTP).EQ.2)) 
     &      ) THEN
                 NDXTP = NDXTP + 1
                 ITP(NDXTP) = IJTP
                 JTP(NDXTP) = IJTP
                 KTP(NDXTP) = KLTP
                 LTP(NDXTP) = KLTP
              END IF
            END DO
          END IF
        END DO
*. Strings differ by single excitation
      ELSE IF( NCREA.EQ.1.AND.NANNI.EQ.1) THEN
*. diagonal excitation plus creation in ICREA1 and annihilation in IANNI1
        DO IDIA = 1, NOBTP
          IF((IDIA.NE.IANNI1.AND.(IR(IDIA).GE.1.OR.IPHGAS(IDIA).EQ.2))
     &   .OR.(IDIA.EQ.IANNI1.AND.(IR(IDIA).GE.2.OR.IPHGAS(IDIA).EQ.2))
     &   )THEN
             NDXTP = NDXTP + 1
             ITP(NDXTP) = MAX(ICREA1,IDIA)
             KTP(NDXTP) = MIN(ICREA1,IDIA)
             JTP(NDXTP) = MAX(IANNI1,IDIA)
             LTP(NDXTP) = MIN(IANNI1,IDIA)
          END IF
        END DO
      ELSE IF(NCREA.EQ.2.AND.NANNI.EQ.2) THEN
*. Strings differ by double excitation
        NDXTP = 1
        ITP(1) = ICREA2
        KTP(1) = ICREA1
        JTP(1) = IANNI2
        LTP(1) = IANNI1
      END IF
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,I4)')
     &  ' Number of connecting double excitations ', NDXTP
        IF(NDXTP.NE.0) THEN
          WRITE(6,*) '  ITYP KTYP LTYP JTYP '
          WRITE(6,*) '  ===================='
          DO  IDX = 1,NDXTP
            WRITE(6,'(1H ,4I5)')ITP(IDX),KTP(IDX),LTP(IDX),JTP(IDX)
          END DO
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,
     &                  ISPGRPCON,IPRNT,iganodeid)
*
* FInd connection matrix for string types
*
* ISPGRPCON(ISPGP,JSPGRP) = 0 => spgrps are identical
*                         = 1 => spgrps are connected by single excitation
*      .                  = 2 => spgrps are connected by double excitation
*       .              . ge.3 => spgrps are connected by triple or
*        .                       higher excitation
*
* Jeppe Olsen, September 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IELFSPGRP(MXPNGAS,*)
*. output
      DIMENSION ISPGRPCON(NSPGRP,NSPGRP)
*
      NTEST = 000
cVOG  NTEST = MAX(NTEST,IPRNT)
*
      DO ISPGRP = 1, NSPGRP
        ISPGRPA = IOFSPGRP-1+ISPGRP
        DO JSPGRP = 1, ISPGRP
          JSPGRPA = IOFSPGRP-1+JSPGRP
          IDIF = 0
          DO IGAS = 1, NGAS
            IDIF = IDIF 
     &    + ABS(IELFSPGRP(IGAS,ISPGRPA)-IELFSPGRP(IGAS,JSPGRPA))
          END DO
          NEXC = IDIF/2
          ISPGRPCON(ISPGRP,JSPGRP) = NEXC
          ISPGRPCON(JSPGRP,ISPGRP) = NEXC
        END DO
      END DO
*
      IF(NTEST.GE.100.and.iganodeid.eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) '==================== '
        WRITE(6,*) 'output from SPGRPCON '
        WRITE(6,*) '==================== '
        WRITE(6,*)
        NEXC1 = 0
        NEXC2 = 0
        DO ISPGRP=1, NSPGRP
          DO JSPGRP = 1, NSPGRP
            IF(ISPGRPCON(ISPGRP,JSPGRP).EQ.1) THEN
              NEXC1 = NEXC1 + 1
            ELSE IF(ISPGRPCON(ISPGRP,JSPGRP).EQ.2) THEN
              NEXC2 = NEXC2 + 1
            END IF
          END DO
        END DO
*
        WRITE(6,*) 
     &  ' single excitation interactions',NEXC1,
     &   '( ',float(NEXC1)*100.0/NSPGRP**2,' % ) '
        WRITE(6,*) 
     &  ' double excitation interactions',NEXC2,
     &   '( ',float(NEXC2)*100.0/NSPGRP**2,' % ) '
*
      END IF
*
      IF(NTEST.GE.1000.and.iganodeid.eq.0) THEN
         WRITE(6,*) ' Supergroup connection matrix '
         CALL IWRTMA(ISPGRPCON,NSPGRP,NSPGRP,NSPGRP,NSPGRP)
      END IF
*
      RETURN
      END 
      SUBROUTINE ADSTN_GAS(IOBSM,IOBTP,ISPGP,ISPGPSM,ISPGPTP,
     &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC)
*
*
* Obtain mappings
* a+IORB !KSTR> = +/-!ISTR> for orbitals of symmetry IOBSM and type IOBTP
* and I strings belonging to supergroup ISPGP wih symmetry ISPGPSM
* and type ISPGPTP(=1=>alpha,=2=>beta)
* 
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR> 
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version 
*               October 96   : Improved version
*
* ======
*. Input
* ======
*
*./BIGGY
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8*8*8)
      DIMENSION IOFFI(LOFFI)
      PARAMETER(MXLNGAS=7)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension 
* (NKSTR,*), Where NKSTR is the number of K-strings of 
*  correct symmetry . Nk is provided by this routine.
*
      CALL LUCIAQENTER('ADSTN ')
*
      IF(NGAS.GT.MXLNGAS) THEN
        WRITE(6,*) ' Ad hoc programming in ADSTN (IOFFI)'
        WRITE(6,*) ' Must be changed - or redimensioned '
        STOP'ADST : IOFFI problem '
      END IF
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADSTN_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) '  ISPGP ISPGPSM ISPGPTP :  ',
     &                ISPGP,ISPGPSM,ISPGPTP
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Problemo : ADSTN_GAS'
C?      WRITE(6,*) ' SCLFAC .ne. 1 '
C?    END IF
*
*. Supergroup and symmetry of K strings
*
      ISPGRPABS = IBSPGPFTP(ISPGPTP)-1+ISPGP
      CALL NEWTYP(ISPGRPABS,1,IOBTP,1,KSPGRPABS)
      CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
      NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
      IF(NTEST.GE.200) WRITE(6,*) 
     & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*
      NORBTS= NOBPTS(IOBTP,IOBSM)
      ZERO =0.0D0
      CALL SETVEC(XI1S,ZERO,NORBTS*NKSTR)
      IZERO = 0    
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
*. First orbital of given GASSpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1+NINOB
*. First orbital of given GASSPace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM) 


*
*. Information about I strings
* =============================
*
*. structure of group of strings defining I strings 
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Number of electrons before active type
      NELB = 0
      DO IGAS = 1, IOBTP -1
        NELB = NELB + NELFGS(IGAS)
      END DO
*. Number of electrons in active space 
      NACGSOB = NOBPT(IOBTP)
   
*. Number of strings per symmetry for each symmetry
      DO IGAS = 1, NGAS
        CALL ICOPVE2(int_mb(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
      END DO
*. Offset and dimension for active group in I strings
      CALL ICOPVE2(int_mb(KISTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               IACIST)
      CALL ICOPVE2(int_mb(KNSTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               NACIST)
C?     WRITE(6,*) ' IACIST and NACIST arrays '
C?     CALL IWRTMA(IACIST,1,NSMST,1,NSMST)
C?     CALL IWRTMA(NACIST,1,NSMST,1,NSMST)
*
*. Generate offsets for I strings with given symmetry in
*  each space
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 2, NGASL 
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3
     &   (ISMFGS(2),NGASL-1,MINVAL(2),MAXVAL(2),NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of Each of NGASL-1 spaces given, symmetry of full space
        ISTSMM1 = 1
        DO IGAS = 2, NGASL 
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE 1
        CALL SYMCOM(2,1,ISTSMM1,ISMGS1,ISPGPSM)
        ISMFGS(1) = ISMGS1
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
*
        IOFFI(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' IOFF, IOFFI(IOFF) NSTRII ',
     &                 IOFF, IOFFI(IOFF),NSTRII
        END IF
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE


*
*. Supergroup and symmetry of K strings
*
CM    CALL NEWTYP(ISPGRPABS,1,IOBTP,1,KSPGRPABS)
CM    CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
CM    NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
CM    IF(NTEST.GE.200) WRITE(6,*) 
CM   & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
*
*. Gas structure of K strings
*
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,KSPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Active group of K-strings
      KACGRP = ITPFGS(IOBTP)
*. Number of strings per symmetry distribution      
      DO IGAS = 1, NGAS
        CALL ICOPVE2(int_mb(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(int_mb(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
 1000 CONTINUE
        IF(KFIRST .EQ. 1 ) THEN
          DO IGAS = 2, NGASL 
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in 2- NGASL 
         CALL NXTNUM3(ISMFGS(2),NGASL-1,MINVAL(2),MAXVAL(2),NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        KFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of spaces 2- NGASL spaces '
          CALL IWRTMA(ISMFGS(2),NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of each of NGASL spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 2, NGASL 
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE 1
        CALL SYMCOM(2,1,ISTSMM1,ISMGS1,KSM)
        ISMFGS(1) = ISMGS1
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IOBTP)
        CALL  SYMCOM(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        ISMFGS(IOBTP) = IACSM
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGAS
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
        ISMFGS(IOBTP) = ISAVE
        IBSTRINI = IOFFI(IOFF)
C?      WRITE(6,*) ' IOFF IBSTRINI ', IOFF,IBSTRINI
*. Number of strings before active GAS space
        NSTB = 1
        DO IGAS = 1, IOBTP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings before active GAS space
        NSTA = 1
        DO IGAS =  IOBTP+1, NGAS
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group 
C?      write(6,*) ' IACSM = ', IACSM
        NIAC  = NACIST(IACSM)
        IIAC =  IACIST(IACSM)
*
        NKAC = NNSTSGP(ISMFGS(IOBTP),IOBTP)
        IKAC = IISTSGP(ISMFGS(IOBTP),IOBTP)
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
C?      write(6,*) ' nstb nsta niac nkac ',
C?   &               nstb,nsta,niac,nkac
*. Obtain annihilation n mapping for all strings of this type
*
        NORBTS= NOBPTS(IOBTP,IOBSM)
*
        NKACT = NSTFGP(KACGRP)
C?      write(6,*) ' KACGRP ', KACGRP
        CALL ADSTN_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,   
     &          int_mb(KSTSTM(KACGRP,1)),int_mb(KSTSTM(KACGRP,2)),
     &          IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &          NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC)
        KSTRBS = KSTRBS + NKSD     
        IF(NGASL-1.GT.0) GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADSTN_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMAT(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
      CALL LUCIAQEXIT('ADSTN ')
      RETURN
      END
      SUBROUTINE ADSTN_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,
     &              IOFFKSP,ICREORB,ICRESTR,
     &              IORBTSF,IORBTF,NORBTS,NSTAK,NSTAKT,NSTAI,
     &              NSTAKTS,ISTAKTS,NELB,NACGSOB,
     *              ISTMAP,SGNMAP,SCLFAC)
*
* Creation mappings from K-strings of given sym in each gasspace
*
* Input 
* NSTB : Number of strings before active gasspace
* NSTA : Number of strings after accive gasspace
* IOFFK : Offset for K group of strings in active gasspace, i.e. start of
*         this symmetry of active K group strings
* IOFFI : Offset for I group of strings in active gasspace, i.e. start of
*         this symmetry of active I group strings
* IOFFISP: Offset for this symmetrydistribution of active I supergroup strings
* IOFFKSP: Offset for this symmetrydistribution of active K supergroup strings
* ICREORB : Orbital part of creation map for active K groupstrings
* ICRESTR : String  part of creation map for active K groupstrings
* IORBTSF   : First active orbital ( first orbital in in active GASspace
*           with required sym)
* IORBTF   : First orbital in active gas space, (can have any sym)
* NORBTS  : Number of orbitals of given symmetry and type            
* NSTAK : Number of K groupstrings with given correct symmetry      
* NSTAKT: Total Number of K groupstrings in active group (all symmetries)
* NSTAKTS: Total Number of K supergroup strings with correct symmetry
* ISTAKTS: Offset for K supergroup strings with hiven symmetrydistribution
* NSTAI : Number of I groupstrings in active gasspace
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION ICREORB(NACGSOB,*), ICRESTR(NACGSOB,*)
*. Output
      DIMENSION ISTMAP(NSTAKTS,*),SGNMAP(NSTAKTS,*)
* 
C?    WRITE(6,*) ' ADSTN_GASSM : NSTA, NSTB, NSTAK',NSTA,NSTB,NSTAK
C?    WRITE(6,*) ' IOFFISP,IOFFKSP', IOFFISP, IOFFKSP
C?    WRITE(6,*) ' IORBTSF IORBTF ', IORBTSF,IORBTF
      IMULTK = NSTAK*NSTB
      IMULTI = NSTAI*NSTB
C?    WRITE(6,*) ' NSTAKT ', NSTAKT 
*
      SIGN0 = (-1)**NELB*SCLFAC
C?    WRITE(6,*) ' NELB sign0 = ', NELB, SIGN0
      DO KSTR = IOFFK, NSTAK+IOFFK-1
        DO IORB = IORBTSF, IORBTSF-1+NORBTS
*. Relative to Type-symmetry start
          IORBRTS = IORB-IORBTSF+1
*. Relative to type start
          IORBRT = IORB-IORBTF+1
C?         write(6,*) 'IORB IORBRT KSTR ', IORB,IORBRT, KSTR
C?         WRITE(6,*) 'ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)',
C?   &                 ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)
          IF(ICREORB(IORBRT,KSTR) .GT. 0 ) THEN
*. Excitation is open, corresponding active I string
            IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
              SIGN = SIGN0
              ISTR = ICRESTR(IORBRT,KSTR)
            ELSE
              SIGN = -SIGN0
              ISTR = -ICRESTR(IORBRT,KSTR)
            END IF
* Relative to start of given symmetry for this group
            ISTR = ISTR - IOFFI+ 1
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
            IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
            IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
C?          write(6,*) ' ISTR IADRK0 IADRI0 = ', ISTR, IADRK0,IADRI0 
*
            NSTAINSTA = NSTAI*NSTA
            NSTAKNSTA = NSTAK*NSTA
            DO IB = 1, NSTB
              DO IA = 1, NSTA
C               IBKA = IADRI0 + (IB-1)*NSTAI*NSTA+IA
C               KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
C?              write(6,*) ' IBKA, KBKA ',IBKA,KBKA
                ISTMAP(IADRK0+IA,IORBRTS) = IADRI0 + IA
                SGNMAP(IADRK0+IA,IORBRTS) = SIGN
              END DO
              IADRI0 = IADRI0 +  NSTAINSTA
              IADRK0 = IADRK0 +  NSTAKNSTA
            END DO
C         ELSE
C            SIGN = 0.0D0
C            ISTR = 0
*. This Creation is inactive for all choices of strings in supergroup
*. before and after the active type. 
C           IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
C?          write(6,*) ' ISTR IADRK0 = ', ISTR, IADRK0
*
C           DO IB = 1, NSTB
C             DO IA = 1, NSTA
C               KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
C?              write(6,*) ' IBKA, KBKA ',IBKA,KBKA
C               IF(ISTMAP(KBKA,IORBRTS).NE.999) THEN
C                 WRITE(6,*) ' overwriting ??? '
C                 WRITE(6,*) ' Element ', (IORBRTS-1)*NSTAKTS+KBKA     
C                 STOP
C               END IF
C               ISTMAP(KBKA,IORBRTS) = ISTR 
C               SGNMAP(KBKA,IORBRTS) = SIGN
C             
C             END DO
C           END DO
          END IF
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
COLD        IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
COLD        IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
C?          write(6,*) ' ISTR IADRK0 IADRI0 = ', ISTR, IADRK0,IADRI0 
*
COLD        DO IB = 1, NSTB
COLD          DO IA = 1, NSTA
COLD            IBKA = IADRI0 + (IB-1)*NSTAI*NSTA+IA
COLD            KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
C?              write(6,*) ' IBKA, KBKA ',IBKA,KBKA
COLD            ISTMAP(KBKA,IORBRTS) = IBKA
COLD            SGNMAP(KBKA,IORBRTS) = SIGN
COLD          END DO
COLD        END DO
C         END IF
*
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADSTN_GASSM '
        WRITE(6,*) ' ======================== '
        NK = NSTB*NSTAK*NSTA
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO IORB = IORBTSF,IORBTSF + NORBTS  - 1 
            IORBR = IORB-IORBTSF+1
            WRITE(6,*) ' Update Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(ISTMAP(1,IORBR),1,NK,1,NK)
            CALL WRTMAT(SGNMAP(1,IORBR),1,NK,1,NK)
          END DO
        END IF
      END IF

      RETURN
      END
      SUBROUTINE COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
*
* Two lists of excitations/annihilations/creations are given.
* COmpress to common nonvanishing entries
*
* Jeppe Olsen, November 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION I1(NKIN,N1),XI1(NKIN,N1)
      DIMENSION I2(NKIN,N2),XI2(NKIN,N2)
*
      NKOUT = 0
      DO K = 1, NKIN
        I1ACT  = 0
        DO I = 1, N1
          IF(I1(K,I).NE.0) I1ACT = 1
        END DO
        I2ACT = 0
        DO I = 1, N2
          IF(I2(K,I).NE.0) I2ACT = 1
        END DO
        IF(I1ACT.EQ.1.AND.I2ACT.EQ.1) THEN
          NKOUT = NKOUT + 1
          IF(NKOUT.NE.K) THEN
            DO I = 1, N1
               I1(NKOUT,I) = I1(K,I)
              XI1(NKOUT,I) =XI1(K,I)
            END DO
            DO I = 1, N2
               I2(NKOUT,I) = I2(K,I)
              XI2(NKOUT,I) =XI2(K,I)
            END DO
          END IF
        END IF
      END DO
*
      RETURN
      END
      subroutine tlucinp
*
* is LUCINP overwritten ??
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      write(6,*) 'nactel nsmob ',nactel,nsmob
*
      return
      end
      SUBROUTINE ABTOR2(SKII,CKJJ,NKA,NIB,NJB,NKB,RHO2B,
     &                  NI,NJ,NK,NL,MAXK,
     &                  KBIB,XKBIB,KBJB,XKBJB,IKORD)
*
* Obtain contributions alpha-beta contributions to two-particle
* density matrix 
*
* Rho2b(ij,kl)  = RHo2b(ij,kl) 
*               + sum(Ka) Skii(Ka,i,Ib)<Ib!Eb(kl)!Jb> Ckjj(Ka,j,Jb)
*
*
* Jeppe Olsen, Fall of 96   
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKJJ(*),SKII(*)
      DIMENSION RHO2B(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*
      IF(IKORD.NE.0) THEN
        WRITE(6,*) ' ABTOR2 : IKORD .NE. 0 '
        WRITE(6,*) ' I am not ready for this '
        STOP     ' ABTOR2 : IKORD .NE. 0 '
      END IF
*
C     CALL LUCIAQENTER('ABTOR')
*. Excitations <Ib!Eb(kl)!Jb>
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
                   KLOFF= ((L-1)*NK + K - 1 )*NI*NJ + 1
                   IMAX = NI
*
C                  IF(IKORD.NE.0) THEN
*. Restrict so (ij) .le. (kl)
C                    IMAX  = K
C                    JKINTOF = INTOF + (K-1)*NJ
C                    DO J = L,NL
C                      XIJILS(J) = XIJKL(JKINTOF-1+J)  
C                    END DO
C                    XIJKL(JKINTOF-1+L) = 0.5D0*XIJKL(JKINTOF-1+L)
C                    DO J = L+1, NL
C                     XIJKL(JKINTOF-1+J) = 0.0D0
C                    END DO
C                  END IF
                   ONE = 1.0D0
                   CALL MATML7(RHO2B(KLOFF),SKII(ISOFF),CKJJ(ICOFF),
     &                         NI,NJ,NKA,IMAX,NKA,NJ,
     &                         ONE,FACTOR ,1)
C                  IF(IKORD.NE.0) THEN
C                     DO J = L,NL
C                       XIJKL(JKINTOF-1+J) =  XIJILS(J) 
C                     END DO
C                  END IF
*
                 END IF
               END DO
*
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
*
C     CALL LUCIAQEXIT('ABTOR')
      RETURN
      END
      SUBROUTINE ADTOR2(RHO2,RHO2T,ITYPE,
     &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
*
* Add contributions to two electron density matrix RHO2
* output density matrix is in the form Rho2(ij,kl),(ij).ge.(kl)
*
*
* Jeppe Olsen, Fall of 96
*
*
* Itype = 1 => alpha-alpha or beta-beta loop
*              input is in form Rho2t(ik,jl), i.ge.k, j.ge.l
* Itype = 2 => alpha-beta loop
*              input is in form Rho2t(ij,kl)
*               
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION RHO2T(*)
*. Input and output
      DIMENSION RHO2(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to ADTOR2 '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' NI NJ NK NL = ', NI,NJ,NK,NL
        WRITE(6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
        WRITE(6,*) ' ITYPE = ',ITYPE
        IF(NTEST.GE.2000) THEN
          WRITE(6,*) ' Initial two body density matrix '
          CALL PRSYM(RHO2,NORB**2)
        END IF
        WRITE(6,*) ' RHO2T : '
        IF(ITYPE.EQ.1) THEN
          IF(IOFF.EQ.KOFF) THEN
            NROW = NI*(NI+1)/2
          ELSE
            NROW = NI*NK
          END IF
          IF(JOFF.EQ.LOFF) THEN
            NCOL = NJ*(NJ+1)/2
          ELSE
            NCOL = NJ*NL
          END IF
        ELSE IF (ITYPE.EQ.2) THEN
          NROW = NI*NJ
          NCOL = NK*NL
        END IF
        CALL WRTMAT(RHO2T,NROW,NCOL,NROW,NCOL)
      END IF
C?    WRITE(6,*) ' Enforced return in ADTOR2 '
C?    RETURN
      NELMNT = NORB**2*(NORB**2+1)/2
*
      NKK = -2810
      NII = -2810
      NJJ = -2810
      NLL = -2810
      LLOFF = -2810
      KKOFF = -2810
      NII = -2810
      IACTIVE = -2810
      SIGN = -.2810
      IF(ITYPE.EQ.1) THEN
*
* =======================================
*     Alpha-alpha or beta-beta term
* =======================================
*
*. Four permutations
      DO IPERM = 1, 4
        IF(IPERM.EQ.1) THEN
          NII = NI
          IIOFF = IOFF
          NJJ = NJ
          JJOFF = JOFF
          NKK = NK
          KKOFF = KOFF
          NLL = NL
          LLOFF = LOFF
          SIGN = 1.0D0
          IACTIVE = 1
        ELSE IF(IPERM.EQ.2) THEN
          IF(IOFF.NE.KOFF) THEN
            NII = NK
            IIOFF = KOFF
            NKK = NI
            KKOFF = IOFF
            NJJ = NJ
            JJOFF = JOFF
            NLL = NL
            LLOFF = LOFF
            IACTIVE = 1
          ELSE 
            IACTIVE = 0
          END IF
          SIGN = -1.0D0
        ELSE IF(IPERM.EQ.3) THEN
          IF(JOFF.NE.LOFF) THEN
            NII = NI
            IIOFF = IOFF
            NKK = NK
            KKOFF = KOFF
            NJJ = NL
            JJOFF = LOFF
            NLL = NJ
            LLOFF = JOFF
            SIGN = -1.0D0
            IACTIVE = 1
          ELSE
            IACTIVE = 0
          END IF
        ELSE IF(IPERM.EQ.4) THEN
          IF(IOFF.NE.KOFF.AND.JOFF.NE.LOFF) THEN
            NKK = NI
            KKOFF = IOFF
            NII = NK
            IIOFF = KOFF
            NJJ = NL
            JJOFF = LOFF
            NLL = NJ
            LLOFF = JOFF
            SIGN = 1.0D0
            IACTIVE = 1
          ELSE
            IACTIVE = 0
          END IF
        END IF
*
        IJOFF = (JJOFF-1)*NORB+IIOFF
        KLOFF = (LLOFF-1)*NORB+KKOFF
C       IF(IACTIVE.EQ.1.AND.IJOFF.GE.KLOFF) THEN
        IF(IACTIVE.EQ.1) THEN
          IJOFF = (JJOFF-1)*NORB+IIOFF
          KLOFF = (LLOFF-1)*NORB+LLOFF
            DO II = 1, NII
              DO JJ = 1, NJJ
                DO KK = 1, NKK
                  DO LL = 1, NLL
                    IJ = (JJ+JJOFF-2)*NORB + II+IIOFF - 1
                    KL = (LL+LLOFF-2)*NORB + KK+KKOFF - 1
                    IF(IJ.GE.KL) THEN
                      IJKL = IJ*(IJ-1)/2+KL
                      IF(IPERM.EQ.1) THEN
                        I = II
                        K = KK
                        J = JJ
                        L = LL
                      ELSE IF(IPERM.EQ.2) THEN
                        I = KK
                        K = II
                        J = JJ
                        L = LL
                      ELSE IF(IPERM.EQ.3) THEN
                        I = II
                        K = KK
                        J = LL
                        L = JJ
                      ELSE IF(IPERM.EQ.4) THEN
                        I = KK
                        K = II
                        J = LL
                        L = JJ
                      END IF
                      IF(IOFF.NE.KOFF) THEN
                        IKIND = (K-1)*NI+I
                        NIK = NI*NK
                        SIGNIK = 1.0D0
                      ELSE
                        IKIND = MAX(I,K)*(MAX(I,K)-1)/2+MIN(I,K)
                        NIK = NI*(NI+1)/2
                        IF(I.EQ.MAX(I,K)) THEN
                          SIGNIK = 1.0D0
                        ELSE
                          SIGNIK = -1.0D0
                        END IF
                      END IF
                      IF(JOFF.NE.LOFF) THEN
                        JLIND = (L-1)*NJ+J
                        SIGNJL = 1.0D0
                      ELSE
                        JLIND = MAX(J,L)*(MAX(J,L)-1)/2+MIN(J,L)
                        IF(J.EQ.MAX(J,L)) THEN
                          SIGNJL = 1.0D0
                        ELSE
                          SIGNJL = -1.0D0
                        END IF
                      END IF
                      IKJLT = (JLIND-1)*NIK+IKIND
                      IF(IJKL.GT.NELMNT) THEN
                         WRITE(6,*) ' Problemo 1 : IJKL .gt. NELMNT'
                         WRITE(6,*) ' IJKL, NELMNT',IJKL,NELMNT
                         WRITE(6,*) ' IJ, KL', IJ,KL
                         WRITE(6,*) ' JJ JJOFF ', JJ,JJOFF
                         WRITE(6,*) ' II IIOFF ', II,IIOFF
                         WRITE(6,*) ' IPERM = ', IPERM
                      END IF
                      RHO2(IJKL) = RHO2(IJKL) 
     &                           - SIGN*SIGNJL*SIGNIK*RHO2T(IKJLT)
*. The minus : Rho2t comes as <a+i a+k aj al>, but we want 
* <a+ia+k al aj>
                    END IF
                  END DO
                END DO
              END DO
            END DO
*. End of active/inactive if
        END IF
*. End of loop over permutations
      END DO
      ELSE IF(ITYPE.EQ.2) THEN
*
* =======================================
*     Alpha-alpha or beta-beta term
* =======================================
*
      DO I = 1, NI
       DO J = 1, NJ
         DO K = 1, NK
           DO L = 1, NL
             IJ = (J+JOFF-2)*NORB + I+IOFF - 1
             KL = (L+LOFF-2)*NORB + K+KOFF - 1
             IF(IJ.EQ.KL) THEN
               FACTOR = 2.0D0
             ELSE 
               FACTOR= 1.0D0
             END IF
             IJKL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
             IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI
     &             + (J-1)*NI + I
                      IF(IJKL.GT.NELMNT.OR.IJKL.LT.0) THEN
                         WRITE(6,*) ' Problemo 2 : IJKL .gt. NELMNT'
                         WRITE(6,*) ' IJKL, NELMNT',IJKL,NELMNT
                         WRITE(6,*) ' I,J,K,L = ', I,J,K,L
                         WRITE(6,*) ' NI,NJ,NK,NL=',NI,NJ,NK,NL
                         WRITE(6,*) ' IOFF, JOFF, KOFF, LOFF = ',
     &                   IOFF,JOFF,KOFF,LOFF
                         WRITE(6,*) ' NORB = ', NORB
                         STOP 'ADTOR2'
                      END IF
             RHO2(IJKL) = RHO2(IJKL) + FACTOR*RHO2T(IJKLT)
C?           WRITE(6,*) ' IJKL, IJKLT ', IJKL, IJKLT
            END DO
          END DO
        END DO
      END DO
*
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Updated two-body density matrix '
         CALL PRSYM(RHO2,NORB**2)
      END IF
*
      RETURN
      END
      SUBROUTINE FOCK_MAT_AB(F,I12L,NSPIN)
*
* Construct Fock matrix, spin-orbital version
*
* F(i,j) = SUM(K) H(i,K) * RHO1(j,K)
*          + SUM(M,K,L) I  (i M K L ) * RHO2( j M K L )
*
*     Andreas, based on Jeppes original below
*
* Unless I12L = 2, only one-electron part is calculated
*. Input
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*
      INCLUDE 'cintfo.inc'
      INCLUDE 'oper.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION F(*)
*
      NTEST = 00
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'FOO   ')
*
      ONE = 1.0D0
      ZERO = 0.0D0

*. for the moment we put F(beta) here:
      IFBETOFF = NTOOB**2
*. Largest set of orbitals with given symmetry and type
      MXTSOB = 0
      DO ISM = 1, NSMOB
      DO IGAS = 1, NGAS
        MXTSOB = MAX(MXTSOB,NOBPTS(IGAS,ISM))
      END DO
      END DO
C?    WRITE(6,*) 'MXTSOB = ', MXTSOB
*. Allocate scratch space for 2-electron integrals and 
*. two-electron densities
      MX4IBLK = MXTSOB ** 4
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT ') !done
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN ') !done
*. And a block of F
      MX2IBLK = MXTSOB ** 2
      CALL MEMMAN(KLFBLK,MX2IBLK,'ADDL  ',2,'KLFBL ')
*. 
*
      ONE = 1.0D0
      II = -2303
      IJ = -2303
      DO IJSPIN = 1, NSPIN
      DO IJSM = 1, NSMOB
        ISM = IJSM
        JSM = IJSM
        NIJS = NOCOBS(IJSM)
*
        IF(IJSM.EQ.1) THEN
         IFOFF = 1
        ELSE
         IFOFF = IFOFF+NOCOBS(IJSM-1)**2
        END IF
*
        DO JGAS = 1, NGAS
          IF(JGAS.EQ.1) THEN
            IJ = 1
          ELSE 
            IJ = IJ + NOBPTS(JGAS-1,JSM)
          END IF
          NJ = NOBPTS(JGAS,IJSM)
          DO IGAS = 1, NGAS
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) 
     &        ' Fock matrix for ISM IGAS JGAS',ISM,IGAS,JGAS
            END IF
            NI = NOBPTS(IGAS,ISM)
            IF(IGAS.EQ.1) THEN
              II = 1
            ELSE 
              II = II + NOBPTS(IGAS-1,ISM)
            END IF
*
*  =======================
*. block F(ijsm,igas,jgas)
*  =======================
*
            CALL SETVEC(dbl_mb(KLFBLK),ZERO,NI*NJ)
* 1 : One-electron part 
            DO KGAS = 1, NGAS
              KSM = IJSM
              NK = NOBPTS(KGAS,KSM)
*. blocks of one-electron integrals and one-electron density
              CALL GETD1(dbl_mb(KLDEN),JSM,JGAS,KSM,KGAS,IJSPIN)
              ISPCAS = IJSPIN
              CALL GETH1(dbl_mb(KLINT),ISM,IGAS,KSM,KGAS)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) 
     &          ' 1-e ints for ISM IGAS KGAS ',ISM,IGAS,KGAS
                CALL WRTMAT(dbl_mb(KLINT),NI,NK,NI,NK)
                WRITE(6,*) 
     &          ' 1-e densi for ISM JGAS KGAS ',ISM,JGAS,KGAS
                CALL WRTMAT(dbl_mb(KLDEN),NJ,NK,NJ,NK)
              END IF
*. And then a matrix multiply( they are pretty much in fashion 
*. these days )
              CALL MATML7(dbl_mb(KLFBLK),dbl_mb(KLINT),dbl_mb(KLDEN),
     &                    NI,NJ,NI,NK,NJ,NK,ONE,ONE,2)
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated block '
                 CALL WRTMAT(dbl_mb(KLFBLK),NI,NJ,NI,NJ)
               END IF
 
            END DO
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One-electron contributions'
              WRITE(6,*) ' =========================='
              CALL WRTMAT(dbl_mb(KLFBLK),NI,NJ,NI,NJ)
            END IF
            IF(I12L.EQ.2) THEN
*. 2 : Two-electron part
            DO KLSPIN = 1, NSPIN
            DO KSM = 1, NSMOB
            DO LSM = 1, NSMOB
*. Obtain MSM
              IF (IJSPIN.EQ.1.AND.KLSPIN.EQ.1) ISPC = 1
              IF (IJSPIN.EQ.2.AND.KLSPIN.EQ.2) ISPC = 2
              IF (IJSPIN.EQ.1.AND.KLSPIN.EQ.2) ISPC = 4
              IF (IJSPIN.EQ.2.AND.KLSPIN.EQ.1) ISPC = 3
              CALL  SYMCOM(3,1,KSM,LSM,KLSM)
              CALL  SYMCOM(3,1,KLSM,ISM,IKLSM)
              IMKLSM = 1
              CALL  SYMCOM(2,1,IKLSM,MSM,IMKLSM)
*
              DO MGAS = 1, NGAS
              DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
                NM = NOBPTS(MGAS,MSM)
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
               
*. Blocks of density matrix and integrals : (K L ! I M),D2(K L, J M)
                IXCHNG = 0
                ICOUL  = 1
                ISPCAS = ISPC
                CALL GETH2(dbl_mb(KLINT),
     &               KSM,KGAS,LSM,LGAS,ISM,IGAS,MSM,MGAS,ISPC)
                
                CALL GETD2 (dbl_mb(KLDEN),
     &               KSM,KGAS,LSM,LGAS,JSM,JGAS,MSM,MGAS,ISPC)
                NKL = NK*NL
                DO M = 1, NM
                  IIOFF = KLINT + (M-1)*NKL*NI
                  IDOFF = KLDEN + (M-1)*NKL*NJ
                  CALL MATML7(dbl_mb(KLFBLK),dbl_mb(IIOFF),
     &                  dbl_mb(IDOFF),NI,NJ,NKL,NI,NKL,NJ,ONE,ONE,1)
                END DO
              END DO
              END DO
              END DO
            END DO
            END DO
            END DO ! KLSPIN
            END IF
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One- + two-electron contributions'
              WRITE(6,*) ' ================================='
              WRITE(6,*) ' IGAS, JGAS = ', IGAS, JGAS
              IF (NSPIN.EQ.2.AND.IJSPIN.EQ.1)
     &             WRITE(6,*) '  (alpha spin)'
              IF (NSPIN.EQ.2.AND.IJSPIN.EQ.2)
     &             WRITE(6,*) '  (beta spin)'
              CALL WRTMAT(dbl_mb(KLFBLK),NI,NJ,NI,NJ)
            END IF
*. Block has been constructed , transfer to -complete- 
*. symmetry blocked Fock matrix
            ISPOFF = (IJSPIN-1)*IFBETOFF
            DO J = 1, NJ
              DO I = 1, NI
C?              WRITE(6,*) 'IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1',
C?   &                      IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1
                F(ISPOFF+IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1 ) = 
     &          dbl_mb(KLFBLK-1+(J-1)*NI+I)
              END DO
            END DO
*
          END DO
        END DO
      END DO
      END DO ! IJSPIN
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from FOCK_MAT '
        WRITE(6,*) ' ====================='
        DO ISPIN = 1, NSPIN
          IF (NSPIN.EQ.2.AND.ISPIN.EQ.1)
     &         WRITE(6,*) '  Alpha part:'
          IF (NSPIN.EQ.2.AND.ISPIN.EQ.2)
     &         WRITE(6,*) '  Beta  part:'
          ISPOFF = (ISPIN-1)*IFBETOFF+1          
          CALL APRBLM2(F(ISPOFF),NOCOBS,NOCOBS,NSMOB,0)
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'FOO   ')
      RETURN
      END
      SUBROUTINE FOCK_MAT(F,I12)
*
* Construct Fock matrix
*
* F(i,j) = SUM(K) H(i,K) * RHO1(j,K)
*          + SUM(M,K,L) I  (i M K L ) * RHO2( j M K L )
*
* Helsingfors, december 11 (1996)
* (after the EFG Winter School)
*
* Unless I12 = 2, only one-electron part is calculated
c      IMPLICIT REAL*8(A-H,O-Z)
*. Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION F(*)
*
      NTEST = 000
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'FOO   ')
*
      ONE = 1.0D0
      ZERO = 0.0D0
*. Largest set of orbitals with given symmetry and type
      MXTSOB = 0
      DO ISM = 1, NSMOB
      DO IGAS = 1, NGAS
        MXTSOB = MAX(MXTSOB,NOBPTS(IGAS,ISM))
      END DO
      END DO
C?    WRITE(6,*) 'MXTSOB = ', MXTSOB
*. Allocate scratch space for 2-electron integrals and 
*. two-electron densities
      MX4IBLK = MXTSOB ** 4
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT ') !done
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN ') !done
*. And a block of F
      MX2IBLK = MXTSOB ** 2
      CALL MEMMAN(KLFBLK,MX2IBLK,'ADDL  ',2,'KLFBL ') !done
*. 
*
      ONE = 1.0D0
      DO IJSM = 1, NSMOB
        ISM = IJSM
        JSM = IJSM
        NIJS = NOCOBS(IJSM)
*
        IF(IJSM.EQ.1) THEN
         IFOFF = 1
        ELSE
         IFOFF = IFOFF+NOCOBS(IJSM-1)**2
        END IF
*
        II = -2810
        DO JGAS = 1, NGAS
          IF(JGAS.EQ.1) THEN
            IJ = 1
          ELSE 
            IJ = IJ + NOBPTS(JGAS-1,JSM)
          END IF
          NJ = NOBPTS(JGAS,IJSM)
          DO IGAS = 1, NGAS
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) 
     &        ' Fock matrix for ISM IGAS JGAS',ISM,IGAS,JGAS
            END IF
            NI = NOBPTS(IGAS,ISM)
            IF(IGAS.EQ.1) THEN
              II = 1
            ELSE 
              II = II + NOBPTS(IGAS-1,ISM)
            END IF
*
*  =======================
*. block F(ijsm,igas,jgas)
*  =======================
*
            CALL SETVEC(dbl_mb(KLFBLK),ZERO,NI*NJ)
* 1 : One-electron part 
            DO KGAS = 1, NGAS
              KSM = IJSM
              NK = NOBPTS(KGAS,KSM)
*. blocks of one-electron integrals and one-electron density
              CALL GETD1(dbl_mb(KLDEN),JSM,JGAS,KSM,KGAS,1)
              CALL GETH1(dbl_mb(KLINT),ISM,IGAS,KSM,KGAS)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) 
     &          ' 1-e ints for ISM IGAS KGAS ',ISM,IGAS,KGAS
                CALL WRTMAT(dbl_mb(KLINT),NI,NK,NI,NK)
                WRITE(6,*) 
     &          ' 1-e densi for ISM JGAS KGAS ',ISM,JGAS,KGAS
                CALL WRTMAT(dbl_mb(KLDEN),NJ,NK,NJ,NK)
              END IF
*. And then a matrix multiply( they are pretty much in fashion 
*. these days )
              CALL MATML7(dbl_mb(KLFBLK),dbl_mb(KLINT),dbl_mb(KLDEN),
     &                    NI,NJ,NI,NK,NJ,NK,ONE,ONE,2)
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated block '
                 CALL WRTMAT(dbl_mb(KLFBLK),NI,NJ,NI,NJ)
               END IF
 
            END DO
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One-electron contributions'
              WRITE(6,*) ' =========================='
              CALL WRTMAT(dbl_mb(KLFBLK),NI,NJ,NI,NJ)
            END IF
            IF(I12.EQ.2) THEN
*. 2 : Two-electron part
            DO KSM = 1, NSMOB
            DO LSM = 1, NSMOB
*. Obtain MSM
              CALL  SYMCOM(3,1,KSM,LSM,KLSM)
              CALL  SYMCOM(3,1,KLSM,ISM,IKLSM)
              IMKLSM = 1
              CALL  SYMCOM(2,1,IKLSM,MSM,IMKLSM)
*
              DO MGAS = 1, NGAS
              DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
                NM = NOBPTS(MGAS,MSM)
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
               
*. Blocks of density matrix and integrals : (K L ! I M),D2(K L, J M)
                IXCHNG = 0
                ICOUL  = 1
                ONE = 1.0D0
                CALL LUCIAGETINT(dbl_mb(KLINT),
     &               KGAS,KSM,LGAS,LSM,IGAS,ISM,MGAS,MSM,
     &               IXCHNG,0,0,ICOUL,ONE,ONE)
                
                CALL GETD2 (dbl_mb(KLDEN),
     &               KSM,KGAS,LSM,LGAS,JSM,JGAS,MSM,MGAS,1)
C               GETINT(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,
C    &                     LTYP,LSM,IXCHNG,0,0,ICOUL)
                NKL = NK*NL
                DO M = 1, NM
                  IIOFF = KLINT + (M-1)*NKL*NI
                  IDOFF = KLDEN + (M-1)*NKL*NJ
                  CALL MATML7(dbl_mb(KLFBLK),dbl_mb(IIOFF),
     &                   dbl_mb(IDOFF),NI,NJ,NKL,NI,NKL,NJ,ONE,ONE,1)
                END DO
              END DO
              END DO
              END DO
            END DO
            END DO
            END IF
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One- + two-electron contributions'
              WRITE(6,*) ' ================================='
              CALL WRTMAT(dbl_mb(KLFBLK),NI,NJ,NI,NJ)
            END IF
*. Block has been constructed , transfer to -complete- 
*. symmetry blocked Fock matrix
            DO J = 1, NJ
              DO I = 1, NI
C?              WRITE(6,*) 'IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1',
C?   &                      IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1
                F(IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1 ) = 
     &          dbl_mb(KLFBLK-1+(J-1)*NI+I)
              END DO
            END DO
*
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from FOCK_MAT '
        WRITE(6,*) ' ====================='
        CALL APRBLM2(F,NOCOBS,NOCOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'FOO   ')
      RETURN
      END
      SUBROUTINE EN_FROM_DENS(E,I12,IMODE)
*
* Obtain contributions from active orbitals to 
* energy from densities and integrals - and add core-energy
*
*
* E = SUM(i,j) H(i,j) * RHO1(i,j)
*          + 0.5*SUM(i,j,K,L) (I J K L ) * RHO2( I J K L )
*
*     if IMODE.NE.0
*
* E = SUM(i,j) H(i,j) * RHO1(i,j)
*          + 0.5*SUM(i,j,K,L) (I J K L ) * 
*                (LAM2( I J K L ) + RHO1(I,J)RHO1(K,L)-RHO1(I,L)RHO1(K,J))
*
* Jeppe Olsen, Early 1997
*              Sept. 98    : I12 added
* Andreas      Dec. 2004   : IMODE added
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*. Input
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
*
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      call luciaqenter('enden')
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Energy under construction '
        WRITE(6,*) ' =========================='
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'EN_FRM')
*
      E1 = 0.0D0
      E2 = 0.0D0
*. Largest set of orbitals with given symmetry and type
      MXTSOB = 0
      DO ISM = 1, NSMOB
      DO IGAS = 1, NGAS
        MXTSOB = MAX(MXTSOB,NOBPTS(IGAS,ISM))
      END DO
      END DO
*. Allocate scratch space for 2-electron integrals and 
*. two-electron densities
      MX4IBLK = MXTSOB ** 4
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT') !done
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN') !done
      ONE = 1.0D0
      DO ISM = 1, NSMOB
       DO JSM = 1, NSMOB
        CALL  SYMCOM(3,1,ISM,JSM,IJSM)
        DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
           NI = NOBPTS(IGAS,ISM)
           NJ = NOBPTS(JGAS,JSM)
           II = IOBPTS(IGAS,ISM)
           IJ = IOBPTS(JGAS,JSM)
           IF(ISM.EQ.JSM) THEN
*
* One-electron part 
* =================
*
*. blocks of one-electron integrals and one-electron density
         
             CALL GETD1(dbl_mb(KLDEN),ISM,IGAS,ISM,JGAS,1)
             CALL GETH1(dbl_mb(KLINT),ISM,IGAS,ISM,JGAS)
             IF(NTEST.GE.100) THEN
               WRITE(6,*) ' Block of 1e integrals ISM,IGAS,JGAS',
     &                    ISM,IGAS,JGAS
               CALL WRTMAT(dbl_mb(KLINT),NI,NJ,NI,NJ)
               WRITE(6,*) ' Block of 1e density ISM,IGAS,JGAS',
     &                    ISM,IGAS,JGAS
               CALL WRTMAT(dbl_mb(KLDEN),NI,NJ,NI,NJ)
             END IF
             E1 = E1 + INPROD(dbl_mb(KLDEN),dbl_mb(KLINT),NI*NJ)
           END IF
*
* Two-electron part 
* =================
*
           IF(I12.EQ.2) THEN
           DO KSM = 1, NSMOB
*. Obtain LSM
             CALL  SYMCOM(3,1,IJSM,KSM,IJKSM)
             IJKLSM = 1
             CALL  SYMCOM(2,1,IJKSM,LSM,IJKLSM)
C?           WRITE(6,*) ' IJSM IJKSM LSM ',IJSM,IJKSM,IJKLSM
*
             DO KGAS = 1, NGAS
             DO LGAS = 1, NGAS
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
*. Blocks of density matrix and integrals 
                IXCHNG = 0
                ICOUL  = 1
                ONE = 1.0D0
                CALL LUCIAGETINT(dbl_mb(KLINT),
     &               IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
     &               IXCHNG,0,0,ICOUL,ONE,ONE)
                CALL GETD2 (dbl_mb(KLDEN),
     &               ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,1)
                IF (IMODE.NE.0) THEN
                  CALL GETD2RED(dbl_mb(KLDEN),
     &                 ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,1)
                END IF
C?              write(6,*) ' Ism Jsm Ksm Lsm' , Ism,Jsm,Ksm,Lsm
C?              write(6,*)
C?   &          ' Igas Jgas Kgas Lgas' , Igas,Jgas,Kgas,Lgas
C?              WRITE(6,*) ' Integral block'
C?              CALL WRTMAT(WORK(KLINT),NI*NJ,NK*NL,NI*NJ,NK*NL)
C?              WRITE(6,*) ' Density block '
C?              CALL WRTMAT(WORK(KLDEN),NI*NJ,NK*NL,NI*NJ,NK*NL)
                NIJKL = NI*NJ*NK*NL
                E2 = E2+0.5D0*INPROD(dbl_mb(KLDEN),dbl_mb(KLINT),NIJKL)
C?              write(6,*) ' Updated 2e-energy ', E2
             END DO
             END DO
           END DO
           END IF
*
          END DO
         END DO
       END DO
      END DO
*
      E = E1 + E2 + ECORE
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from EN_FROM_DEN' 
        WRITE(6,*)
        WRITE(6,*) ' One-electron energy ', E1
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-electron energy ', E2
        END IF
        WRITE(6,*)
        WRITE(6,*) ' Total energy : ', E
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'EN_FRM')
      call luciaqexit('enden')
*
      RETURN
      END
      SUBROUTINE GETD1(RHO1B,ISM,IGAS,JSM,JGAS,ISPIN)
*
* Extract TS block of one-elecvtron density matrix
*
c      IMPLICIT REAL*8(A,H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. output
      DIMENSION RHO1B(*)
*
      IF(IGAS.EQ.0.OR.IGAS.GT.NGAS.OR.
     &   JGAS.EQ.0.OR.JGAS.GT.NGAS    ) THEN
        WRITE(6,*) ' GETD1 called for inactive or secondary space'
        WRITE(6,'(A,2I4)') ' IGAS, JGAS = ', IGAS, JGAS
        STOP       ' GETD1 called for inactive or secondary space'
      END IF
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
*
      II = IOBPTS_AC(IGAS,ISM)
      IJ = IOBPTS_AC(JGAS,JSM)
    
*
      IOFF = (ISPIN-1)*NACOB**2
*
      DO I = 1, NI
        DO J = 1, NJ
          IABS = I-1+II
          JABS = J-1+IJ
          IJABS = (JABS-1)*NACOB + IABS
          RHO1B((J-1)*NI+I) = dbl_mb(KRHO1-1+IOFF+IJABS)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Block of one-electron density matrix'
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
        WRITE(6,*) 'IGAS,ISM,JGAS,JSM',IGAS,ISM,JGAS,JSM
        CALL WRTMAT(RHO1B,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END 
      FUNCTION GETD1E(IORB,ITP,ISM,JORB,JTP,JSM)
*
* One-electron density matrix for active
* orbitals (IORB,ITP,ISM),(JORB,JTP,JSM)
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'oper.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      ISPIN = ISPCAS
      IF (ISPCAS.EQ.0) ISPIN=1
      IF (ISPCAS.GT.2) THEN
        STOP 'GETD1E'
      END IF
      IOFF = (ISPIN-1)*NTOOB**2
*
      II = IOBPTS(ITP,ISM)
      JJ = IOBPTS(JTP,JSM)

      IABS = IORB-1+II
      JABS = JORB-1+JJ
      IJABS = (JABS-1)*NTOOB + IABS
      GETD1E = dbl_mb(KRHO1-1+IOFF+IJABS)

      RETURN
      END
      SUBROUTINE GETD2(RHO2B,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,
     &                 ISPC)
*. Extract given TS block from the 2e-density matrix
*.
*. Jeppe Olsen, Some day in Hfors CITY, winter 1996
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. Initial implementation, KISS to the MAX !!!
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc' 
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION RHO2B(*)
*. output is in the form i,j,k,l
*
      IF(IGAS.EQ.0.OR.IGAS.GT.NGAS.OR.
     &   JGAS.EQ.0.OR.JGAS.GT.NGAS.OR.
     &   KGAS.EQ.0.OR.KGAS.GT.NGAS.OR.
     &   LGAS.EQ.0.OR.LGAS.GT.NGAS    ) THEN
        WRITE(6,*) ' GETD2 called for inactive or secondary space'
        WRITE(6,'(A,4I4)') ' IGAS, JGAS,KGAS, LGAS = ', 
     &                       IGAS, JGAS,KGAS, LGAS
        STOP       ' GETD2 called for inactive or secondary space'
      END IF
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
*
      IELMNT = 0
      DO L = 1, NL
        DO K = 1, NK
          DO J = 1, NJ
            DO I = 1, NI
              IELMNT = IELMNT + 1
              RHO2B(IELMNT) = GETD2E(I,IGAS,ISM,J,JGAS,JSM,
     &                               K,KGAS,KSM,L,LGAS,LSM,ISPC)
            END DO
          END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Block of two-electron density' 
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Type and symmetry of the 4 orbitals (i j k l )'
        WRITE(6,'(1H ,8(1X,I4) )')
     &  IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
*
        NIJ = NI*NJ
        NKL = NK*NL
        CALL WRTMAT(RHO2B,NIJ,NKL,NIJ,NKL)
      END IF 
*
      RETURN
      END
      SUBROUTINE GETD2RED(RHO2B,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,
     &                 ISPC)
*. Update given TS block of the irreducible 2e-density matrix
*. with the reducible part from the 1e-density
*.
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION RHO2B(*)
*. output is in the form i,j,k,l
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
*
      IF(ISPC.EQ.0) FAC=0.5D0  ! closed shell case
      IF(ISPC.NE.0) FAC=0.5D0  ! spin-orbital case
*
      IELMNT = 0
      DO L = 1, NL
        DO K = 1, NK
          DO J = 1, NJ
            DO I = 1, NI
              IELMNT = IELMNT + 1
              RHO2B(IELMNT) = RHO2B(IELMNT) +
     &             GETD1E(I,IGAS,ISM,J,JGAS,JSM)*
     &                             GETD1E(K,KGAS,KSM,L,LGAS,LSM,ISPC)
     &             -FAC*
     &             GETD1E(I,IGAS,ISM,L,LGAS,LSM)*
     &                             GETD1E(K,KGAS,KSM,J,JGAS,JSM,ISPC)
            END DO
          END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' BLock of two-electron density' 
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Type and symmetry of the 4 orbitals (i j k l )'
        WRITE(6,'(1H ,8(1X,I4) )')
     &  IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
*
        NIJ = NI*NJ
        NKL = NK*NL
        CALL WRTMAT(RHO2B,NIJ,NKL,NIJ,NKL)
      END IF 
*
      RETURN
      END
      SUBROUTINE GETD2_A(RHO2B,
     &     IORB,ISM,IGAS,JORB,JSM,JGAS,
     &     KORB,KSM,KGAS,LORB,LSM,LGAS,
     &     ISPC)
*. Extract given TS block from the 2e-density matrix
*.
*. Jeppe Olsen, Some day in Hfors CITY, winter 1996
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. Initial implementation, KISS to the MAX !!!
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION RHO2B(*)
*. output is in the form i,j,k,l
*
      IF (IORB.EQ.0) THEN
        NI = NOBPTS(IGAS,ISM)
        II = 1
      ELSE
        NI = MIN(1,NOBPTS(IGAS,ISM))
        II = IORB
      END IF
      IF (JORB.EQ.0) THEN
        NJ = NOBPTS(JGAS,JSM)
        JJ = 1
      ELSE
        NJ = MIN(1,NOBPTS(JGAS,JSM))
        JJ = JORB
      END IF
      IF (KORB.EQ.0) THEN
        NK = NOBPTS(KGAS,KSM)
        KK = 1
      ELSE
        NK = MIN(1,NOBPTS(KGAS,KSM))
        KK = KORB
      END IF
      IF (LORB.EQ.0) THEN
        NL = NOBPTS(LGAS,LSM)
        LL = 1
      ELSE
        NL = MIN(1,NOBPTS(LGAS,LSM))
        LL = LORB
      END IF
*
      IELMNT = 0
      DO L = LL, (LL-1)+NL
        DO K = KK, (KK-1)+NK
          DO J = JJ, (JJ-1)+NJ
            DO I = II, (II-1)+NI
              IELMNT = IELMNT + 1
              RHO2B(IELMNT) = GETD2E(I,IGAS,ISM,J,JGAS,JSM,
     &                               K,KGAS,KSM,L,LGAS,LSM,ISPC)
            END DO
          END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' BLock of two-electron density' 
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Type and symmetry of the 4 orbitals (i j k l )'
        WRITE(6,'(1H ,8(1X,I4) )')
     &  IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
*
        NIJ = NI*NJ
        NKL = NK*NL
        CALL WRTMAT(RHO2B,NIJ,NKL,NIJ,NKL)
      END IF 
*
      RETURN
      END
      FUNCTION GETD2E(I,IGAS,ISM,J,JGAS,JSM,K,KGAS,KSM,L,LGAS,LSM,ISPC)
*
* Obtain element of two-electron density matrix
* Currently stored without symmetry
*
*. 2-electron density is assumed stored in work(krho2)
* modified for density only over active orbitals, July 2010
*      
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Compiler warnings..
      IJKL = -2303
*
      IOFF = 0
      IF (ISPC.GE.2) IOFF = IOFF + NACOB**2 * (NACOB**2+1) / 2
      IF (ISPC.GE.3) IOFF = IOFF + NACOB**2 * (NACOB**2+1) / 2
*
      IABS = I + IOBPTS_AC(IGAS,ISM)-1
      JABS = J + IOBPTS_AC(JGAS,JSM)-1
      KABS = K + IOBPTS_AC(KGAS,KSM)-1
      LABS = L + IOBPTS_AC(LGAS,LSM)-1
*
      IJ = (JABS-1)*NACOB+IABS
      KL = (LABS-1)*NACOB+KABS
      IF (ISPC.LT.3) THEN
        IF(IJ.GE.KL) THEN
          IJKL = IJ*(IJ-1)/2+KL
        ELSE
          IJKL = KL*(KL-1)/2+IJ
        END IF
      ELSE IF (ISPC.EQ.3) THEN
        IJKL = (KL-1)*NACOB*NACOB+IJ
      ELSE IF (ISPC.EQ.4) THEN
        IJKL = (IJ-1)*NACOB*NACOB+KL        
      END IF
*
      X = dbl_mb(KRHO2-1+IOFF+IJKL)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Element of two electron density matrix'
        WRITE(6,*) 'I,IGAS,ISM,J,JGAS,JSM,K,KGAS,KSM,L,LGAS,LSM'
        WRITE(6,'( 1H ,(12I4) )')
     &  I,IGAS,ISM,J,JGAS,JSM,K,KGAS,KSM,L,LGAS,LSM
        WRITE(6,*) ' IJ, KL ', IJ,KL
        WRITE(6,*) ' IABS JABS KABS LABS',IABS,JABS,KABS,LABS
        WRITE(6,*) 'Address and value', IJKL, X
      END IF
*
      GETD2E = X
*
      RETURN
      END  
      SUBROUTINE GETH2(XINT,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,
     &                 ISPC)
*. Extract given TS block from the 2e-integral list
*. KISS variant of GETINT
*. AK from JO's GETD2 --- 2004
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'oper.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION XINT(*)
*. output is in the form i,j,k,l
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
*
      ISPCAS = ISPC
      IELMNT = 0
      DO L = 1, NL
        DO K = 1, NK
          DO J = 1, NJ
            DO I = 1, NI
              IELMNT = IELMNT + 1
              XINT(IELMNT) = GTIJKL(IOBPTS(IGAS,ISM)+I-1,
     &                              IOBPTS(JGAS,JSM)+J-1,
     &                              IOBPTS(KGAS,KSM)+K-1,
     &                              IOBPTS(LGAS,LSM)+L-1)
            END DO
          END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' BLock of two-electron integrals' 
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        WRITE(6,*) ' Type and symmetry of the 4 orbitals (i j k l )'
        WRITE(6,'(1H ,8(1X,I4) )')
     &  IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
*
        NIJ = NI*NJ
        NKL = NK*NL
        CALL WRTMAT(XINT,NIJ,NKL,NIJ,NKL)
      END IF 
*
      RETURN
      END
      SUBROUTINE GETH2_A(XINT,
     &     IORB,ISM,IGAS,JORB,JSM,JGAS,
     &     KORB,KSM,KGAS,LORB,LSM,LGAS,
     &     ISPC)
*. Extract given TS block from the 2e-integral list
*. KISS variant of GETINT
*. AK from JO's GETD2 --- 2004
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'oper.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION XINT(*)
*. output is in the form i,j,k,l
*
      IF (IORB.EQ.0) THEN
        NI = NOBPTS(IGAS,ISM)
        II = 1
      ELSE
        NI = MIN(1,NOBPTS(IGAS,ISM))
        II = IORB
      END IF
      IF (JORB.EQ.0) THEN
        NJ = NOBPTS(JGAS,JSM)
        JJ = 1
      ELSE
        NJ = MIN(1,NOBPTS(JGAS,JSM))
        JJ = JORB
      END IF
      IF (KORB.EQ.0) THEN
        NK = NOBPTS(KGAS,KSM)
        KK = 1
      ELSE
        NK = MIN(1,NOBPTS(KGAS,KSM))
        KK = KORB
      END IF
      IF (LORB.EQ.0) THEN
        NL = NOBPTS(LGAS,LSM)
        LL = 1
      ELSE
        NL = MIN(1,NOBPTS(LGAS,LSM))
        LL = LORB
      END IF
*
      ISPCAS = ISPC
      IELMNT = 0
      DO L = LL, (LL-1)+NL
        DO K = KK, (KK-1)+NK
          DO J = JJ, (JJ-1)+NJ
            DO I = II, (II-1)+NI
              IELMNT = IELMNT + 1
              XINT(IELMNT) = GTIJKL(IOBPTS(IGAS,ISM)+I-1,
     &                              IOBPTS(JGAS,JSM)+J-1,
     &                              IOBPTS(KGAS,KSM)+K-1,
     &                              IOBPTS(LGAS,LSM)+L-1)
            END DO
          END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' BLock of two-electron integrals' 
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        WRITE(6,*) ' Type and symmetry of the 4 orbitals (i j k l )'
        WRITE(6,'(1H ,8(1X,I4) )')
     &  IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
*
        NIJ = NI*NJ
        NKL = NK*NL
        CALL WRTMAT(XINT,NIJ,NKL,NIJ,NKL)
      END IF 
*
      RETURN
      END
      SUBROUTINE SPSPCLS(ISPSPCLS,ICLS,NCLS)
*
* Obtain mapping a-supergroup X b-supergroup => class
*
* Classes are specified by ICLS
*
* Jeppe Olsen, Jan 97 
*
*. Modified Oct. 2004 to improve efficiency for very many types
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      INTEGER ICLS(*)
*. OUtput
      INTEGER ISPSPCLS(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SPSPCL')
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
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
*. construct arrays ordering the alphastring types according 
*. to accumulated occupations after two orbital spaces
*
      CALL MEMMAN(KLNSTRFSPLT ,(NAEL+1)*(NAEL+1),'ADDL  ',1,'NSFSPL')!d
      CALL MEMMAN(KLIBSTRFSPLT,(NAEL+1)*(NAEL+1),'ADDL  ',1,'BSFSPL')!d
      CALL MEMMAN(KLISTRFSPLT ,NOCTPA   ,'ADDL  ',1,'ISFSPL')!d
      CALL GROUP_STRTP(NELFSPGP(1,IOCTPA),NOCTPA,NGAS,NAEL,ISPLIT1,
     &     ISPLIT2,int_mb(KLNSTRFSPLT),int_mb(KLIBSTRFSPLT),
     &     int_mb(KLISTRFSPLT),MXPNGAS)
C     GROUP_STRTP(ISTRTP,NSTRTP,NGAS,NEL,ISPLIT1,ISPLIT2,
C    &           NSTRFSPLT,IBSTRFSPLT,ISTRFSPLT)
      CALL SPSPCLS_GAS(NOCTPA,NOCTPB,
     &            ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,NGAS,ISPSPCLS,ICLS,NCLS,IPRDIA,ISPLIT1,
     &            ISPLIT2,NAEL,int_mb(KLNSTRFSPLT),int_mb(KLIBSTRFSPLT),
     &            int_mb(KLISTRFSPLT))
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SPSPCL')
      RETURN
      END
*
      SUBROUTINE SPSPCLS_GAS(NOCTPA,NOCTPB,IOCA,IOCB,NELFGP,
     &           NGAS,ISPSPCLS,ICLS,NCLS,IPRNT,
     &           ISPLIT1,ISPLIT2,NAEL,NSTRFSPLT,IBSTRFSPLT,
     &           ISTRFSPLT)
*
* Obtain mapping a-supergroup X b-supergroup => class
*
*. Jeppe Olsen, modified by including info that speeds process
*. up for large expansions, oct. 2004
*
* =====
*.Input
* =====
*
* NOCTPA : Number of alpha types
* NOCTPB : Number of beta types
*
* IOCA(IGAS,ISTR) occupation of AS IGAS for alpha string type ISTR
* IOCB(IGAS,ISTR) occupation of AS IGAS for beta  string type ISTR
*
* MXPNGAS : Largest allowed number of gas spaces
* NGAS    : Actual number of gas spaces

*
* ======
*.Output
* ======
*
* ISPSPCLS(IATP,IBTP) => Class of this block of determinants
*                        =0 indicates unallowed(class less) combination
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*.Input
      INTEGER IOCA(MXPNGAS,NOCTPA),IOCB(MXPNGAS,NOCTPB)
      INTEGER NELFGP(*)
      INTEGER ICLS(NGAS,NCLS)
      INTEGER NSTRFSPLT(NAEL+1,NAEL+1), IBSTRFSPLT(NAEL+1,NAEL+1)
      INTEGER ISTRFSPLT(*)
*.Output
      INTEGER ISPSPCLS(NOCTPA,NOCTPB)
*. Local scratch 
       INTEGER IIA(MXPNGAS)

*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ISPSPCLS_GAS entered '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' IOCA and IOCB (groups)'
        CALL IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NGAS)
        WRITE(6,*)
        CALL IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NGAS)
        WRITE(6,*) 
        WRITE(6,*) ' ICLS '
        CALL IWRTMA(ICLS,NGAS,NCLS,NGAS,NCLS)
      END IF
*
      I_NEW_OR_OLD = 1
*
      IF(I_NEW_OR_OLD.EQ.2) THEN
*. Old way, scales cubically
        DO 100 IATP = 1, NOCTPA
          DO 90 IBTP = 1, NOCTPB
            IICLS = 0
            DO KCLS = 1, NCLS
              IAMOKAY = 1
              DO IGAS = 1, NGAS
                IEL = NELFGP(IOCA(IGAS,IATP))+NELFGP(IOCB(IGAS,IBTP))
                IF(IEL.NE.ICLS(IGAS,KCLS)) IAMOKAY = 0
              END DO
              IF(IAMOKAY.EQ.1) IICLS=KCLS
            END DO
            ISPSPCLS(IATP,IBTP) = IICLS
   90     CONTINUE
  100   CONTINUE
      ELSE 
*. New reduced way, scales less than cubically
        IZERO = 0
        CALL ISETVC(ISPSPCLS,IZERO,NOCTPA*NOCTPB)
        DO IBTP = 1, NOCTPB
          DO KCLS = 1, NCLS
*. Find the required occupation of A
            DO JGAS = 1, NGAS
              IIA(JGAS) = ICLS(JGAS,KCLS)-NELFGP(IOCB(JGAS,IBTP))
            END DO
C?          WRITE(6,*) ' IBTP, KCLS = ', IBTP, KCLS
C?          WRITE(6,*) ' IIA '
C?          CALL IWRTMA(IIA,1,NGAS,1,NGAS)
*. Check that no occs are negative 
            INEG = 0
            DO IGAS = 1, NGAS
              IF(IIA(IGAS).LT.0) INEG = 1
            END DO
            IF(INEG.EQ.0) THEN
              IEL1 = IELSUM(IIA,ISPLIT1)
              IEL2 = IELSUM(IIA,ISPLIT2)
C?            WRITE(6,*) ' IEL1, IEL2 = ', IEL1, IEL2
              IB_A = IBSTRFSPLT(IEL1+1,IEL2+1)
              N_A  = NSTRFSPLT (IEL1+1,IEL2+1)
C?            WRITE(6,*) ' IB_A, N_A ', IB_A, N_A
              DO IIATP = IB_A, IB_A+N_A-1
                IATP = ISTRFSPLT(IIATP)
                IAMOKAY = 1
                DO JGAS = 1, NGAS
                  IF(NELFGP(IOCA(JGAS,IATP)).NE.IIA(JGAS))IAMOKAY = 0
                END DO
                IF(IAMOKAY.EQ.1) ISPSPCLS(IATP,IBTP) = KCLS
              END DO
CMOVED        ISPSPCLS(IATP,IBTP) = IICLS
COLD        ELSE
COLD          ISPSPCLS(IATP,IBTP) = 0    
            END IF
          END DO
        END DO
      END IF
*
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Matrix giving classes for alpha-beta supergroups'
        WRITE(6,*)
        CALL IWRTMA(ISPSPCLS,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
C     BLKCLS(WORK(KLCIBT),NBLOCKS,WORK(KLBLKCLS),WORK(KLSPSPCL),
C    &            NOCTPA,NOCTPB)
      SUBROUTINE BLKCLS(IBLKS,NBLKS,IBLKCLS,ISPSPCL,
     &                  NCLS,LCLS,NOCTPA,NOCTPB,RLCLS)
*
* Class of each block, and dimension of each class
*
* Jeppe Olsen
*
*. Last modification; Nov. 2, 2012; Jeppe Olsen; Test output
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLKS(8,NBLKS)
      INTEGER ISPSPCL(NOCTPA,NOCTPB)
*. Output
      INTEGER IBLKCLS(NBLKS),LCLS(NCLS)
      DIMENSION RLCLS(NCLS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from BLKCLS '
       WRITE(6,*) ' ================='
       WRITE(6,*)
       WRITE(6,*) ' ISPSPCL'
       CALL IWRTMA(ISPSPCL,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
       WRITE(6,*) ' NCLS =', NCLS
      END IF
*  
      IZERO = 0
      CALL ISETVC(LCLS,IZERO,NCLS)
      ZERO = 0.0D0
      CALL SETVEC(RLCLS,ZERO,NCLS)
      DO JBLK = 1, NBLKS
        IICLS = ISPSPCL(IBLKS(1,JBLK),IBLKS(2,JBLK))
        IBLKCLS(JBLK) = IICLS
        LCLS(IICLS) = LCLS(IICLS) + IBLKS(8,JBLK)
        RLCLS(IICLS) = RLCLS(IICLS) + FLOAT(IBLKS(8,JBLK))
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' BLKCLS Speaking '
        WRITE(6,*) ' ==============='
        WRITE(6,*)
        WRITE(6,*) ' Dimension of each class (as integer)'
        CALL IWRTMA(LCLS,1,NCLS,1,NCLS)
        WRITE(6,*) ' Dimension of each class (as Real)'
        CALL WRTMAT(RLCLS,1,NCLS,1,NCLS)
        WRITE(6,*)
        WRITE(6,*) ' Class of each block : '
        CALL IWRTMA(IBLKCLS,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END
      SUBROUTINE DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,
     &                     NOBTP,IL,IR)
*
* Obtain types of I,J,K,l so
* <L!a+I a+K a L a J!R> is nonvanishing
* only combinations with type(I) .ge. type(K) and type(J).ge.type(L)
* are included
*
      INTEGER IL(NOBTP),IR(NOBTP)
      INTEGER ITP(*),JTP(*),KTP(*),LTP(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DXTYP_GAS in action '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' Occupation of left string '
        CALL IWRTMA(IL,1,NOBTP,1,NOBTP)
        WRITE(6,*) ' Occupation of right string '
        CALL IWRTMA(IR,1,NOBTP,1,NOBTP)
      END IF
*
*. Number of differing occupations
      NANNI = 0
      NCREA = 0
      NDIFT = 0
*
      ICREA1 = 0
      ICREA2 = 0
      IANNI1 = 0
      IANNI2 = 0
      DO IOBTP = 1, NOBTP
        NDIFT = NDIFT + ABS(IL(IOBTP)-IR(IOBTP))
        NDIF = IL(IOBTP)-IR(IOBTP)
        IF(NDIF.EQ.2) THEN
*. two electrons of type IOBTP must be created
          ICREA1 = IOBTP
          ICREA2 = IOBTP
          NCREA = NCREA + 2
        ELSE IF (NDIF .EQ. -2 ) THEN
*. Two electrons of type IOBTP must be annihilated
          IANNI1 = IOBTP
          IANNI2 = IOBTP
          NANNI = NANNI + 2
        ELSE IF (NDIF.EQ.1) THEN
*. one electron of type IOBTP must be created
          IF(NCREA.EQ.0) THEN
            ICREA1 = IOBTP
          ELSE
            ICREA2 = IOBTP
          END IF
          NCREA = NCREA + 1
        ELSE IF (NDIF.EQ.-1) THEN
* One electron of type IOBTP must be annihilated
          IF(NANNI.EQ.0) THEN
            IANNI1 = IOBTP
          ELSE
            IANNI2 = IOBTP
          END IF
          NANNI = NANNI + 1
        END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NCREA, NANNI ', NCREA, NANNI
        WRITE(6,*)  ' ICREA2, IANNI2 ', ICREA2,IANNI2
        WRITE(6,*)  ' ICREA11,IANNI11 ', ICREA11,IANNI11
        WRITE(6,*)  ' ICREA21,IANNI21 ', ICREA21,IANNI21
      END IF
*
      NDXTP = 0
      IF(NDIFT.GT.4) THEN
        NDXTP = 0
      ELSE
      IF(NCREA.EQ.0.AND.NANNI.EQ.0) THEN
*. strings identical, include diagonal excitions  itp = jtp, ktp=ltp 
        DO IJTP = 1, NOBTP
          IF(IR(IJTP).GE.1) THEN
            DO KLTP = 1, IJTP 
              IF((IJTP.NE.KLTP.AND.IR(KLTP).GE.1).OR.
     &           (IJTP.EQ.KLTP.AND.IR(KLTP).GE.2)) THEN
                 NDXTP = NDXTP + 1
                 ITP(NDXTP) = IJTP
                 JTP(NDXTP) = IJTP
                 KTP(NDXTP) = KLTP
                 LTP(NDXTP) = KLTP
              END IF
            END DO
          END IF
        END DO
*. Strings differ by single excitation
      ELSE IF( NCREA.EQ.1.AND.NANNI.EQ.1) THEN
*. diagonal excitation plus creation in ICREA1 and annihilation in IANNI1
        DO IDIA = 1, NOBTP
          IF((IDIA.NE.IANNI1.AND.IR(IDIA).GE.1).OR.
     &       (IDIA.EQ.IANNI1.AND.IR(IDIA).GE.2)) THEN
             NDXTP = NDXTP + 1
             ITP(NDXTP) = MAX(ICREA1,IDIA)
             KTP(NDXTP) = MIN(ICREA1,IDIA)
             JTP(NDXTP) = MAX(IANNI1,IDIA)
             LTP(NDXTP) = MIN(IANNI1,IDIA)
          END IF
        END DO
      ELSE IF(NCREA.EQ.2.AND.NANNI.EQ.2) THEN
*. Strings differ by double excitation
        NDXTP = 1
        ITP(1) = ICREA2
        KTP(1) = ICREA1
        JTP(1) = IANNI2
        LTP(1) = IANNI1
      END IF
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,I4)')
     &  ' Number of connecting double excitations ', NDXTP
        IF(NDXTP.NE.0) THEN
          WRITE(6,*) '  ITYP KTYP LTYP JTYP '
          WRITE(6,*) '  ===================='
          DO  IDX = 1,NDXTP
            WRITE(6,'(1H ,4I5)')ITP(IDX),KTP(IDX),LTP(IDX),JTP(IDX)
          END DO
        END IF
      END IF
*
      RETURN
      END
      FUNCTION IBASSPC_FOR_CLS(ICLS)
* 
*. Obtain base space for occupation class ICLS
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER ICLS(NGAS)
*
      IBASE = 0
      NEL = -2810
      DO ISPC = 1, NCMBSPC
        DO JJSPC = 1, LCMBSPC(ISPC)
          JSPC = ICMBSPC(JJSPC,ISPC)
*. Test for occupation constraints in CI space JSPC
          I_AM_OKAY = 1
          DO IGAS = 1, NGAS
            IF(IGAS.EQ.1) THEN
              NEL = ICLS(IGAS)
            ELSE
              NEL = NEL + ICLS(IGAS)
            END IF
*
            IF(NEL.LT.IGSOCCX(IGAS,1,JSPC).OR.
     &         NEL.GT.IGSOCCX(IGAS,2,JSPC)    ) THEN
                I_AM_OKAY = 0
            END IF
          END DO
*         ^ End of loop over gasspaces for given cispace
*
*. Perhaps occupation constraints in ensemble gaspace
          IM_IN = 1
          IF(I_CHECK_ENSGS.EQ.1) 
     &    CALL  CHECK_IS_OCC_IN_ENGSOCC(ICLS,JSPC,IM_IN)
C               CHECK_IS_OCC_IN_ENGSOCC(IGSOCCL,ISPC,IM_IN)
          IF(IM_IN.EQ.1.AND.I_AM_OKAY.EQ.1.AND.IBASE.EQ.0) THEN
            IBASE = ISPC
          END IF
*
        END DO
*       ^ End of loop over cisspaces for given combination space
      END DO
*     ^ End of loop over combinations apaces
*
      IBASSPC_FOR_CLS = IBASE    
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation class and its basespace '
        CALL IWRTMA(ICLS,1,NGAS,1,NGAS)
        WRITE(6,*) IBASE
      END IF 
*
      RETURN
      END 
CTOBE       SUBROUTINE CLS_TO_SPSP 
CTOBE      &(ICLS,NCLS,I_CLS_TO_SPSP,N_CLS_TO_SPSP,IB_CLS_TO_SPSP,
CTOBE      & ISPSPCLS,NOCTPA,NOCTPB)
CTOBE * 
CTOBE * Combination of supergroups belonging to given class
CTOBE *
CTOBE * Jeppe Olsen, Nov 99 
CTOBE *
CTOBE * Input :
CTOBE *
CTOBE * ISPSPCLS : Class for given pair of supergroups
CTOBE *
CTOBE *. Output 
CTOBE *
CTOBE * N_CLS_TO_SPSP  : Number of supergroup combinations for class
CTOBE * IB_CLS_TO_SPSP : Base  for supergroup combinations of given class 
CTOBE * I_CLS_TO_SPSP  : supergroup combinations of given class
CTOBE *
CTOBE       INCLUDE 'implicit.inc'
CTOBE       INCLUDE 'mxpdim.inc'
CTOBE       INCLUDE 'gasstr.inc'
CTOBE       INCLUDE 'cgas.inc'
CTOBE       INCLUDE 'wrkspc.inc'
CTOBE *. Input 
CTOBE       INTEGER ISPSPCLS(NOCTPA,NOCTPB)
CTOBE *. Output
CTOBE       INTEGER I_CLS_TO_SPSP(2,*),N_CLS_TO_SPSP(NCLS)
CTOBE       INTEGER IB_CLS_TO_SPSP(NCLS)
CTOBE *
CTOBE       IZERO = 0
CTOBE       CALL ISETVC(N_CLS_TO_SPSP,IZERO,NCLS)
*
CTOBE       DO IOCTPA = 1, NOCTPA                  
CTOBE        DO IOCTPB = 1, NOCTPB                    
CTOBE          ICLS = ISPSPCLS(IOCTPA,IOCTPB)
CTOBE          N_CLS_TO_SPSP(ICLS) = N_CLS_TO_SPSP(ICLS) + 1
CTOBE        END DO
CTOBE       END DO
CTOBE *
CTOBE       IB_CLS_TO_SPSP(1) = 1
CTOBE       DO ICLS = 2, NCLS
CTOBE         IB_CLS_TO_SPSP(ICLS) =  
CTOBE      &  IB_CLS_TO_SPSP(ICLS-1) + N_CLS_TO_SPSP(ICLS-1)
CTOBE       END DO
CTOBE *
CTOBE       CALL ISETVC(N_CLS_TO_SPSP,IZERO,NCLS)
CTOBE       DO IOCTPA = 1, NOCTPA                  
CTOBE        DO IOCTPB = 1, NOCTPB                    
CTOBE          ICLS = ISPSPCLS(IOCTPA,IOCTPB)
CTOBE          I_CLS_TO_SPSP(1,
CTOBE          N_CLS_TO_SPSP(ICLS) = N_CLS_TO_SPSP(ICLS) + 1
CTOBE        END DO
CTOBE       END DO
CTOBE 
      SUBROUTINE ADTOR2S(RHO2S,RHO2TS,ITYPE,
     &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
*
* Add contributions to two electron density matrix spin-density matrix
* RHO2SS ( i.e. RHO2AA or RHO2B) or RHO2AB
*
* Jeppe Olsen, Adapted from ADTOR2, Sept. 2004 
*
*
* Itype = 1 => alpha-alpha or beta-beta loop
*              input is in form Rho2ts(ik,jl), i.ge.k, j.ge.l
* Itype = 2 => alpha-beta loop
*              input is in form Rho2ts(ij,kl)
*               
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION RHO2TS(*)
*. Input and output
      DIMENSION RHO2S(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to ADTOR2S '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' NI NJ NK NL = ', NI,NJ,NK,NL
        WRITE(6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
        WRITE(6,*) ' ITYPE = ',ITYPE
*
        IF(ITYPE.EQ.1) THEN
          IF(IOFF.EQ.KOFF) THEN
            NROW = NI*(NI+1)/2
          ELSE
            NROW = NI*NK
          END IF
          IF(JOFF.EQ.LOFF) THEN
            NCOL = NJ*(NJ+1)/2
          ELSE
            NCOL = NJ*NL
          END IF
        ELSE IF (ITYPE.EQ.2) THEN
          NROW = NI*NJ
          NCOL = NK*NL
        END IF  
        IF(ITYPE.EQ.1) THEN
          WRITE(6,*) ' Input matrix in the form (IK,JL) '
        ELSE
          WRITE(6,*) ' Input matrix in the form (IJ,KL) '
        END IF
        CALL WRTMAT(RHO2TS,NROW,NCOL,NROW,NCOL)
      END IF
*
      IF(ITYPE.EQ.1) THEN
*
* =======================================
*     Alpha-alpha or beta-beta term
* =======================================
*
*
       IJOFF = (JOFF-1)*NORB+IOFF
       KLOFF = (LOFF-1)*NORB+LOFF
        DO I = 1, NI
          DO J = 1, NJ
            IF(IOFF.EQ.KOFF) THEN 
              KMAX = I
            ELSE 
              KMAX = NK
            END IF
            DO K = 1, KMAX
              IF(LOFF.EQ.JOFF) THEN
               LMAX = J
              ELSE
               LMAX = NL
              END IF
              DO L = 1, LMAX
                IF(IOFF.NE.KOFF) THEN
                  IKIND = (K-1)*NI+I
                  NIK = NI*NK
                ELSE
                  IKIND = I*(I-1)/2 + K
                  NIK = NI*(NI+1)/2
                END IF
                IK_OUT = (I+IOFF-1)*(I+IOFF-2)/2 + K + KOFF - 1
                IF(JOFF.NE.LOFF) THEN
                  JLIND = (L-1)*NJ+J
                ELSE
                  JLIND = J*(J-1)/2+L
                END IF
                JL_OUT = (J+JOFF-1)*(J+JOFF-2)/2 + L + LOFF - 1
                IKJL_OUT = (JL_OUT-1)*NORB*(NORB+1)/2 + IK_OUT
                IKJL_IN = (JLIND-1)*NIK+IKIND
                RHO2S(IKJL_OUT) = RHO2S(IKJL_OUT)
     &                      + RHO2TS(IKJL_IN)
              END DO
            END DO
          END DO
        END DO
      ELSE IF(ITYPE.EQ.2) THEN
*
* ====================
*     Alpha-beta term
* ====================
*
      DO I = 1, NI
       DO J = 1, NJ
         DO K = 1, NK
           DO L = 1, NL
             IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI + (J-1)*NI + I
             IKLJ_OUT = (J+JOFF-2)*NORB**3+(L+LOFF-2)*NORB**2
     &                + (K+KOFF-2)*NORB   +(I+IOFF-1)
             RHO2S(IKLJ_OUT) = RHO2S(IKLJ_OUT) + RHO2TS(IJKLT)
            END DO
          END DO
        END DO
      END DO
*
      END IF
*
      IF(NTEST.GE.1000) THEN
        IF(ITYPE.EQ.1)  THEN
          WRITE(6,*) ' Updated two-body spin-SS density matrix '
          NDIM = NORB*(NORB+1)/2
        ELSE 
          WRITE(6,*) ' Updated two-body spin AB density matrix '
          NDIM = NORB**2
        END IF
        CALL WRTMAT(RHO2S,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE GROUP_STRTP(ISTRTP,NSTRTP,NGAS,NEL,ISPLIT1,ISPLIT2,
     &           NSTRFSPLT,IBSTRFSPLT,ISTRFSPLT,MXPNGAS)
*
* A number of string types ( supergroups...) with NEL electrons are given. 
* Order them according to the accumulatice number 
* of electrons after orbitals ISPLIT1 and ISPLIT2 have been 
* occupied (ISPLIT1 and ISPLIT2 are also output)
*
*. Jeppe Olsen, Oct. 2004, to get some routines running
*. with very many stringtypes
*
      INCLUDE 'implicit.inc'
*. input
      INTEGER ISTRTP(MXPNGAS,NSTRTP)
*. output 
      INTEGER NSTRFSPLT(NEL+1,NEL+1),IBSTRFSPLT(NEL+1,NEL+1)
      INTEGER ISTRFSPLT(NSTRTP)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Number of elecs per strtyp '
        CALL IWRTMA(ISTRTP,NGAS,NSTRTP,MXPNGAS,NSTRTP)
      END IF
*. Find the two split orbitals to be used for splitting
*. Could be refined by finding lowest and highesr gasspace 
*. with variable occupations, and only dividing the
*. orbital spaces with variable occupation
*
*. Find the two gasspaces used for splitting the types
      ISPLIT1 = (NGAS+1)/3
      ISPLIT2 = (2*NGAS+2)/3
      WRITE(6,*) ' GAS spaces used for splitting ', ISPLIT1,ISPLIT2
      IF(ISPLIT1.EQ.0) ISPLIT1 = 1
*. Number of strings with given accumulated occ
      IZERO = 0
      CALL ISETVC(NSTRFSPLT,IZERO,(NEL+1)*(NEL+1))
      DO ITP = 1, NSTRTP
       IEL1 = IELSUM(ISTRTP(1,ITP),ISPLIT1)
       IEL2 = IELSUM(ISTRTP(1,ITP),ISPLIT2)
C?     WRITE(6,*) ' ITP, IEL1, IEL2 ', ITP, IEL1, IEL2
       NSTRFSPLT(IEL1+1,IEL2+1) =  NSTRFSPLT(IEL1+1,IEL2+1) + 1
      END DO
*. And offsets
      IOFF  = 1
      DO IEL1 = 0, NEL
        DO IEL2 = 0, NEL
        IBSTRFSPLT(IEL1+1,IEL2+1) = IOFF
        IOFF = IOFF +  NSTRFSPLT(IEL1+1,IEL2+1)
        END DO
      END DO
*. And do the mappings
      CALL ISETVC(NSTRFSPLT,IZERO,(NEL+1)*(NEL+1))
      DO ITP = 1, NSTRTP
       IEL1 = IELSUM(ISTRTP(1,ITP),ISPLIT1)
       IEL2 = IELSUM(ISTRTP(1,ITP),ISPLIT2)
       NSTRFSPLT(IEL1+1,IEL2+1) =  NSTRFSPLT(IEL1+1,IEL2+1) + 1
       IADR = IBSTRFSPLT(IEL1+1,IEL2+1) +  NSTRFSPLT(IEL1+1,IEL2+1) -1
       ISTRFSPLT(IADR) = ITP
      END DO
*
      I_CHECK = 1
      IF(I_CHECK.EQ.1) THEN
*. Some debugging
        IZEROES = 0
        ISUM = 0
        DO I = 1, NSTRTP
          IF(ISTRFSPLT(I).EQ.0) IZEROES = 1
          ISUM = ISUM + ISTRFSPLT(I)
        END DO
        IF(IZEROES.EQ.1) THEN
          WRITE(6,*) ' Some types are missing '
          WRITE(6,*) ' The stringtypes split-reordered '
          CALL IWRTMA(ISTRFSPLT,1,NSTRTP,1,NSTRTP)
          STOP ' Some types are missing'
        END IF
        IF(ISUM.NE.NSTRTP*(NSTRTP+1)/2) THEN
          WRITE(6,*) ' Not all stringtypes are listed once '
          WRITE(6,*) ' Expected and actual sum ', 
     &    ISUM,NSTRTP*(NSTRTP+1)/2
          WRITE(6,*) ' The stringtypes split-reordered '
          CALL IWRTMA(ISTRFSPLT,1,NSTRTP,1,NSTRTP)
        END IF
      END IF
*     ^ End of checking should be performed
          
      
   
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' The string numbers split according to '
        WRITE(6,*) ' occupation of two orbital numbers '
        WRITE(6,*) ' (Orbital numbers : )', ISPLIT1,ISPLIT2
        CALL IWRTMA(ISTRFSPLT,1,NSTRTP,1,NSTRTP)
        WRITE(6,*) ' The NSTRFSPLT array '
        CALL IWRTMA(NSTRFSPLT,NEL+1,NEL+1,NEL+1,NEL+1)
      END IF
*
      RETURN
      END 
      SUBROUTINE FOCK_MAT_STANDARD(F,I12,FI,FA)
*
* Construct Fock matrix using the standard definition 
*
* F(j,i) =   SUM(K)     H(i,K)     * RHO1(j,K)
*          + SUM(M,K,L) (i M K L ) * RHO2(j M K L)
*
* Which becomes
*
* j: inactive: F(j,i) = 2(FI(JI) + 2 FA(JA))
* j: active  : F(j,i) = sum(k:active) D(j,k) FI(k,i)
*                     + sum(klm:active) d(jklm) (ik!lm)
* j: secondary:F(j,i) = 0
*
*
* Modified from FOCK_MAT   
*
* Jeppe Olsen
* Updated June 2010 with inclusion of inactive orbitals...
*
* Unless I12 = 2, only one-electron part is calculated
c      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION FI(*), FA(*)
      integer FI,FA
*. Output
CNW   DIMENSION F(*)
      integer F
*
      call luciaqenter('fmats')
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) 
        WRITE(6,*) ' Output from FOCK_MAT_STANDARD '
        WRITE(6,*) ' ------------------------------'
        WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        if (ga_nodeid().eq.0)
     &  WRITE(6,*) ' Input inactive and active Fock matrices '
CNW     CALL APRBLM2(FI,NTOOBS,NTOOBS,NSMOB,1)
        call ga_print(FI)
        if (ga_nodeid().eq.0) WRITE(6,*)
CNW     CALL APRBLM2(FA,NTOOBS,NTOOBS,NSMOB,1)
        call ga_print(FA)
      END IF
*
      ONE = 1.0D0
      ZERO = 0.0D0
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'FOO   ')
*
CNW   CALL COPVEC(FI,WORK(KINT1),NINT1)
      if (FI.ne.KINT1) call ga_copy(FI,KINT1)
      LEN_F =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
CNW   CALL SETVEC(F,ZERO,LEN_F)
      call ga_zero(F)
*. Largest set of orbitals with given symmetry and type
      MXTSOB_AC = 0
      MXTSOB    = 0
      DO ISM = 1, NSMOB
        DO IGAS = 1, NGAS
          MXTSOB_AC = MAX(MXTSOB_AC,NOBPTS(IGAS,ISM))
        END DO
        MXTSOB = MAX(MXTSOB,MXTSOB_AC,NINOBS(ISM),NSCOBS(ISM))
      END DO
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' MXTSOB, MXTSOB_AC = ', MXTSOB, MXTSOB_AC
      END IF
*. Allocate scratch space for 2-electron integrals and 
*. two-electron densities
      CALL MEMMAN(KLINT,MXTSOB_AC**3*MXTSOB,'ADDL  ',2,'KLINT ')!d
      CALL MEMMAN(KLDEN,MXTSOB_AC**4,'ADDL  ',2,'KLDEN ') !done
*. And a block of F
      MX2IBLK = MXTSOB ** 2
      CALL MEMMAN(KLFBLK,MX2IBLK,'ADDL  ',2,'KLFBL ') !done
*
      ONE = 1.0D0
      II = -2303
      IJ = -2303
      DO IJSM = 1, NSMOB
        ISM = IJSM
        JSM = IJSM
        NIJS = NTOOBS(IJSM)
*
        IF(IJSM.EQ.1) THEN
         IFOFF = 1
         IFOFFS= 1
        ELSE
         IFOFF = IFOFF+NTOOBS(IJSM-1)**2
         IFOFFS= IFOFFS+NTOOBS(IJSM-1)*(NTOOBS(IJSM-1)+1)/2
        END IF
*
        DO JGAS = 0, NGAS
         IF(JGAS.EQ.0) THEN
           IJ = 1
         ELSE IF(JGAS.EQ.1) THEN
           IJ = NINOBS(JSM)+1
         ELSE
           IJ = IJ + NOBPTS(JGAS-1,JSM)
         END IF
         IF(JGAS.EQ.0) THEN
          NJ = NINOBS(IJSM)
         ELSE
          NJ = NOBPTS(JGAS,IJSM)
         END IF
          DO IGAS = 0, NGAS+1
           IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) 
     &       ' Info for ISM IGAS JGAS',ISM,IGAS,JGAS
           END IF
*
           IF(IGAS.EQ.0) THEN
            NI = NINOBS(ISM)
           ELSE IF(IGAS.LE.NGAS) THEN
            NI = NOBPTS(IGAS,ISM)
           ELSE
            NI = NSCOBS(ISM)
           END IF
*
           IF(IGAS.EQ.0) THEN
             II = 1
           ELSE IF(IGAS.EQ.1) THEN 
             II = 1 + NINOBS(ISM)
           ELSE 
             II = II + NOBPTS(IGAS-1,ISM)
           END IF
           IF(NI*NJ.NE.0) THEN
*
*  =======================
*. block F(ijsm,jgas,igas)
*  =======================
*
            CALL SETVEC(dbl_mb(KLFBLK),ZERO,NI*NJ)
            IF(JGAS.EQ.0) THEN
*
*. Inactive part
*
            call ga_distribution(FI,ga_nodeid(),ilo,ihi,idum,idum)
            if (ihi.gt.0) then
            call ga_access(FI, ilo, ihi, 1, 1, ikfi, idum)
            call ga_access(FA, ilo, ihi, 1, 1, ikfa, idum)
             DO I = 1, NI
              DO J = 1, NJ
               IJMAX = MAX(I+II-1,J+IJ-1)
               IJMIN = MIN(I+II-1,J+IJ-1)
               IJSYM = IJMAX*(IJMAX-1)/2 + IJMIN + IFOFFS - 1
CNW            dbl_mb(KLFBLK-1+(I-1)*NJ+J) = 2.0D0*
CNW  &         (FI(IFOFFS-1+IJSYM)+FA(IFOFFS-1+IJSYM))
               if (IJSYM.ge.ilo.and.IJSYM.le.ihi) then 
               dbl_mb(KLFBLK-1+(I-1)*NJ+J) = 
     &           2.0D0*(dbl_mb(ikfi+IJSYM-ilo)+dbl_mb(ikfa+IJSYM-ilo))
c              WRITE(6,*) ' IJSM, I, J, FI(IJ), FA(IJ), F(IJ) =',
c    &         IJSM,I,J,dbl_mb(ikfi+IJSYM-ilo),dbl_mb(ikfa+IJSYM-ilo),
c    &         dbl_mb(KLFBLK-1+(I-1)*NJ+J)
C?   &         IJSM,I,J,FI(IFOFFS-1+IJSYM), FA(IFOFFS-1+IJSYM),
C?   &         WORK(KLFBLK-1+(I-1)*NJ+J)
               endif
              END DO
             END DO
            call ga_release(FI, ilo, ihi, 1, 1)
            call ga_release(FA, ilo, ihi, 1, 1)
            endif
            call ga_sync()
            call ga_dgop(1,dbl_mb(KLFBLK),ni*nj,'+')
            ELSE IF(JGAS.LE.NGAS) THEN
*
*. Active part
*
* 1 : One-electron part 
             DO KGAS = 1, NGAS
              KSM = IJSM
              NK = NOBPTS(KGAS,KSM)
*. blocks of one-electron integrals and one-electron density
              CALL GETD1(dbl_mb(KLDEN),JSM,JGAS,KSM,KGAS,1)
              CALL GETH1(dbl_mb(KLINT),ISM,IGAS,KSM,KGAS)
              IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
                WRITE(6,*) 
     &          ' 1-e ints for ISM IGAS KGAS ',ISM,IGAS,KGAS
                CALL WRTMAT(dbl_mb(KLINT),NI,NK,NI,NK)
                WRITE(6,*) 
     &          ' 1-e densi for ISM JGAS KGAS ',ISM,JGAS,KGAS
                CALL WRTMAT(dbl_mb(KLDEN),NJ,NK,NJ,NK)
              END IF
*. And then a matrix multiply( they are pretty much in fashion 
*. these days )
              CALL MATML7(dbl_mb(KLFBLK),dbl_mb(KLDEN),dbl_mb(KLINT),
     &                    NJ,NI,NJ,NK,NI,NK,ONE,ONE,2)
               IF(NTEST.GE.1000.and. ga_nodeid().eq.0) THEN
                 WRITE(6,*) ' Updated block '
                 CALL WRTMAT(dbl_mb(KLFBLK),NJ,NI,NJ,NI)
               END IF
 
             END DO
             IF(NTEST.GE.1000.and. ga_nodeid().eq.0) THEN
              WRITE(6,*) ' One-electron contributions'
              WRITE(6,*) ' =========================='
              CALL WRTMAT(dbl_mb(KLFBLK),NJ,NI,NJ,NI)
             END IF
             IF(I12.EQ.2) THEN
*. 2 : Two-electron part
             DO KSM = 1, NSMOB
             DO LSM = 1, NSMOB
*. Obtain MSM
              CALL  SYMCOM(3,1,KSM,LSM,KLSM)
              CALL  SYMCOM(3,1,KLSM,ISM,IKLSM)
              IMKLSM = 1
              CALL  SYMCOM(2,1,IKLSM,MSM,IMKLSM)
*
              DO MGAS = 1, NGAS
              DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
               ISKIP = 0
               XFACTOR = 1.0D0
*. Check that K,L,J are occupied and ISKIP = 0
                NM = NOBPTS(MGAS,MSM)
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
                IF(NM*NK*NL.NE.0) THEN
               
*. Blocks of density matrix and integrals: D2(K L, J M) (K L ! I M)
*. Modern form where (oo!og) must be accessed with fourth index being
*. general (why does reorganizations lead to less flexibility ???
                  IXCHNG = 0
                  ICOUL  = 1
                  ONE = 1.0D0
* Obtain ( L K ! M I)
                  CALL LUCIAGETINT(dbl_mb(KLINT),
     &                 LGAS,LSM,KGAS,KSM,MGAS,MSM,IGAS,ISM,
     &                 IXCHNG,0,0,ICOUL,ONE,ONE)
*. Obtain Rho2(L K, M J)
                  CALL GETD2 (dbl_mb(KLDEN),
     &                 LSM,LGAS,KSM,KGAS,MSM,MGAS,JSM,JGAS,1)
                  NKLM = NK*NL*NM
*. And bring it home with a matrix mutiply, Sum(KLM) Rho2(lkm,j)Int2(lkm,i)
                  CALL MATML7(dbl_mb(KLFBLK),dbl_mb(KLDEN),
     &                        dbl_mb(KLINT),
     &                        NJ,NI,NKLM,NJ,NKLM,NI,ONE,XFACTOR,1)
                END IF !End if nonvanishing block
              END DO
              END DO
              END DO
             END DO
             END DO
             END IF
             IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
              WRITE(6,*) ' One- + two-electron contributions'
              WRITE(6,*) ' ================================='
              CALL WRTMAT(dbl_mb(KLFBLK),NJ,NI,NJ,NI)
             END IF
            END IF
*           ^ End of inactive/active switch
*. Block has been constructed , transfer to -complete- 
*. symmetry blocked Fock matrix
*
*
            call ga_distribution(F,ga_nodeid(),ilo,ihi,idum,idum)
            DO J = 1, NJ
              DO I = 1, NI
C?              WRITE(6,*) 'IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1',
C?   &                      IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1
                IINDX = IFOFF-1+(I+II-1-1)*NIJS + J+IJ-1
                if (iindx.ge.ilo.and.iindx.le.ihi) then
                call ga_put(F,IINDX,IINDX,1,1,
     &                      dbl_mb(KLFBLK-1+(I-1)*NJ+J),1)
                endif
CNW             F(IFOFF-1+(I+II-1-1)*NIJS + J+IJ-1 ) = 
CNW  &          dbl_mb(KLFBLK-1+(I-1)*NJ+J)
              END DO
            END DO
*
           END IF
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        if (ga_nodeid().eq.0) then
        WRITE(6,*)
        WRITE(6,*) ' Final F from FOCK_MAT_STANDARD '
        WRITE(6,*) ' =============================='
        endif
CNW     CALL APRBLM2(F,NTOOBS,NTOOBS,NSMOB,0)
        call ga_print(F)
      END IF
      call luciaqexit('fmats')
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'FOO   ')
      RETURN
      END
      SUBROUTINE STR_TO_PHFORM(ISTR, NEL, IPHSTR, IPHGAS,
     &                         ITFT,NINOB,NACOB,NPHEL)
*
*. Convert a string from standard form- specifying occupied electrons
*  to ph form - specifying unoccpied holes and occupied particles, 
*  starting from ninob + 1
*
* In connection with changing to ph form of operator, a sign change is
* required for holes but not particles. This sign is also contained in IPHSTR
*
*. Jeppe Olsen
*. General input
      INTEGER IPHGAS(*), ITFT(NACOB+NINOB)
*. Specific input
      INTEGER ISTR(NEL)
*. Output
      INTEGER IPHSTR(*)
*
      IEL = 1
      NPHEL = 0
C?    WRITE(6,*) ' NINOB, NACOB = ', NINOB, NACOB
      DO IORB = NINOB+1, NINOB + NACOB
C?       WRITE(6,*) ' IORB, ITFT, IPHGAS = ', 
C?   &               IORB, ITFT(IORB), IPHGAS(ITFT(IORB))
*. Is orbital occupied?
       I_AM_OCCUPIED = 0
       IF(IEL.LE.NEL) THEN
         IF(ISTR(IEL).EQ.IORB) I_AM_OCCUPIED = 1
       END IF
       IF(I_AM_OCCUPIED.EQ.1) THEN
*. Orbital is in string
         IF(IPHGAS(ITFT(IORB)).EQ.1) THEN
           NPHEL = NPHEL + 1
           IPHSTR(NPHEL) = IORB
           IEL = IEL + 1
C?         WRITE(6,*) ' particle electron added, IEL, IORB, NPHEL = ',
C?   &                  IEL, IORB, NPHEL
         ELSE
           IEL = IEL + 1
         END IF
       ELSE 
*. Orbital is not in string
         IF(IPHGAS(ITFT(IORB)).EQ.2) THEN
           NPHEL = NPHEL + 1
           IPHSTR(NPHEL) = -IORB
         END IF
       END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Standard and PH form(with sign) of string'
        CALL IWRTMA(ISTR,1,NEL,1,NEL)
        CALL IWRTMA(IPHSTR,1,NPHEL,1,NPHEL)
      END IF
*
      RETURN
      END
      SUBROUTINE GASDIAS_PH(NAEL,IASTR,NBEL,IBSTR,
     &           NACOB,DIAG,NSMST,H,XA,XB,SCR,RJ,RK,
     &           NSSOA,NSSOB,LUDIA,ECORE,
     &           PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,I12,
     &           IBLTP,NBLOCK,IBLKFO,IPHGAS,ITFT,IPHSTR,NPHELFSPGP,
     &           IBSPGPA,IBSPGPB,NINOB,ISCR)
*
* Calculate determinant diagonal using PH formalism
* Turbo-ras version
*
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, June 2010
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
C     REAL * 8  INPROD
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H(NTOOB)
      DIMENSION IPHGAS(*), ITFT(NTOOB), NPHELFSPGP(*)
*.
      INCLUDE 'cprnt.inc'
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NTOOB),XB(NTOOB),SCR(2*NTOOB),ISCR(NTOOB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IPHSTR(*)
      DIMENSION RJKAA(*), buffer(100)
*. Output
cnw   DIMENSION DIAG(*)
      integer diag, ilo, ihi
*
      integer nxtask, next, itask, icount, nproc
      external nxtask
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
      NOCOB = NINOB + NACOB
*
      if (ga_nodeid().eq.0) THEN 
      IF( NTEST .GE. 20.OR.IPRINTEGRAL.GE.100) THEN 
        WRITE(6,*) ' GASDIAS_PH in action:'
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Core energy ', ECORE
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NOCOB,NOCOB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NOCOB,NOCOB,NTOOB,NTOOB)
        END IF
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK               
          WRITE(6,'(10X,8I3,2I8)') (IBLKFO(II,IBLOCK),II=1,8)
cVOG      WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' IBLTP: '
        CALL IWRTMA(IBLTP,1,NSMST,1,NSMST)
*
        WRITE(6,*) ' I12 = ',I12
      END IF
      end if
*
*  Diagonal elements according to Handys formulae
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
      IF(LUDIA.GT.0) CALL REWINO(LUDIA)
      call ga_distribution(DIAG,ga_nodeid(),ilo,ihi,idum,idum)
*
      DO IBLK = 1, NBLOCK
*
       if (((IDET+1).ge.ilo).and.((IDET+1).le.ihi)) then
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
*. Number of ph operators in alpha and beta
        NPHELA = NPHELFSPGP(IATP-1+IBSPGPA)
        NPHELB = NPHELFSPGP(IBTP-1+IBSPGPB)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NOCOB,0,IDUM,IDUM)
        IOFF =  1
        DO IA = 1, NSSOA(IASM,IATP)
*. Convert to PH form and save in IPHSTR
C       STR_TO_PHFORM(ISTR, NEL, IPHSTR, IPHGAS,
C    &                ITFT,NINOB,NACOB,NPHEL)
          CALL STR_TO_PHFORM(IASTR(1,IA),NAEL,IPHSTR((IA-1)*NPHELA+1),
     &                       IPHGAS,ITFT,NINOB,NACOB,NPHELA_2)
          EAA = 0.0D0
          DO IEL = 1, NPHELA
            IF(IPHSTR((IA-1)*NPHELA+IEL).GT.0) THEN
              IAEL = IPHSTR((IA-1)*NPHELA+IEL) 
              SIGNI = 1.0D0
            ELSE
              IAEL = -IPHSTR((IA-1)*NPHELA+IEL) 
              SIGNI = -1.0D0
            END IF
            EAA = EAA + H(IAEL)*SIGNI
            IF(I12.EQ.2) THEN
              DO JEL = 1, NPHELA
                IF(IPHSTR((IA-1)*NPHELA+JEL).GT.0) THEN
                  JAEL = IPHSTR((IA-1)*NPHELA+JEL) 
                  SIGNJ = 1.0D0
                ELSE
                  JAEL = -IPHSTR((IA-1)*NPHELA+JEL) 
                  SIGNJ = -1.0D0
                END IF
                EAA =   EAA + 0.5D0*RK(JAEL,IAEL )*SIGNI*SIGNJ
              END DO   
            END IF
          END DO
          RJKAA(IA-IOFF+1) = EAA 
        END DO
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NOCOB,0,IDUM,IDUM)
        IBSTRT = 1                
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
*. Obtain in iscr ph-form of string
          CALL STR_TO_PHFORM(IBSTR(1,IB),NBEL,ISCR(1),
     &                       IPHGAS,ITFT,NINOB,NACOB,NPHELB_2)
*
*. Terms depending only on IB
*
          HB = 0.0D0
          RJBB = 0.0D0
          CALL SETVEC(XB,0.0D0,NOCOB)
*
          DO IEL = 1, NPHELB
            IF(ISCR(IEL).GT.0) THEN
              IBEL = ISCR(IEL)
              SIGNIB = 1.0D0
            ELSE 
              IBEL = -ISCR(IEL)
              SIGNIB = -1.0D0
            END IF
            HB = HB + H(IBEL)*SIGNIB
*
            IF(I12.EQ.2) THEN
              DO JEL = 1, NPHELB
                IF(ISCR(JEL).GT.0) THEN
                  JBEL = ISCR(JEL)
                  SIGNJB = 1.0D0
                ELSE
                  JBEL = -ISCR(JEL)
                  SIGNJB = -1.0D0
                END IF
                RJBB = RJBB + RK(JBEL,IBEL )*SIGNIB*SIGNJB
              END DO
*
              DO IORB = NINOB+1, NOCOB
                XB(IORB) = XB(IORB) + RJ(IORB,IBEL)*SIGNIB
              END DO 
            END IF
          END DO
          EB = HB + 0.5D0*RJBB + ECORE
*
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1                 
          END IF
          IASTOP = NSSOA(IASM,IATP) 
*
cVOG      if (itask.eq.next) then
          icount = 0
          call SETVEC(buffer,0.0d0,IASTOP)
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            ITDET = ITDET + 1
            icount = icount + 1
            buffer(icount) = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NPHELA
              IF(IPHSTR((IA-1)*NPHELA+IEL).GT.0) THEN
                buffer(icount) = buffer(icount) 
     &                           + XB(IPHSTR((IA-1)*NPHELA+IEL)) 
              ELSE
                buffer(icount) = buffer(icount) 
     &                           - XB(-IPHSTR((IA-1)*NPHELA+IEL)) 
              END IF
            END DO
CNW         DIAG(IDET) = X
CNW         IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
CNW  &         IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
            IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &         IB.EQ.IA) buffer(icount) = buffer(icount) + XADD
cVOG        if (icount.eq.100.OR.IA.eq.IASTOP) then
            if (IA.eq.IASTOP) then
               call ga_put(DIAG,IDET-icount+1,IDET,1,1,buffer,1)
               call SETVEC(buffer,0.0d0,IASTOP)
               icount = 0
            endif
          END DO
*         ^ End of loop over alpha strings|
        END DO
*       ^ End of loop over betastrings
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
       else
        IDET = IDET + IABS(IBLKFO(8,IBLK))
       end if
      END DO
*        ^ End of loop over blocks
cVOG  next = nxtask(-nproc, 1)
      call ga_sync()
      IF(NTEST.GE.10) WRITE(6,*)
     &' Number of diagonal elements generated ',ITDET
*
      IF(NTEST .GE.100 .AND.ICISTR.LE.1 ) THEN
        if (ga_nodeid().eq.0) WRITE(6,*) ' CIDIAGONAL ',DIAG
CNW     CALL WRTMAT(DIAG(1),1,IDET,1,IDET)
        call ga_print(diag)
      END IF
*
      IF ( ICISTR.GE.2 ) CALL ITODS(-1,1,-1,LUDIA)
*
      RETURN
      END
*
      SUBROUTINE FIND_INTARR_IN_ORD_INTARR_NLIST(INTARR,INTARR_LIST,
     &           NELMNT,NELMNT_MAX,NARR,IREO,IMET,INUM)
* An array of integers INTARR(NELMNT) are given and an ordered list of such 
* arrays INTARR_ARR(NELMNT,NARR) are given. Find the address of INTARR in the 
* ordered list. The arrays is assumed ordered according to last differing 
* digit.
* If the array is not in the list INUM = 0 is returned
*
* IMET = 1: Search through all elements in list
* IMET = 2: Search by bisection 
*
*. Jeppe Olsen, for speeding up the SPGP_AC code
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER INTARR(NELMNT), INTARR_LIST(NELMNT_MAX,NARR)
*. IREO is included to allow for a reordering (to be programmed, perhaps)
*
      INUM = 0
      NTEST = 0000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Wellcome to FIND_INTARR_IN_ORD... '
        WRITE(6,*) ' ================================== '
      END IF
*
      IPROBLEMO = 0
*
      IF(IMET.EQ.1) THEN
*
*. Just go through the table and check
*
       DO IARR = 1, NARR
         IM_FOUND = 1
         DO IELMNT = 1, NELMNT
          IF(INTARR(IELMNT).NE.INTARR_LIST(IELMNT,IARR)) 
     &    IM_FOUND = 0
         END DO! of loop over elements
         IF(IM_FOUND.EQ.1) THEN
           INUM = IARR
           GOTO 2303
         END IF
       END DO !of loop over arrays in list
 2303  CONTINUE
      END IF
*
      IF(IMET.EQ.2) THEN
*
*. Check through a binary search
*
        IMIN = 1
        IMAX = NARR
        MAXLOOP = 30
        NLOOP = 0
*. Loop over bisections
 1001   CONTINUE
          NLOOP = NLOOP + 1
          IF(NLOOP.GT.MAXLOOP) THEN
*. Seems to  be trapped in another of Jeppes errors
            IPROBLEMO = 1
            WRITE(6,*) ' Problem .. '
            WRITE(6,*) ' IMIN, IMAX, IELMNT(OLD) ',
     &      IMIN,IMAX, IELMNT
            GOTO 1002
          END IF
*
          IELMNT = (IMAX+IMIN)/2
          IF(NTEST.GE.1000) WRITE(6,'(A,4(2X,I6))')
     &    ' NLOOP, IMIN, IMAX, IELMNT = ', NLOOP, IMIN, IMAX, IELMNT
          
*. Is the array IELMNT before or after the array in question
* Last digit is used for comparison
          IFL = 2
          CALL COMP_TWO_INTARR(INTARR,INTARR_LIST(1,IELMNT),NELMNT,IFL,
     &                         IFIRST)
          IF(IFIRST.EQ.0) THEN
* element has been located 
            INUM = IELMNT
            GOTO 1002
          END IF
          IF(IMIN.GE.IMAX) THEN 
* Min and max are identical or inverted and we have not identified element
            GOTO 1002
          END IF
          IF(IFIRST.EQ.1) THEN
* Element IELMNT is above target so reduce MAX
            IMAX = IELMNT - 1
          ELSE
* Element IELMNT is below target so increase MIN
            IMIN = IELMNT +1
          END IF
          GOTO 1001
 1002   CONTINUE ! End of loop over bisections
      END IF ! switch of method
*. 
      IF(NTEST.GE.100.OR.IPROBLEMO.NE.0) THEN
        WRITE(6,*) ' Input array '
        CALL IWRTMA(INTARR,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' List of possible output groups '
        CALL IWRTMA(INTARR_LIST,NELMNT,NARR,NELMNT_MAX,NARR)
        WRITE(6,*) ' Input in list = ', INUM
      END IF
*
      IF(IPROBLEMO.EQ.1) THEN
       WRITE(6,*) ' Trapped in loop in FINT_INTARR... '
       STOP       ' Trapped in loop in FINT_INTARR... '
      END IF
*
      RETURN
      END
       
      

      SUBROUTINE COMP_TWO_INTARR(I1,I2,NELMNT,IFL,IFIRST)
*
* Compare two integer arrays I1 and I2 and decide which of thise 
* comes first in a ordered list of integer arrays. 
* IFL = 1 => The ordering is based on first differing digit
*      =2 => The ordering is bases in last differing digit
*
* Output: IFIRST
* IFIRST = 1 => I1 comes first
* IFIRST = 2 => I2 comes first
* IFIRST = 0 => The two arrays are identical
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER I1(NELMNT),I2(NELMNT)
*
* Decide the order in which array is searched
*
      IF(IFL.EQ.1) THEN
        ISTART = 1
        IEND = NELMNT
        ISTEP = 1
      ELSE
        ISTART = NELMNT
        IEND = 1
        ISTEP = -1
      END IF
*
*. Find first differing element and compare
*
      IFIRST = 0
      DO IELMNT = ISTART, IEND, ISTEP
        IF(I1(IELMNT).NE.I2(IELMNT)) THEN
          IF(I1(IELMNT).LT.I2(IELMNT)) THEN
           IFIRST = 1
          ELSE 
           IFIRST = 2
          END IF
          GOTO 2303
        END IF
      END DO
 2303 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Two integer arrays to be compared '
        CALL IWRTMA(I1,1,NELMNT,1,NELMNT)
        CALL IWRTMA(I2,1,NELMNT,1,NELMNT)
        WRITE(6,'(A, I1)') ' Ifirst = ', IFIRST
      END IF
*
      RETURN
      END
      SUBROUTINE CHECK_IS_OCC_IN_ENGSOCC(IGSOCCL,ISPC,IM_IN)
*
* An accution IGSOCC is given. Check to see if this allowed
* according to the restrictions in ensemble I of gaspaces 
*
* IF ISPC > 0, check against the occupations in CISPACE ISPC
* IF ISPC =-1, check against the occupations in the compound space
*
*. Jeppe Olsen, Geneva, Feb. 15, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IGSOCCL(NGAS)
*
      NTEST = 00
*. Number of electrons in the ensemble 
      NELEC_ENS = 0
      DO I = 1, LENSGS(1)
        NELEC_ENS = NELEC_ENS + IGSOCCL(IENSGS(I,1))
      END DO
*
      IM_IN = 0
      IF(ISPC.GE.1) THEN
        DO I = 1, NELVAL_IN_ENSGS(ISPC)
          IF(NELEC_ENS.EQ.IEL_IN_ENSGS(I,ISPC))IM_IN = 1
        END DO
      ELSE IF (ISPC.EQ.-1) THEN
        DO I = 1, NELVAL_IN_ENSGS_CMP
          IF(NELEC_ENS.EQ.IEL_IN_ENSGS_CMP(I))IM_IN = 1
        END DO
      ELSE 
        WRITE(6,*) 
     &  ' Illegal value of ISPC in CHECK_IS_OCC_IN_ENGSOCC ', ISPC
        STOP ' Illegal value of ISPC in CHECK_IS_OCC_IN_ENGSOCC '
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Occupation '
       CALL IWRTMA(IGSOCCL,1,NGAS,1,NGAS)
       WRITE(6,*) ' CI-space in question ', ISPC
       WRITE(6,'(A,I3)') ' Number of electrons in ensemble 1 ',NELEC_ENS
       IF(IM_IN.EQ.1) THEN
        WRITE(6,*) ' Allowed '
       ELSE
        WRITE(6,*) ' Not allowed '
       END IF
      END IF
*
      RETURN
      END
      SUBROUTINE Z_ENSOCC_CMP
*
* Obtain the allowed occupation of the ensemble gas pace (1) for the 
* compound space
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Find the largest number of electrons included 
      NEL_MAX = 0
      DO ISPC = 1, NCISPC
        DO I = 1, NELVAL_IN_ENSGS(ISPC)
          NEL_MAX = MAX(NEL_MAX,IEL_IN_ENSGS(I,ISPC))
         END DO
      END DO
*
      IEL_VAL = 0 
      DO IEL = 0, NEL_MAX
        IEL_IS_IN = 0
        DO ISPC = 1, NCISPC
          DO I = 1, NELVAL_IN_ENSGS(ISPC)
            IF(IEL.EQ.IEL_IN_ENSGS(I,ISPC)) IEL_IS_IN = 1
          END DO
        END DO
        IF(IEL_IS_IN.EQ.1) THEN
          IEL_VAL = IEL_VAL + 1
          IEL_IN_ENSGS_CMP(IEL_VAL) = IEL
        END IF
      END DO
      NELVAL_IN_ENSGS_CMP = IEL_VAL
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,I3)') 
     &  ' Number of occupations of ensemble gas in compound state',
     &   NELVAL_IN_ENSGS_CMP
        N = NELVAL_IN_ENSGS_CMP
        WRITE(6,*) ' And the allowed occupations '
        CALL IWRTMA(IEL_IN_ENSGS_CMP,1,N,1,N)
      END IF
*
      RETURN
      END
      SUBROUTINE EXPCIV_CSF(ISM,ISPCIN,LUIN,ISPCUT,LUUT,LBLK,
     &                  NROOT,ICOPY,NTESTG)
*
* Expand CI vector in CI space ISPCIN to CI vector in ISPCUT
* Input vector is supposed to be on LUIN
* Output vector will be placed on unit LUUT
*. If ICOPY .ne. 0 the output vectors will be copied to LUIN
*
* Storage form is defined by ICISTR 
*
* CSF version, march 2012
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'spinfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      CALL LUCIAQENTER('EXPCV')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPCIV')
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from EXPCIV_CSF'
        WRITE(6,*) ' LUIN, LUUT= ', LUIN,LUUT
        WRITE(6,*) ' NROOT, NTESTG = ', NROOT, NTESTG
      END IF
*
* The occupation classes of the input and output vectors
*
      CALL MEMMAN(KLOCCLSIN,NOCCLS_MAX,'ADDL  ',1,'OCCLIN') !d
      CALL MEMMAN(KLOCCLSUT,NOCCLS_MAX,'ADDL  ',1,'OCCLUT') !d
*
C     OCCLS_IN_CISPACE(NOCCLS_ACT,IOCCLS_ACT,
C    &           NOCCLS,IOCCLS_OCC,NGAS,
C    &           NMNMX_SPC,IMNMX_SPC,MNMX_OCC,ISPC)
      CALL OCCLS_IN_CISPACE(NOCCLSIN,int_mb(KLOCCLSIN),
     &       NOCCLS_MAX,int_mb(KIOCCLS),NGAS,
     &       LCMBSPC(ISPCIN),ICMBSPC(1,ISPCIN),IGSOCCX,ISPCIN)
      CALL OCCLS_IN_CISPACE(NOCCLSUT,int_mb(KLOCCLSUT),
     &       NOCCLS_MAX,int_mb(KIOCCLS),NGAS,
     &       LCMBSPC(ISPCUT),ICMBSPC(1,ISPCUT),IGSOCCX,ISPCUT)

*
*. Vectors for holding input and output expansions
*
      IF(ICISTR.EQ.1) THEN
*. Complete vector is one record
        NCSFIN = NCSF_PER_SYM_GN(ISM, ISPCIN)
        NCSFUT = NCSF_PER_SYM_GN(ISM, ISPCUT)
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' NCSCFIN, NCSFUT = ',NCSFIN, NCSFUT
        LENGTHIN = NCSFIN
        LENGTHUT = NCSFUT
      ELSE
        LENGTHIN = NCS_FOR_OCCLS_MAX
        LENGTHUT = NCS_FOR_OCCLS_MAX
      END IF
*. 
      CALL MEMMAN(KLVEC,LENGTHIN,'ADDL  ',2,'VCCSFI') !done
      CALL MEMMAN(KLVECUT,LENGTHUT,'ADDL  ',2,'VCCSFU') !done
*. Array for holding dimension as function of occupation class
      CALL MEMMAN(KLLOCCLS_SM,NOCCLS_MAX,'ADDL  ',1,'OC_LSM') !d
*
*. Initial vectors in initial files
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial vectors in EXPCIV '
        WRITE(6,*) ' =========================='
*
        CALL REWINO(LUIN)
        DO IROOT = 1, NROOT
          WRITE(6,*) ' Root number ', IROOT 
          IF(ICISTR.EQ.1) THEN
            CALL FRMDSC(dbl_mb(KLVEC),NCSFIN,-1,LUIN,IMZERO,IAMPACK)
            WRITE(6,*) ' NCSFIN = ', NCSFIN
            CALL WRTMAT(dbl_mb(KLVEC),1,NCSFIN,1,NCSFIN)
          ELSE
            CALL WRTVCD(dbl_mb(KLVEC),LUIN,0,-1)
          END IF
        END DO
      END IF
*
*. And then the work
*
      CALL REWINO(LUIN)
      CALL REWINO(LUUT)
*
*  Zero block before any blocks have been read in.
*  Use IDIAG to decide
      IF(IDIAG.EQ.1) THEN
        IAMPACK = 0
      ELSE
        IAMPACK = 1
      END IF
C     WRITE(6,*) ' IAMPACK in EXPCIV ', IAMPACK
*
* Well, in the current implementation, there is a (big) difference
* between the way the CSF's are stored for ICNFBAT/ICISTR = 1 and
* ICNFBAT/ICISTR = 2 or more
*
* For ICNFBAT = 1 the CSF are stored as
* =====================================
* Loop over number of open orbitals
*  Loop over occupation classes
*   Loop over configurations
*    Loop over CSF's of this config
*    End loop
*   End loop
*  End loop
* End loop
*
* For ICNFBAT > 1 the CSF are stored as
* =====================================
* Loop over occupation classes
*  Loop over number of open orbitals
*   Loop over configurations
*    Loop over CSF's of this config
*    End loop
*   End loop
*  End loop
* End loop
* 
* Due to this (design flaw), the two cases must be handled seperately
*
 
      IF(ICNFBAT.GT.1) THEN
*
*. Occupation classes are stored separately
*  =========================================
*
* Number of CSF's per occ class for given symmetry
*
       IF(NTEST.GE.1000) WRITE(6,*) ' ISM, NIRREP = ', ISM, NIRREP
       CALL EXTRROW(int_mb(KNCS_FOR_OCCLS),ISM,NIRREP,NOCCLS_MAX,
     &              int_mb(KLLOCCLS_SM))
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Number of CSF of right sym per occls'
         CALL IWRTMA(int_mb(KLLOCCLS_SM),1,NOCCLS_MAX,1,NOCCLS_MAX)
       END IF
*
       DO IROOT = 1, NROOT
C        EXP_BLKVEC(LU_IN,NBLK_IN, IBLK_IN,
C    &              LU_OUT,NBLK_OUT,IBLK_OUT,
C    &              LBLK,ITASK,VEC,VEC_OUT,IREW,ICISTR)
         INCORE = 0
         CALL EXP_BLKVEC(LUIN,NOCCLSIN,int_mb(KLOCCLSIN),
     &        LUUT,NOCCLSUT,int_mb(KLOCCLSUT),
     &        int_mb(KLLOCCLS_SM),1,dbl_mb(KLVEC),dbl_mb(KLVECUT),
     &        0,ICISTR,INCORE)
*. The file is positioned before EOV vector mark, read this to allow 
*. correct read of next root
         CALL IFRMDS(LBL,1,-1,LUIN)
       END DO
*
      ELSE
*
* ================================================================
* Complete vector stored in one record, but ordered according to 
* number of open orbitals
* ================================================================
*
       DO IROOT = 1, NROOT
         CALL FRMDSC(dbl_mb(KLVEC),NCSFIN,-1,LUIN,IMZERO,IAMPACK)
         IB_IN = 1
         IB_UT = 1
         DO IOPEN = 0, MAXOP
         IF(NPCSCNF(IOPEN+1).NE.0) THEN
*. Number of configurations per occupation class for this IOPEN
C   NCN_PER_OP_SM(MAXOP+1,NIRREP,NOCCLS_MAX)
           IROW = (ISM-1)*(MAXOP+1) + IOPEN + 1
           NROW = NIRREP*(MAXOP+1)
           CALL EXTRROW(int_mb(KNCN_PER_OP_SM),IROW,NROW,NOCCLS_MAX,
     &                int_mb(KLLOCCLS_SM))
*. Number of CSFs per occupation class for this IOPEN
           NCSF_PER_CONF = NPCSCNF(IOPEN+1)
           CALL ISCLVEC(int_mb(KLLOCCLS_SM),NCSF_PER_CONF,NOCCLS_MAX)
C                   IELSUM_IND(IACT,NACT,IVEC)
           INCORE = 1
           CALL EXP_BLKVEC(LUIN,NOCCLSIN,int_mb(KLOCCLSIN),
     &          LUUT,NOCCLSUT,int_mb(KLOCCLSUT),
     &          int_mb(KLLOCCLS_SM),1,
     &          dbl_mb(KLVEC-1+IB_IN),dbl_mb(KLVECUT-1+IB_UT),
     &          0,ICISTR,INCORE)
           NCSFINL = IELSUM_IND(int_mb(KLOCCLSIN),NOCCLSIN,
     &              int_mb(KLLOCCLS_SM))
           NCSFUTL = IELSUM_IND(int_mb(KLOCCLSUT),NOCCLSUT,
     &              int_mb(KLLOCCLS_SM))
           IB_IN = IB_IN + NCSFINL
           IB_UT = IB_UT + NCSFUTL
*
         END IF! NPCSCNF .NE. 0
         END DO! IOPEN
         CALL TODSC(dbl_mb(KLVECUT),NCSFUT,-1,LUUT)
       END DO !IROOT
      END IF !ICNBAT switch
*
      IF(ICOPY.NE.0) THEN
*. Copy expanded vectors to LUIN
        CALL REWINO(LUIN)
        CALL REWINO(LUUT)
        DO IROOT = 1, NROOT
         IF(ICISTR.EQ.1) THEN
           CALL FRMDSC(dbl_mb(KLVECUT),NCSFUT,-1,LUUT,IMZERO,IAMPACK)
           CALL  TODSC(dbl_mb(KLVECUT),NCSFUT,-1,LUIN)
         ELSE
           CALL COPVCD(LUUT,LUIN,dbl_mb(KLVECUT),0,-1)
         END IF
        END DO
      END IF! ICOPY
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output  vectors in EXPCIV '
        WRITE(6,*) ' =========================='
        IF(ICOPY.EQ.0) THEN
          LLUUT = LUUT
        ELSE
          LLUUT = LUIN
        END IF
*
        CALL REWINO(LLUUT)
        DO IROOT = 1, NROOT
         WRITE(6,*) ' Root number ', IROOT 
         IF(ICISTR.EQ.1) THEN
          CALL FRMDSC(dbl_mb(KLVECUT),NCSFUT,-1,LLUUT,IMZERO,IAMPACK)
          CALL WRTMAT(dbl_mb(KLVECUT),1,NCSFUT,1,NCOMBUT)
         ELSE
          CALL WRTVCD(dbl_mb(KLVECUT),LLUUT,0,-1)
         END IF
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXPCIV')
C?    STOP ' Enforced stop at end of EXPCIV_CSF '
      CALL LUCIAQEXIT('EXPCV')
*
      RETURN
      END 
      SUBROUTINE ABEXPQ(A,B,IASM,IBSM,AB)
*
* Contribution from active orbitals to
* expectation value of product of two one-electron operators
*
* Jeppe Olsen, May 2012
*
* <0!AB!0> = sym(ijkl) A(ij)B(kl) d(ijkl) + sum(ij) rho1(ij) (AB)(ij)
* Obtain contributions from active orbitals to 
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      REAL*8 INPROD
*
*. Input, A and B are required to be symmetrypacked
*. complete  form
      DIMENSION A(*),B(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ABEXP '
        WRITE(6,*) ' ================ '
      END IF
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Energy under construction '
        WRITE(6,*) ' =========================='
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'EN_FRM')
*
      E1 = 0.0D0
      E2 = 0.0D0
*. Largest set of orbitals with given symmetry and type
      MXTSOB = 0
      DO ISM = 1, NSMOB
      DO IGAS = 1, NGAS
        MXTSOB = MAX(MXTSOB,NOBPTS(IGAS,ISM))
      END DO
      END DO
*. Allocate scratch space for 2-electron integrals and 
*. two-electron densities
      MX4IBLK = MXTSOB ** 4
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT') !done
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN') !done
      ONE = 1.0D0
      DO ISM = 1, NSMOB
       DO JSM = 1, NSMOB
        CALL  SYMCOM(3,1,ISM,JSM,IJSM)
        DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
           NI = NOBPTS(IGAS,ISM)
           NJ = NOBPTS(JGAS,JSM)
           II = IOBPTS(IGAS,ISM)
           IJ = IOBPTS(JGAS,JSM)
           IF(ISM.EQ.JSM) THEN
*
* One-electron part 
* =================
*
*. blocks of one-electron integrals and one-electron density
         
             CALL GETD1(dbl_mb(KLDEN),ISM,IGAS,ISM,JGAS,1)
             CALL GETH1(dbl_mb(KLINT),ISM,IGAS,ISM,JGAS)
             IF(NTEST.GE.100) THEN
               WRITE(6,*) ' Block of 1e integrals ISM,IGAS,JGAS',
     &                    ISM,IGAS,JGAS
               CALL WRTMAT(dbl_mb(KLINT),NI,NJ,NI,NJ)
               WRITE(6,*) ' Block of 1e density ISM,IGAS,JGAS',
     &                    ISM,IGAS,JGAS
               CALL WRTMAT(dbl_mb(KLDEN),NI,NJ,NI,NJ)
             END IF
             E1 = E1 + INPROD(dbl_mb(KLDEN),dbl_mb(KLINT),NI*NJ)
           END IF
*
* Two-electron part 
* =================
*
           IF(I12.EQ.2) THEN
           DO KSM = 1, NSMOB
*. Obtain LSM
             CALL  SYMCOM(3,1,IJSM,KSM,IJKSM)
             IJKLSM = 1
             CALL  SYMCOM(2,1,IJKSM,LSM,IJKLSM)
C?           WRITE(6,*) ' IJSM IJKSM LSM ',IJSM,IJKSM,IJKLSM
*
             DO KGAS = 1, NGAS
             DO LGAS = 1, NGAS
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
*. Blocks of density matrix and integrals 
                IXCHNG = 0
                ICOUL  = 1
                ONE = 1.0D0
                CALL LUCIAGETINT(dbl_mb(KLINT),
     &               IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
     &               IXCHNG,0,0,ICOUL,ONE,ONE)
                CALL GETD2 (dbl_mb(KLDEN),
     &               ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,1)
                IF (IMODE.NE.0) THEN
                  CALL GETD2RED(dbl_mb(KLDEN),
     &                 ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,1)
                END IF
C?              write(6,*) ' Ism Jsm Ksm Lsm' , Ism,Jsm,Ksm,Lsm
C?              write(6,*)
C?   &          ' Igas Jgas Kgas Lgas' , Igas,Jgas,Kgas,Lgas
C?              WRITE(6,*) ' Integral block'
C?              CALL WRTMAT(WORK(KLINT),NI*NJ,NK*NL,NI*NJ,NK*NL)
C?              WRITE(6,*) ' Density block '
C?              CALL WRTMAT(WORK(KLDEN),NI*NJ,NK*NL,NI*NJ,NK*NL)
                NIJKL = NI*NJ*NK*NL
                E2 = E2+0.5D0*INPROD(dbl_mb(KLDEN),dbl_mb(KLINT),NIJKL)
C?              write(6,*) ' Updated 2e-energy ', E2
             END DO
             END DO
           END DO
           END IF
*
          END DO
         END DO
       END DO
      END DO
*
      E = E1 + E2 + ECORE
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from EN_FROM_DEN' 
        WRITE(6,*)
        WRITE(6,*) ' One-electron energy ', E1
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-electron energy ', E2
        END IF
        WRITE(6,*)
        WRITE(6,*) ' Total energy : ', E
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'EN_FRM')
*
      RETURN
      END
      SUBROUTINE FOCK_MAT_NORT(F1,F2,I12,FI,FA)
*
* Construct the two Fock matrices needed to compute orbital gradient
* in nonorthogonal theory:
*
* F1(j,i) =   SUM(K)     H(i,K)     * RHO1(j,K)
*          + SUM(K,L,M) (i K L M ) * RHO2(j K L M)
*
* F2(j,i) =   SUM(K)     H(K,i)     * RHO1(K,j)
*          + SUM(M,K,L) (K i L M ) * RHO2(K j L M)
*
*
* The matrices F1 and F2 are on output in the actual basis, so transformations
* from the biobase are performed
*
* The matrices are first calculated in the mixed basis
* (creation operators are in the bio-base, annihilation in the
*  actual basis).
*
*
* F may be written as
*
* j: inactive: FX(j,i) = 2(FI(JI) + 2 FA(JA)), X = 1, 2
*
* j: active  : F1(j,i) = sum(k:active) D(j,k) FI(i,k)
*                      + sum(klm:active) d(jklm) (ik!lm)
* j: active  : F2(j,i) = sum(k:active) D(k,j) FI(k,i)
*                      + sum(klm:active) d(kjlm) (ki!lm)
* j: secondary:FX(j,i) = 0
*
*
* Modified from FOCK_MAT   
*
* Jeppe Olsen, June 2012 in Zurich
*
* Unless I12 = 2, only one-electron part is calculated
c      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      DIMENSION FI(*), FA(*)
*. Output
      DIMENSION F1(*), F2(*)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
        WRITE(6,*) ' Output from FOCK_MAT_NORT'
        WRITE(6,*) ' ------------------------------'
        WRITE(6,*)
        WRITE(6,*) ' Input inactive and active Fock matrices '
        CALL APRBLM2(FI,NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*)
        CALL APRBLM2(FA,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      ONE = 1.0D0
      ZERO = 0.0D0
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'FOONOR')
*
* The one-body density matrix is on input given in the actual basis,
* transform, so the first index is in the biorthonormal basis
*
      LEN_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      LEN_R = NACOB**2
     
      CALL MEMMAN(KLRHO1,LEN_R,'ADDL  ',2,'RHO1L ') !done
      CALL MEMMAN(KLRHO1B,LEN_R,'ADDL  ',2,'RHO1S ') !done
      CALL MEMMAN(KLRHO1_SAVE,LEN_R,'ADDL  ',2,'RHO1S ') !done
      CALL MEMMAN(KLCBIOA,LEN_R,'ADDL  ',2,'CBIOAC') !done
*
      CALL COPVEC(dbl_mb(KRHO1),dbl_mb(KLRHO1_SAVE),LEN_R)
*.  Obtain, in KLCBIOA, CBIO over active orbitals only
      CALL EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT
     &     (WORK(KCBIO),dbl_mb(KLCBIOA),1)
*. Obtain rho1 in symmetry block form
C            REORHO1(RHO1I,RHO1O,IRHO1SM)
      CALL REORHO1(dbl_mb(KRHO1),dbl_mb(KLRHO1),1,1)
*. transform RHO1 to bio-actual MO basis
      CALL TR_BIOMAT(dbl_mb(KLRHO1),dbl_mb(KLRHO1B),dbl_mb(KLCBIOA),
     &               NACOBS,1,2,1,1)
*. Transfer back to full matrix over active orbitals
      CALL REORHO1(dbl_mb(KRHO1),dbl_mb(KLRHO1B),1,2)
*
      CALL COPVEC(FI,WORK(KINT1),LEN_F)
*
      CALL SETVEC(F1,ZERO,LEN_F)
      CALL SETVEC(F2,ZERO,LEN_F)
*
*. Largest set of orbitals with given symmetry and type
      MXTSOB_AC = 0
      MXTSOB    = 0
      DO ISM = 1, NSMOB
        DO IGAS = 1, NGAS
          MXTSOB_AC = MAX(MXTSOB_AC,NOBPTS(IGAS,ISM))
        END DO
        MXTSOB = MAX(MXTSOB,MXTSOB_AC,NINOBS(ISM),NSCOBS(ISM))
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' MXTSOB, MXTSOB_AC = ', MXTSOB, MXTSOB_AC
      END IF
*. Allocate scratch space for 2-electron integrals and 
*. two-electron densities
      CALL MEMMAN(KLINT,MXTSOB_AC**3*MXTSOB,'ADDL  ',2,'KLINT ') !done
      CALL MEMMAN(KLDEN,MXTSOB_AC**4,'ADDL  ',2,'KLDEN ') !done
*. And  block of F1, F2
      MX2IBLK = MXTSOB ** 2
      CALL MEMMAN(KLF1BLK,LEN_F,'ADDL  ',2,'F1BL  ') !done
      CALL MEMMAN(KLF2BLK,MX2IBLK,'ADDL  ',2,'F2BL  ') !done
*
      ONE = 1.0D0
      II = -2303
      IJ = -2303
      DO IJSM = 1, NSMOB
        ISM = IJSM
        JSM = IJSM
        NIJS = NTOOBS(IJSM)
*
        IF(IJSM.EQ.1) THEN
         IFOFF = 1
        ELSE
         IFOFF = IFOFF+NTOOBS(IJSM-1)**2
        END IF
*
        DO JGAS = 0, NGAS
         IF(JGAS.EQ.0) THEN
           IJ = 1
         ELSE IF(JGAS.EQ.1) THEN
           IJ = NINOBS(JSM)+1
         ELSE
           IJ = IJ + NOBPTS(JGAS-1,JSM)
         END IF
         IF(JGAS.EQ.0) THEN
          NJ = NINOBS(IJSM)
         ELSE
          NJ = NOBPTS(JGAS,IJSM)
         END IF
          DO IGAS = 0, NGAS+1
           IF(NTEST.GE.1000) THEN
             WRITE(6,*) 
     &       ' Info for ISM IGAS JGAS',ISM,IGAS,JGAS
           END IF
*
           IF(IGAS.EQ.0) THEN
            NI = NINOBS(ISM)
           ELSE IF(IGAS.LE.NGAS) THEN
            NI = NOBPTS(IGAS,ISM)
           ELSE
            NI = NSCOBS(ISM)
           END IF
*
           IF(IGAS.EQ.0) THEN
             II = 1
           ELSE IF(IGAS.EQ.1) THEN 
             II = 1 + NINOBS(ISM)
           ELSE 
             II = II + NOBPTS(IGAS-1,ISM)
           END IF
           IF(NI*NJ.NE.0) THEN
*
*  =======================
*. block F(ijsm,jgas,igas)
*  =======================
*
            CALL SETVEC(dbl_mb(KLF1BLK),ZERO,NI*NJ)
            CALL SETVEC(dbl_mb(KLF2BLK),ZERO,NI*NJ)
            IF(JGAS.EQ.0) THEN
*
*. Inactive part
*
             DO I = 1, NI
              DO J = 1, NJ
*. Addres in sym-block
               IJS = (IJ-1+J-1)*NIJS + I+II-1
               JIS = (II-1+I-1)*NIJS + J+IJ-1
*. Address in sym-gas block
               IJT = (J-1)*NI + I
               JIT = (I-1)*NJ + J
  
               dbl_mb(KLF1BLK-1+JIT) = 2.0D0*
     &         (FI(IFOFF-1+IJS)+FA(IFOFF-1+IJS))
               dbl_mb(KLF2BLK-1+JIT) = 2.0D0*
     &         (FI(IFOFF-1+JIS)+FA(IFOFF-1+JIS))
               IF(NTEST.GE.1000) THEN
                 WRITE(6,'(A,3I3,3(1X,F13.5))') 
     &           ' IJSM, I, J, FI(IJ), FA(IJ), F(JI) =',
     &             IJSM, I, J, FI(IFOFF-1+IJ), FA(IFOFF-1+IJ),
     &           dbl_mb(KLF1BLK-1+JI)
              END IF
              END DO
             END DO
            ELSE IF(JGAS.LE.NGAS) THEN
*
*. Active part
*
* 1 : One-electron part 
             DO KGAS = 1, NGAS
              KSM = IJSM
              NK = NOBPTS(KGAS,KSM)
*
*. For F1:
*
              CALL GETD1(dbl_mb(KLDEN),JSM,JGAS,KSM,KGAS,1)
              CALL GETH1(dbl_mb(KLINT),ISM,IGAS,KSM,KGAS)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) 
     &          ' 1-e ints for ISM IGAS KGAS ',ISM,IGAS,KGAS
                CALL WRTMAT(dbl_mb(KLINT),NI,NK,NI,NK)
                WRITE(6,*) 
     &          ' 1-e densi for ISM JGAS KGAS ',ISM,JGAS,KGAS
                CALL WRTMAT(dbl_mb(KLDEN),NJ,NK,NJ,NK)
              END IF
              CALL MATML7(dbl_mb(KLF1BLK),dbl_mb(KLDEN),dbl_mb(KLINT),
     &                    NJ,NI,NJ,NK,NI,NK,ONE,ONE,2)
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated F1 block '
                 CALL WRTMAT(dbl_mb(KLF1BLK),NJ,NI,NJ,NI)
               END IF
*
*. For F2:
*
              CALL GETD1(dbl_mb(KLDEN),KSM,KGAS,JSM,JGAS,1)
              CALL GETH1(dbl_mb(KLINT),KSM,KGAS,ISM,IGAS)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) 
     &          ' 1-e ints for KSM KGAS IGAS ',KSM,KGAS,IGAS
                CALL WRTMAT(dbl_mb(KLINT),NK,NI,NK,NI)
                WRITE(6,*) 
     &          ' 1-e densi for KSM KGAS JGAS ',KSM,KGAS,JGAS
                CALL WRTMAT(dbl_mb(KLDEN),NK,NJ,NK,NJ)
              END IF
              CALL MATML7(dbl_mb(KLF2BLK),dbl_mb(KLDEN),dbl_mb(KLINT),
     &                    NJ,NI,NK,NJ,NK,NI,ONE,ONE,1)
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated F2 block '
                 CALL WRTMAT(dbl_mb(KLF2BLK),NJ,NI,NJ,NI)
               END IF
             END DO
*
             IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One-electron contributions to F1 and F2 '
              WRITE(6,*) ' ========================================'
              CALL WRTMAT(dbl_mb(KLF1BLK),NJ,NI,NJ,NI)
              WRITE(6,*)
              CALL WRTMAT(dbl_mb(KLF2BLK),NJ,NI,NJ,NI)
             END IF
             IF(I12.EQ.2) THEN
*. 2 : Two-electron part
             DO KSM = 1, NSMOB
             DO LSM = 1, NSMOB
*. Obtain MSM
              CALL  SYMCOM(3,1,KSM,LSM,KLSM)
              CALL  SYMCOM(3,1,KLSM,ISM,IKLSM)
              IMKLSM = 1
              CALL  SYMCOM(2,1,IKLSM,MSM,IMKLSM)
*
              DO MGAS = 1, NGAS
              DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
               ISKIP = 0
               XFACTOR = 1.0D0
*. Check that K,L,J are occupied and ISKIP = 0
                NM = NOBPTS(MGAS,MSM)
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
                IF(NM*NK*NL.NE.0) THEN
* 
* For F1: + sum(klm:active) d(lmjk) (lm!ik)
*
                  IF(NTEST.GE.1000) THEN
                   WRITE(6,*) ' For F1: '
                  END IF
*. Blocks of density matrix and integrals: D2(L M J K) (L M ! I K)
                  IXCHNG = 0
                  ICOUL  = 1
                  ONE = 1.0D0
* Obtain ( L M ! I K)
                  CALL LUCIAGETINT(dbl_mb(KLINT),
     &                 LGAS,LSM,MGAS,MSM,IGAS,ISM,KGAS,KSM,
     &                 IXCHNG,0,0,ICOUL,ONE,ONE)
*. Obtain Rho2(L M, J K)
                  CALL GETD2 (dbl_mb(KLDEN),
     &                 LSM,LGAS,MSM,MGAS,JSM,JGAS,KSM,KGAS,1)
                  NLM = NL*NM
                  DO K = 1, NK
                    IIOFF = (K-1)*NL*NM*NI + 1
                    IDOFF = (K-1)*NL*NM*NJ + 1
                    CALL MATML7(dbl_mb(KLF1BLK),
     &                   dbl_mb(KLDEN-1+IDOFF),dbl_mb(KLINT-1+IIOFF),
     &              NJ,NI,NLM,NJ,NLM,NI,ONE,ONE,1)
                  END DO
                  IF(NTEST.GE.1000) THEN
                    WRITE(6,*) ' Updated F1 block '
                    CALL WRTMAT(dbl_mb(KLF1BLK),NJ,NI,NJ,NI)
                  END IF
* 
* For F2: + sum(klm:active) d(lmkj) (lm!ki)
*
                  IF(NTEST.GE.1000) THEN
                   WRITE(6,*) ' For F2: '
                  END IF
*. Blocks of density matrix and integrals: D2(L M K J) (L M ! K I)
                  IXCHNG = 0
                  ICOUL  = 1
                  ONE = 1.0D0
* Obtain ( L M ! K I)
                  CALL LUCIAGETINT(dbl_mb(KLINT),
     &                 LGAS,LSM,MGAS,MSM,KGAS,KSM,IGAS,ISM,
     &                 IXCHNG,0,0,ICOUL,ONE,ONE)
*. Obtain Rho2(L M, K J)
                  CALL GETD2 (dbl_mb(KLDEN),
     &                 LSM,LGAS,MSM,MGAS,KSM,KGAS,JSM,JGAS,1)
                  NKLM = NK*NL*NM
                  CALL MATML7(dbl_mb(KLF2BLK),dbl_mb(KLDEN),
     &            dbl_mb(KLINT),NJ,NI,NKLM,NJ,NKLM,NI,ONE,ONE,1)
                END IF !End if nonvanishing block
              END DO
              END DO
              END DO
             END DO
             END DO
             END IF
             IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One- + two-electron contributions'
              WRITE(6,*) ' ================================='
              CALL WRTMAT(dbl_mb(KLF1BLK),NJ,NI,NJ,NI)
              WRITE(6,*)
              CALL WRTMAT(dbl_mb(KLF2BLK),NJ,NI,NJ,NI)
             END IF
            END IF
*           ^ End of inactive/active switch
*. Blocks has been constructed , transfer to -complete- 
*. symmetry blocked Fock matrix
*
*
            DO J = 1, NJ
              DO I = 1, NI
C?              WRITE(6,*) 'IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1',
C?   &                      IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1
                F1(IFOFF-1+(I+II-1-1)*NIJS + J+IJ-1 ) = 
     &          dbl_mb(KLF1BLK-1+(I-1)*NJ+J)
                F2(IFOFF-1+(I+II-1-1)*NIJS + J+IJ-1 ) = 
     &          dbl_mb(KLF2BLK-1+(I-1)*NJ+J)
              END DO
            END DO
*
           END IF
          END DO
        END DO
      END DO
*. Clean up
      CALL COPVEC(dbl_mb(KLRHO1_SAVE),dbl_mb(KRHO1),LEN_R)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' F1 and F2 before BIOTR '
        WRITE(6,*) ' ======================='
        CALL APRBLM2(F1,NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*)
        WRITE(6,*)
        CALL APRBLM2(F2,NTOOBS,NTOOBS,NSMOB,0)
      END IF
* The matrices F1 and F2 have now been obtained in the basis
* where the integrals have been calculated. For F1, the first
* index is in biobase, second index is in original basis. For F2 it is 
* the other way around...
*. Transform so both indeces are in the original basis
*
C     TR_BIOMAT(XIN,XOUT,CBIO,NORB_PSM,
C    &            INB_IN,INB_OUT,JNB_IN,JNB_OUT)
      CALL TR_BIOMAT(F1,dbl_mb(KLF1BLK),WORK(KCBIO),NTOOBS,
     &     2,1,1,1)
      CALL COPVEC(dbl_mb(KLF1BLK),F1,LEN_F)
*
      CALL TR_BIOMAT(F2,dbl_mb(KLF1BLK),WORK(KCBIO),NTOOBS,
     &     1,1,2,1)
      CALL COPVEC(dbl_mb(KLF1BLK),F2,LEN_F)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Final F1 and F2 from FOCK_MAT_STANDARD '
        WRITE(6,*) ' ======================================='
        CALL APRBLM2(F1,NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*)
        WRITE(6,*)
        CALL APRBLM2(F2,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'FOONOR')
      RETURN
      END
      SUBROUTINE GETSTR_ALLSM_GNSPGP
     &           (ISTRTP,NGRPA,IGRPA,NSTRPSM,IBSTRPSM,NEL,ISTR)
*
* Obtain all strings obtained as a product of NGAS groups
*
* =====
* Input 
* =====
*
* ISTRTP  : Type of of superstrings ( alpha => 1, beta => 2 )
* NGRPA  :  Number of Orbital spaces in general group
* IGRPA   : The groups of the NGASA spaces
* NEL    : Number of electrons in general supergroup
*
* ======
* Output 
* ======
*
* NSTRPSM : Number of superstrings generated
* IBSTRPSM : Occupation of superstring
* ISTR: The actual strings
*
*
* Jeppe Olsen, April 2013
*
*. Input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INTEGER IGRPA(NGRPA)
*. output
      INTEGER ISTR(NEL,*), NSTRPSM(NSMOB), IBSTRPSM(NSMOB)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      NTEST = 00
      CALL LUCIAQENTER('GTSTAS')
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_ALLSM_GNSGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,*) ' Number of orbital groups in general group', 
     &  NGRPA
        WRITE(6,*) ' And the groups '
        CALL IWRTMA3(IGRPA, 1, NGRPA, 1, NGRPA)
        WRITE(6,*)
        WRITE(6,*) ' NEL = ', NEL
      END IF
*
* Number of strings per sym and offset for given general supergroup
*
C     GET_DIM_GNSPGP(NGRPA,IGRPA,NSTPSM,IBSTPSM)
      CALL GET_DIM_GNSPGP(NGRPA,IGRPA,NSTRPSM,IBSTRPSM)
* We now have offset to address of strings of given sym. 
*.Zero NSTRPSM as this also will be used to pointer for start
* of strings relattive to offset (did you get that Jeppe?)
      IZERO = 0
      CALL ISETVC(NSTRPSM,IZERO,NSMST)
*
*
*. Largest and lowest symmetries active in each GAS space
*
      NEL = 0
      DO IGRP = 1, NGRPA
        NELFGS(IGRP) = NELFGP(IGRPA(IGRP))          
        NEL = NEL + NELFGS(IGRP)
        IF(NELFGS(IGRP).GT.0) NGRPL = IGRP
      END DO
      IF(NGRPL.EQ.0) NGRPL = 1
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Number of elecs per Group '
         CALL IWRTMA(NELFGS,1,NGRPA,1,NGRPA)
      END IF 
*
*. Number of strings per Group and offsets for strings of given sym
*
      DO IGRP = 1, NGRPA
C?      WRITE(6,*) ' (IGRPA(IGRP)-1)*NSMST+1 = ',
C?   &               (IGRPA(IGRP)-1)*NSMST+1
        CALL ICOPVE2(int_mb(KNSTSGP(1)),(IGRPA(IGRP)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGRP))
        CALL ICOPVE2(int_mb(KISTSGP(1)),(IGRPA(IGRP)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGRP))
      END DO
*
* Largest and lowest active symmetries for each GAS space
*
      DO IGRP = 1, NGRPA
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGRP).GT.0) MAXVAL(IGRP) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGRP).GT.0) MINVAL(IGRP) = ISMST
        END DO
      END DO
      IF(NTEST.GE.1000) THEN
      WRITE(6,*) 'The MINVAL and MAXVAL arrays '
        CALL IWRTMA3(MINVAL,1,NGRPA,1,NGRPA)
        CALL IWRTMA3(MAXVAL,1,NGRPA,1,NGRPA)
      END IF
*
*. Loop over symmetries of each GRP
*
      IFIRST = 1
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGRP = 1, NGRPA 
            ISMFGS(IGRP) = MINVAL(IGRP)
          END DO
        ELSE
*. Next distribution of symmetries in NGRPA
         CALL NXTNUM3(ISMFGS,NGRPA,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGRPA spaces '
          CALL IWRTMA(ISMFGS,NGRPA,1,NGRPA,1)
        END IF
*. Symmetry of NGRPA spaces given, symmetry of total space
        ISTSM = 1
        DO IGRP = 1, NGRPA
          CALL  SYMCOM(3,1,ISTSM,ISMFGS(IGRP),JSTSM)
          ISTSM = JSTSM
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGRPA,1,NGRPA)
          WRITE(6,*) ' Total symmetry = ', ISTSM
        END IF
*. Obtain all strings of this symmetry distribution
CT      CALL LUCIAQENTER('GASSM')
        IB = IBSTRPSM(ISTSM)+NSTRPSM(ISTSM)
        CALL GETSTRN_GASSM_GNSPGP(ISMFGS,IGRPA,ISTR(1,IB),NSTR,NEL,
     &                          NNSTSGP,IISTSGP,NGRPA)
        NSTRPSM(ISTSM) = NSTRPSM(ISTSM) + NSTR
C     GETSTRN_GASSM_GNSPGP(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
C    &                          NNSTSGP,IISTSGP,NGRPL)
CT      CALL LUCIAQEXIT('GASSM')
*. ready for next symmetry distribution 
        GOTO 1000
 1001 CONTINUE
*. End of loop over symmetry distributions
*
      IF(NTEST.GE.100) THEN
       DO ISM = 1, NSMST
        IB = IBSTRPSM(ISM)
        N  = NSTRPSM(ISM)
        WRITE(6,*) 
     &  ' Symmetry and number of strings generated ', ISM,N
        WRITE(6,*)
        WRITE(6,*) ' Strings: '
        WRITE(6,*)
        CALL PRTSTR(ISTR(1,IB),NEL,N)
       END DO
      END IF
*
      CALL LUCIAQEXIT('GTSTAS')
      RETURN
      END 
      SUBROUTINE GET_DIM_GNSPGP(NGRPA,IGRPA,NSTPSM,IBSTPSM)
*
* A general supergroup is specified by NGRPA groups IGRPA
* Find number of strings in this supergroup per symmetry and 
* offsets to given symmetry
*
*. Jeppe Olsen, April 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'strbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*. Input
      INTEGER IGRPA(NGRPA)
*. Output
      INTEGER NSTPSM(NSMST), IBSTPSM(NSMST)
*. Local scratch
      INTEGER NSTFTP(MXPNSMST), NSTFTP2(MXPNSMST)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_DIM_GNSPGP'
      END IF
*
      IZERO = 0
*
      CALL ISETVC(NSTPSM,IZERO,NSMST)
      NSTPSM(1) = 1
*
      DO IGRP = 1, NGRPA
        CALL ICOPVE2(int_mb(KNSTSGP(1)),(IGRPA(IGRP)-1)*NSMST+1,NSMST,
     &               NSTFTP)
        CALL ICOPVE(NSTPSM,NSTFTP2,NSMST)
        CALL DIM_PROD_TWO_STGRPS(NSTPSM,NSTFTP2,NSTFTP,NSMST)
      END DO
*
      IBSTPSM(1) = 1
      DO ISM = 2, NSMST
        IBSTPSM(ISM) = IBSTPSM(ISM-1)+NSTPSM(ISM-1)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Obtained number of strings per sym:'
        CALL IWRTMA(NSTPSM,1,NSMST,1,NSMST)
        WRITE(6,*) ' Offset to strings with given sym: '
        CALL IWRTMA(IBSTPSM,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_PROD_TWO_STGRPS(N12PSM,N1PSM,N2PSM,NSM)
*
* The number of strings per symmetry for two groups are given in
* N1PSM, N2PSM
* Find number of strings per sym of product of the spaces
*
* D2H hardwired pt
*
*. Jeppe Olsen, April 29, 2013
*
      INCLUDE 'implicit.inc' 
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER N1PSM(NSM),N2PSM(NSM)
*. Output
      INTEGER N12PSM(NSM)
*
      DO I12SM = 1, NSM
       N12PSM(I12SM) = 0
       DO I1SM = 1, NSM
         I2SM = MULTD2H(I12SM,I1SM)
         N12PSM(I12SM) = N12PSM(I12SM) + N1PSM(I1SM)*N2PSM(I2SM)
       END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GETSTRN_GASSM_GNSPGP(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
     &                          NNSTSGP,IISTSGP,NGRPL)
*
* Obtain all superstrings containing  strings of given sym and type 
* for a general supergroup
*
* General?: A general number of groups, NGRPL may be used
*
* ( Superstring :contains electrons belonging to all gasspaces  
*        string :contains electrons belonging to a given GAS space
* A super string is thus a product of NGAS strings )
*
* Jeppe Olsen, April 29, 2013, simple modification of GETSTRN_GASSM_SPGP
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input 
      INTEGER ITPFGS(NGRPL), ISMFGS(NGRPL)
      INTEGER NNSTSGP(MXPNSMST,NGRPL), IISTSGP(MXPNSMST,NGRPL)
*. Local scratch 
      INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS)
*. Output 
      INTEGER ISTROC(NEL,*)
*. Number of strings per GAS space
      DO IGRP = 1, NGRPL
        NSTFGS(IGRP)  = NNSTSGP(ISMFGS(IGRP),IGRP)
        IBSTFGS(IGRP) = IISTSGP(ISMFGS(IGRP),IGRP)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GETSTR_GASSM_GNSPGP speaking '
        WRITE(6,*) '  =========================== '
        WRITE(6,*) ' ISMFGS,ITPFGS (input) '
        CALL IWRTMA(ISMFGS,1,NGRPL,1,NGRPL)
        CALL IWRTMA(ITPFGS,1,NGRPL,1,NGRPL)
        WRITE(6,*)
        WRITE(6,*) ' NSTFGS, IBSTFGS ( intermediate results ) '
        CALL IWRTMA(NSTFGS,1,NGRPL,1,NGRPL)
        CALL IWRTMA(IBSTFGS,1,NGRPL,1,NGRPL)
      END IF
*
      NSTRTOT = 1
      DO IGRP = 1, NGRPL
        NSTRTOT = NSTRTOT*NSTFGS(IGRP)
      END DO
C     WRITE(6,*) ' NSTRTOT = ', NSTRTOT
      IF(NGRPL.EQ.0) GOTO 2810
*
*
      IF(NSTRTOT.EQ.0) GOTO 1001
*. Loop over GAS spaces
      DO IGRP = 1, NGRPL
*. Number of electrons in GRP = 1, IGRP - 1
        IF(IGRP.EQ.1) THEN
          NELB = 0
        ELSE
          NELB = NELB +  NELFGP(ITPFGS(IGRP-1))
        END IF
*. Number of electron in IGRP
        NELI = NELFGP(ITPFGS(IGRP))
C?      WRITE(6,*) ' IGRP, NELI and NELB ', IGRP, NELI,NELB
        IF(NELI.GT.0) THEN
         
*. The order of strings corresponds to a matrix A(I(after),Igas,I(before))
*. where I(after) loops over strings in IGAS+1 - IGASL and
*  I(before) loop over strings in 1 - IGAS -1
          NSTA = 1
          DO JGRP = IGRP+1, NGRPL
            NSTA = NSTA * NSTFGS(JGRP)
          END DO
*
          NSTB =  1
          DO JGRP = 1, IGRP-1
            NSTB = NSTB * NSTFGS(JGRP)
          END DO
*
          NSTI = NSTFGS(IGRP)
         
C?        write(6,*) ' before call to add_str_group '
          IF(NTEST.GE.200) THEN
            WRITE(6,*) ' NSTI,NSTB,NSTA,NELB,NELI,NEL ',
     &                   NSTI,NSTB,NSTA,NELB,NELI,NEL
            WRITE(6,*) ' IBSTFGS(IGRP),KOC()',
     &                   IBSTFGS(IGRP),KOCSTR(ITPFGS(IGRP))
          END IF
*
          CALL ADD_STR_GROUP(NSTI,
     &          IBSTFGS(IGRP),
     &          int_mb(KOCSTR(ITPFGS(IGRP))),
     &          NSTB,NSTA,ISTROC,NELB+1,NELI,NEL)
*. Loop over strings in IGRP 
        END IF
      END DO
 1001 CONTINUE
 2810 CONTINUE
      NSTR = NSTRTOT  
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from  GETSTR_GASSM_SPGP ' 
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and type strings : '
        WRITE(6,*)
        WRITE(6,*) '   GP    Sym  Type '
        WRITE(6,*) ' =================='
        DO IGRP = 1, NGRPL
          WRITE(6,'(3I6)') IGRP,ISMFGS(IGRP),ITPFGS(IGRP)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Number of strings generated : ', NSTR
        WRITE(6,*) ' Strings generated '
        CALL PRTSTR(ISTROC,NEL,NSTR)
      END IF
*
      RETURN
      END 
