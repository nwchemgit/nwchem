      SUBROUTINE GET_SUPSYM_INFO
*
* Obtain supersymmetry info
*
*. Jeppe Olsen, May 2012
*               Revised July 2012
*  Last revision; Jeppe Olsen; March 5, 2013; call to GET_SUPSYM_FOR_BASIS changed to include
*                                             irrep info
*
* Last modified July 8, 2012 (Jeppe)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTSPSM')
*. Parity for standard symmetry 
      IF(CSUPSYM.EQ.'LINEAR'.AND.INVCNT.EQ.1) THEN
        CALL SET_PARITY_FOR_STASYM(IPA_FOR_STASYM,NIRREP)
      ELSE
        IZERO = 0
        CALL ISETVC(IPA_FOR_STASYM,IZERO,NIRREP)
      END IF
* Labels for the basis functions
      CALL GET_SUPSYM_LABELS_FOR_ORBITALS
     &(int_mb(KCSUPSYM_FOR_ORB),int_mb(KLVAL_FOR_ORB),
     & int_mb(KMLVAL_FOR_ORB),int_mb(KPA_FOR_ORB))
*. Relation between supersymmetry and irreps
      CALL SYM_AND_IRREP_FOR_SUPSYM(
     & int_mb(KL_FOR_SUPSYM),int_mb(KML_FOR_SUPSYM),
     & int_mb(KPA_FOR_SUPSYM),int_mb(KIRREP_FOR_SUPSYM),
     & int_mb(KNSUPSYM_FOR_IRREP), int_mb(KIBSUPSYM_FOR_IRREP),
     & int_mb(KISUPSYM_FOR_IRREP) )
C     SYM_AND_IRREP_FOR_SUPSYM
C    &           (L_FOR_SUPSYM,ML_FOR_SUPSYM,IPA_FOR_SUPSYM,IRREP_FOR_SUPSYM,
C    &           NSUPSYM_FOR_IRREP, IBSUPSYM_FOR_IRREP,ISUPSYM_FOR_IRREP)
*
*. Info on the supersymmetry of the basis set
*
      CALL GET_SUPSYM_FOR_BASIS(int_mb(KISUPSYM_FOR_BAS),
     &     int_mb(KNBAS_FOR_SUP_STA_SYM),int_mb(KIBBAS_FOR_SUP_STA_SYM),
     &     int_mb(KIBAS_FOR_SUP_STA_SYM),int_mb(KISHELL_FOR_BAS), 
     &     int_mb(KNBAS_FOR_SHELL),int_mb(KIBBAS_FOR_SHELL),
     &     int_mb(KIBAS_FOR_SHELL),
     &     int_mb(KNSUPSYM_FOR_IRREP),int_mb(KIBSUPSYM_FOR_IRREP),
     &     int_mb(KISUPSYM_FOR_IRREP))
*
* Obtain mappings from standard to super symmetry components
* simply by reading info in NBAS_FOR_SUP_STA_SYM
      CALL SUP_TO_STASYM(int_mb(KNBAS_FOR_SUP_STA_SYM))
*. Number, offsets and actual numbers of orbitals of 
*. given symmetry
*
* The GENSMOB arrays allowing a symmetric treatment of standard 
* and super symmetry
*. Just let the supersymmetry components be general symmetry
      NGENSMOB= N_SUPSYM
      CALL ICOPVE(NBAS_SUPSYM,NBAS_GENSMOB,NGENSMOB)
      CALL ICOPVE(IBBAS_SUPSYM,IBBAS_GENSMOB,NGENSMOB)
      CALL ICOPVE(IBAS_SUPSYM,ISTA_TO_GENSM_REO,NTOOB)
      CALL ICOPVE(ISTASM_FOR_SUPSYM,ISTASM_FOR_GENSM,NGENSMOB)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTSPSM')
      RETURN
      END
      SUBROUTINE GET_MAX_SUPSYM_IRREP
*
* Obtain Max supersymmetry Irrep from read of labels of orbitals
*. (saved in ORBINP)
*
* Jeppe Olsen, May 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
*
*. General input
      CHARACTER*4 AO_CENT, AO_TYPE, CHAR4
      COMMON/AOLABELS/AO_CENT(MXPORB),AO_TYPE(MXPORB)
*. A label is in general of the form nlm, where l is the 
*. standard spectroscopic term
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*)  ' =============================='
        WRITE(6,*)  ' Info from GET_MAX_SUPSYM_IRREP'
        WRITE(6,*)  ' =============================='
        WRITE(6,*)
      END IF
*. Obtained largest L-value
      LMAX = 0
      DO IORB = 1, NTOOB
       CHAR4 = AO_TYPE(IORB)
       IF(NTEST.GE.1000) WRITE(6,'(A,A)') ' AO_TYPE = ', AO_TYPE(IORB)
       IF(CHAR4(2:2).EQ.'s') THEN
         LMAX = MAX(0,LMAX)
       ELSE IF(CHAR4(2:2).EQ.'p') THEN
         LMAX = MAX(1,LMAX)
       ELSE IF(CHAR4(2:2).EQ.'d') THEN
         LMAX = MAX(2,LMAX)
       ELSE IF(CHAR4(2:2).EQ.'f') THEN
         LMAX = MAX(3,LMAX)
       ELSE IF(CHAR4(2:2).EQ.'g') THEN
         LMAX = MAX(4,LMAX)
       ELSE IF(CHAR4(2:2).EQ.'h') THEN
         LMAX = MAX(5,LMAX)
       ELSE IF(CHAR4(2:2).EQ.'i') THEN
         LMAX = MAX(6,LMAX)
       ELSE
         WRITE(6,*) ' Unknown form of orbital ', CHAR4
         STOP  ' Unknown form of orbital '
       END IF
      END DO
*
C?    IF(NTEST.GE.100) WRITE(6,*) ' Lmax = ', LMAX
*
      IF(CSUPSYM.EQ.'ATOMIC') THEN
        N_SUPSYM = LMAX + 1 + LMAX*(LMAX+1)
        N_SUPSYM_IRREP = LMAX +1
      ELSE IF (CSUPSYM.EQ.'LINEAR') THEN
        N_SUPSYM_IRREP = (INVCNT+1)*(LMAX + 1)
        N_SUPSYM = (INVCNT+1)*(2*LMAX + 1)
      END IF
      LMAX_ORB = LMAX
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' N_SUPSYM_IRREP = ', N_SUPSYM_IRREP
        WRITE(6,*) ' N_SUPSYM       = ', N_SUPSYM
        WRITE(6,*) ' Lmax           = ', LMAX
        WRITE(6,*) ' INVCNT         = ', INVCNT
      END IF
*
      RETURN
      END
      SUBROUTINE GET_SUPSYM_LABELS_FOR_ORBITALS
     &           (CSUPSYM_FOR_ORB,LVAL_FOR_ORB,MLVAL_FOR_ORB,
     &            IPA_FOR_ORB)
*
* Obtain labels for symmetry for the various basis functions
*
*. Jeppe Olsen, early morning May 23, 2012 (when I should do other things...)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
*. Output
      CHARACTER*3 CSUPSYM_FOR_ORB(NTOOB)
      INTEGER LVAL_FOR_ORB(NTOOB)
      INTEGER MLVAL_FOR_ORB(NTOOB)
      INTEGER IPA_FOR_ORB(NTOOB)
*
*
      CHARACTER*4 AO_CENT, AO_TYPE, CHAR4, CHAR4B
      COMMON/AOLABELS/AO_CENT(MXPORB),AO_TYPE(MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ========================================'
        WRITE(6,*) ' Info from GET_SUPSYM_LABELS_FOR_ORBITALS'
        WRITE(6,*) ' ========================================'
        WRITE(6,*) 
      END IF
* 
      DO IORB = 1, NTOOB
        CHAR4 = AO_TYPE(IORB)
        IF(CSUPSYM.EQ.'ATOMIC') THEN
          CSUPSYM_FOR_ORB(IORB) = CHAR4(1:2)
        ELSE IF (CSUPSYM.EQ.'LINEAR') THEN
          IF(CHAR4(1:2).EQ.'1s') THEN
            CSUPSYM_FOR_ORB(IORB) = '0 '
            CHAR4B = '0    '
          ELSE
           CHAR4B = '    '
           CHAR4B(1:2) = CHAR4(3:4)
           CSUPSYM_FOR_ORB(IORB) = CHAR4B(1:2)
          END IF
*. And parity labels
          IF(INVCNT.EQ.1) THEN
            IF(IPA_FOR_STASYM(ISMFSO(IORB)).EQ.1.) THEN
             CHAR4B(3:3) = 'g'
             IPA_FOR_ORB(IORB) = 1
            ELSE
             CHAR4B(3:3) = 'u'
             IPA_FOR_ORB(IORB) =-1
            ENDIF
            CSUPSYM_FOR_ORB(IORB) = CHAR4B(1:3)
          END IF !parity is present
        END IF
        CALL ORB_LABEL_TO_LML
     &  (AO_TYPE(IORB),LVAL_FOR_ORB(IORB),MLVAL_FOR_ORB(IORB))
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Labels for basis functions '
        WRITE(6,*) ' ==========================='
        DO IORB = 1, NTOOB
          WRITE(6,'(2X,I4,2X,A3)') IORB, CSUPSYM_FOR_ORB(IORB)
        END DO
*
        IF(INVCNT.EQ.0) THEN
         WRITE(6,*) ' L and ML-values for orbitals '
         WRITE(6,*) ' ============================='
         WRITE(6,*)
         WRITE(6,*) ' Orbital L   Ml '
         WRITE(6,*) ' ==============='
         DO IORB = 1, NTOOB
           WRITE(6,'(2X,I3,4X,I2,2X,I2)')
     &     IORB,LVAL_FOR_ORB(IORB),MLVAL_FOR_ORB(IORB)
         END DO
        ELSE
         WRITE(6,*) ' L, ML, and parity for orbitals '
         WRITE(6,*) ' ============================='
         WRITE(6,*)
         WRITE(6,*) ' Orbital L   Ml Parity'
         WRITE(6,*) ' ====================='
         DO IORB = 1, NTOOB
           WRITE(6,'(2X,I3,4X,I2,2X,I2,3X,I3)')
     &     IORB,LVAL_FOR_ORB(IORB),MLVAL_FOR_ORB(IORB),
     &     IPA_FOR_ORB(IORB)
         END DO
        ENDIF! explicit center of inversion 
*
      END IF
*
      RETURN
      END
      SUBROUTINE ORB_LABEL_TO_LML(ORB_LABEL,LVAL,MLVAL)
*
* A char*4 lavel ORB_LABEL is given, obtain 
* corresponding values of L and ML
*
* Very unelegant routine, but gets the work done (I hope)
*
* Jeppe Olsen, May 2012
*
*. Input
      CHARACTER*4 ORB_LABEL
*
      NTEST = 00
*
*. Lvalue
*
      IF(ORB_LABEL(2:2).EQ.'s') THEN
        LVAL = 0
      ELSE IF( ORB_LABEL(2:2).EQ.'p') THEN
        LVAL = 1
      ELSE IF( ORB_LABEL(2:2).EQ.'d') THEN
        LVAL = 2
      ELSE IF(ORB_LABEL(2:2).EQ.'f') THEN
        LVAL = 3
      ELSE IF(ORB_LABEL(2:2).EQ.'g') THEN
        LVAL = 4
      ELSE IF(ORB_LABEL(2:2).EQ.'h') THEN
        LVAL = 5
      ELSE IF(ORB_LABEL(2:2).EQ.'i') THEN
        LVAL = 6
      ELSE 
        WRITE(6,*) ' Unknown ORB_LABEL: ', ORB_LABEL
        STOP       ' Unknown ORB_LABEL '
      END IF
*
*. Ml-value
*
      IF(ORB_LABEL(2:2).EQ.'s') THEN
        MLVAL = 0
      ELSE IF(ORB_LABEL(2:3).EQ.'px') THEN
        MLVAL = 1
      ELSE IF(ORB_LABEL(2:3).EQ.'py') THEN
        MLVAL = -1
      ELSE IF(ORB_LABEL(2:3).EQ.'pz') THEN
        MLVAL = 0
      ELSE IF(ORB_LABEL(3:3).EQ.'0') THEN
        MLVAL = 0
      ELSE IF(ORB_LABEL(3:4).EQ.'1+') THEN
        MLVAL = 1
      ELSE IF(ORB_LABEL(3:4).EQ.'1-') THEN
        MLVAL = -1
      ELSE IF(ORB_LABEL(3:4).EQ.'2+') THEN
        MLVAL = 2
      ELSE IF(ORB_LABEL(3:4).EQ.'2-') THEN
        MLVAL = -2
      ELSE IF(ORB_LABEL(3:4).EQ.'3+') THEN
        MLVAL = 3
      ELSE IF(ORB_LABEL(3:4).EQ.'3-') THEN
        MLVAL = -3
      ELSE IF(ORB_LABEL(3:4).EQ.'4+') THEN
        MLVAL = 4
      ELSE IF(ORB_LABEL(3:4).EQ.'4-') THEN
        MLVAL = -4
      ELSE IF(ORB_LABEL(3:4).EQ.'5+') THEN
        MLVAL = 5
      ELSE IF(ORB_LABEL(3:4).EQ.'5-') THEN
        MLVAL = -5
      ELSE IF(ORB_LABEL(3:4).EQ.'6+') THEN
        MLVAL = 6
      ELSE IF(ORB_LABEL(3:4).EQ.'6-') THEN
        MLVAL = -6
      ELSE 
        WRITE(6,*) ' Unknown ORB_LABEL: ', ORB_LABEL
        STOP       ' Unknown ORB_LABEL: '
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,A,4X,I2,2X,I2)')
     &  ' ORB_LABEL, LVAL, MLVAL ', ORB_LABEL,LVAL,MLVAL
      END IF
*
      RETURN
      END
      SUBROUTINE SYM_AND_IRREP_FOR_SUPSYM(
     &           L_FOR_SUPSYM,ML_FOR_SUPSYM,
     &           IPA_FOR_SUPSYM,IRREP_FOR_SUPSYM,
     &           NSUPSYM_FOR_IRREP, IBSUPSYM_FOR_IRREP,
     &           ISUPSYM_FOR_IRREP)
*. Obtain tables L_FOR_SUPSYM,ML_FOR_SUPSYM,IPA_FOR_SUPSYM
*  giving L, ML and parity for symmetry
*. for linear molecules, L is set to zero
*  Obtain IRREP_FOR_SUPSYM giving irrep for symmetry 
*  Obtain NSUPSYM_FOR_IRREP, IBSUPSYM_FOR_IRREP,ISUPSYM_FOR_IRREP giving symmetry for irrep
*  
* Irrep: well irrep, say s,p,d for ATOMIC case, !Lambda! for LINEAR,
*        |Lambda| g for linear with inversion
* symmetry: given component of irrep, say p+1 for atomic, Pi+1 for linear
*
*. Jeppe Olsen, May 23, 2012, inversion added June11
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
*. Output
      INTEGER L_FOR_SUPSYM(*), ML_FOR_SUPSYM(*),IPA_FOR_SUPSYM(*)
      INTEGER IRREP_FOR_SUPSYM(*)
      INTEGER NSUPSYM_FOR_IRREP(*), IBSUPSYM_FOR_IRREP(*) 
      INTEGER ISUPSYM_FOR_IRREP(*)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================================='
        WRITE(6,*) ' Info from SYM_AND_IRREP_FOR_SUPSYM '
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
        WRITE(6,*) ' LMAX_ORB = ', LMAX_ORB
      END IF
*
      IF(CSUPSYM.EQ.'ATOMIC') THEN
        ISYM = 0
        IRREP = 0
        DO L = 0,  LMAX_ORB
          IRREP = IRREP + 1
          NDEG = 2*L + 1
          NSUPSYM_FOR_IRREP(IRREP) = NDEG
          IBSUPSYM_FOR_IRREP(IRREP) = ISYM + 1
          DO ML = -L,L
            ISYM = ISYM + 1
            L_FOR_SUPSYM(ISYM) = L
            ML_FOR_SUPSYM(ISYM) = ML
            ISUPSYM_FOR_IRREP(ISYM) = ISYM
            IRREP_FOR_SUPSYM(ISYM) = IRREP
          END DO
        END DO
      ELSE IF (CSUPSYM.EQ.'LINEAR') THEN
        ISYM = 0
        IRREP = 0
        IF(INVCNT.EQ.0) THEN
         NPARITY = 1
        ELSE
         NPARITY = 2
        END IF
        DO LAMBDA = 0, LMAX_ORB
         DO IPARITY = 1, NPARITY
          IRREP = IRREP + 1
          IF(LAMBDA.EQ.0) THEN
            NDEG = 1
          ELSE 
            NDEG = 2
          END IF
          NSUPSYM_FOR_IRREP(IRREP) = NDEG
          IBSUPSYM_FOR_IRREP(IRREP) = ISYM + 1
          DO ICOMP = 1, NDEG
            ISYM = ISYM + 1
            IRREP_FOR_SUPSYM(ISYM) = IRREP
            ISUPSYM_FOR_IRREP(ISYM) = ISYM
            L_FOR_SUPSYM(ISYM) = 0
            IF(NDEG.EQ.1) THEN
             ML_FOR_SUPSYM(ISYM) = 0
            ELSE
             IF(ICOMP.EQ.1) THEN
               ML_FOR_SUPSYM(ISYM) = -LAMBDA
             ELSE
               ML_FOR_SUPSYM(ISYM) = LAMBDA
             END IF
            END IF
            IF(INVCNT.EQ.0) THEN
             IPA_FOR_SUPSYM(ISYM) = 0
            ELSE
             IF(IPARITY.EQ.1) THEN
              IPA_FOR_SUPSYM(ISYM) = 1
             ELSE
              IPA_FOR_SUPSYM(ISYM) =-1
             END IF
            END IF! Parity present
          END DO ! loop over ICOMP
         END DO! Loop over parity
        END DO !loop over lambda
      END IF! Switch of SUPSYM
*
      IF(NTEST.GE.100) THEN
       IF(CSUPSYM.EQ.'ATOMIC') THEN
         WRITE(6,*) ' L and ML values for symmetries '
         WRITE(6,*) ' ============================== '
         WRITE(6,*)
         WRITE(6,*) ' Symmetry   L    Ml   '
         WRITE(6,*) ' ====================='
         DO ISYM = 1, N_SUPSYM
           WRITE(6,'(4X,I2,6X,I2,4X,I2)') 
     &     ISYM, L_FOR_SUPSYM(ISYM), ML_FOR_SUPSYM(ISYM)
         END DO
       ELSE IF(CSUPSYM.EQ.'LINEAR') THEN
         IF(INVCNT.EQ.0) THEN
          WRITE(6,*) ' ML values for symmetries '
          WRITE(6,*) ' ============================== '
          WRITE(6,*)
          WRITE(6,*) ' Symmetry   Ml   '
          WRITE(6,*) ' ====================='
          DO ISYM = 1, N_SUPSYM
            WRITE(6,'(3X,I2,7X,I2)') 
     &      ISYM, ML_FOR_SUPSYM(ISYM)
          END DO
         ELSE
          WRITE(6,*) ' Ml and parity for symmetries '
          WRITE(6,*) ' ============================== '
          WRITE(6,*)
          WRITE(6,*) ' Symmetry   Ml, parity   '
          WRITE(6,*) ' ====================='
          DO ISYM = 1, N_SUPSYM
            WRITE(6,'(3X,I2,7X,I2,5X,I2)') 
     &      ISYM, ML_FOR_SUPSYM(ISYM),IPA_FOR_SUPSYM(ISYM)
          END DO
         END IF ! parity is considered
       END IF ! form of supersymmetry
*
       WRITE(6,*) ' Irrep for the various symmetries '
       WRITE(6,*) ' ================================ '
       WRITE(6,*) 
       WRITE(6,*) ' Symmetry   Irrep  '
       WRITE(6,*) ' ================= '
       DO ISYM = 1, N_SUPSYM
         WRITE(6,'(5X,I2,7X,I2)')
     &   ISYM, IRREP_FOR_SUPSYM(ISYM)
       END DO
*
       WRITE(6,*) ' Symmetry for the various irreps '
       WRITE(6,*) ' ================================'
       WRITE(6,*)
       DO IRREP = 1, N_SUPSYM_IRREP
         WRITE(6,*) ' Supersymmetries for IRREP = ', IRREP
         IB = IBSUPSYM_FOR_IRREP(IRREP)
         NDEG = NSUPSYM_FOR_IRREP(IRREP)
         CALL IWRTMA(ISUPSYM_FOR_IRREP(IB),1,NDEG,1,NDEG)
       END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE SYM_AND_IRREP_FOR_LMLPA(L,ML,IPA,ISYM,IRREP)
*
* An orbital has label L and ML, and parity IPA
* obtain corresponding symmetry and irrep number
*
*. Jeppe Olsen, May 2012
*               June 2012: parity added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      NTEST = 00
*
      CALL SYM_AND_IRREP_FOR_LMLPA_S(L,ML,IPA,ISYM,IRREP,
     &     CSUPSYM,int_mb(KL_FOR_SUPSYM),int_mb(KML_FOR_SUPSYM),
     &     int_mb(KPA_FOR_SUPSYM),
     &     int_mb(KIRREP_FOR_SUPSYM),N_SUPSYM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from  SYM_AND_IRREP_FOR_LML '
        WRITE(6,*) ' Input: L and ML ', L, ML
        WRITE(6,*) ' Output: ISYM, IRREP ', ISYM, IRREP
      END IF
*
      RETURN
      END
C     SYM_AND_IRREP_FOR_LMLPA_S(L,ML,IPA,ISYM,IRREP,
C    &     CSUPSYM,WORK(KL_FOR_SUPSYM),WORK(KML_FOR_SUPSYM),
C    &     WORK(KIRREP_FOR_SUPSYM))
      SUBROUTINE SYM_AND_IRREP_FOR_LMLPA_S(L,ML,IPA,ISYM,IRREP,
     &     CSUPSYM,L_FOR_SUPSYM,ML_FOR_SUPSYM,
     &     IPA_FOR_SUPSYM,
     &     IRREP_FOR_SUPSYM,N_SUPSYM)
*
* Obtain supersymmetry and irrep for orbital with given L and ML,
* and  parity
*
* Jeppe Olsen, May 23, 2012 
*
      INCLUDE 'implicit.inc'
      CHARACTER*6 CSUPSYM
*. General input
      INTEGER L_FOR_SUPSYM(N_SUPSYM),ML_FOR_SUPSYM(N_SUPSYM)
      INTEGER IPA_FOR_SUPSYM(N_SUPSYM)
      INTEGER IRREP_FOR_SUPSYM(N_SUPSYM)
*
      NTEST = 00
      ISYM = 0   
      IF(CSUPSYM.EQ.'ATOMIC') THEN
       DO ISUPSYM = 1, N_SUPSYM
        IF(L_FOR_SUPSYM(ISUPSYM).EQ.L.AND.
     &     ML_FOR_SUPSYM(ISUPSYM).EQ.ML ) THEN
           ISYM = ISUPSYM
        END IF
       END DO
      ELSE IF(CSUPSYM.EQ.'LINEAR') THEN
       DO ISUPSYM = 1, N_SUPSYM
        IF(ML_FOR_SUPSYM(ISUPSYM).EQ.ML.AND.
     &     IPA_FOR_SUPSYM(ISUPSYM).EQ.IPA) THEN
           ISYM = ISUPSYM
        END IF
       END DO
      ELSE
        WRITE(6,*) ' Unknown type of supersymmetry ', CSUPSYM
        STOP       ' Unknown type of supersymmetry '
      END IF
*
      IF(ISYM.EQ.0) THEN
        WRITE(6,*) 
     &  ' Supersymmetry was not found for L, ML, IPA = ', 
     &   L, ML, IPA
        STOP ' Supersymmetry was not found for L, ML, IPA '
      END IF
*. And irrep
      IRREP = IRREP_FOR_SUPSYM(ISYM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,5I3)') ' L, ML, IPA (in) ISYM, IRREP(out) =  ',
     &               L, ML, IPA, ISYM, IRREP
      END IF
*
      RETURN
      END
C     GET_SUPSYM_FOR_BASIS(WORK(KISUPSYM_FOR_BAS),
C    &      WORK(KNBAS_FOR_SUP_STA_SYM),WORK(KIBBAS_FOR_SUP_STA_SYM),
C    &      WORK(KIBAS_FOR_SUP_STA_SYM),WORK(KIRREP_FOR_BAS), 
C    &      WORK(KNBAS_FOR_IRREP),WORK(KIBBAS_FOR_IRREP),
C    &      WORK(KIBAS_FOR_IRREP))
      SUBROUTINE GET_SUPSYM_FOR_BASIS(
     &           ISUPSYM_FOR_BAS, NBAS_FOR_SUP_STA_SYM, 
     &           IBBAS_FOR_SUP_STA_SYM, IBAS_FOR_SUP_STA_SYM,
     &           ISHELL_FOR_BAS,NBAS_FOR_SHELL,IBBAS_FOR_SHELL,
     &           IBAS_FOR_SHELL, NSUPSYM_FOR_IRREP,IBSUPSYM_FOR_IRREP,
     &           ISUPSYM_FOR_IRREP)
*
* Obtain supersymmetry info for basis set from labels 
* of basis functions
*
*. Jeppe Olsen, May 23, 2012
* Last modification; Jeppe Olsen; March 6, 2013; Info on irrep <=> basis added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      
*. General input
      CHARACTER*4 AO_CENT, AO_TYPE, CHAR4
      COMMON/AOLABELS/AO_CENT(MXPORB),AO_TYPE(MXPORB)
      INTEGER NSUPSYM_FOR_IRREP(*), IBSUPSYM_FOR_IRREP(*)
      INTEGER ISUPSYM_FOR_IRREP(*)
*. Output
      INTEGER ISUPSYM_FOR_BAS(*)
      INTEGER NBAS_FOR_SUP_STA_SYM(N_SUPSYM,NSMOB), 
     &        IBBAS_FOR_SUP_STA_SYM(N_SUPSYM,NSMOB),
     &        IBAS_FOR_SUP_STA_SYM(*)
      INTEGER ISHELL_FOR_BAS(NTOOB),NBAS_FOR_SHELL(NTOOB)
      INTEGER IBBAS_FOR_SHELL(NTOOB),IBAS_FOR_SHELL(NTOOB)
*. Output is also the arrays NBAS_SUPSYM, IBBAS_SUPSYM, IBAS_SUPSYM in /ORBINP/
*
      NTEST = 100
*
      CALL GET_SUPSYM_FOR_BASIS_S(ISUPSYM_FOR_BAS,
     &     int_mb(KLVAL_FOR_ORB),int_mb(KMLVAL_FOR_ORB),
     &     int_mb(KPA_FOR_ORB),NTOOB)
*
*
      IZERO = 0
      CALL ISETVC(NBAS_FOR_SUP_STA_SYM,IZERO,NSMOB*N_SUPSYM)
      CALL ISETVC(NBAS_SUPSYM,IZERO,N_SUPSYM)
      IBAS= 0
      DO ISTASYM = 1, NSMOB
        NBASS = NTOOBS(ISTASYM) 
        DO IIBAS = 1, NBASS
          IBAS = IBAS + 1
          ISUPSYM = ISUPSYM_FOR_BAS(IBAS)
          NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM) =
     &    NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM) + 1
          NBAS_SUPSYM(ISUPSYM) = NBAS_SUPSYM(ISUPSYM) + 1
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
     &  ' Number of basis functions per super- and standard-sym'
        WRITE(6,*) 
     &  ' ====================================================='
        CALL IWRTMA3(NBAS_FOR_SUP_STA_SYM,N_SUPSYM,NSMOB,N_SUPSYM,NSMOB)
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Number of basis functions per super-sym'
        WRITE(6,*) 
     &  ' ======================================='
        CALL IWRTMA3(NBAS_SUPSYM,1,N_SUPSYM,1,N_SUPSYM)
      END IF
*
*. And then the basisfunctions of a given super and standard symmetry
*
      IB = 1
      DO ISTASYM = 1, NSMOB
       DO ISUPSYM = 1, N_SUPSYM
        IBBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM) = IB
        NBAS = NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM)
        IB = IB + NBAS
       END DO
      END DO
*
      CALL ISETVC(NBAS_FOR_SUP_STA_SYM,IZERO,NSMOB*N_SUPSYM)
      IBAS = 0
      DO ISTASYM = 1, NSMOB
        NBAS = NTOOBS(ISTASYM)
        DO IIBAS = 1, NBAS
          IBAS = IBAS + 1
          ISUPSYM = ISUPSYM_FOR_BAS(IBAS)
C?        WRITE(6,*) ' TEST: ISTASYM, ISUPSYM = ', 
C?   &                       ISTASYM, ISUPSYM 
          NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM) = 
     &    NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM) + 1
          N = NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM)
          IB = IBBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM)
          IBAS_FOR_SUP_STA_SYM(IB-1+N) = IBAS
        END DO
      END DO
*
* and the orbitals of a given supersymmetry
*
      IB = 1
      DO ISUPSYM = 1, N_SUPSYM
        IBBAS_SUPSYM(ISUPSYM) = IB
        NBAS = NBAS_SUPSYM(ISUPSYM)
        IB = IB + NBAS
      END DO
*
      CALL ISETVC(NBAS_SUPSYM,IZERO,N_SUPSYM)
      DO IBAS = 1, NTOOB
        ISUPSYM = ISUPSYM_FOR_BAS(IBAS)
        NBAS_SUPSYM(ISUPSYM) = 
     &  NBAS_SUPSYM(ISUPSYM) + 1
        N = NBAS_SUPSYM(ISUPSYM)
        IB = IBBAS_SUPSYM(ISUPSYM)
        IBAS_SUPSYM(IB-1+N) = IBAS
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Basis functions of given super and standard sym '
        WRITE(6,*) ' ================================================'
        DO ISTASYM = 1, NSMOB
          DO ISUPSYM = 1, N_SUPSYM
            N = NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM)
            IB = IBBAS_FOR_SUP_STA_SYM(ISUPSYM,ISTASYM)
            IF(N.NE.0) THEN
             WRITE(6,'(A,2(1X,I3))') 
     &       ' Standard and supersymmetry ', ISTASYM, ISUPSYM
             CALL IWRTMA3(IBAS_FOR_SUP_STA_SYM(IB),1,N,1,N)
            END IF
          END DO
        END DO
*
        WRITE(6,*) ' Basis functions of given supersym '
        WRITE(6,*) ' =================================='
        DO ISUPSYM = 1, N_SUPSYM
          N = NBAS_SUPSYM(ISUPSYM)
          IB = IBBAS_SUPSYM(ISUPSYM)
          IF(N.NE.0) THEN
           WRITE(6,'(A,1(1X,I3))') 
     &     ' Supersymmetry ',  ISUPSYM
           CALL IWRTMA3(IBAS_SUPSYM(IB),1,N,1,N)
          END IF
        END DO
      END IF
*
*. Relation between basis functions /MO's in symmetry order and
*. Shells in symmetry order
*
*
*. a) The number of the shell corresponding to each
*     basis functions
*. A bit of inefficient coding (but I do not want to create a scratch array today...)
*
      IB_SHELL = 1
COLD  DO IRREP = 1, NIRREP
      DO IRREP = 1, N_SUPSYM_IRREP
*. Supersymmmetries of this irrep
        IB = IBSUPSYM_FOR_IRREP(IRREP)
        N  = NSUPSYM_FOR_IRREP(IRREP)
        IF(NTEST.GE.10) WRITE(6,*) ' IRREP, IB, N = ', IRREP, IB, N
        DO IISUPSYM = IB, IB-1+N
         ISUPSYM = ISUPSYM_FOR_IRREP(IISUPSYM)
         IBO = IBBAS_SUPSYM(ISUPSYM)
         NO = NBAS_SUPSYM(ISUPSYM)
         IBS = IB_SHELL
         IF(NTEST.GE.10) THEN
           WRITE(6,*) ' ISUPSYM, IBO, NO, IBS = ', 
     &                  ISUPSYM, IBO, NO, IBS
         END IF
         DO IIORB = IBO, IBO-1+NO
           IORB = IBAS_SUPSYM(IIORB)
           ISHELL_FOR_BAS(IORB) = IBS-1+IIORB-IBO+1
           IF(NTEST.GE.10) THEN
             WRITE(6,*) ' IIORB, IORB, RHS = ', 
     &                    IIORB, IORB, IBS-1+IIORB-IBO+1
           END IF
         END DO
        END DO
        IB_SHELL = IB_SHELL + NO
      END DO ! loop over IRREPS
      NTOSH = IB_SHELL - 1
      WRITE(6,*) ' Total number of shells = ', NTOSH
*
      IF(NTEST.GE.0) THEN
        WRITE(6,*) ' Shell number for orbitals in sym.order'
        WRITE(6,*) ' ======================================'
        CALL IWRTMA3(ISHELL_FOR_BAS,1,NTOOB,1,NTOOB)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_SUPSYM_FOR_BASIS_S(ISUPSYM_FOR_BAS,
     &           LVAL_FOR_ORB,MLVAL_FOR_ORB,IPA_FOR_ORB,
     &           NTOOB)
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LVAL_FOR_ORB(NTOOB),MLVAL_FOR_ORB(NTOOB)
      INTEGER IPA_FOR_ORB(NTOOB)
*. Output
      INTEGER ISUPSYM_FOR_BAS(NTOOB)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info for GET_SUPSYM_FOR_BASIS_S'
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
      END IF
     
*
      DO IORB = 1, NTOOB
        L  =  LVAL_FOR_ORB(IORB)
        ML = MLVAL_FOR_ORB(IORB)
        IPA = IPA_FOR_ORB(IORB)
        CALL SYM_AND_IRREP_FOR_LMLPA(L,ML,IPA,ISYM,IRREP)
C            SYM_AND_IRREP_FOR_LMLPA(L,ML,IPA,ISYM,IRREP)
        ISUPSYM_FOR_BAS(IORB) = ISYM
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Basis functions and their supersymmetry '
       WRITE(6,*) ' ========================================'
       WRITE(6,*)
       WRITE(6,*) ' Orbital   Supersymmetry '
       WRITE(6,*) ' ======================= '
       DO IORB = 1, NTOOB
         WRITE(6,'(2X,I3,5X,I2)') IORB, ISUPSYM_FOR_BAS(IORB)
       END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE GEN_GENSMOB(NGENSMOB,NBAS_GENSMOB,IBBAS_GENSMOB,
     &          ISTA_TO_GENSM_REO,
     &          NBAS_FOR_SUP_STA_SYM,IBBAS_FOR_SUP_STA_SYM,
     &          IBAS_FOR_SUP_STA_SYM,N_SUPSYM,NSMOB,ICHECK)
*
* Set the GENSMOB arrays for supersymmetry case
*
*. Jeppe Olsen, May 23, 2012
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NBAS_FOR_SUP_STA_SYM(N_SUPSYM,NSMOB),
     &        IBBAS_FOR_SUP_STA_SYM(N_SUPSYM,NSMOB),
     &        IBAS_FOR_SUP_STA_SYM(*)
*. Output
      INTEGER NBAS_GENSMOB(N_SUPSYM*NSMOB),
     &        IBBAS_GENSMOB(N_SUPSYM*NSMOB),
     &        ISTA_TO_GENSM_REO(*)
*
      NTEST = 100
      IF(NTEST.GE.100)  THEN
        WRITE(6,*) ' GEN_GENSMOB reporting'
        WRITE(6,*) ' ====================='
        WRITE(6,*) 
        WRITE(6,*) ' ICHECK = ', ICHECK
        WRITE(6,*) ' N_SUPSYM, NSMOB = ', N_SUPSYM, NSMOB
      END IF
*
      NGENSMOB = N_SUPSYM*NSMOB
      IGENSMOB = 0
      IB_OUT = 1
*
      DO ISMOB = 1, NSMOB
       DO ISUPSYM = 1, N_SUPSYM
        IF(NTEST.GE.1000) WRITE(6,'(A,2(2X,I2))')
     &  ' ISMOB, ISUPSUM = ', ISMOB, ISUPSUM
        IGENSMOB = IGENSMOB + 1
        N = NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISMOB)
        IB_IN = IBBAS_FOR_SUP_STA_SYM(ISUPSYM,ISMOB)
        NBAS_GENSMOB(IGENSMOB) = N
        IBBAS_GENSMOB(IGENSMOB) = IB_OUT
        WRITE(6,*) ' TEST: IB, IB_OUT, N = ',
     &                     IB, IB_OUT, N
        CALL ICOPVE(
     &  IBAS_FOR_SUP_STA_SYM(IB),ISTA_TO_GENSM_REO(IB_OUT),N)
        IB_OUT = IB_OUT + N
       END DO
      END DO
      NORB = IB_OUT - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Information on GENSMOB arrays '
        WRITE(6,*) ' =============================='
        WRITE(6,*)
        WRITE(6,*) ' Number of general orbital symmetries = ', NGENSMOB
        WRITE(6,*) ' Number of orbitals per general orbital symmetries'
        CALL IWRTMA3(NBAS_GENSMOB,1,NGENSMOB,1,NGENSMOB)
        WRITE(6,*) ' Offsets of orbitals for general orbital symmetry'
        CALL IWRTMA3(IBBAS_GENSMOB,1,NGENSMOB,1,NGENSMOB)
        WRITE(6,*) ' Reorder array, standard to general symmetry '
        CALL IWRTMA3(ISTA_TO_GENSM_REO,1,NORB,1,NSMOB)
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_CMO_STA_GEN(CMO_STA,CMO_GEN,
     &           IDO_REORDER,IREO,IWAY)
*
* Reform CMO  matrix between standard and general symmetry forms
*
* IWAY = 1: Standard => general symmetry
* IWAY = 2: general symmetry => standard
*
* IF IDO_REORDER = 0, then the orbitals are ordered as
*                         the basis-functions
*                    = 1, the reorder array IREO from the above order 
*                          is used
*. Jeppe Olsen, May 24, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
*. Input and output
      DIMENSION CMO_STA(*), CMO_GEN(*)
      INTEGER IREO(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from REFORM_CMO_STA_GEN'
        WRITE(6,*) ' ============================'
      END IF
*
      IOFF_C_GEN = 1
      DO IGENSM = 1, NGENSMOB
        NGEN = NBAS_GENSMOB(IGENSM)
        IBGEN = IBBAS_GENSMOB(IGENSM)
        ISYM = ISTASM_FOR_GENSM(IGENSM)
        NSTA = NTOOBS(ISYM)
        IF(NTEST.GE.1000) WRITE(6,'(A,4I4)')
     &  'IGENSM, ISYM, NGEN, NSTA = ', IGENSM, ISYM, NGEN, NSTA
*. Start of symmetry block in C of standard expansion
        IBSTA_C = 1
        DO JSYM = 1, ISYM-1
          IBSTA_C = IBSTA_C + NTOOBS(JSYM)**2
        END DO
*. Start of orbitals of symmetry ISYM
        IBSTA = ITOOBS(ISYM)
*
        DO IGEN = 1, NGEN
          ISTA = ISTA_TO_GENSM_REO(IBGEN-1+IGEN)
          IF(IDO_REORDER.EQ.1) THEN
            ISTA = IREO(ISTA)
          END IF
          IOFF_C_STA = IBSTA_C-1 + (ISTA-IBSTA)*NSTA + 1
C         REFORM_SINGLE_CMO_STA_GEN( CMO_STA,CMO_GEN,ISYM,IGENSM,IWAY)
          IF(NTEST.GE.1000) WRITE(6,'(A,2I6)')
     &    ' IOFF_C_STA, IOFF_C_GEN = ', IOFF_C_STA, IOFF_C_GEN
          CALL REFORM_SINGLE_CMO_STA_GEN(
     &    CMO_STA(IOFF_C_STA),CMO_GEN(IOFF_C_GEN),ISYM,IGENSM,IWAY)
          IOFF_C_GEN = IOFF_C_GEN + NGEN
        END DO !loop over orbitals of given general symmetry
      END DO! Loop over general symmetries
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CMO coefficients in general form '
        CALL APRBLM2(CMO_GEN,NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,0)
        WRITE(6,*) ' CMO coefficients in standard form '
        CALL APRBLM2(CMO_STA,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_SINGLE_CMO_STA_GEN
     &           (CMO_STA,CMO_GEN,ISYM,IGENSM,IWAY)
*
* Reform between standard and general symmetry form of expansions of
* a single orbital of symmetry ISYM and general symmetry IGENSM
*
* IWAY = 1: Standard => general symmetry
* IWAY = 2: general symmetry => standard
*
*. Jeppe Olsen, May 23 - FCN has just become Danish Champions in soccer-
*               FCK is down to silver...
*
      INCLUDE 'implicit.inc'
*. Input and output: Coefs of a single MO
      DIMENSION CMO_STA(*), CMO_GEN(*)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*
      NTEST = 000
      IF(NTEST.GE.100) WRITE(6,*) ' Info from REFORM_CMO... '
*
      NGEN = NBAS_GENSMOB(IGENSM)
      IBGEN = IBBAS_GENSMOB(IGENSM)
      NSTA = NTOOBS(ISYM)
      IF(NTEST.GE.1000) WRITE(6,'(A,4I4)')
     &' IGENSM, ISYM, NGEN, NSTA = ', IGENSM, ISYM, NGEN, NSTA
*. Start for given symmetry
      IB_SYM = ITOOBS(ISYM)
*
      IF(IWAY.EQ.1) THEN
*. Standard => supersymmetry
        DO IGEN = 1, NGEN
         CMO_GEN(IGEN) = 
     &   CMO_STA(ISTA_TO_GENSM_REO(IBGEN-1+IGEN)-IB_SYM+1)
        END DO
      ELSE
*. Supersymmetry to standard
        ZERO = 0.0D0
        CALL SETVEC(CMO_STA,ZERO,NSTA)
        DO IGEN = 1, NGEN
         CMO_STA(ISTA_TO_GENSM_REO(IBGEN-1+IGEN)-IB_SYM+1) =
     &   CMO_GEN(IGEN) 
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Standard => supersymmetry '
        ELSE 
          WRITE(6,*) ' Supersymmetry => standard '
        END IF
        WRITE(6,*) ' Single MO in standard form '
        CALL WRTMAT(CMO_STA,1,NSTA,1,NSTA)
        WRITE(6,*) ' Single MO in general symmetry form '
        CALL WRTMAT(CMO_GEN,1,NGEN,1,NGEN)
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_MAT_STA_SUP(ASTA,ASUP,IPACK,IWAY)
*
* Reform a matrix between standard and supersymmetry order
*
* IWAY = 1: Standard => supersymmetry
* IWAY = 2: Supersymmetry => standard
*
* IPACK = 0: Full blocked matrix
*       = 1: Lower half blocked matrix
*
*. Jeppe Olsen, May 23
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input and output
      DIMENSION ASTA(*),ASUP(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from REFORM_MAT_STA_SUP '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        WRITE(6,*) ' IWAY, IPACK = ', IWAY, IPACK 
      END IF
*
      IB_STA = 1
      IB_SUP = 1
      IJ_SUP = 0
      IJB_SYM = 1
      DO ISYM = 1, NSMOB
       IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' Info for ISYM = ', ISYM
         WRITE(6,*) ' ======================= '
         WRITE(6,*)
       END IF
       NSTA = NTOOBS(ISYM)
       IBSTA = ITOOBS(ISYM)
       IF(IWAY.EQ.2) THEN
*. Zero symmetry block
        IF(IPACK.EQ.0) THEN
          NSMBLK = NSTA*NSTA
        ELSE
          NSMBLK = NSTA*(NSTA+1)/2
        END IF
        ZERO = 0.0D0
C?      WRITE(6,*) ' IJB_SYM, NSMBLK = ', IJB_SYM, NSMBLK
        CALL SETVEC(ASTA(IJB_SYM),ZERO,NSMBLK)
       END IF
       NSUPSYM = NSUP_FOR_STA_SYM(ISYM) 
       IBSUPSYM = IBSUP_FOR_STA_SYM(ISYM)
       DO IISUPSYM = IBSUPSYM, IBSUPSYM - 1 + NSUPSYM
        ISUPSYM = ISUP_FOR_STA_SYM(IISUPSYM)
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' IISUPSYM, ISUPSYM= ', IISUPSYM, ISUPSYM
        NSUP = NBAS_GENSMOB(ISUPSYM)
        IBSUP = IBBAS_GENSMOB(ISUPSYM)
*. Start of general symmetry block ISUPSYM
        IB_SUP_MAT = 1
        DO JSUPSYM = 1, ISUPSYM-1
          NS = NBAS_GENSMOB(JSUPSYM)
          IF(IPACK.EQ.0) THEN
            IB_SUP_MAT = IB_SUP_MAT + NS**2
          ELSE
            IB_SUP_MAT = IB_SUP_MAT + NS*(NS+1)/2
          END IF
        END DO
        IJ_SUP = IB_SUP_MAT - 1
*
        DO ISUP = 1, NSUP
          IF(IPACK.EQ.0) THEN
           JSUP_MAX = NSUP
          ELSE 
           JSUP_MAX = ISUP
          END IF
          DO JSUP = 1, JSUP_MAX
            IJ_SUP = IJ_SUP + 1
*. Indeces of orbitals and matrix  in standard form
            I_STA = ISTA_TO_GENSM_REO(IBSUP-1+ISUP)-IBSTA+1
            J_STA = ISTA_TO_GENSM_REO(IBSUP-1+JSUP)-IBSTA+1
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' I_STA, J_STA =', I_STA, J_STA
            IF(IPACK.EQ.0) THEN
              IJ_STA = (J_STA-1)*NSTA+I_STA
            ELSE
              IJ_STA = I_STA*(I_STA-1)/2 + J_STA
            END IF
            IF(NTEST.GE.1000) WRITE(6,*)
     &      ' ISUP, JSUP, IJ_STA = ', ISUP, JSUP, IJ_STA
            IF(NTEST.GE.1000) WRITE(6,*)
     &      'IJ_SUP, IJB_SYM-1+IJ_STA ', IJ_SUP, IJB_SYM-1+IJ_STA
            IF(IWAY.EQ.1) THEN
              ASUP(IJ_SUP) = ASTA(IJB_SYM-1+IJ_STA)
            ELSE
              ASTA(IJB_SYM-1+IJ_STA) = ASUP(IJ_SUP)
            END IF
          END DO! Loop over JSUP
        END DO ! Loop over ISUP
       END DO! Loop over super symmetry
*. Update pointer to start of standard 
       IF(IPACK.EQ.0) THEN
         IJB_SYM = IJB_SYM + NSTA**2
       ELSE
         IJB_SYM = IJB_SYM + NSTA*(NSTA+1)/2
       END IF
C?     WRITE(6,*) ' NSTA and updated IJB_SYM = ', NSTA, IJB_SYM
      END DO ! End of loop over Symmetries
*
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*)  ' Standard => supersymmetry '
        ELSE
          WRITE(6,*)  ' Supersymmetry => Standard'
        END IF
        WRITE(6,*) ' Matrix in standard symmetry form '
        CALL APRBLM2(ASTA,NTOOBS,NTOOBS,NSMOB,IPACK)
        WRITE(6,*) ' Matrix in Super symmetry form '
        CALL APRBLM2(ASUP,NBAS_GENSMOB, NBAS_GENSMOB,NGENSMOB,IPACK)
      END IF
*
      RETURN
      END
      SUBROUTINE N_SUPSYM_IRREP_TO_SUPSYM(N_PER_IRREP,N_PER_SUPSYM)
*
* Obtain number of orbitals per super symmetry from 
* number of shells per supersymmetry irrep
*
*. Jeppe Olsen, May 23, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER N_PER_IRREP(*)
* Output
      INTEGER N_PER_SUPSYM(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from N_SUPSYM_IRREP_TO_SUPSYM '
        WRITE(6,*) ' ================================== '
        WRITE(6,*)
      END IF
*
      CALL N_SUPSYM_IRREP_TO_SUPSYM_S(
     &     N_PER_IRREP,N_PER_SUPSYM,
     &     N_SUPSYM_IRREP, N_SUPSYM,
     &     int_mb(KNSUPSYM_FOR_IRREP),
     &     int_mb(KIBSUPSYM_FOR_IRREP),
     &     int_mb(KISUPSYM_FOR_IRREP)  )
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of shells per supsym irrep '
        CALL IWRTMA3(N_PER_IRREP,1,N_SUPSYM_IRREP,1,N_SUPSYM_IRREP)
        WRITE(6,*) ' Number of orbitals per supersymmetry '
        CALL IWRTMA3(N_PER_SUPSYM,1,N_SUPSYM,1,N_SUPSYM)
      END IF
*
      RETURN
      END
      SUBROUTINE N_SUPSYM_IRREP_TO_SUPSYM_S(
     &     N_PER_IRREP,N_PER_SUPSYM,
     &     N_SUPSYM_IRREP, N_SUPSYM,
     &     NSUPSYM_FOR_IRREP,
     &     IBSUPSYM_FOR_IRREP,
     &     ISUPSYM_FOR_IRREP)
*
* Obtain number of orbitals per super symmetry from 
* number of shells per supersymmetry irrep - slave routine
*
*. Jeppe Olsen, May 23, 2012
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER N_PER_IRREP(N_SUPSYM_IRREP)
*. Output
      INTEGER N_PER_SUPSYM(N_SUPSYM)
*. General info
      INTEGER  NSUPSYM_FOR_IRREP(*),
     &         IBSUPSYM_FOR_IRREP(*),
     &         ISUPSYM_FOR_IRREP(*)
*
      NTEST = 00
      IF(NTEST.GE.100) WRITE(6,*) 'From N_SUPSYM_IRREP_TO_SUPSYM_S'
      IZERO = 0
      CALL ISETVC(N_PER_SUPSYM,IZERO,N_SUPSYM)
*
      DO IRREP = 1, N_SUPSYM_IRREP
        NSHL = N_PER_IRREP(IRREP)
        NCOMP = NSUPSYM_FOR_IRREP(IRREP)
        IBCOMP = IBSUPSYM_FOR_IRREP(IRREP)
        IF(NTEST.GE.1000)
     &  WRITE(6,*) ' IRREP, NCOMP, IBCOMP = ', IRREP, NCOMP, IBCOMP
        DO IICOMP = IBCOMP, IBCOMP-1+NCOMP
          ICOMP =  ISUPSYM_FOR_IRREP(IICOMP)
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IICOMP, ICOMP = ', IICOMP, ICOMP
          N_PER_SUPSYM(ICOMP) = NSHL
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE SUP_TO_STASYM(NBAS_FOR_SUP_STA_SYM)
*
* Super-symmetry components for each standard symmetry
*
* Output is in /ORBIBP/: NSUP_FOR_STA_SYM,IBSUP_FOR_STA_SYM,ISUP_FOR_STA_SYM,
*                        ISTASM_FOR_SUPSYM
*. Jeppe Olsen, May 24, 2012
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Specific input
      INTEGER NBAS_FOR_SUP_STA_SYM(N_SUPSYM,NSMOB)
*
      DO ISYM = 1, NSMOB
       NCOMP = 0
       DO ISUPSYM = 1, N_SUPSYM
        IF(NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISYM).NE.0) NCOMP = NCOMP + 1
       END DO
       NSUP_FOR_STA_SYM(ISYM) = NCOMP
      END DO
*
      IB = 1
      DO ISYM = 1, NSMOB
       IBSUP_FOR_STA_SYM(ISYM) = IB
       IB = IB + NSUP_FOR_STA_SYM(ISYM)
      END DO
*
      DO ISYM = 1, NSMOB
       IB = IBSUP_FOR_STA_SYM(ISYM)
       NCOMP = 0
       DO ISUPSYM = 1, N_SUPSYM
        IF(NBAS_FOR_SUP_STA_SYM(ISUPSYM,ISYM).NE.0) THEN
          NCOMP = NCOMP + 1
          ISUP_FOR_STA_SYM(IB-1+NCOMP) = ISUPSYM
          ISTASM_FOR_SUPSYM(ISUPSYM) = ISYM
        END IF
       END DO
      END DO
*
      NTEST = 100
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Supersymmetries for standard symmetry '
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        DO ISYM = 1, NSMOB
          WRITE(6,*)
          WRITE(6,*) ' Supersymmetries for standard symmetry', ISYM
          IB = IBSUP_FOR_STA_SYM(ISYM)
          N  = NSUP_FOR_STA_SYM(ISYM)
          CALL IWRTMA3(ISUP_FOR_STA_SYM(IB),1,N,1,N)
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Standard symmetry of supersymmetries '
        WRITE(6,*) ' ==================================== '
        WRITE(6,*)
        WRITE(6,*) ' Super    Standard '
        WRITE(6,*) ' =================='
        DO ISUP = 1, N_SUPSYM
          WRITE(6,'(2X,I2,8X,I2)') ISUP, ISTASM_FOR_SUPSYM(ISUP)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_MAT_STA_GEN(ASTA,AGEN,IPACK,IWAY)
*
* Reform a matrix between standard and general symmetry order
* as defined by the *GENSM* arrays - may be super or standard
* symmetry
*
* IWAY = 1: Standard => supersymmetry
* IWAY = 2: Supersymmetry => standard
*
* IPACK = 0: Full blocked matrix
*       = 1: Lower half blocked matrix
*
*. Jeppe Olsen, May 23, rewritten May 24
*. Last modification: Oct. 1, 2012: Jeppe Olsen: Removed bug for nonsymmetric matrices
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input and output
      DIMENSION ASTA(*),AGEN(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from REFORM_MAT_STA_GEN '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        WRITE(6,*) ' ISTA_TO_GENSM_REO '
        CALL IWRTMA3(ISTA_TO_GENSM_REO,1,NTOOB,1,NTOOB)
      END IF
*
      IF(IWAY.EQ.2) THEN
*. Zero matrix
       LEN_STA = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,IPACK)
       ZERO = 0.0D0
       CALL SETVEC(ASTA,ZERO,LEN_STA)
      END IF
*. Loop over blocks of general symmetry
      IB_GEN_MAT = 1
      IJ_GEN = 0
      DO IGENSM = 1, NGENSMOB
*. Standard symmetry
        ISTASM = ISTASM_FOR_GENSM(IGENSM)
        IB_STA_ORB = ITOOBS(ISTASM)
        NSTA = NTOOBS(ISTASM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IGENSM, ISTASM = ', IGENSM, ISTASM
        END IF
        IB_STA_MAT = 1
        DO JSTASM = 1, ISTASM-1
          N = NTOOBS(JSTASM)
          IF(IPACK.EQ.0) THEN
            IB_STA_MAT = IB_STA_MAT + N**2
          ELSE
            IB_STA_MAT = IB_STA_MAT + N*(N+1)/2
          END IF
        END DO
        IF(NTEST.GE.1000) WRITE(6,*) ' IB_STA_MAT = ', IB_STA_MAT
*. Loop over pairs of orbitals in general symmetry block
        NGEN = NBAS_GENSMOB(IGENSM)
        IBGEN = IBBAS_GENSMOB(IGENSM)
        IJ_GEN0 = IJ_GEN
        DO IBASGN = 1, NGEN
         IF(IPACK.EQ.0) THEN
           JBASGN_MAX = NGEN
         ELSE
           JBASGN_MAX = IBASGN
         END IF
         DO JBASGN = 1,  JBASGN_MAX 
*. The corresponding orbitals in standard order
           IBAS = ISTA_TO_GENSM_REO(IBGEN-1+IBASGN)-IB_STA_ORB + 1
           JBAS = ISTA_TO_GENSM_REO(IBGEN-1+JBASGN)-IB_STA_ORB + 1
           IF(IPACK.EQ.0) THEN
            IJ_STA = IB_STA_MAT - 1 + (JBAS-1)*NSTA + IBAS
           ELSE
            IJ_STA = IB_STA_MAT - 1 + IBAS*(IBAS-1)/2 + JBAS
           END IF
           IF(NTEST.GE.1000)
     &     WRITE(6,'(A,4I3)') ' IBASGN, IBAS, JBASGN, JBAS = ',
     &     IBASGN, IBAS, JBASGN, JBAS
*
           IF(IPACK.EQ.0) THEN
             IJ_GEN = IJ_GEN0 + (JBASGN-1)*NGEN + IBASGN
           ELSE
             IJ_GEN = IJ_GEN + 1
           END IF
           IF(NTEST.GE.1000) 
     &     WRITE(6,'(A,2I5)') ' IJ_GEN, IJ_STA = ', IJ_GEN, IJ_STA
           IF(IWAY.EQ.1) THEN
            AGEN(IJ_GEN) = ASTA(IJ_STA)
           ELSE
            ASTA(IJ_STA) = AGEN(IJ_GEN)
           END IF
         END DO! loop over IBASGN
        END DO! loop over JBASGN
      END DO! loop over IGENSM
*
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*)  ' Standard => supersymmetry '
        ELSE
          WRITE(6,*)  ' Supersymmetry => Standard'
        END IF
        WRITE(6,*) ' Matrix in standard symmetry form '
        CALL APRBLM2(ASTA,NTOOBS,NTOOBS,NSMOB,IPACK)
        WRITE(6,*) ' Matrix in general symmetry form '
        CALL APRBLM2(AGEN,NBAS_GENSMOB, NBAS_GENSMOB,NGENSMOB,IPACK)
      END IF
*
      RETURN
      END
      SUBROUTINE N_SUPSYM_TO_STASYM(N_PER_SUPSYM,N_PER_STASYM)
* A list of integers for supersymmetries are given
* Reform to list of integers over standard symmetries
*
*. Jeppe Olsen, May 24, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER N_PER_SUPSYM(*)
*. Output
      INTEGER N_PER_STASYM(*)
*
      NTEST = 00
*
      DO ISTASYM = 1, NSMOB
        NSUP = NSUP_FOR_STA_SYM(ISTASYM)
        IBSUP = IBSUP_FOR_STA_SYM(ISTASYM)
        N = 0
        DO IISUPSYM = IBSUP, IBSUP + NSUP - 1
         ISUPSYM = ISUP_FOR_STA_SYM(IISUPSYM)
         IF(NTEST.GE.1000)
     &   WRITE(6,'(A,3I4)') ' ISTASYM, IISUPSYM, ISUPSYM = ',
     &                        ISTASYM, IISUPSYM, ISUPSYM
         N = N + N_PER_SUPSYM(ISUPSYM)
        END DO
        N_PER_STASYM(ISTASYM) = N
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from N_SUPSYM_TO_STASYM:'
        WRITE(6,*) ' ==============================='
        WRITE(6,*) ' Integer list over supersymmetries (input) '
        CALL IWRTMA3(N_PER_SUPSYM,1,N_SUPSYM,1,N_SUPSYM)
        WRITE(6,*) ' Integer list over standard symmetries (output) '
        CALL IWRTMA3(N_PER_STASYM,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
      SUBROUTINE ORDER_SUPSYM_ORBITALS(NSPC,ISPC,MO_SUPSYM,IREO,
     &           ISUPSYM_FOR_BAS)
*
* A set of orbital spaces ISPC defined in terms of 
* supersymmetries are given. Order the orbitals accordingly.
*
*. Well, the deal is that the basis functions within a 
*. given standard symmetry is defined from input, whereas the 
*. ordering of the molecular orbitals in a given standard 
*. symmetry is not specified. LUCIA does as a standard 
*. use a canonical order where the MO's are ordered in the 
*. same way as the basis functions.
*. However, we typically define say the occupied orbitals 
*  or the CAS orbitals of a given symmetry as the lowest in 
*  a given stadard symmetry. To accomplish this, it is 
*  useful to introduce a reordering of the orbitals 
*  in a given symmetry. THis is what this routine is about...
*
* Jeppe Olsen, May 24, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER ISPC(MXP_NSUPSYM,NSPC), ISUPSYM_FOR_BAS(*)
*. Output: IREO: New order index from standard order index
      INTEGER MO_SUPSYM(NTOOB),IREO(*)
*. Local scratch: dim of number of irreps
      INTEGER ISCR1(1000),ISCR2(1000),ISCR3(1000)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from ORDER_SUPSYM_ORBITALS '
       WRITE(6,*) ' The requested division of orbitals'
       CALL IWRTMA3(ISPC,N_SUPSYM,NSPC,MXP_NSUPSYM,NSPC)
      END IF
*
* Obtain first MO_SUPSYM
*
      IORB = 0
      DO ISYM = 1, NSMOB
        NSUP = NSUP_FOR_STA_SYM(ISYM)
        IBSUP = IBSUP_FOR_STA_SYM(ISYM)
*. Loop over the spac division for this symmetry
        DO JSPC = 1, NSPC
          DO IISUPSYM = IBSUP, IBSUP + NSUP -1
            ISUPSYM = ISUP_FOR_STA_SYM(IISUPSYM)
            NSUPSPC = ISPC(ISUPSYM,JSPC)
            NORB_SUPSPC = ISPC(ISUPSYM,JSPC)
            DO IIORB = 1, NORB_SUPSPC
              IORB = IORB + 1
              MO_SUPSYM(IORB) = ISUPSYM
              IF(NTEST.GE.10000)
     &        WRITE(6,*) ' IORB, ISUPSYM ', IORB, ISUPSYM
            END DO! loop over orbital IORB
          END DO! Loop over supersymmetries 
        END DO ! Loop over orbital spaces
      END DO ! Loop over standard symmetries
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Supersymmetry of occ-ordered orbitals'
        CALL IWRTMA3(MO_SUPSYM,1,NTOOB,1,NTOOB)
      END IF
*
*. Counting index in ISCR2 for orbital with given supersymmetry
*
      IZERO = 0
      CALL ISETVC(ISCR1,IZERO,N_SUPSYM)
      DO IORB = 1, NTOOB
        ISUPSYM = MO_SUPSYM(IORB)
        ICOUNT = ISCR1(ISUPSYM)+1
        ISCR2(IORB) = ICOUNT
        ISCR1(ISUPSYM) = ISCR1(ISUPSYM) + 1
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Counting index of reordered orbitals'
        CALL IWRTMA3(ISCR2,1,NTOOB,1,NTOOB)
      END IF
*
*. Loop over orbitals in standard order and find the corresponding new 
*. order. Not elegant(quadratic scaling!!), but I'm to tired for elegance
*
      IZERO = 0
      CALL ISETVC(ISCR1,IZERO, N_SUPSYM)
      MONE = -1
      CALL ISETVC(IREO,MONE,NTOOB)
      DO ISYM = 1, NSMOB
       NORB = NTOOBS(ISYM)
       IBORB = ITOOBS(ISYM)
       DO IORB_STA = IBORB, IBORB + NORB -1
         ISUPSYM = ISUPSYM_FOR_BAS(IORB_STA)
C?       WRITE(6,*) ' ISUPSYM_FOR_BAS(1) = ', ISUPSYM_FOR_BAS(1)
         IF(NTEST.GE.1000)
     &   WRITE(6,*) ' IORB_STA, ISUPSYM = ', IORB_STA, ISUPSYM
         ISCR1(ISUPSYM) = ISCR1(ISUPSYM) + 1
         ICOUNT = ISCR1(ISUPSYM)
         IF(NTEST.GE.1000)
     &   WRITE(6,'(A,3(1X,I2))') ' IORB_STA, ISUPSUM, ICOUNT = ',
     &                             IORB_STA, ISUPSYM, ICOUNT
*. Find orbital with this supersymmetry and count number
         DO IORB = 1, NTOOB
           IF(ISUPSYM.EQ.MO_SUPSYM(IORB).AND.
     &        ICOUNT.EQ.ISCR2(IORB)) THEN
*. Match
              IREO(IORB_STA) = IORB
           END IF
         END DO! Loop over orbitals in new order
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reordering of orbitals (standard => occ/gas)'
        CALL IWRTMA3(IREO,1,NTOOB,1,NTOOB)
      END IF
*
      RETURN
      END
      SUBROUTINE REO_CMOAO(CIN,COUT,IREO,ICOPY,IWAY)
*
* Reorder coefficients of MOAO matrix
*
* IWAY = 1: COUT(IREO(I)) = CIN(I)
* IWAY = 2: COUT(I) = CIN(IREO(I)
*
*. If ICOPY = 1, COUT is copied over CIN
*
* Jeppe Olsen, May 25, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input (and perhaps output)
      DIMENSION CIN(*), IREO(NTOOB)
*. Output 
      DIMENSION COUT(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from REO_CMOAO'
        WRITE(6,*) ' ==================='
        WRITE(6,*) 
        WRITE(6,*) ' Input matrix '
C       CALL APRBLM2(CIN,NTOOBS,NTOOBS,NSMOB,0)
        CALL APRBLM_F7(CIN,NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*) ' Reorder array '
        CALL IWRTMA3(IREO,1,NTOOB,1,NTOOB)
      END IF
*
      IB_SYM_MAT = 1
      DO ISYM = 1, NSMOB
        NORB = NTOOBS(ISYM)
        IBORB = ITOOBS(ISYM)
        DO IORB_IN = IBORB, IBORB + NORB - 1
          IORB_REO = IREO(IORB_IN)
          IOFF_IN = IB_SYM_MAT + (IORB_IN-IBORB)*NORB
          IOFF_REO = IB_SYM_MAT + (IORB_REO-IBORB)*NORB
          IF(NTEST.GE.1000) WRITE(6,'(A,4(1X,I3))')
     &    'IORB_IN, IORB_REO, IOFF_IN,IOFF_REO = ',
     &     IORB_IN, IORB_REO, IOFF_IN,IOFF_REO
          IF(IWAY.EQ.1) THEN
            CALL COPVEC(CIN(IOFF_IN),COUT(IOFF_REO),NORB)
          ELSE
            CALL COPVEC(CIN(IOFF_REO),COUT(IOFF_IN),NORB)
          END IF
        END DO ! Loop over IORB_IN
        IB_SYM_MAT = IB_SYM_MAT + NORB*NORB
      END DO! Loop over ISYM
*
      LEN_C = IB_SYM_MAT
      IF(ICOPY.EQ.1) CALL COPVEC(COUT,CIN,LEN_C)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Reordered MO matrix '
        WRITE(6,*) ' ===================='
        WRITE(6,*)
        IF(ICOPY.EQ.1) THEN
C         CALL APRBLM2(CIN,NTOOBS,NTOOBS,NSMOB,0)
          CALL APRBLM_F7(CIN,NTOOBS,NTOOBS,NSMOB,0)
        ELSE
C         CALL APRBLM2(COUT,NTOOBS,NTOOBS,NSMOB,0)
          CALL APRBLM_F7(COUT,NTOOBS,NTOOBS,NSMOB,0)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE SET_HF_DIST_SUPSYM
*
* Define the matrix HF_DSV_SUPSYM,HF_DSV_STASYM,HF_DSV_GNSYM 
* giving number of orbitals per supsym, stasym, gnsym
* in Doubly, Singly and virtual spaces 
*
* Input and output are all in ORBINP
*
* Jeppe Olsen, May 26, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from SET_HF_DIST_SUPSYM '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' NHFD_SUPSYM, NHFS_SUPSYM, NBAS_SUPSYM( input ) '
        CALL IWRTMA3(NHFD_SUPSYM,1,N_SUPSYM,1,N_SUPSYM)
        CALL IWRTMA3(NHFS_SUPSYM,1,N_SUPSYM,1,N_SUPSYM)
        CALL IWRTMA3(NBAS_SUPSYM,1,N_SUPSYM,1,N_SUPSYM)
*
        WRITE(6,*) ' NHFD_STASYM, NHFS_STASYM, ( input ) '
        CALL IWRTMA3(NHFD_STASYM,1,NSMOB,1,MXPOBS)
        CALL IWRTMA3(NHFS_STASYM,1,NSMOB,1,MXPOBS)
      END IF
*
*. Just copy the info on doubly and singly occupied orbitals
*
      CALL ICOPVE(NHFD_SUPSYM,NHF_DSV_SUPSYM(1,1),N_SUPSYM)
      CALL ICOPVE(NHFS_SUPSYM,NHF_DSV_SUPSYM(1,2),N_SUPSYM)
*
      CALL ICOPVE(NHFD_STASYM,NHF_DSV_STASYM(1,1),N_SUPSYM)
      CALL ICOPVE(NHFS_STASYM,NHF_DSV_STASYM(1,2),N_SUPSYM)
*. And the virtual spaces as the NBAS - NDOUBLE - NSINGLE
      IONE = 1
      IMONE = -1
C          IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
      CALL IVCSUM(NHF_DSV_SUPSYM(1,3),NBAS_SUPSYM,NHF_DSV_SUPSYM(1,1),
     &     IONE,IMONE,N_SUPSYM)
      CALL IVCSUM(NHF_DSV_SUPSYM(1,3),NHF_DSV_SUPSYM(1,3),
     &     NHF_DSV_SUPSYM(1,2),IONE,IMONE,N_SUPSYM)
*
      CALL IVCSUM(NHF_DSV_STASYM(1,3),NTOOBS,NHF_DSV_STASYM(1,1),
     &     IONE,IMONE,NSMOB)
      CALL IVCSUM(NHF_DSV_STASYM(1,3),NHF_DSV_STASYM(1,3),
     &     NHF_DSV_STASYM(1,2),IONE,IMONE,NSMOB)
*
      CALL ICOPVE(NHF_DSV_SUPSYM(1,1),NHF_DSV_GNSYM(1,1),N_SUPSYM)
      CALL ICOPVE(NHF_DSV_SUPSYM(1,2),NHF_DSV_GNSYM(1,2),N_SUPSYM)
      CALL ICOPVE(NHF_DSV_SUPSYM(1,3),NHF_DSV_GNSYM(1,3),N_SUPSYM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Doubly, singly, and unoccupied orbitals per supersymmetry'
        WRITE(6,*) 
     &  ' ========================================================='
        CALL IWRTMA3(NHF_DSV_SUPSYM,N_SUPSYM,3,MXP_NSUPSYM,3)
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Doubly, singly, and unoccupied orbitals per standard symmetry'
        WRITE(6,*) 
     &  ' ============================================================='
        CALL IWRTMA3(NHF_DSV_STASYM,NSMOB,3,MXPOBS,3)
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Doubly, singly, and unoccupied orbitals per gensymmetry'
        WRITE(6,*) 
     &  ' ========================================================='
        CALL IWRTMA3(NHF_DSV_GNSYM,N_SUPSYM,3,MXP_NSUPSYM,3)
      END IF
*
      RETURN
      END
      SUBROUTINE SET_HF_DIST_STASYM
*
* Define the matrix HF_DSV_STASYM,HF_DSV_GNSYM 
* giving number of orbitals per supsym, stasym, gnsym
* in Doubly, Singly and virtual spaces 
*
* Input and output are all in ORBINP. 
* The number of doubly and singly occupied orbitals were given in 
* NHFD_IRREP_SUPSYM, NHFS_IRREP_SUPSYM - not logical, I admit
*
* Jeppe Olsen, May 26, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from SET_HF_DIST_STASYM '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
      END IF
*
*. Just copy the info on doubly and singly occupied orbitals
*
      CALL ICOPVE(NHFD_IRREP_SUPSYM, NHFD_SUPSYM,NSMOB)
      CALL ICOPVE(NHFS_IRREP_SUPSYM, NHFS_SUPSYM,NSMOB)
*
      CALL ICOPVE(NHFD_SUPSYM,NHF_DSV_STASYM(1,1),NSMOB)
      WRITE(6,*) ' NHF_DSV_STASYM(*,1): '
      CALL IWRTMA3(NHF_DSV_STASYM(1,1),1,NSMOB,1,NSMOB)
      CALL ICOPVE(NHFS_SUPSYM,NHF_DSV_STASYM(1,2),NSMOB)
*. And the virtual spaces as the NBAS - NDOUBLE - NSINGLE
      IONE = 1
      IMONE = -1
C          IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
*
      CALL IVCSUM(NHF_DSV_STASYM(1,3),NTOOBS,NHF_DSV_STASYM(1,1),
     &     IONE,IMONE,NSMOB)
      CALL IVCSUM(NHF_DSV_STASYM(1,3),NHF_DSV_STASYM(1,3),
     &     NHF_DSV_STASYM(1,2),IONE,IMONE,NSMOB)
*
      CALL ICOPVE(NHF_DSV_STASYM(1,1),NHF_DSV_GNSYM(1,1),NSMOB)
      CALL ICOPVE(NHF_DSV_STASYM(1,2),NHF_DSV_GNSYM(1,2),NSMOB)
      CALL ICOPVE(NHF_DSV_STASYM(1,3),NHF_DSV_GNSYM(1,3),NSMOB)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Doubly, singly, and unoccupied orbitals per standard symmetry'
        WRITE(6,*) 
     &  ' ============================================================='
        CALL IWRTMA3(NHF_DSV_STASYM,NSMOB,3,MXPOBS,3)
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Doubly, singly, and unoccupied orbitals per gensymmetry'
        WRITE(6,*) 
     &  ' ========================================================='
        CALL IWRTMA3(NHF_DSV_GNSYM,NSMOB,3,MXP_NSUPSYM,3)
      END IF
*
      RETURN
      END
      SUBROUTINE ORDER_GAS_SUPSYM_ORBITALS
*
* Obtain the order of the orbitals according to the
* symmetry specified by NGAS_IRREP_SUPSYM
*
* Jeppe Olsen, May 26
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from ORDER_GAS_SUPSYM_ORBITALS '
        WRITE(6,*) ' ===================================== '
        WRITE(6,*)
      END IF
*. Reform from irrep to super-symmetry and standard symmetry
      DO IGAS = 0, NGAS + 1
        CALL N_SUPSYM_IRREP_TO_SUPSYM(
     &       NGAS_IRREP_SUPSYM(1,IGAS),NGAS_SUPSYM(1,IGAS))
        CALL N_SUPSYM_TO_STASYM(
     &       NGAS_SUPSYM(1,IGAS),NGAS_STASYM(1,IGAS))
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GA-spaces over supersymmetries '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL IWRTMA3(NGAS_SUPSYM,N_SUPSYM,NGAS+2,MXP_NSUPSYM,NGAS+2)
        WRITE(6,*) ' GA-spaces over standard symmetries '
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
        CALL IWRTMA3(NGAS_STASYM,NSMOB,NGAS+2,MXPOBS,NGAS+2)
      END IF
*
* Obtain the reordering for this set of orbitals
*
      CALL ORDER_SUPSYM_ORBITALS(NGAS+2,NGAS_SUPSYM,
     &     int_mb(KMO_SUPSYM),int_mb(KMO_STA_TO_ACT_REO),
     &     int_mb(KISUPSYM_FOR_BAS)                )
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Super-symmetry reordering of orbitals, GAS '
        CALL IWRTMA3(int_mb(KMO_STA_TO_ACT_REO),1,NTOOB,1,NTOOB)
      END IF
*
      RETURN
      END
      SUBROUTINE SET_PARITY_FOR_STASYM(IPA_FOR_STASYM,NIRREP)
*
* Define parity (+1/-1) for the various standard symmetries
* Currently assuming D2H in the standard MOLECULE order
*
* Jeppe Olsen, June 2012
*
      INCLUDE 'implicit.inc'
*. Output
      INTEGER IPA_FOR_STASYM(*)
*
      NTEST = 100
*
      IF(NIRREP.EQ.8) THEN
*. the gerade symmetries
        IPA_FOR_STASYM(1) = 1
        IPA_FOR_STASYM(4) = 1
        IPA_FOR_STASYM(6) = 1
        IPA_FOR_STASYM(7) = 1
*. and the ungerade symmetries
        IPA_FOR_STASYM(2) =-1
        IPA_FOR_STASYM(3) =-1
        IPA_FOR_STASYM(5) =-1
        IPA_FOR_STASYM(8) =-1
      ELSE 
        WRITE(6,*) ' SET_PARITY_FOR_STASYM not set for NIRREP = ',NIRREP
        STOP       ' SET_PARITY_FOR_STASYM not set for NIRREP  '
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Parity of standard symmetries '
        CALL IWRTMA3(IPA_FOR_STASYM,1,NIRREP,1,NIRREP)
      END IF
*
      RETURN
      END
      SUBROUTINE SUPSYM_FROM_CMOAO(CMOAO,ISUPSYM_FOR_BAS,
     &           ISUPSYM_FOR_MOS)
*
* A CMOAO matrix is given in CMOAO and the supersymmetry of the basis functions is given in
* ISUPSYM_FOR_BAS. Obtain the supersymmetry of the MO's
*
*. Jeppe Olsen, July 3, 2012
*. Last modification; Oct. 3, 2012; Jeppe Olsen; Removed bug
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER ISUPSYM_FOR_BAS(*)
      DIMENSION CMOAO(*)
*. Output
      INTEGER ISUPSYM_FOR_MOS(NTOOB)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from SUPSYM_FROM_CMOAO '
        WRITE(6,*) ' ============================'
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input CMO '
        WRITE(6,*)
        CALL PRINT_CMOAO(CMOAO)
      END IF
*
      IJOFF = 1
      IORB  = 0
      DO ISM = 1, NSMOB
        NI = NTOOBS(ISM)
        DO I = 1, NI
          IORB = IORB + 1
          ICOFF = IJOFF-1+(I-1)*NI + 1
          ISUPSYM = ISUPSYM_FOR_MO(CMOAO(ICOFF),ISM,ISUPSYM_FOR_BAS)
          ISUPSYM_FOR_MOS(IORB) = ISUPSYM
        END DO
        IJOFF = IJOFF + NI**2
      END DO! loop over ISM
*. Check that all orbitals were assigned supesymmetry
      NZERO = 0
      DO IOB = 1, NTOOB
        IF(ISUPSYM_FOR_MOS(IOB).EQ.0) NZERO = NZERO + 1
      END DO
*
      IF(NTEST.GE.100.OR.NZERO.NE.0) THEN
        WRITE(6,*) ' Supersymmetries determined from CMOAO matrix '
        CALL IWRTMA3(ISUPSYM_FOR_MOS,1,NTOOB,1,NTOOB)
      END IF
*
      IF(NZERO.NE.0) THEN
        WRITE(6,*) ' Not all orbitals have well-defined super-symmetry'
        WRITE(6,*) ' Number of orbitals with unassigned super-sym ', 
     &             NZERO
        STOP       ' Not all orbitals have well-defined super-symmetry'
      END IF
*
      RETURN
      END 
      FUNCTION ISUPSYM_FOR_MO(CMO,ISM, ISUPSYM_FOR_BAS)
*
* The MO-coefficients of an MO of standard symmetry ISM is given in CMO.
*. Obtain the super symmetry of this orbital. If the orbital does not have 
* a well-defined super-symmetry, a zero is returned.
*
*. Jeppe Olsen, July 3, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Input
       INTEGER ISUPSYM_FOR_BAS(*)
       DIMENSION CMO(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ISUPSYM_FOR_MO '
        WRITE(6,*) ' =========================='
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Standard symmetry of CMO ', ISM
        WRITE(6,*) ' Expansion CMO: '
        NI = NTOOBS(ISM)
        CALL WRTMAT_F7(CMO,1,NI,1,NI)
      END IF
*
*. Offset to basis functions of standard symmetry ISM
      I_OFF = 1
      DO IISM = 1, ISM-1
       I_OFF = I_OFF + NTOOBS(IISM)
      END DO
      NI = NTOOBS(ISM)
*. Find supersymmetry of basis function with largest expansion coefficient
      ISUPSYM_MAX = 0
      IMAX_BAS = 0
      C_MAX = 0.0D0
      DO IBAS = 1, NI
        IF(ABS(CMO(IBAS)).GT.ABS(C_MAX))THEN
          C_MAX = CMO(IBAS)
          IMAX_BAS = IBAS
        END IF
      END DO
*
      IF(C_MAX.EQ.0.0D0) THEN
*. vanishing expansion, supersymmetry is set to zero
       WRITE(6,*) ' Warning: vanishing CMO in ISUPSYM_FOR_MO'
       ISUPSYM = 0
      ELSE
       ISUPSYM = ISUPSYM_FOR_BAS(I_OFF-1+IMAX_BAS)
       IF(NTEST.GE.1000) 
     & WRITE(6,*) ' Super-sym of bf with max coef ', ISUPSYM
*. Make sure that all non-vanishing coefficients have the same supersymmetry
*.(could be relaxed by using a non-vanishing threshold)
*. Note, I am pt using a rather high threshold, so some contamination
*. is allowed. Done in connection with some MCVB calculations
       THRES = 1.0D-6
       NMISS = 0
       CSUM = 0.0D0
       DO IBAS = 1, NI
         IF(ABS(CMO(IBAS)).GT.THRES) THEN
          JSUPSYM = ISUPSYM_FOR_BAS(I_OFF-1+IBAS)
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IBAS, JSUPSYM = ', IBAS,JSUPSYM
          IF(JSUPSYM.NE.ISUPSYM) THEN 
            NMISS = NMISS + 1
            CSUM = CSUM + ABS(CMO(IBAS))
            ILAST = IBAS
          END IF
         END IF
       END DO
       IF(NMISS.NE.0) THEN
*. Not well defined supersymmetry, set to zero
        WRITE(6,*) ' Problem, ISM, NI, ISUPSYM = ', ISM, NI, ISUPSYM
        WRITE(6,*) 
     &  ' Number of coefficients with deviating supersymmetry ', NMISS
        WRITE(6,*) 
     &  ' Sum of coefficients with deviating supersymmetry ', CSUM
        WRITE(6,*) 
     &  ' Last included basis function with wrong supersymmetry ', 
     &    ILAST
        WRITE(6,*) ' Expansion CMO: '
        CALL WRTMAT_F7(CMO,1,NI,1,NI)
        ISUPSYM = 0
       END IF
      END IF ! CMO was nonvanishing
*
      ISUPSYM_FOR_MO = ISUPSYM
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Supersymmetry of MO: ', ISUPSYM
      END IF
*
      RETURN
      END
      SUBROUTINE EXTR_CP_GASBLKS_FROM_GENSYM_MAT
     &           (AS,ASG,IEORC,IGAS_F,IGAS_L,IPAK)
*
* Reform between two diagonal blockings of orbital matrices
* S: Blocks over general symmetry, all blocks between IGAS_F and 
*       IGAS_L
* SG: Blocks over general symmetry and GAS
*
*. Jeppe Olsen, July 2012
*. Last modification, Sept. 24 2012, Jeppe Olsen, Debugged..
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Input and output
      DIMENSION AS(*),ASG(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from EXTR_CP_GASBLKS_FROM_GENSYM_MAT'
       WRITE(6,*) ' ========================================='
      END IF
*
      IJOFF_S = 1
      IJOFF_SG = 1
      IOFF_S = 1
      IOFF_SG = 1
*
      DO IGENSM = 1, NGENSMOB 
       IF(NTEST.GE.1000) WRITE(6,*) ' IGENSM = ', IGENSM
*. Number of orbitals in AFULL for this supersymmetry
       NORB_S = 0
       DO IGAS = IGAS_F, IGAS_L
        NORB_S = NORB_S + NGAS_GNSYM(IGENSM,IGAS)
       END DO
       IOFF_S = 1
       IOFF_SG = 1
       DO IGAS = IGAS_F, IGAS_L
         IF(NTEST.GE.1000) WRITE(6,*)  ' IGAS = ', IGAS
         NORB_SG = NGAS_GNSYM(IGENSM,IGAS)
*. Loop over pairs of orbitals in this symmetry-gas block
         DO IORB_SG = 1,NORB_SG
          IORB_S = IOFF_SG - 1 + IORB_SG  
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IORB_S, IORB_SG = ', IORB_S, IORB_SG
          IF(IPAK.EQ.1) THEN
           JORB_MX = IORB_SG
          ELSE
           JORB_MX = NORB_SG
          END IF
          DO JORB_SG = 1, JORB_MX
            JORB_S = IOFF_SG - 1 + JORB_SG  
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' JORB_S, JORB_SG = ', JORB_S, JORB_SG
            IF(IPAK.EQ.0) THEN
             IJ_SG = (JORB_SG-1)*NORB_SG + IORB_SG
             IJ_S =  (JORB_S-1)*NORB_S + IORB_S
            ELSE
             IJ_SG = IORB_SG*(IORB_SG-1)/2 + JORB_SG
             IJ_S = IORB_S*(IORB_S-1)/2 + JORB_S
            END IF
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' IJ_SG, IJ_S = ', IJ_SG, IJ_S
            IJ_SG_ABS = IJOFF_SG - 1 + IJ_SG
            IJ_S_ABS = IJOFF_S - 1 + IJ_S
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' IJ_SG_ABS, IJ_S_ABS = ', IJ_SG_ABS, IJ_S_ABS
*
            IF(IEORC.EQ.1) THEN
             AS(IJ_S_ABS) = ASG(IJ_SG_ABS)
            ELSE
             ASG(IJ_SG_ABS) = AS(IJ_S_ABS)
            END IF
          END DO ! Loop over JORB_SG
         END DO ! Loop over IORB_SG
*
         IOFF_SG = IOFF_SG + NORB_SG
         IF(IPAK.EQ.0) THEN
           IJOFF_SG = IJOFF_SG + NORB_SG**2
         ELSE
           IJOFF_SG = IJOFF_SG + NORB_SG*(NORB_SG+1)/2
         END IF
         IF(NTEST.GE.1000)
     &   WRITE(6,*) ' Updated IJOFF_SG = ', IJOFF_SG
       END DO ! loop over IGAS
*
       IF(IPAK.EQ.0) THEN
         IJOFF_S = IJOFF_S + NORB_S**2
       ELSE
         IJOFF_S = IJOFF_S + NORB_S*(NORB_S+1)/2
       END IF
       IF(NTEST.GE.1000)
     & WRITE(6,*) ' Updated IJOFF_S = ', IJOFF_S
      END DO !loop over IGENSM
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Matrix over general symmetries '
        CALL WRT_SG_MAT(AS,1,IGAS_F,IGAS_L,IPAK,1)
        WRITE(6,*) ' Matrix over general symmetries and Gasblocks'
        CALL WRT_SG_MAT(ASG,2,IGAS_F,IGAS_L,IPAK,1)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_SG_MAT(A,IS_OR_SG,IGAS_F,IGAS_L,IPAK,IEXT)
*
* Matrix A consists of diagonal general symmetry-blocks (IS_OR_SG = 1)
* or general symmetry-gas blocks(IS_OR_SG=2), both for gas-paces between
* IGAS_F and IGAS_L. Print this
* IEXT = 0 => Compact (F7)
* IEXT = 1 => Standard output
*
*. Jeppe Olsen, July 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION A(*)
*
      IOFF = 1
      DO IGENSM = 1, NGENSMOB
        NI_S = 0
        DO IGAS = IGAS_F, IGAS_L
          NI_S = NI_S +  NGAS_GNSYM(IGENSM, IGAS)
        END DO
        WRITE(6,'(A, I2)') ' Symmetry block ', IGENSM
        IF(IS_OR_SG.EQ.1) THEN
          IF(IPAK.EQ.0) THEN
            IF(IEXT.EQ.0) THEN
              CALL WRTMAT_F7(A(IOFF),NI_S,NI_S,NI_S,NI_S)
            ELSE
              CALL WRTMAT(A(IOFF),NI_S,NI_S,NI_S,NI_S)
            END IF
          ELSE
            IF(IEXT.EQ.0) THEN
              CALL PRSYM_F7(A(IOFF),NI_S)
            ELSE
              CALL PRSYM(A(IOFF),NI_S)
            END IF
          END IF
        ELSE IF(IS_OR_SG .EQ. 2) THEN
          DO IGAS = IGAS_F, IGAS_L
            WRITE(6,*) ' Diagonal block with GAS = ', IGAS
            NI_SG = NGAS_GNSYM(IGENSM,IGAS)
C?          WRITE(6,*) ' IGENSM, IGAS, NI_SG = ',
C?   &                   IGENSM, IGAS, NI_SG 
            IF(IPAK.EQ.0) THEN
              IF(IEXT.EQ.0) THEN
                CALL WRTMAT_F7(A(IOFF),NI_SG,NI_SG,NI_SG,NI_SG)
              ELSE
                CALL WRTMAT(A(IOFF),NI_SG,NI_SG,NI_SG,NI_SG)
              END IF
              IOFF = IOFF + NI_SG**2
            ELSE
              IF(IEXT.EQ.0) THEN
                CALL PRSYM_F7(A(IOFF),NI_SG)
              ELSE
                CALL PRSYM(A(IOFF),NI_SG)
              END IF
              IOFF = IOFF + NI_SG*(NI_SG+1)/2
            END IF
          END DO
        END IF !IS_OR_SG switch
*
        IF(IS_OR_SG.EQ.1) THEN
          IF(IPAK.EQ.0) THEN
            IOFF = IOFF + NI_S**2
          ELSE
            IOFF = IOFF + NI_S*(NI_S+1)/2
          END IF
        END IF
*
      END DO ! Loop over general symmetry blocks
*
      RETURN
      END
      SUBROUTINE REFORM_RHO1_TO_GNSM(
     &           RHO1_ST,RHO1_GNSM_ST,IWAY,IREO_GNSYM_TO_TS)
*
* Reform between standard and general symmetry-order
* of total symmetric density matrix
*
*. Jeppe Olsen, July 2012
*
* Last modified, July 8, 2012 (Jeppe)
*
* IWAY = 1 => Standard form to general symmetry blocked
* IWAY = 2 => General symmetry blockes to standard from 
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IREO_GNSYM_TO_TS(*)
*. Input and output
      DIMENSION RHO1_ST(*), RHO1_GNSM_ST(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from REFORM_RHO1_TO_GNSM '
        WRITE(6,*) ' IWAY = ', IWAY
        WRITE(6,*) ' IREO_GNSYM_TO_TS: '
        CALL IWRTMA(IREO_GNSYM_TO_TS,1,NACOB,1,NACOB)
      END IF
*
      IJOFF_S = 1
      IOFF_S = 1
      DO IGENSM = 1, NGENSMOB
        NACT_S = 0  
        DO IGAS = 1, NGAS
         NACT_S = NACT_S + NGAS_GNSYM(IGENSM,IGAS)
        END DO
        DO IORB_S = 1, NACT_S
        DO JORB_S = 1, NACT_S
*
         IORB = IREO_GNSYM_TO_TS(IOFF_S-1+IORB_S)
         JORB = IREO_GNSYM_TO_TS(IOFF_S-1+JORB_S)
C?       WRITE(6,*) ' IGENSM, IORB_S, JORB_S, IORB, JORB = ',
C?   &                IGENSM, IORB_S, JORB_S, IORB, JORB
*
         IF(IWAY.EQ.1) THEN
           RHO1_GNSM_ST(IJOFF_S-1+(JORB_S-1)*NACT_S+IORB_S) =
     &     RHO1_ST((JORB-1)*NACOB + IORB)
         ELSE
           RHO1_ST((JORB-1)*NACOB + IORB) = 
     &     RHO1_GNSM_ST(IJOFF_S-1+(JORB_S-1)*NACT_S+IORB_S) 
         END IF ! IWAY switch
        END DO
        END DO! End of loops over orbitals IORB_S, JORB_S
        IOFF_S = IOFF_S + NACT_S
        IJOFF_S = IJOFF_S + NACT_S**2
      END DO ! Loop over symmetries
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' 1-body density as NACOB x NACOB matrix'
        CALL WRTMAT(RHO1_ST,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*)
        WRITE(6,*) ' 1-body density as blocks over general symmetry'
C       WRT_SG_MAT(A,IS_OR_SG,IGAS_F,IGAS_L,IPAK,IEXT)
        CALL WRT_SG_MAT(RHO1_GNSM_ST,1,1,NGAS,0,1)
      END IF
*
      RETURN
      END
      SUBROUTINE REO_ACT_ORB_TO_GNSM(IMO_GNSYM, IREO_GNSYM_TO_TS)
*
* Reorder array for active orbitals from General symmetry to
* standard order for active orbitals.
*
*. Jeppe Olsen, July 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IMO_GNSYM(NTOOB)
*. Output
      INTEGER IREO_GNSYM_TO_TS(NACOB)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from REO_ACT_ORB_TO_GNSM'
        WRITE(6,*) ' Input: IMO_GNSYM '
        CALL IWRTMA3(IMO_GNSYM,1,NTOOB,1,NTOOB)
      END IF
*
      IOFF_GNSYM = 1
      IACOB = 0
      DO IGNSYM = 1, NGENSMOB
C?      WRITE(6,*) ' IGNSYM = ', IGNSYM
        NINOB_S = NGAS_GNSYM(IGNSYM,0)
        NACOB_S = 0
        DO IGAS = 1, NGAS
          NACOB_S = NACOB_S + NGAS_GNSYM(IGNSYM,IGAS)
        END DO
C?      WRITE(6,*) ' NINOB_S, NACOB_S = ', NINOB_S, NACOB_S
        IOB_S = 0
*. Find orbitals NINOB_S+1, NINOB_S + NACOB_S of sym IGNSYM
        DO IORB = 1, NTOOB
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' TESTJ: IORB, GNSYM = ', IORB, IMO_GNSYM(IORB)
          IF(IMO_GNSYM(IORB).EQ.IGNSYM) THEN
            IOB_S = IOB_S + 1
            IF(NINOB_S.LT.IOB_S.AND.IOB_S.LE.NINOB_S+NACOB_S) THEN
* Orbital IORB is active orbital IACOB in general sym order
              IACOB = IACOB + 1
*. Orbital IORB address in standard type order 
              IACOB_STA = IREOST(IORB)-NINOB
C?            WRITE(6,*) ' IACOB, IACOB_STA '
              IREO_GNSYM_TO_TS(IACOB) = IACOB_STA
            END IF
          END IF
        END DO! loop over IORB
      END DO ! Loop over IGNSYM
* 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Reorder array for active orbs: General => Type order '
        CALL IWRTMA3(IREO_GNSYM_TO_TS,1,NACOB,1,NACOB)
      END IF
*
      RETURN
      END
      SUBROUTINE LEN_GAS_GS_BLOCKS(LEN_GAS_GS,N_GAS_GS,IGAS_F,IGAS_L)
*
* Obtain in LEN_GAS_GS the dimension of each symmetry-gas block 
* for general symmetry
*
*. Jeppe Olsen, July 2012
*
*. Last revision: July 8, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Output
      INTEGER LEN_GAS_GS((NGAS+2)*NGENSMOB)
*
      NTEST = 000
*
      NBLK = 0
      DO IGENSM = 1, NGENSMOB
       DO IGAS = IGAS_F, IGAS_L
         NBLK = NBLK + 1
         LEN_GAS_GS(NBLK) = NGAS_GNSYM(IGENSM,IGAS)
       END DO
      END DO
      N_GAS_GS = NBLK
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of General symmetry gas-blocks ', NBLK
        WRITE(6,*) ' Length of general-symmetry gas-block '
        CALL IWRTMA3(LEN_GAS_GS,1,NBLK,1,NLBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE INVERT_REO(IREO,IREO_INV, NDIM)
* Obtain inverse mapping of reordering IREO(I)
* If K = IREO(I) then IREO_INV(K) = I
* 
* IREO_INV(IREO(I)) = I
*
*. Jeppe Olsen, Oct. 2, 2012
*. Last modification; Oct. 2, 2012; Jeppe Olsen, original version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IREO(NDIM)
*. Output
      INTEGER IREO_INV(NDIM)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from INVERT_REO '
      END IF
*
      DO I = 1, NDIM
        K = IREO(I)
        IREO_INV(K) = I
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' IREO and IREO_INV '
       CALL IWRTMA3(IREO,1,NDIM,1,NDIM)
       WRITE(6,*)
       CALL IWRTMA3(IREO_INV,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE COMB_TWO_REO(IREO3,IREO2,IREO1,NDIM)
*
* Two reorder arrays IREO1, IREO2 are given, combine these
* Obtain IREO3(I) = IREO2(IREO1(I)) 
*
*. Jeppe Olsen, Oct. 2, 2012
*. Last modification; Oct. 2, 2012; Jeppe Olsen, original version
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IREO1(NDIM), IREO2(NDIM)
*. Output
      INTEGER IREO3(NDIM)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from COMB_TWO_REO '
      END IF
      
      DO I = 1, NDIM
        IREO3(I) = IREO2(IREO1(I))
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input: IREO1, IREO2 '
        CALL IWRTMA3(IREO1,1,NDIM,1,NDIM)
        WRITE(6,*)
        CALL IWRTMA3(IREO2,1,NDIM,1,NDIM)
        WRITE(6,*) ' Output: IREO3 '
        CALL IWRTMA3(IREO3,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_IACT_TO_GENSM_REO(IACT_TO_GENSM_REO,
     &           ISTA_TO_GENSM_REO, MO_STA_TO_ACT_REO, NTOOB)
*
* Obtain reorder array going from orbitals in super-symmetry order to 
* actual supersymmetry order
*
*. Jeppe Olsen, Oct.2, 2012
*. Last revision; Oct. 3, 2012; Jeppe Olsen; debugged
*
      INCLUDE 'implicit.inc'
COLD  INCLUDE 'mxpdim.inc'
*. Input
      INTEGER ISTA_TO_GENSM_REO(NTOOB),MO_STA_TO_ACT_REO(NTOOB)
*. Output
      INTEGER IACT_TO_GENSM_REO(NTOOB)
*. local scratch
COLD  INTEGER IREO(MXPORB)
*. Array from actual order to standard order
COLD  INVERT_REO(IREO,IREO_INV, NDIM)
COLD  CALL INVERT_REO(MO_STA_TO_ACT_REO,IREO,NTOOB)
*. Combine IACT_TO_GENSM_REO(I) = MO_STA_TO_ACT_REO(ISTA_TO_GENSM_REO((I))
C     COMB_TWO_REO(IREO3,IREO2,IREO1,NDIM)
      CALL COMB_TWO_REO(IACT_TO_GENSM_REO, 
     &     MO_STA_TO_ACT_REO, ISTA_TO_GENSM_REO,NTOOB)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IACT_TO_GENSM_REO array '
        WRITE(6,*) ' ======================= '
        CALL IWRTMA3(IACT_TO_GENSM_REO,1,NTOOB,1,NTOOB)
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_CMO_SUP_SHL(CMO_SUP,CMO_SHL,IWAY)
*
* Reform CMO  matrix between supersymmetry order and shell order.
*. Outer routine
*
*. Jeppe Olsen, March 5, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
* IWAY = 1: Super-symmetry to shell order
* IWAY = 2: Shell-order to super-symmetry order
*
      CALL REFORM_CMO_SUP_SHL_IN(CMO_SUP,CMO_SHL,IWAY,
     &     int_mb(KNSUPSYM_FOR_IRREP),int_mb(KIBSUPSYM_FOR_IRREP),
     &     int_mb(KISUPSYM_FOR_IRREP))
*
      RETURN
      END
      SUBROUTINE REFORM_CMO_SUP_SHL_IN(CMO_SUP,CMO_SHL,IWAY,
     &     NSUPSYM_FOR_IRREP,IBSUPSYM_FOR_IRREP,
     &     ISUPSYM_FOR_IRREP)

*
* In shell-order the orbitals are arranged as
*
* Loop over irreducible representations
*  Loop over shells of this supersymmetry
*   Loop over subshells for a given shell 
*   End of loop over subshell of a given shell
*  End of loop over shells
* End of loop over irreps
*
*. Jeppe Olsen, March 5, 2013
*
* IWAY = 1: Super-symmetry to shell order
* IWAY = 2: Shell-order to super-symmetry order
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      INTEGER IBSUPSYM_FOR_IRREP(*),NSUPSYM_FOR_IRREP(*)
      INTEGER ISUPSYM_FOR_IRREP(*)
*
*. Input and output
      DIMENSION CMO_SUP(*), CMO_SHL(*)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from REFORM_CMO_SUP_SHL'
        WRITE(6,*) ' ============================'
        WRITE(6,*) 
        WRITE(6,*) ' IWAY = ', IWAY
C?      WRITE(6,*) ' First 10 elements of CMO_SUP '
C?      CALL WRTMAT(CMO_SUP,1,10,1,10)
      END IF
*
      IB_SHL = 1
      DO IRREP = 1, N_SUPSYM_IRREP
*. Number of shells for this irrep
        NSHELL = NBAS_SUPSYM(IBSUPSYM_FOR_IRREP(IRREP))
*
        IB = IBSUPSYM_FOR_IRREP(IRREP)
        NDEG = NSUPSYM_FOR_IRREP(IRREP)
        DO ISHELL = 1, NSHELL
        DO IISUPSYM = IB, IB + NDEG-1
         ISUPSYM = ISUPSYM_FOR_IRREP(IISUPSYM)
*. Offset to shell ISHELL of supersymmetry ISUPSYM
          IB_SUP = IB_CMOSUP_ORB(ISUPSYM,ISHELL)
          IF(IWAY.EQ.1) THEN
            CALL COPVEC(CMO_SUP(IB_SUP),CMO_SHL(IB_SHL),NSHELL)
C?          WRITE(6,*) ' Elements copied '
C?          CALL WRTMAT(CMO_SHL(IB_SHL),1,NSHELL,1,NSHELL)
          ELSE
            CALL COPVEC(CMO_SHL(IB_SHL),CMO_SUP(IB_SUP),NSHELL)
          END IF
C?        WRITE(6,'(A,3I4)') ' ISHELL IISUPSYM, ISHELL = ',
C?   &                 ISHELL IISUPSYM, ISHELL 
C?        WRITE(6,'(A,2I4)') ' IB_SHL, IB_SUP = ', IB_SHL, IB_SUP
          IB_SHL = IB_SHL + NSHELL
        END DO ! loop over supersymmetries
       END DO ! End of loop over shells
      END DO ! End of loop over irreps
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' C arranged according to shells '
        WRITE(6,*) ' ================================'
        WRITE(6,*)
        CALL PRINT_CSHELL(CMO_SHL)
      END IF
*
      RETURN
      END
      FUNCTION IB_CMOSUP_ORB(ISUPSYM,ISHELL)
*
* Determine offset of orbital with given supersymmetry and shell number
* in supersymmetry-ordered matrix
*
*. Jeppe Olsen, March 5, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*
      NTEST = 00
*. Offset to supersymmetry block
      IOFF = 1
      DO JSUPSYM = 1, ISUPSYM-1
        NSHELL = NBAS_SUPSYM(JSUPSYM)
        IOFF = IOFF + NSHELL**2
      END DO
*. And to the given shell in the block
      IOFF = IOFF + (ISHELL-1)*NBAS_SUPSYM(ISUPSYM)
*
      IB_CMOSUP_ORB = IOFF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from IB_CMOSUP_ORB: '
        WRITE(6,*) ' ISUPSYM, ISHELL, IB_CMOSUP_ORB = ',
     &               ISUPSYM, ISHELL, IB_CMOSUP_ORB
      END IF
*
      RETURN
      END
      SUBROUTINE PRINT_CSHELL(CSHELL)
*
* A MO-AO expansion CSHELL is given in SHELL ordered form. 
* Print it! 
*. Outer part
*. Jeppe Olsen, March 5, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION CSHELL(*)
*
      CALL PRINT_CSHELLIN(CSHELL,int_mb(KNSUPSYM_FOR_IRREP),
     &                    int_mb(KIBSUPSYM_FOR_IRREP))
*
      RETURN
      END
      SUBROUTINE PRINT_CSHELLIN(CSHELL,NSUPSYM_FOR_IRREP,
     &                          IBSUPSYM_FOR_IRREP)
*
* A MO-AO expansion CSHELL is given in SHELL ordered form. 
* Print it! 
*. Outer part
*. Jeppe Olsen, March 5, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
*
      DIMENSION  NSUPSYM_FOR_IRREP(*)
      DIMENSION IBSUPSYM_FOR_IRREP(*)
*
*. Jeppe Olsen, March 5, 2013
*
*. CMO to be printed
      DIMENSION CSHELL(*)
*
      IB = 1
      DO IRREP = 1, N_SUPSYM_IRREP
        WRITE(6,*) ' Irrep number ', IRREP
        WRITE(6,*) ' ====================='
        NDEG =  NSUPSYM_FOR_IRREP(IRREP)
        NSHELL = NBAS_SUPSYM(IBSUPSYM_FOR_IRREP(IRREP))
        DO ISHELL = 1, NSHELL
          WRITE(6,*) ' Subshells for shell ', ISHELL
          CALL WRTMAT(CSHELL(IB),NSHELL,NDEG,NSHELL,NDEG)
          IB = IB + NDEG*NSHELL
C?        WRITE(6,*) ' TESTY, NDEG, NSHELL, IB = ',
C?   &                        NDEG, NSHELL, IB
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE EXC_OO_TO_SS(NOOEX,IOOEX,NSSEX,ISSEX,
     &           NOOFSSX,IBOOFSSX,IOOFSSX)
*
* A set of NOOEXC orbital excitations are given by IOOEXC. Obtain the 
* corresponding shell excitations and obtain the mappings between the
* Orbital and shell excitations
*
*. Jeppe Olsen, March 5/6 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
*. Input
      RETURN
      END
      SUBROUTINE NONRED_SS_EXC(NOOEX,IOOEXC,NSSEX)
*
* A set of orbital excitation is given by IOOEXC
* Obtain the corresponding set of shell-excitations
*
*. Input:
*    NOOEX: Number of orbital excitations
*    IOOEX: Orbital excitations
*.Output: (mainly as pointers)
*    NSSEX: Number of shell excitations
*    KISSEXC: THe shell excitations in compact form
*    KNIOOFSS: Number of orbital excitations for a gvien shell excitation
*    KIBOOFSS: The offset for orbitial excitations
*    KIOOFSS: The actual orbital excitations for a given shell excitation
*   
* The information returned is pointers defined in the routine
*  
*
*. Jeppe Olsen, March 6, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOOEXC(2,NOOEX)
*
      IDUM = 0
*. No mark as allocated memory should be conserved
*
      CALL MEMMAN(KLSSEXE,NTOSH**2,'ADDL  ',1,'SS_EXE')
      CALL MEMMAN(KLACT_TO_STA,NTOOB,'ADDL  ',1,'REACST')
*. Obtain mapping from actual to standard mapping of orbitals from inverse
C     INVERT_MAP(MAP,MAPINV,NELMNT)
      CALL INVERT_MAP(int_mb(KMO_STA_TO_ACT_REO),int_mb(KLACT_TO_STA),
     &                NTOOB)

* 
*. Obtain the number of shell-shell excitations
*
C     GET_SSEXC(NOOEX,IOOEXC,NSSEX,ISSEXE,ISSEXC,
C    &           NOOFSS,IBOOFSS,IOOFSS,IFLAG,IACT_TO_STA,ISHL_FOR_STA)
      CALL GET_SSEXC(NOOEX,IOOEXC,NSSEX,int_mb(KLSSEXE),
     &                IDUM,IDUM,IDUM,IDUM,1,
     &                int_mb(KLACT_TO_STA),int_mb(KISHELL_FOR_BAS))
*. Allocate space for the shell excitations
      CALL MEMMAN(KISSEXC,2*NSSEX,'ADDL  ',1,'ISSEXC')
      CALL MEMMAN(KNOOFSS,NSSEX,'ADDL  ',1,'NOOFSS')
      CALL MEMMAN(KIBOOFSS,NSSEX,'ADDL  ',1,'BOOFSS')
      CALL MEMMAN(KIOOFSS,NOOEX,'ADDL  ',1,'IOOFSS')
*
*. And the actual shell-shell excitations
*
C     GET_SSEXC(NOOEX,IOOEXC,NSSEX,ISSEXE,ISSEXC,
C    &           NOOFSS,IBOOFSS,IOOFSS,IFLAG,IACT_TO_STA,ISHL_FOR_STA)
      CALL GET_SSEXC(NOOEX,IOOEXC,NSSEX,int_mb(KLSSEXE),
     &                int_mb(KISSEXC),int_mb(KNOOFSS),int_mb(KIBOOFSS),
     &                int_mb(KIOOFSS),2,
     &                int_mb(KLACT_TO_STA),int_mb(KISHELL_FOR_BAS))
*
      RETURN
      END
      SUBROUTINE GET_SSEXC(NOOEX,IOOEXC,NSSEX,ISSEXE,ISSEXC,
     &           NOOFSS,IBOOFSS,IOOFSS,IFLAG,IACT_TO_STA,ISHL_FOR_STA)
*
* IFLAG = 1: Obtain the allowed shell-shell excitation in expanded form in ISSEXE 
* IFLAG = 2: Use the ISSEXE array to obtain the actual components of the 
*            shell-excitations in compact form
*
*. Jeppe Olsen, March 6, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
*. Input
      INTEGER IOOEXC(2,NOOEX)
*. Mapping from actual(symmetry-type) order to standard order
      INTEGER IACT_TO_STA(*)
*. Shell number of given basis function/standard ordered MO
      INTEGER ISHL_FOR_STA(*)
*. Output
      INTEGER ISSEXE(NTOSH,NTOSH)
      INTEGER ISSEXC(2,*),NOOFSS(*),IBOOFSS(*),IOOFSS(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GET_SSEXC '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,*) ' IFLAG = ', IFLAG
        WRITE(6,*) ' NTOSH = ', NTOSH
        WRITE(6,*) ' NOOEX = ', NOOEX
      END IF
*
*. Set up array of shell-shell excitations
*
      IF(IFLAG.EQ.1) THEN
*
* Obtain the array ISSEXE(ISHL,JSHL) giving the number of 
* orbital excitations between these shells.
*
*. INITIALIZE
       IZERO = 0
       CALL ISETVC(ISSEXE,IZERO,NTOSH**2)
       DO JOOEX = 1, NOOEX
*. Orbitals in symmetry order
        IIORB = IREOTS(IOOEXC(1,JOOEX))
        JJORB = IREOTS(IOOEXC(2,JOOEX))
C?      WRITE(6,*) ' JOOEX, IIORB, JJORB = ', JOOEX, IIORB, JJORB
*. The same numbers in standard order
        IORB = IACT_TO_STA(IIORB)
        JORB = IACT_TO_STA(JJORB)
*. Shell numbers of these orbitals
        ISHL = ISHL_FOR_STA(IORB)
C?      WRITE(6,*) ' IORB, ISHL = ', IORB, ISHL
        JSHL = ISHL_FOR_STA(JORB)
C?      WRITE(6,*) ' JORB, JSHL = ', JORB, JSHL
*. Enroll
        ISSEXE(ISHL,JSHL) = ISSEXE(ISHL,JSHL) + 1
       END DO
*. Total number of shell-shell excitations
       NSSEX = 0
       DO ISHL = 1, NTOSH
        DO JSHL = 1, NTOSH
         IF(ISSEXE(ISHL,JSHL).NE.0) THEN
           NSSEX = NSSEX + 1
         END IF
        END DO
       END DO
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Shell-Shell excitation array '
         CALL IWRTMA3(ISSEXE,NTOSH,NTOSH,NTOSH,NTOSH)
       END IF
*
      END IF ! IFLAG = 1
*
      IF(IFLAG.EQ.2) THEN
*
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Shell-Shell excitation array(input now) '
         CALL IWRTMA3(ISSEXE,NTOSH,NTOSH,NTOSH,NTOSH)
        END IF
*. Obtain the shell-shell-excitations in compact form and pointer to start of OO for given SS
        ISSEX = 0
        DO ISHL = 1, NTOSH
         DO JSHL = 1, NTOSH
           IF(ISSEXE(ISHL,JSHL).NE.0) THEN
             ISSEX = ISSEX + 1
C?           WRITE(6,*) ' ISHL, JSHL, ISSEX = ', ISHL, JSHL, ISSEX
             ISSEXC(1,ISSEX) = ISHL
             ISSEXC(2,ISSEX) = JSHL
             NOOFSS(ISSEX) = ISSEXE(ISHL,JSHL)
           END IF
         END DO
        END DO
C?      WRITE(6,*) ' NOOFSS: '
C?      CALL IWRTMA3(NOOFSS,1,ISSEX,1,ISSEX)
*. Pointers to start of Orbital excitations for given shell excitation
        IBS = 1
        DO ISSEX = 1, NSSEX
          IBOOFSS(ISSEX) = IBS
          IBS = IBS + NOOFSS(ISSEX)
*. And clear for later use
          NOOFSS(ISSEX) = 0
        END DO
C?      WRITE(6,*) ' IBOOFSS: '
C?      CALL IWRTMA3(IBOOFSS,1,ISSEX,1,ISSEX)
*. Change ISSEXE to give the number of a given shell excitation
        ISSEX = 0
        DO ISHL = 1, NTOSH
         DO JSHL = 1, NTOSH
          IF(ISSEXE(ISHL,JSHL).NE.0) THEN
            ISSEX = ISSEX + 1
            ISSEXE(ISHL,JSHL) = ISSEX
          END IF
         END DO
        END DO
*. And the orbital excitations of a given shell excitation
        DO IOOEX = 1, NOOEX
*. Orbitals in symmetry-type order
         IIORB = IREOTS(IOOEXC(1,IOOEX))
         JJORB = IREOTS(IOOEXC(2,IOOEX))
*. The same numbers in standard order
         IORB = IACT_TO_STA(IIORB)
         JORB = IACT_TO_STA(JJORB)
*. Shell numbers of these orbitals
         ISHL = ISHL_FOR_STA(IORB)
         JSHL = ISHL_FOR_STA(JORB)
         IJEXC = ISSEXE(ISHL,JSHL)
         NOOFSS(IJEXC) = NOOFSS(IJEXC) + 1
         IB = IBOOFSS(IJEXC)
         IOOFSS(IB + NOOFSS(IJEXC)-1) = IOOEX
        END DO
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' The orbital excitations of shell excitations'
          WRITE(6,*) ' ============================================'
          DO ISSEX = 1, NSSEX
            WRITE(6,*)
            WRITE(6,*) ' Shell excitation ', ISSEX
            IB = IBOOFSS(ISSEX)
            N  = NOOFSS(ISSEX)
            CALL IWRTMA3(IOOFSS(IB),1,N,1,N) 
          END DO
        END IF ! NTEST large 
      END IF ! IFLAG = 2
*
   
*
      RETURN
      END
      SUBROUTINE INVERT_MAP(MAP,MAPINV,NELMNT)
*
* A mappping MAP(I) is given. Obtain inverse mapping
*
*. Jeppe Olsen, March 6, 2013
      INCLUDE 'implicit.inc'
*. Input
      INTEGER MAP(NELMNT)
*. Output
      INTEGER MAPINV(NELMNT)
*
      DO I = 1, NELMNT
       MAPINV(MAP(I)) = I
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Map (input) and inverted map(output)'
        WRITE(6,*) ' ====================================='
        WRITE(6,*)
        CALL IWRTMA3(MAP,1,NELMNT,1,NELMNT)
        CALL IWRTMA3(MAPINV,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
      SUBROUTINE SHELL_AVERAGE_ORBEXC(VECIN,NSSEX,NOOFSS,IBOOFSS,
     &           IOOFSS,VECUT,NOOEX,ICOPY)
*
*  A vector V over orbital excitation is given.
* Average over orbital excitations belonging to identical shell excitation.
*
* If ICOPY = 1, the output vector is copied to the input vector
*
*. Jeppe Olsen, March 7, 2013
*
      INCLUDE 'implicit.inc'
*. General input
      INTEGER NOOFSS(NSSEX),IBOOFSS(NSSEX),IOOFSS(*)
*. Specific input 
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECUT(*)
*
      NTEST = 00
*
      DO ISSEX = 1, NSSEX
*. Obtain average value for this shell-excitation
       IB = IBOOFSS(ISSEX)
       N  = NOOFSS(ISSEX)
       AVE = 0.0D0
       DO I = 1, N
        AVE = AVE + VECIN(IOOFSS(IB-1+I))
       END DO
       AVE = AVE/FLOAT(N) 
*. And spread out
       DO I = 1, N
         VECUT(IOOFSS(IB-1+I)) = AVE
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Input and output from averaging over shell-components '
        WRITE(6,*) 
     &  ' ======================================================'
        CALL WRT_2VEC(VECIN,VECUT,NOOEX)
      END IF
*
      IF(ICOPY.NE.0) THEN
        CALL COPVEC(VECUT,VECIN,NOOEX)
      END IF
*
      RETURN
      END
      SUBROUTINE AVE_DENS_OVER_SUBSHELLS(RHO1,RHO1AVE)
* 
* Average the one-particle density over subshells belonging to a given shell
* Outer routine
*
*. Jeppe Olsen, March 8, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input 
      DIMENSION RHO1(NACOB,NACOB)
*. Output
      DIMENSION RHO1AVE(NACOB,NACOB)
*
      CALL AVE_DENS_OVER_SUBSHELLS_IN(RHO1,RHO1VE,NSHTO,NACOB,
     &     int_mb(KNBAS_FOR_SHELL),int_mb(KIBBAS_FOR_SHELL),
     &     int_mb(KIBAS_FOR_SHELL),ITPFSO,NGAS)
*
      RETURN
      END
      SUBROUTINE AVE_DENS_OVER_SUBSHELLS_IN(RHO1,RHO1AVE,NSHTO,NACOB,
     &           NBAS_FOR_SHELL,IBBAS_FOR_SHELL,IBAS_FOR_SHELL,
     &           ITPFSO,NGAS)
*
* Average the one-particle density over subshells belonging to a given shell
*
*. We do not have a list of active shells. What we instead will do is to 
* connect a shell with the first subshell of this shell, and then we can check
* the gas space of this subshell
*
*
      INCLUDE 'implicit.inc'
      INTEGER NBAS_FOR_SHELL(*), IBBAS_FOR_SHELL(*)
      INTEGER IBAS_FOR_SHELL(*), ITPFSO(*)
*. Input
      DIMENSION RHO1(NACOB,NACOB)
*. Output
      DIMENSION RHO1AVE(NACOB,NACOB)
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from AVE_DENS_OVER_SUBSHELLS_IN '
      END IF
      WRITE(6,*) ' THIS ROUTINE HAS NEVER  BEEN TESTED OR USED '
*
      DO ISHELL = 1, NSHTO
       IB = IBBAS_FOR_SHELL(ISHELL)
       IORB = IBAS_FOR_SHELL(IB)
       IF(0.LT.ITPFSO(IORB).AND.ITPFSO(IORB).LE.NGAS) THEN
*. Shell is active
         DO JSHELL = 1, NSHTO
          JB = IBBAS_FOR_SHELL(JSHELL)
          JORB = IBAS_FOR_SHELL(JB)
          IF(0.LT.ITPFSO(JORB).AND.ITPFSO(JORB).LE.NGAS) THEN
*
* Obtain average value for these shells, non-diagonal and diagonal
*
            AVE = 0.0D0
            AVED = 0.0D0
*
            NI = NBAS_FOR_SHELL(ISHELL)
            NJ = NBAS_FOR_SHELL(JSHELL)
            DO IIORB = IB, IB - 1 + NI
             DO JJORB = JB, JB - 1 + NI
              IORB = IBAS_FOR_SHELL(IIORB)
              JORB = IBAS_FOR_SHELL(JJORB) 
*
              IF(IORB.NE.JORB) THEN
                AVE = AVE + RHO1(IORB,JORB)
              ELSE
                AVED = AVED + RHO1(IORB,JORB)
              END IF
             END DO
            END DO
*
            IF(ISHELL.EQ.JSHELL) THEN   
              AVED = AVED/FLOAT(NI)
              AVE = AVE/(FLOAT(NI)*(FLOAT(NI)-1))
            ELSE
              AVE = AVE/(FLOAT(NI)*(FLOAT(NI)))
            END IF
*
* Scatter average values out
*
            DO IIORB = IB, IB - 1 + NI
             DO JJORB = JB, JB - 1 + NI
              IORB = IBAS_FOR_SHELL(IIORB)
              JORB = IBAS_FOR_SHELL(JJORB) 
              RHO1AVE(IORB,JORB) = AVE
             END DO
            END DO
*
            IF(ISHELL.EQ.JSHELL) THEN
             DO IIORB = IB, IB-1+NI
              IORB = IBAS_FOR_SHELL(IIORB)
              RHO1AVE(IORB,IORB) = AVED 
             END DO
            END IF
*
          END IF
         END DO ! For jshell
       END IF
      END DO ! For Ishell
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Original and averaged density matrix '
        WRITE(6,*) ' ====================================='
        WRITE(6,*)
        CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*)
        CALL WRTMAT(RHO1AVE,NACOB,NACOB,NACOB,NACOB)
      END IF
*
      RETURN 
      END
      SUBROUTINE AVE_SUPSYM_MAT(ASUP,NOBPSPSM,IPACK)
*
* A matrix ASUP  over blocks of supersymmetry is given
* Average over blocks belonging to the same irrep.
*
*. Outer part
*
*. Jeppe Olsen, March 9, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'  
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*. Input
      DIMENSION NOBPSPSM(*)
*. Input and output
      DIMENSION ASUP(*)
*
      CALL AVE_SUPSYM_MAT_IN(ASUP,NOBPSPSM,IPACK,
     &     N_SUPSYM_IRREP,N_SUPSYM,int_mb(KNSUPSYM_FOR_IRREP),
     &     int_mb(KIBSUPSYM_FOR_IRREP),int_mb(KISUPSYM_FOR_IRREP))
*
      RETURN
      END
      SUBROUTINE AVE_SUPSYM_MAT_IN(ASUP,NOBPSPSM,IPACK,
     &     N_SUPSYM_IRREP,N_SUPSYM,NSUPSYM_FOR_IRREP,
     &     IBSUPSYM_FOR_IRREP,ISUPSYM_FOR_IRREP)
*
* A matrix ASUP over supersymmmetries is given with NOBPSPSM orbitals per supersymmetry
* Average over supersymmetries belonging to the same super-symmetry irrep
*
*. Jeppe Olsen, March 9, 2013
*
*. General input
      INCLUDE 'implicit.inc'
      INTEGER NSUPSYM_FOR_IRREP(*), IBSUPSYM_FOR_IRREP(*)
      INTEGER ISUPSYM_FOR_IRREP(*)
*.
      INTEGER NOBPSPSM(*)
*. Input and output
      DIMENSION ASUP(*)
*
      NTEST = 100
      ONE = 1.0D0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from AVE_SUPSYM_MAT_IN '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' N_SUPSYM = ', N_SUPSYM
        WRITE(6,*) ' NOBPSPSM: '
        CALL IWRTMA(NOBPSPSM,1,N_SUPSYM,1,N_SUPSYM)
        WRITE(6,*) ' Input matrix over supersymmetries '
        WRITE(6,*)
        CALL APRBLM2(ASUP,NOBPSPSM,NOBPSPSM,N_SUPSYM,IPACK)
      END IF

*
      DO IRREP = 1, N_SUPSYM_IRREP
       NSPSM = NSUPSYM_FOR_IRREP(IRREP)
       IF(NTEST.GE.1000) WRITE(6,*) ' Info for IRREP = ', IRREP
       IF(NSPSM.GT.1) THEN
*
*. Average over the various components in the first supersymmetry of the  given irrep
*
        IB = IBSUPSYM_FOR_IRREP(IRREP)
        N  = NSUPSYM_FOR_IRREP(IRREP)
        ISUPSYM1 = ISUPSYM_FOR_IRREP(IB)
        L1 = NOBPSPSM(ISUPSYM1)
        IF(IPACK.EQ.0) THEN
          LBLK = L1*L1
        ELSE
          LBLK = L1*(L1+1)/2
        END IF
        IOFF1 = IOFF_BLCK(ISUPSYM1,NOBPSPSM,NOBPSPSM,IPAK)
        WRITE(6,*) ' ISUPSYM1, IOFF1 = ', ISUPSYM1, IOFF1
*. And terms from the remaining supersymmetries
        DO IISUPSYM = 2, N
          ISUPSYM = ISUPSYM_FOR_IRREP(IB-1+IISUPSYM)
          IOFF = IOFF_BLCK(ISUPSYM,NOBPSPSM,NOBPSPSM,IPAK)
          CALL VECSUM(ASUP(IOFF1),ASUP(IOFF1),ASUP(IOFF),ONE,ONE,LBLK)
        END DO
        FACTOR = 1.0D0/FLOAT(N)
C?      WRITE(6,*) ' FACTOR = ', FACTOR
        CALL SCALVE(ASUP(IOFF1),FACTOR,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Averaged block '
          CALL APRBLM2(ASUP(IOFF1),L1,L1,1,IPACK)
        END IF
*
*. And copy the average to the remaining blocks
*
        DO IISUPSYM = 2, N
          ISUPSYM = ISUPSYM_FOR_IRREP(IB-1+IISUPSYM)
          IOFF = IOFF_BLCK(ISUPSYM,NOBPSPSM,NOBPSPSM,IPAK)
          CALL COPVEC(ASUP(IOFF1),ASUP(IOFF),LBLK)
        END DO
*
       END IF ! Irrep was degenerate
      END DO ! Loop over irreps
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Averaged matrix '
        CALL APRBLM2(ASUP,NOBPSPSM,NOBPSPSM,N_SUPSYM,IPACK)
      END IF
*
      RETURN
      END
      FUNCTION IOFF_BLCK(IBLK,LR,LC,IPAK)
*
* Offset to block IBLK in matrix with LR/LC row/colomn elements per block
*
* Jeppe Olsen, March 9, (did it really take me about 25 years to write this function)
*
      INCLUDE 'implicit.inc'
*
      INTEGER LR(*), LC(*)
*
      IOFF = 1
      DO JBLK = 1, IBLK - 1
       IF(IPAK.EQ.0) THEN
         LBLK = LR(JBLK)*LC(JBLK)
       ELSE
         LBLK = LR(JBLK)*(LR(JBLK)+1)/2
       END IF
       IOFF = IOFF + LBLK
      END DO
*
      IOFF_BLCK = IOFF
*
      RETURN
      END
      SUBROUTINE PRINT_CMO_AS_SHELLS(CMO,IFORM)
*
* A CMOA matrix CMO is given in form defined by CMO. 
* Print as subshells of a given shell
*
* IFORM = 1: Input CMO is in standard order
* IFORM = 2: Input CMO is in actual(gas) order
* IFORM = 3: Input CMO is in supersymmetry-order
* IFORM = 4: Input CMO is in Shell-order
*
*. Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION CMO(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PR_SHL')
*
      LCMO = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0) 
      CALL MEMMAN(KLMO1,LCMO,'ADDL  ',2,'LMO1  ')
*
*. Reform to shell format 
*
      CALL REFORM_CMO(CMO,IFORM,dbl_mb(KLMO1),4)
*
*. And print
*
      CALL PRINT_CSHELL(dbl_mb(KLMO1))
*

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PR_SHL')
*
      RETURN
      END
      SUBROUTINE REFORM_CMO(C_IN,IFORM_IN, C_OUT, IFORM_OUT)
*
* Reform matrix from FORM IFORM_IN to FORM IFORM_OUT
*
* IFORM = 1 =>  standard form
* IFORM = 2 =>  gas-ordered form
* IFORM = 3 =>  super-symmetry form
* IFORM = 4 =>  shell form
*
*. Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION C_IN(*)
*. Output
      DIMENSION C_OUT(*)
*
      NTEST = 000
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from REFORM_CMO '
        WRITE(6,*) ' ==================== '
        WRITE(6,*) 
        WRITE(6,*) ' IFORM_IN,IFORM_OUT = ', IFORM_IN, IFORM_OUT
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'REFCMO')
      LCMO = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0) 
      CALL MEMMAN(KLMO1,LCMO,'ADDL  ',2,'LMO1  ')
*
*. Reform from input form to standard form
*
C REFORM_CMO_TO_STANDARD(CMO,CMOST,IFORM,IWAY)
      CALL REFORM_CMO_TO_STANDARD(C_IN,dbl_mb(KLMO1),IFORM_IN,1)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Intermediate CMO in standard form '
        CALL APRBLM2(dbl_mb(KLMO1),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
*. Reform from standard to output form
*
      CALL REFORM_CMO_TO_STANDARD(C_OUT,dbl_mb(KLMO1),IFORM_OUT,2)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input CMO to REFORM_CMO: '
        CALL PRINT_CMO_ARBFORM(C_IN,IFORM_IN)
        WRITE(6,*) ' Output CMO from REFORM_CMO: '
        CALL PRINT_CMO_ARBFORM(C_OUT,IFORM_OUT)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'REFCMO')
*
      RETURN
      END
      SUBROUTINE REFORM_CMO_TO_STANDARD(CMO,CMOST,IFORM,IWAY)
*
* Reform between various super-symmetry forms of a CMO matrix.
*
* A CMO matrix, CMO is given in a from specified by IFORM. Reform this to 
* standard form and save in CMOST
*
* IFORM = 1 => CMO is in standard form
* IFORM = 2 => CMO is in gas-ordered form
* IFORM = 3 => CMO is in super-symmetry form
* IFORM = 4 => CMO is in shell form
*
* IWAY = 1: From the general form in CMO to CMOST
* IWAY = 2: From the standard form in CMOST to CMO
*
*. Jeppe Olsen, March 9, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input or output
      DIMENSION CMO(*), CMOST(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from REFORM_CMO_TO_STANDARD' 
        WRITE(6,*) ' ================================'
        WRITE(6,*)
        WRITE(6,*) ' Iway and Iform = ', IWAY, IFORM
        WRITE(6,*) ' Input matrix: '
        IF(IWAY.EQ.1) THEN
C              PRINT_CMO_ARBFORM(CMO,IFORM)
          CALL PRINT_CMO_ARBFORM(CMO,IFORM)
        ELSE
          CALL PRINT_CMO_ARBFORM(CMOST,1)
       END IF
*
      END IF !Ntest is large
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'RFSPCM')
      LCMO = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0) 
      CALL MEMMAN(KLMO1,LCMO,'ADDL  ',2,'LMO1  ')
* 
      IF(IWAY.EQ.1) THEN
*
* CMO in some form => CMOST in standard form
*
*
        IF(IFORM.EQ.4) THEN
*. Shell in CMO => supersymmetry form in KLMO1
C              REFORM_CMO_SUP_SHL(CMO_SUP,CMO_SHL,IWAY)
          CALL REFORM_CMO_SUP_SHL(dbl_mb(KLMO1),CMO,2)
        END IF
*
        IF (IFORM.GE.3) THEN
          IF(IFORM.EQ.3) CALL COPVEC(CMO,dbl_mb(KLMO1),LCMO)
*. Supersymmetry in KLMO1 => standard in CMOST
C                REFORM_CMO_STA_GEN(CMO_STA,CMO_GEN,IDO_REORDER,IREO,IWAY)
           CALL  REFORM_CMO_STA_GEN(CMOST,dbl_mb(KLMO1),0,0,2)
        ELSE IF (IFORM.EQ.2) THEN
*. Actual => standard
C               REO_CMOAO(CIN,COUT,IREO,ICOPY,IWAY)
           CALL REO_CMOAO(CMO,CMOST,
     &          int_mb(KMO_STA_TO_ACT_REO),0,2)
        ELSE IF(IFORM.EQ.1) THEN
*. Easy living, just copy
          CALL COPVEC(CMO,CMOST,LCMO)
        END IF
      ELSE IF (IWAY.EQ.2) THEN
*
* CMOST in standard form => CMO in some form 
*
         IF(IFORM.GE.3) THEN
* Standard in CMOST to supersymmetry in CMO
           CALL REFORM_CMO_STA_GEN(CMOST,CMO,0,0,1)
         END IF
         IF(IFORM.EQ.4) THEN
           CALL COPVEC(CMO,dbl_mb(KLMO1),LCMO)
*. Supersymmetry to shell
C               REFORM_CMO_SUP_SHL(CMO_SUP,CMO_SHL,IWAY)
           CALL REFORM_CMO_SUP_SHL(dbl_mb(KLMO1),CMO,1)
         END IF
*
         IF(IFORM.EQ.2) THEN
*. Standard => Actual/gas ordered 
C               REO_CMOAO(CIN,COUT,IREO,ICOPY,IWAY)
           CALL REO_CMOAO(CMOST,CMO,
     &          int_mb(KMO_STA_TO_ACT_REO),0,1)
         ELSE IF(IFORM.EQ.1) THEN
*. Standard => standard
           CALL COPVEC(CMOST,CMO,LCMO)
         END IF
      END IF ! switch IWAY
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output matrix from REFORM_CMO_TO_STANDARD: '
        IF(IWAY.EQ.2) THEN
          WRITE(6,*) ' Output matrix is of general type'
C              PRINT_CMO_ARBFORM(CMO,IFORM)
          CALL PRINT_CMO_ARBFORM(CMO,IFORM)
        ELSE
          WRITE(6,*) ' Output matrix is STANDARD type '
          CALL PRINT_CMO_ARBFORM(CMOST,1)
       END IF
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'RFSPCM')
*
      RETURN
      END
      SUBROUTINE PRINT_CMO_ARBFORM(CMO,IFORM)
*
* Print a CMO matrix with form defined by IFORM
*
*. IFORM = 1: Standard form
*. IFORM = 2: Actual/gas ordered
*. IFORM = 3: Supersymmetry ordered
*. IFORM = 4: Shell ordered
*
*. Jeppe Olsen, March 10, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      DIMENSION CMO(*)
*
      IF(IFORM.EQ.1.OR.IFORM.EQ.2) THEN
       CALL PRINT_CMOAO(CMO)
      ELSE IF(IFORM.EQ.3) THEN
       CALL APRBLM2(CMO,NBAS_SUPSYM,NBAS_SUPSYM,N_SUPSYM)
      ELSE IF(IFORM.EQ.4) THEN
       CALL PRINT_CSHELL(CMO)
      ELSE
       WRITE(6,*) ' PRINT_CMO_ARBFORM: Unknown IFORM = ', IFORM
       STOP       ' PRINT_CMO_ARBFORM: Unknown IFORM  '
      END IF
*
      RETURN
      END
      SUBROUTINE REO_2SUPSYM_ORDERS(ISUPSYM1,ISUPSYM2,IREO12)
* Two sequences of orbitals defined by their supersymmetries are given
*
* Obtain reordering array IREO12(I1) = I2
*
* Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER ISUPSYM1(NTOOB),ISUPSYM2(NTOOB)
*. Output 
      INTEGER IREO12(NTOOB)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output REO_2SUPSYM_ORDERS' 
       WRITE(6,*) ' Supersymmetry lists 1 and 2 '
       CALL IWRTMA3(ISUPSYM1,1,NTOOB,1,NTOOB)
       CALL IWRTMA3(ISUPSYM2,1,NTOOB,1,NTOOB)
      END IF
*
      DO ISUPSYM = 1, N_SUPSYM
        NORB1 = 0
        DO IORB1 = 1, NTOOB
         IF(ISUPSYM1(IORB1).EQ.ISUPSYM) THEN
           NORB1 = NORB1 + 1
*. Find orbital NORB1 of symmetry ISUPSYM in ISUPSYM2
           NORB2 = 0
           DO IIORB2 = 1, NTOOB
             IF(ISUPSYM2(IIORB2).EQ.ISUPSYM) THEN
               NORB2 = NORB2 + 1
               IF(NORB2.EQ.NORB1) THEN
                 IORB2 = IIORB2
               END IF
             END IF
           END DO
           IREO12(IORB1) = IORB2
         END IF  
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Obtained redorder array '
        CALL IWRTMA3(IREO12,1,NTOOB,1,NTOOB)
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_OCC_ORDER_SUPSYM(IMO_OCCORD_SUPSYM)
*
* Obtain the super-symmetry for the required order of orbitals 
* in ocupation/gas order
*
*. Jeppe Olsen, March 2013
*
*
*.
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Output
      INTEGER IMO_OCCORD_SUPSYM(NTOOB)
*
      NTEST = 100 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GET_OCC_ORDER_SUPSYM '
        WRITE(6,*) ' ================================='
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NGAS_SUPSYM: '
        CALL IWRTMA3(NGAS_SUPSYM,N_SUPSYM,NGAS+2,MXP_NSUPSYM,NGAS+2)
      END IF
*
      IZERO = 0
      CALL ISETVC(IMO_OCCORD_SUPSYM,IZERO,NTOOB)
*
      IORB = 0
      DO ISYM = 1, NSMOB
        NSUP = NSUP_FOR_STA_SYM(ISYM)
        IBSUP = IBSUP_FOR_STA_SYM(ISYM)
        IF(NTEST.GE.1000) WRITE(6,*) ' ISYM, NSUP = ', ISYM, NSUP
*. Loop over the spac division for this symmetry
        DO JGAS= 0, NGAS + 1
          DO IISUPSYM = IBSUP, IBSUP + NSUP -1
            ISUPSYM = ISUP_FOR_STA_SYM(IISUPSYM)
            NORB = NGAS_SUPSYM(ISUPSYM,JGAS)
            IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' ISUPSYM, JGAS, NORB = ',
     &                   ISUPSYM, JGAS, NORB 
            END IF
            DO IIORB = 1, NORB
              IORB = IORB + 1
              IMO_OCCORD_SUPSYM(IORB) = ISUPSYM
              IF(NTEST.GE.10000)
     &        WRITE(6,*) ' IORB, ISUPSYM ', IORB, ISUPSYM
            END DO! loop over orbital IIORB
          END DO! Loop over supersymmetries  for given standard
        END DO ! Loop over orbital spaces
      END DO ! Loop over standard symmetries
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' GET_OCC_ORDER: Supersymmetry of occ-ordered orbitals'
        CALL IWRTMA3(IMO_OCCORD_SUPSYM,1,NTOOB,1,NTOOB)
      END IF
*
      RETURN
      END 
      SUBROUTINE ANA_SUBSHELLS_CMO(CMO,IFORM,XMAX,MAXIRR,MAXSHL,
     &                              IALIGN)
*
* A CMOA matrix CMO is given in form defined by CMO. 
* Analyize the matrix for differences between shells belonginh
* to the same shell, and align if required.
*
* IFORM = 1: Input CMO is in standard order
* IFORM = 2: Input CMO is in actual(gas) order
* IFORM = 3: Input CMO is in supersymmetry-order
* IFORM = 4: Input CMO is in Shell-order
*
*. Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION CMO(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'AN_SHL')
*
      LCMO = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0) 
      CALL MEMMAN(KLMO1,LCMO,'ADDL  ',2,'LMO1  ')
*
*. Reform to shell format 
*
      CALL REFORM_CMO(CMO,IFORM,dbl_mb(KLMO1),4)
*
*. And analyze /align
*
      CALL ANA_SUBSHELLS_CMO_IN(dbl_mb(KLMO1),
     &     int_mb(KNSUPSYM_FOR_IRREP),int_mb(KIBSUPSYM_FOR_IRREP),
     &     XMAX,MAXIRR,MAXSHL,IALIGN)
*
      IF(IALIGN.EQ.1) THEN
*. Pump aligned MOs back in CMO
        CALL REFORM_CMO(dbl_mb(KLMO1),4,CMO,IFORM)
      END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'AN_SHL')
*
      RETURN
      END
      SUBROUTINE ANA_SUBSHELLS_CMO_IN(CSHELL,NSUPSYM_FOR_IRREP,
     &           IBSUPSYM_FOR_IRREP,XMAX,MAXIRR,MAXSHL,IALIGN)
*
* A MO-AO expansion CSHELL is given in SHELL ordered form. 
*
* Find largest deviation between equivalent subshells
* align the subshells if IALIGN = 1
*
*. Jeppe Olsen, March 5, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
*
      INTEGER  NSUPSYM_FOR_IRREP(*)
      INTEGER  IBSUPSYM_FOR_IRREP(*)
*
*. Jeppe Olsen, March, 2013
*
*. CMO to be analyzed
      DIMENSION CSHELL(*)
*
      NTEST = 100
*
      IB = 1
      XMAX = -1.0D0
      DO IRREP = 1, N_SUPSYM_IRREP
        NDEG =  NSUPSYM_FOR_IRREP(IRREP)
        NSHELL = NBAS_SUPSYM(IBSUPSYM_FOR_IRREP(IRREP))
C?      WRITE(6,*) ' IRREP, NDEG, NSHELL = ', IRREP, NSHELL, NSHELL
        DO ISHELL = 1, NSHELL
         IBS = IB + (ISHELL-1)*NSHELL*NDEG
*. Check differences to first sub shell
          DO ISUB = 2, NDEG
           DO IBAS = 1, NSHELL
             DIF = CSHELL(IBS-1+(ISUB-1)*NSHELL+IBAS) -
     &             CSHELL(IBS-1+(1   -1)*NSHELL+IBAS)
C?           WRITE(6,*) ' IR,ISH,IB, DIF = ',
C?   &       IRREP,ISHELL,IBAS,DIF
C?           WRITE(6,*) 'CSHELL1, CSHELLI =',
C?   &       CSHELL(IBS-1+(1-1)*NSHELL+IBAS),
C?   &       CSHELL(IBS-1+(ISUB-1)*NSHELL+IBAS)
             IF(ABS(DIF).GT.XMAX) THEN
               XMAX = ABS(DIF)
               MAXIRR = IRREP
               MAXSHL= ISHELL
               MAXBAS = IBAS
             END IF
           END DO
          END DO
        END DO
*
        IF(IALIGN.EQ.1) THEN
*. Align: Copy the first subshell to the remaining
          DO ISHELL = 1, NSHELL
           IBS = IB + (ISHELL-1)*NSHELL*NDEG
           IB1 = IBS
           DO ISUB = 2, NDEG
             IBI =IBS+(ISUB-1)*NSHELL
             CALL COPVEC(CSHELL(IB1),CSHELL(IBI),NSHELL)
           END DO
          END DO
        END IF
*
        IB = IB + NDEG*NSHELL*NSHELL
C?      WRITE(6,*) ' TESTY, NDEG, NSHELL, IB = ',
C?   &                      NDEG, NSHELL, IB
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ANA_SUBSHELLS_ ...' 
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Largest difference of subshells: ',XMAX
        WRITE(6,*) ' Occurs for, IRREP, shell, basis func ', 
     &              MAXIRR,MAXSHL, MAXBAS
      END IF
*
      RETURN
      END


c $Id$
