      SUBROUTINE PERT_SUBSPACE(NPERT,H0,V,S,ECORE)
*
* Perturbation calculation has been performed
* giving a subspace of correction vectors.
*
* Analyze this subspace 
*
* H0, V and S are matrices in subspace
*
*. Jeppe Olsen, July 98
*.              Summer of 99 : Improved stability of orthogonalization
*
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      DIMENSION S(*),H0(*),V(*)
*
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'PERT_S')
*
      NTEST = 10
      IF(NTEST.GE.1) THEN
         WRITE(6,*)
         WRITE(6,*) ' Welcome to pert_subspace '
         WRITE(6,*) ' ======================== '
         WRITE(6,*)
      END IF
*. Dimension of subspace : number of correction vectors + 1
      NDIM = NPERT + 1
*. Dump S, H0 and V on LU98 in format E22.15
      IDUMP_TO98 = 1
      IF(IDUMP_TO98.EQ.1) THEN
        LU98 = 98
        CALL REWINO(LU98)
*.S
        DO IJ = 1, NDIM*(NDIM+1)/2
          WRITE(LU98,'(E22.15)') S(IJ)
        END DO
*.H0
        DO IJ = 1, NDIM*(NDIM+1)/2
          WRITE(LU98,'(E22.15)') H0(IJ)
        END DO
*.V 
        DO IJ = 1, NDIM*(NDIM+1)/2
          WRITE(LU98,'(E22.15)') V(IJ)
        END DO
      END IF
* 
*. A bit of local memory
      LENNY = NDIM ** 2   
      IDUM = 0
      CALL MEMMAN(KLMAT1,LENNY,'ADDL  ',2,'MAT1  ')
      CALL MEMMAN(KLMAT2,LENNY,'ADDL  ',2,'MAT2  ')
      CALL MEMMAN(KLMAT3,3*LENNY,'ADDL  ',2,'MAT3  ')
      CALL MEMMAN(KLMAT4,LENNY,'ADDL  ',2,'MAT4  ')
      CALL MEMMAN(KLMAT5,LENNY,'ADDL  ',2,'MAT5  ')
      CALL MEMMAN(KLMAT6,LENNY,'ADDL  ',2,'MAT6  ')
      CALL MEMMAN(KLMAT7,LENNY,'ADDL  ',2,'MAT6  ')
      CALL MEMMAN(KLVEC1,NDIM ,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,NDIM ,'ADDL  ',2,'VEC2  ')
      CALL MEMMAN(KLVEC3,NDIM ,'ADDL  ',2,'VEC3  ')
*
*. Orthonormalize using symmetric orthogonalization or modified GS
*
*
      I_SYM_OR_MGS = 2
*. Metric in complete matrix form
      WRITE(6,*) ' Input S matrix '
      CALL TRIPAK(dbl_mb(KLMAT5),S,2,NDIM,NDIM)
      CALL WRTMAT(dbl_mb(KLMAT5),NDIM,NDIM,NDIM,NDIM)
      IF( I_SYM_OR_MGS .EQ.1 ) THEN
*. S**(-1/2)
C           SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)      
        CALL SQRTMT(dbl_mb(KLMAT5),NDIM,2,dbl_mb(KLMAT2),
     &             dbl_mb(KLMAT1),dbl_mb(KLMAT3))
        IF(NTEST.GE.5) THEN
          WRITE(6,*) ' S-1/2 matrix '
          CALL WRTMAT(dbl_mb(KLMAT1),NDIM,NDIM,NDIM,NDIM)
        END IF
       ELSE
*. Modified Gram-Schmidt
         CALL MGS3(dbl_mb(KLMAT1),dbl_mb(KLMAT5),NDIM,dbl_mb(KLMAT2))
       END IF
*. Transform H0 and V to orthogonal basis
C          TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
      CALL TRAN_SYM_BLOC_MAT(H0,dbl_mb(KLMAT1),1,NDIM,dbl_mb(KLMAT2),
     &                       dbl_mb(KLMAT3))
      CALL COPVEC(dbl_mb(KLMAT2),H0,NDIM*(NDIM+1)/2)
      CALL TRAN_SYM_BLOC_MAT(V ,dbl_mb(KLMAT1),1,NDIM,dbl_mb(KLMAT2),
     &                       dbl_mb(KLMAT3))
      CALL COPVEC(dbl_mb(KLMAT2),V ,NDIM*(NDIM+1)/2)
*
      WRITE(6,*) ' H0 in orthonormal basis '
      CALL PRSYM (H0,NDIM)
      WRITE(6,*) ' V  in orthonormal basis '
      CALL PRSYM (V ,NDIM)
*. Find Metrix in orthonormal basis to check for inaccuracies
      CALL TRAN_SYM_BLOC_MAT(S ,dbl_mb(KLMAT1),1,NDIM,
     &                       dbl_mb(KLMAT2),dbl_mb(KLMAT3))
      WRITE(6,*) ' S in orthonormal basis '
      CALL PRSYM(dbl_mb(KLMAT2),NDIM)
      
*. Find basis where H0 is diagonal and transform
*. Diagonalize H0, eigenvectors in MAT2
      CALL EIGENL(H0,dbl_mb(KLMAT2),NDIM,0,1)
      CALL COPDIA(H0,dbl_mb(KLMAT3),NDIM,1)
*. And put back
      ZERO = 0.0D0
      CALL SETVEC(H0,ZERO,NDIM*(NDIM+1)/2)
      DO I = 1, NDIM
        H0(I*(I+1)/2) = dbl_mb(KLMAT3-1+I)
      END DO
*. Transform  V to basis that diagonalizes H0
      CALL TRAN_SYM_BLOC_MAT(V ,dbl_mb(KLMAT2),1,NDIM,dbl_mb(KLMAT4),
     &                       dbl_mb(KLMAT3))
      CALL COPVEC(dbl_mb(KLMAT4),V,NDIM*(NDIM+1)/2)
*
      WRITE(6,*) ' H0 in basis of sub space zero order states '
      WRITE(6,*) ' ========================================== '
      WRITE(6,*)
      CALL PRSYM(H0,NDIM)
      WRITE(6,*)
      WRITE(6,*) ' V  in basis of sub space zero order states '
      WRITE(6,*) ' ========================================== '
      CALL PRSYM(V,NDIM)
*. Eigenvalues H = H0+V, eigenvectors in MAT4
      ONE = 1.0D0
      CALL VECSUM(dbl_mb(KLMAT3),H0,V,ONE,ONE,NDIM*(NDIM+1)/2)
      CALL EIGENL(dbl_mb(KLMAT3),dbl_mb(KLMAT4),NDIM,1,1)
      CALL COPDIA(dbl_mb(KLMAT3),dbl_mb(KLVEC1),NDIM,1)
      DO I = 1, NDIM
        dbl_mb(KLVEC1-1+I) = dbl_mb(KLVEC1-1+I) + ECORE 
      END DO
*.
      WRITE(6,*)
      WRITE(6,*) ' Eigenvalues of H (with core-energy)in subspace '
      WRITE(6,*) ' =============================================== '
      WRITE(6,*)
      CALL WRTMAT(dbl_mb(KLVEC1),NDIM,1,NDIM,1)
*
*. Perturbation expansion in subspace 
*
*. Expand H0 and V  to complete matrices
C       TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(dbl_mb(KLMAT5),H0,2,NDIM,NDIM)
      CALL TRIPAK(dbl_mb(KLMAT6),V ,2,NDIM,NDIM)
*
      MAXORD = 500
      LEN = NDIM*(1+MAXORD)
      CALL MEMMAN(KLC,LEN,'ADDL  ',2,'KLC   ')
      CALL MEMMAN(KLEN,MAXORD+1,'ADDL  ',2,'KLC   ')
*. Zero order state
      ZERO = 0.0D0
      CALL SETVEC(dbl_mb(KLC),ZERO,NDIM)
      dbl_mb(KLC) = ONE
C          MATPERT(H0,V,NDIM,NORD,EN,C,VEC1,VEC2,VEC3)
      CALL MATPERT(dbl_mb(KLMAT5),dbl_mb(KLMAT6),NDIM,MAXORD,
     &             dbl_mb(KLEN),dbl_mb(KLC),dbl_mb(KLVEC1),
     &             dbl_mb(KLVEC2),dbl_mb(KLVEC3),ECORE )
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PERT_S')
*
      RETURN
      END
      SUBROUTINE MGS3(X,S,NDIM,SCR1)
*
* Modified Gram-Schmidt procedure by forward orthogonalization
*
*  watch out for zero columns indicating linear dependency
*
* Jeppe Olsen, Summer of 99
*
* S is input overlap matrix, X is output set of orthonormalized vectors
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. input 
      DIMENSION S(NDIM,NDIM)
*. Output
      DIMENSION X(NDIM,NDIM)
*. Scratch : vector of length NDIM
      DIMENSION SCR1(*)
*
*. Initialize X to unit matrix
*
      ZERO = 0.0D0
      ONE = 1.0D0
      CALL SETVEC(X,ZERO,NDIM**2)
      CALL SETDIA(X,ONE,NDIM,0)     
C          SETDIA(MATRIX,VALUE,NDIM,IPACK)
*
      DO IVEC = 1, NDIM
*. Normalize vector IVEC
        CALL MATVCB(S,X(1,IVEC),SCR1,NDIM,NDIM,0)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
*. avoid NaN's by putting norm to at least zero
        XNORM = INPROD(X(1,IVEC),SCR1,NDIM)
        
        IF (XNORM.LE.0D0) THEN
          FACTOR = 0.0D0
        ELSE
          FACTOR = 1.0D0/SQRT(XNORM)
        END IF
        CALL SCALVE(X(1,IVEC), FACTOR, NDIM)
        CALL SCALVE(SCR1,FACTOR,NDIM)
*. Subtract X(1,IVEC) from all remaining vectors
        DO JVEC = IVEC+1,NDIM
          XSX = INPROD(SCR1,X(1,JVEC),NDIM)
          CALL VECSUM(X(1,JVEC),X(1,JVEC),X(1,IVEC),ONE,-XSX,NDIM) 
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthogonalization information '
        WRITE(6,*) ' Overlap matrix '
        CALL WRTMAT(S,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Orthonormalized vectors '
        CALL WRTMAT(X,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE MATPERT(H0,V,NDIM,NORD,EN,C,VEC1,VEC2,VEC3,ECORE)
*
* Perturbation expansion of simple  eigenvalue problem 
*
* Explicit matrix version  
*
*
*. Obtain corrections to energy and wawe functions 
*
*. The normalization condition used is C(K)T  C(0) = 0
*
* The energy corrections are   
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I) 
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
*
*
* Jeppe Summer of 98           
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION H0(NDIM**2),V(NDIM**2)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(0:NORD)
*. Scratch 
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
*
*. Zero order energy
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(H0,C,VEC1,NDIM,NDIM,0)
      E0   = INPROD(VEC1,C,NDIM)
*
      WRITE(6,*) 'E0  = ', E0   
      EN(0) = E0    
*. Save diagonal of H0 - E(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = H0((I-1)*NDIM+I)-E0
      END DO
C?    WRITE(6,*) ' Zero order diagonal '
C?    CALL WRTMAT(VEC3,1,NDIM,1,NDIM)
*. And then start the iterations
      DO IORD = 1, NORD
*
*  =================
*. Energy correction
*  =================
*
* E(n) =  C(0)T V C(N-1) 
        CALL MATVCB(V,C(1,IORD+1-1),VEC1,NDIM,NDIM,0)
        EN(IORD) = INPROD(C,VEC1,NDIM)
C?      WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ==========================
*. Wave function corrections
*  ==========================
*
* C(N) = (H(0)-E(0))-1 (-VC(N-1)
*                           +Sum(K=1,N)E(K)C(N-K))
        CALL MATVCB(V,C(1,IORD+1-1),VEC2,NDIM,NDIM,0)
        ONEM = -1.0D0
        CALL SCALVE(VEC2,ONEM,NDIM)
C?      write(6,*) ' first term to rhs '
C?      CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
*
        ONE = 1.0D0
        DO K = 1, IORD 
          CALL VECSUM(VEC2,VEC2,C(1,IORD+1-K),ONE,EN(K),NDIM)
        END DO
*. Check overlap with zero order state ( should be zero )
        OVLAP = INPROD(C(1,1),VEC2,NDIM)
        FACTOR = -OVLAP
C?      WRITE(6,*) ' OVLAP = ',OVLAP
        CALL VECSUM(VEC2,VEC2,C(1,1),ONE,FACTOR,NDIM)
*. Multiply with (H0(0)-E(0))-1
C            DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
        ZERO = 0.0D0
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
*
        CALL COPVEC(VEC1,C(1,IORD+1),NDIM)
*
C?      WRITE(6,*) ' Eigenfunction correction ', IORD
C?      CALL WRTMAT(C(1,IORD+1),1,NDIM,1,NDIM)
      END DO
* 
      WRITE(6,*) ' Energy corrections : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) '   Order             Correction '
      WRITE(6,*) ' ===================================='
      DO IORD = 1, NORD
        WRITE(6,'(1H ,3X,I3,E20.8)')IORD,EN(IORD)
      END DO
*
      ETOT = E0 + ECORE
      DO IORD = 1, NORD
        ETOT = ETOT + EN(IORD)
      END DO
      WRITE(6,*) ' Zero-order energy ', E0 + ECORE
      WRITE(6,*) ' Sum(K=0,NORD) E(K) ', ETOT 
*
      RETURN
      END
      SUBROUTINE SXSTRS(ISTROF,NIST,KSTROF,
     &                   NI,IOFF,NJ,JOFF,
     &                   IAMAPO,IAMAPS,LAMAP,IAMPFL,
     &                   ICMAPO,ICMAPS,LCMAP,ICMPFL,
     &                   ICMPO,ICMPL,
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
* ICMPFL : Creation map complete ?
* LCMAP  : Row dimension of Creation     map
*
* Jeppe Olsen, March 1994
      IMPLICIT REAL*8(A-H,O-Z)
*
*. Input
*
      INTEGER IAMAPO(LAMAP,*), IAMAPS(LAMAP,*)
C     INTEGER ICMAPO(LCMAP,*), ICMAPS(LCMAP,*)
      INTEGER ICMAPO(*),ICMAPS(*)
      INTEGER ICMPO(*),ICMPL(*)
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
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' SXFSTS : NTEST = ', NTEST
        WRITE(6,*) ' LCMAP LAMAP ', LCMAP,LAMAP
        WRITE(6,*) ' IAMPFL, ICMPFL ',IAMPFL, ICMPFL
        WRITE(6,*) ' MXSXST ', MXSXST 
      END IF
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
          IF(ICMPFL.EQ.1) THEN
            JISTRO = (JISTR-1)*LCMAP+1
            NIEFF = LCMAP
          ELSE
            JISTRO = ICMPO(JISTR)
            NIEFF = ICMPL(JISTR)
          END IF
C?        WRITE(6,*) ' JISTRO NIEFF ', JISTRO,NIEFF
C         DO 1001 IORB = IOFF,IOFF+NI-1
          DO 1001 IIORB = 1, NIEFF 
*
            IF(ICMPFL.EQ.1) THEN
              IORB = IIORB-1+IOFF
              IIIORB = IORB
            ELSE
              IORB = ABS(ICMAPO(JISTRO-1+IIORB))
              IIIORB = IIORB
              IF(IORB.LT.IOFF .OR. IORB .GT. IOFF + NI -1 ) 
     &        GOTO 1001
            END IF
*
C?          WRITE(6,*) ' IORB = ', IORB
            IF(IOFFDG.NE.0 .AND. IORB.EQ.JORB) GOTO 1001
            IJACT = 0
            IF(ICMAPO(JISTRO-1+IIIORB).EQ.+IORB) THEN
              IJACT  = 1
              IF(ICMAPS(JISTRO-1+IIIORB).GT.0) THEN
                IJISTR = ICMAPS(JISTRO-1+IIIORB)
                SIJ = SJ
              ELSE 
                IJISTR = -ICMAPS(JISTRO-1+IIIORB)
                SIJ = -SJ 
              END IF
            END IF
C?          WRITE(6,*) ' IIIORB IJISTR ', IIIORB, IJISTR
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
      SUBROUTINE H1STRS(H,ISTROF,NIST,KSTROF,
     &                   NI,IOFF,NJ,JOFF,
     &                   IAMAPO,IAMAPS,LAMAP,IAMPFL,
     &                   ICMAPO,ICMAPS,ICMPO,ICMPL,LCMAP,ICMPFL,
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
* ICMPFL : Creation     map complete ?
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
C     INTEGER ICMAPO(LCMAP,*), ICMAPS(LCMAP,*)
      INTEGER ICMAPO(*),ICMAPS(*)
      INTEGER ICMPL(*),ICMPO(*)
      DIMENSION H(NI,NJ)
      DIMENSION IPIJKK(NTOOB,NTOOB),XIJKK(NTOOB,*)
      DIMENSION IOCSTR(NEL,*)
*. Output
      INTEGER IEXSTR(MXSXST,*)
      INTEGER NEX(*)
      DIMENSION FACSTR(MXSXST,*)
*. To get rid of annoying and incorrect compiler warnings
      JISTR = 0
      SJ = 0.0D0
      IJISTR = 0
      SIJ = 0.0D0
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' H1STRS : NTEST = ', NTEST
        WRITE(6,*) ' LCMAP LAMAP ', LCMAP,LAMAP
        WRITE(6,*) ' IAMPFL, ICMPFL ',IAMPFL, ICMPFL
        WRITE(6,*) ' MXSXST ', MXSXST 
      END IF
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
          IF(ICMPFL.EQ.1) THEN
            JISTRO = (JISTR-1)*LCMAP
            NIEFF = NI
          ELSE
            JISTRO = ICMPO(JISTR)
            NIEFF =  ICMPL(JISTR)
          END IF
C         DO 1001 IORB = IOFF,IOFF+NI-1
          DO 1001 IIORB = 1, NIEFF
            IF(ICMPFL.EQ.1) THEN 
              IORB = IIORB-1+IOFF
              IIIORB = IORB
            ELSE
              IORB = ABS(ICMAPO(JISTRO-1+IIORB))
              IIIORB = IIORB
              IF(IORB.LT.IOFF.OR.IORB.GT.IOFF+NI-1)
     &        GOTO 1001
            END IF
C?          WRITE(6,*) ' IORB = ', IORB
            IJACT = 0
            IF(ICMAPO(JISTRO-1+IIIORB).EQ.+IORB) THEN
              IJACT  = 1
              IF(ICMAPS(JISTRO-1+IIIORB).GT.0) THEN
                IJISTR = ICMAPS(JISTRO-1+IIIORB)
                SIJ = SJ
              ELSE 
                IJISTR = -ICMAPS(JISTRO-1+IIIORB)
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
      SUBROUTINE DIABLK(IASM,IATP,IBSM,IBTP,IFULL,DIAG)
*
* Obtain diagonal block of determinant block 
* IASM,IATP,IBSM,IBTP
*
*
*. Type of operator is taken from OPER
*. Type of CI space is taken from CANDS
*
c      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'oper.inc'
*
      INCLUDE 'cands.inc'
C     COMMON/OPER/I12,IPERTOP
      INCLUDE 'cintfo.inc'
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL LUCIAQENTER('DIABL')
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRDIA)
*
** Specifications of internal space
*
      ISM = ICSM
      ISPC = ICSPC
*
      IAGP = IASTFI(ISPC)
      IBGP = IBSTFI(ISPC)
      NAEL = NELEC(IAGP)
      NBEL = NELEC(IBGP)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' DIABLK speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IAGP IBGP NAEL NBEL '
        WRITE(6,*)   IAGP,IBGP,NAEL,NBEL
        WRITE(6,*) ' I12,IPERTOP',I12,IPERTOP
      END IF
 
      NOCTPA = NOCTYP(IAGP)
      NOCTPB = NOCTYP(IBGP)
      MNRS1C = MNR1IC(ISPC)
      MXRS3C = MXR3IC(ISPC)
      IF(NTEST.GE.10) THEN
        write(6,*) ' NOCTPA NOCTPB MNRS1C,MXRS3C'
        write(6,*)   NOCTPA,NOCTPB,MNRS1C,MXRS3C
      END IF
*
**. Local memory
*
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUM,'DIABLK')
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NTOOB,   'ADDL  ',2,'KLH1D ')
      CALL MEMMAN(KLSMOS,NSMST,   'ADDL  ',2,'KLSMOS')
      CALL MEMMAN(KLBLTP,NSMST,   'ADDL  ',2,'KLSMOS')
      CALL MEMMAN(KLIOIO,NOCTPA*NOCTPB,   'ADDL  ',2,'KLIOIO')
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KLSVST,NSMST,   'ADDL  ',2,'KLSVST')
      ELSE
        KLSVST = 1
      END IF
      MAXA = IMNMX(int_mb(KNSTSO(IAGP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*
** Info on block structure of internal state
*
      IF(IDC.EQ.3.OR.IDC.EQ.4)
     &CALL SIGVST(int_mb(KLSVST),NSMST)
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,int_mb(KLBLTP),int_mb(KLSVST))
      STOP ' Update call to IAIBCM_GAS'
      CALL IAIBCM_GAS(MNRS1C,MXRS3C,NOCTPA,NOCTPB,
     &     WORK(KEL1(IAGP)),WORK(KEL3(IAGP)),
     &     WORK(KEL1(IBGP)),WORK(KEL3(IBGP)),WORK(KLIOIO),IPRDIA)
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1) 
      CALL GT1DIA(WORK(KLH1D))
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1) 
      ECORES = 0.0D0
      IF(I12.EQ.2)
     &CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*
      CALL DIABLKS(NAEL,int_mb(KOCSTR(IAGP)),NBEL,int_mb(KOCSTR(IBGP)),
     &             NACOB,NSMST,WORK(KLH1D),
     &             ISMOST(1,ISM),WORK(KLBLTP),
     &             WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &             WORK(KLK),int_mb(KNSTSO(IAGP)),int_mb(KNSTSO(IBGP)),
     &             WORK(KLIOIO),NOCTPA,NOCTPB,int_mb(KISTSO(IAGP)),
     &             int_mb(KISTSO(IBGP)),ECORES,
     &             PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &             WORK(KLRJKA),I12,
     &             IASM,IATP,IBSM,IBTP,IFULL,DIAG)
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIABLK')
      CALL LUCIAQEXIT('DIABL')
*
      RETURN
      END
      SUBROUTINE DIABLKS(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,NSMST,H,
     &                  ISMOST,IBLTP,XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,IOCOC,NOCTPA,NOCTPB,
     &                  ISSOA,ISSOB,ECORE,
     &                  PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,I12,
     &                  IASM,IATP,IBSM,IBTP,IFULL,DIAG)
*
*
* Obtain specific block of diagonal
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
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DIABLKS in action '
        WRITE(6,*) ' ================= '
        WRITE(6,*) 
        WRITE(6,*) ' IATP IBTP IASM IBSM '
        WRITE(6,*)   IATP,IBTP,IASM,IBSM
      END IF
*
      IF( NTEST .GE. 2000 ) THEN
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
      ONE = 1.0D0
      ONEM = -1.0D0
      IF(I12.EQ.2) 
     &CALL VECSUM(RK,RK,RJ,ONEM,ONE,NTOOB **2)
*
        IF(IFULL.EQ.0.AND.IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
          IF(IREST1.EQ.1) THEN
            MXBTP = IATP
          ELSE
            MXBTP = NOCTPB
          END IF
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
*
          IBSTRT = ISSOB(IBTP,IBSM)
          IBSTOP = IBSTRT + NSSOB(IBTP,IBSM)-1
*
          IDET = 0
          DO 899 IB = IBSTRT,IBSTOP
            IBREL = IB - IBSTRT + 1
*
*. Array for terms depending only upon IB
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
                DO  JEL = 1, NBEL
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
              IASTRT = ISSOA(IATP,IASM) - 1 + IBREL
            ELSE
              IASTRT = ISSOA(IATP,IASM)
            END IF
            IASTOP = ISSOA(IATP,IASM) + NSSOA(IATP,IASM) - 1
            DO 800 IA = IASTRT,IASTOP
              IDET = IDET + 1
              X = EB + RJKAA(IA-IOFF+1)
              DO 890 IEL = 1, NAEL
                X = X +XB(IASTR(IEL,IA)) 
  890         CONTINUE
              DIAG(IDET) = X
  800       CONTINUE
  899     CONTINUE
*
          IF(Ntest.ge.1000) THEN
             WRITE(6,*) ' Next batch of diagonal elements '
             CALL WRTMAT(DIAG,1,IDET,1,IDET)                
          END IF
*
      RETURN
      END
      SUBROUTINE HTV(VECIN,VECOUT,LUIN,LUOUT)
*
* Full operator times vector
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'oper.inc'
*
* complete operator in action
      IPERTOP = 0
      I12 = 2
*
      CALL MV7(VECIN,VECOUT,LUIN,LUOUT,0,0)
*
      RETURN
      END
      SUBROUTINE SIMPRT(LURF,LUN,LUVN,
     &           EN,SCR,MAXORD,VEC1,VEC2,LU1,LU2,
     &           LBLK,IH0DIA,LUH0,S,V,H0,ECORE,ECORE_H,ECORE_HEX,
     &           EREF,IE0AVEX,LUHI0,EFINAL)
*
* Solve the perturbation equations
*
* E(n) = <0!V!0(n-1)>
*
* !0(n)> = (H-E)-1( E(n)!0> -V!0(n-1)>
*                  + sum(l=1,n-1) E(l)!0(n-l)> )
*
* Simplified version for total symmetric perturbation
*
* Alternative expressions for the energy corrections
* are invoked using the 2n+1 rule 
*
* E(n+1) =  <0!v!n>  
*
*        = <k!v!n-k> - sum(j=0,k-1)sum(m=k-j,n-j)<j+1!n-m-j>E(m)
*
* A note on the perturbation :
*
* The perturbation is in general of the form 
*
* H0 = QH(apr)Q + E0P, P = |0><0>, Q = 1-P
*
* Where Q is some approcimation to the hamiltonian.
*
* In order to calculate ((H0-alpha)** -1 |x>, where <x|0>=0 , one must 
* distinguish between two cases 
*
* =====================================
* 1 |0> is an eigenfunction for H(apr) 
* =====================================
*
* In this case (H0-alpha)** -1 |x> = (H(apr)-alpha)** -1 |x>
*
* =========================================
* 2 |0> is  not an eigenfunction for H(apr)
* =========================================
*
*               (H0-alpha)** -1 |x> = (H(apr)-alpha)** -1 |x>
*             - (H(apr)-alpha)** -1 |0> <x|(H(apr)-alpha)** -1|0>
*                                       -------------------------
*                                       <X|(H(apr)-alpha)** -1|0>
*         
*
* Input
* =====
* LURF : file containing reference vector
* LUN : file number for file to contain perturbation vectors
* LUN : file number for file to contain perturbation vectors
* MAXORD : Order through which the equations should be solved
* VEC1,VEC2 : Scratch vectors ,omplete or blocks of vectorS
* LU1, LU2, : scratch files
* S V : vectors of size MAXORD*(MAXORD-1)/2
* Ecore : Core energy
* Eref  : Exact energy of reference state
* IE0AVEX : choice of zero order energy : 1 => E0 = <0|H0|0>
*                                         2 => E0 = EREF ( as supplied )
*
*
* Output
* ======
* LUN : contains the MAXORD correction vectors
* LUVN : contains the perturbation times the last correction vector
* EN : Contains the energy corrections through order 2*MAXORD+1
* EFINAL: The energy correct through order 2*MAXORD + 1
*
* Internal links
* ===============
*
* Solutions of linear eqs : HINTV                   
* Hamiltonian times vector: MV7
* H0 times vector         : H0TVM
*
*. Please do not go beyond perturbation level 100
* Jeppe Olsen ,  Summer of 94
*                Winter 96 : Nondiagonal H(apr), general H(apr),
*                option for diagonal on disc eliminated 
*                (only direct calculation allowed now)
*
*                Winter of 99 : (H0-1)** -1 |0> on LUHI0 added
      IMPLICIT REAL*8(A-H,O-Z)
      REAL * 8 INPRDD
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION EN(*)
      DIMENSION S(*),V(*),H0(*)
*. For communicating with H0TVM
      COMMON/CENOT/E0
      INCLUDE 'cshift.inc'
*. For communicating with MV7
      INCLUDE 'oper.inc'
*. A bit of  scratch
      DIMENSION SCR(*)
*
      NTEST = 5
*. Use direct diagonal routines
      IDIDIA = 1

      ONE = 1.0D0
      ONEM = -1.0D0
      ZERO = 0.0D0
      CALL SETVEC(V,ZERO,(MAXORD+1)*(MAXORD+1+1)/2)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to SIMPRT '
        WRITE(6,*) ' LBLK = ', LBLK
        WRITE(6,*) ' ECORE_H = ', ECORE_H
        WRITE(6,*) ' ECORE_HEX = ', ECORE_HEX
      END IF
*
      IF(NTEST.GE.600) THEN
        WRITE(6,*) ' initial reference '
        CALL WRTVCD(VEC1,LURF,1,LBLK)
      END IF
*
* ===============================================================
* 1 :                   Initialization 
* ===============================================================
*
*
        IF(IE0AVEX.GE.2) THEN
          E0RF = EREF
        ELSE
*  ===============
*. E0RF = <0!H(apr)!0>
*  ===============
*
          IF(IH0DIA.NE.0) THEN
*. Diagonal H0, simple
            CALL REWINO(LU1)
            CALL REWINO(LURF)
            CALL DIA0TRM_GAS(1,LURF,LU1,VEC1,VEC2,0.0D0)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' LU1 according to DIATRM '
              CALL WRTVCD(VEC1,LU1,1,LBLK)
            END IF
          ELSE
*. multiply with H(apr)            
            E0 = 0.0D0
            IPROJ = 0
            SHIFT = ECORE_H
            IPERTOP = 1
CJAN25      CALL MV7(VEC1,VEC2,LURF,LU1,0,0)
            CALL H0TVM(VEC1,VEC2,LURF,LU1)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' H(apr) times vector '
              CALL WRTVCD(VEC1,LU1,1,LBLK)
            END IF
          END IF
          E0RF = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK) 
        END IF
*
        WRITE(6,*) ' E0RF = ', E0RF 
        ENERGY = E0RF
*. Check of |0> is an eigenfunction for H(apr) ( not H0 ! ) 
* Calculate H(apr)|0> - <0!H(apr)!0> |0>
        IPERTOP = 1
        CALL MV7(VEC1,VEC2,LURF,LU1,0,0)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H(apr) times vector, part II '
          CALL WRTVCD(VEC1,LU1,1,LBLK)
        END IF
        HAPR00 = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK) 
        CALL VECSMD(VEC1,VEC2,ONE,-HAPR00,LU1,LURF,LUVN,1,LBLK)
        XNORM = INPRDD(VEC1,VEC2,LUVN,LUVN,1,LBLK)
*
        IF(ABS(XNORM/HAPR00) .LE. 1.0D-12) THEN
         IHAPREIG = 1
        ELSE
         IHAPREIG = 0
        END IF
*
        IF(IHAPREIG.EQ.0) THEN
*         HAPRM100 =  <0!(H(apr)-E0)**-1 |0>
          CALL DIA0TRM_GAS(2,LURF,LU1,VEC1,VEC2,-E0RF)
          HAPRM100  = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
          WRITE(6,*) ' HAPRM100', HAPRM100
*. Obtain (H0-E0) ** (-1) |0> ( diagonal approx )
          CALL DIA0TRM_GAS(2,LURF,LUHI0,VEC1,VEC2,-E0RF)
C         SHIFT = -(E0RF-ECORE_H)
C         SHIFT_DIA = -E0RF
C         E0 = E0RF
C         IAPR = 1
C         IPERTOP = 1
C         IPROJ = 0
C         WRITE(6,*) ' SHIFT before call to HINTV ', SHIFT
C         CALL COPVCD(LURF,LU2,VEC1,1,LBLK)
C         CALL HINTV(LU2,LUHI0,SHIFT,SHIFT_DIA,VEC1,VEC2,LBLK,0,0 )
        ENDIF


      WRITE(6,*) '  HAPR00,  XNORM, IHAPREIG, HAPRM100 : ',
     &              HAPR00,  XNORM, IHAPREIG, HAPRM100
      
*
*
*. V times initial vector  on  LUVN
* ==================================
*
*. H0+V !0(0)> on LU1 (ECORE_HEX missing)
        IAPR = 0
        IPERTOP = 0
        CALL REWINO(LURF)
        CALL REWINO(LU1)
        CALL HTV(VEC1,VEC2,LURF,LU1)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Output from HTV '
          CALL WRTVCD(VEC1,LU1,1,LBLK)
        END IF
*
* V|0> = (H - H0) !0(0)> = (H' + ECORE_HEX - E0) !0(0)> on LUVN
        E0RFM =  - (E0RF-ECORE_HEX)
        CALL VECSMD(VEC1,VEC2,ONE,E0RFM,LU1,LURF,LUVN,1,LBLK)
*
*.S(11),V(11),EN(1)
*=================
        S(1) = INPRDD(VEC1,VEC2,LURF,LURF,1,LBLK)
        V(1) = INPRDD(VEC1,VEC2,LUVN,LURF,1,LBLK)
        H0(1) = E0RF
        EN(1) = V(1)
        IF(NTEST.GE.1) WRITE(6,*)
     &  ' Energy correction , n and E(n) ', 1,EN(1)
*
*
* =======================================================================
*.2               Loop over orders of correction vectors 
* =======================================================================
*
      DO 1000 IORD = 1, MAXORD
*
* On entrance :  correction vectors 1 - IORD-1 on LUN
*                V!0(IORD-1)>                  on LUVN
*
*
* !0(n)>
* ======
*
*  E(n)!0(0)> -V!0(n-1)> on LU1
        CALL VECSMD(VEC1,VEC2,EN(IORD),ONEM,LURF,LUVN,LU1,1,LBLK)
*.  sum( l = 1,  IORD -1) (E( n- l) !0(l)> on LU2
        IF(IORD.GT.1) THEN
          DO II = 1, IORD -1
            SCR(II) = EN(IORD-II)
          END DO 
          CALL MVCSMD(LUN,SCR,LU2,LUVN,VEC1,VEC2,IORD-1,1,LBLK)
*. add on LUVN
          CALL VECSMD(VEC1,VEC2,ONE,ONE,LU2,LU1,LUVN,1,LBLK)
        ELSE
          CALL COPVCD(LU1,LUVN,VEC1,1,LBLK)
          CALL REWINO(LUN)
        END IF
*. project !0> component out, SAVE on LU2
        OVLAP = INPRDD(VEC1,VEC2,LURF,LUVN,1,LBLK)
        IF(NTEST.GE.2)  write(6,*) ' ovlap1  ', OVLAP
        CALL  VECSMD(VEC1,VEC2,ONE,-OVLAP,LUVN,LURF,LU2,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' RHS of lin.eq ,order ', IORD
          CALL WRTVCD(VEC1,LU2,1,LBLK)
        END IF
*. Multiply with (H0-E0)-1, save result on LUN
        IF(IH0DIA.NE.0) THEN
*. Multiply with inverted Diagonal 
          CALL REWINO(LU2)
          CALL REWINO(LU1)
          CALL DIA0TRM_GAS(2,LU2,LU1,VEC1,VEC2,-E0RF)
*. Increased printlevel for Frank Jensen
          ITOFRANK = 0
          IF(NTEST.GE.1000.OR.ITOFRANK.EQ.1) THEN
            WRITE(6,*) ' new correction vector of order ', IORD      
            CALL WRTVCD(VEC1,LU1,1,LBLK)
          END IF
          IF(IHAPREIG.EQ.0) THEN
*. Orthogonalize with (H(apr)-E0)**-1|0>
            OVLAP = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
            IF(NTEST.GE.2)  write(6,*) ' ovlap2  ', OVLAP
*. Set (H(apr)-E0) ** 1 |0>
            CALL REWINO(LURF)
            CALL REWINO(LUVN)
            CALL DIA0TRM_GAS(2,LURF,LUVN,VEC1,VEC2,-E0RF)
            FACTOR = (-OVLAP)/HAPRM100
            CALL  VECSMD(VEC1,VEC2,ONE,FACTOR,LU1,LUVN,LU2,1,LBLK)
*. Save on LU1
            CALL COPVCD(LU2,LU1,VEC1,1,LBLK)
          END IF
        ELSE 
*. Solve set of linear equations 
          SHIFT = -(E0RF-ECORE_H)
          SHIFT_DIA = -E0RF
          E0 = E0RF
          IAPR = 1
          IPERTOP = 1
          IPROJ = 1
          WRITE(6,*) ' SHIFT before call to HINTV ', SHIFT
          CALL HINTV(LU2,LU1,SHIFT,SHIFT_DIA,VEC1,VEC2,LBLK,LURF,LUHI0)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' new correction vector '            
            CALL WRTVCD(VEC1,LU1,1,LBLK)
          END IF
        END If
        X0N = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
        WRITE(6,*) ' Overlap <0!N> ', X0N
*Save on LUN 
        CALL REWINO(LU1)
        CALL COPVCD(LU1,LUN,VEC1,0,LBLK)
*. V!0(n)> on LUVN = (H - H0 )!0(n)> = (H'(holeform) - (H0 -ECORE_HEX))!0(n)>
* ================
*
*. H0+V !0(n)> on LU2 ( except ECORE_HEX ) |0(n)>
        CALL REWINO(LU1)
        CALL REWINO(LU2)
        CALL HTV(VEC1,VEC2,LU1,LU2)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H times correction vector '
          CALL WRTVCD(VEC1,LU2,1,LBLK)
        END IF
*. Test : E(N+1) = <0! V !n> = <0! H !n>
        ENN = INPRDD(VEC1,VEC2,LURF,LU2,1,LBLK)
        IF(NTEST.GE.1) WRITE(6,*) ' TEST : ENN = ', ENN
*
*
* H0 |0(n)> = Q H apr |0(n)> on LUVN ( and include -ECORE_HEX missi
*
        E0 = E0RF
        SHIFT = ECORE_H-ECORE_HEX
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Input to H0TVM '
          CALL WRTVCD(VEC1,LU1,1,LBLK)
        END IF
        IPERTOP = 1
        IPROJ = 1
        CALL H0TVM(VEC1,VEC2,LU1,LUVN)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H0 times correction vector '
          CALL WRTVCD(VEC1,LUVN,1,LBLK)
        END IF
*. Project |0> out
        OVLAP = INPRDD(VEC1,VEC2,LURF,LUVN,1,LBLK)
        IF(NTEST.GE.2)  write(6,*) ' ovlap3  ', OVLAP
        OVLAPM = -OVLAP
        CALL  VECSMD(VEC1,VEC2,ONE,OVLAPM,LUVN,LURF,LU1,1,LBLK)
        CALL COPVCD(LU1,LUVN,VEC1,1,LBLK)
        OVLAP = INPRDD(VEC1,VEC2,LU1,LURF,1,LBLK)
*
* (H - H0) !0(n)> on LUVN
C       CALL VECSMD(VEC1,VEC2,ONE,ONEM,LU2,LUVN,LU1,1,LBLK)
        CALL VECSMD(VEC1,VEC2,ONE,ONEM,LU2,LU1,LUVN,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' V ! 0(n)> '
          CALL WRTVCD(VEC1,LUVN,1,LBLK)
        END IF
*. E(n+1) = <0(n)!V!0>
* ====================
        EN(IORD+1) = INPRDD(VEC1,VEC2,LUVN,LURF,1,LBLK)
        ENERGY = ENERGY + EN(IORD)
        IF(NTEST.GE.1) WRITE(6,*)
     &  ' Energy correction , n and E(n) ', IORD+1,EN(IORD+1)
*
*. Augment matrices H0, S and V
* ============================
*
*
*  H0(ij) = <0(i-1)!H0!0(j-1)>
*
      CALL REWINO(LUN)
      DO JORD = 0, IORD
        IJ = (IORD+1)*(IORD+1-1)/2 + JORD+1
        IF(JORD.NE.0) THEN
          CALL REWINO(LU1)
          H0(IJ) = INPRDD(VEC1,VEC2,LUN ,LU1,0,LBLK)
        ELSE
          H0(IJ) = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
        END IF
      END DO
*
*  s(ij) = <0(i-1)!0(j-1)>
*
*. Place correction vector !0(n)> on LU1
      CALL SKPVCD(LUN,IORD-1,VEC1,1,LBLK)
      CALL REWINO(LU1)
      CALL COPVCD(LUN,LU1,VEC1,0,LBLK)
*
      CALL REWINO(LUN)
      DO JORD = 0, IORD
        IJ = (IORD+1)*(IORD+1-1)/2 + JORD+1
        IF(JORD.NE.0) THEN
          CALL REWINO(LU1)
          S(IJ) = INPRDD(VEC1,VEC2,LU1,LUN,0,LBLK)
        ELSE
          S(IJ) = 0.0D0
        END IF
      END DO
*
*  v(ij) = <0(i-1)!v!0(j-1)>
*
      CALL REWINO(LUN)
      DO JORD = 0, IORD
        IJ = (IORD+1)*(IORD+1-1)/2 + JORD+1
        IF(JORD.NE.0) THEN
          CALL REWINO(LUVN)
          V(IJ) = INPRDD(VEC1,VEC2,LUN,LUVN,0,LBLK)
        ELSE
          V(IJ) = INPRDD(VEC1,VEC2,LURF,LUVN,1,LBLK)
        END IF
      END DO
C!    IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Updated S matrix '
        CALL PRSYM(S,IORD+1)
        WRITE(6,*) ' Updated V matrix '
        CALL PRSYM(V,IORD+1)
        WRITE(6,*) ' Updated H0 matrix '
        CALL PRSYM(H0,IORD+1)
C!    END IF
*
*. Obtain additional energy expressions by 2n+1 rule.
* ===================================================
* 
        DO N = IORD+1,2*IORD+1
*. E(N) = <K!V!N-K-1> - Sum(j=0,k-1)sum(m=k-j,n-1-j)<j+1!n-1-m-j>E(m)
*. Use K = IORD
           X = V((IORD+1)*(IORD+1-1)/2+N-IORD-1+1)
           DO J = 0,IORD-1
             DO M = IORD-J,N-1-J
               II = MAX(J+1,N-1-M-J)
               JJ = MIN(J+1,N-1-M-J)
               IJ = (II+1)*(II+1-1)/2 + JJ + 1
               X = X - S(IJ)*EN(M)
             END DO
           END DO
           EN(N) = X
        END DO
*
        IF(NTEST.GE.2 ) THEN
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*) ' =========================================== '
          WRITE(6,*) ' Energy corrections obtained in iteration ', IORD
          WRITE(6,*) ' =========================================== '
          WRITE(6,*)
          WRITE(6,*)
     &    '   Order       Energy correction      Total Energy '
          WRITE(6,*)
     &   ' ========================================================='
          ENERGY = E0RF+ECORE-ECORE_HEX
          DO JORD = 1, 2*IORD+1
            ENERGY = ENERGY + EN(JORD) 
            WRITE(6,'(4X,I2,8X,1E18.10,6X,1E18.10)')
     &      JORD,EN(JORD),ENERGY
          END DO
        END IF
*. I can't wait to see the output, so FLUSH 
        LUOUT = 6
C       CALL  FLUSH(LUOUT)
*
 1000 CONTINUE
*
      WRITE(6,*)
      WRITE(6,*) ' Zero order energy : ', E0RF+ECORE-ECORE_HEX
      WRITE(6,*)
      WRITE(6,*) ' =========================================== '
      WRITE(6,*) ' Energy corrections obtained as <0!V!0(n-1)> '
      WRITE(6,*) ' =========================================== '
      WRITE(6,*)
      WRITE(6,*)
     &'   Order         Energy correction          Total Energy '
      WRITE(6,*)
     &' ========================================================='
      ENERGY = E0RF+ECORE-ECORE_HEX
      DO IORD = 1, MAXORD
        ENERGY = ENERGY + EN(IORD) 
        WRITE(6,'(4X,I2,8X,1E20.12,6X,1E22.14)')
     &  IORD,EN(IORD),ENERGY
      END DO
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) ' =========================================== '
      WRITE(6,*) ' Energy corrections obtained from 2n+1 rule  '
      WRITE(6,*) ' =========================================== '
      WRITE(6,*)
*
      WRITE(6,*)
     &  '   Order         Energy correction        Total Energy '
      WRITE(6,*)
     & ' ========================================================='
      DO JORD =MAXORD+1, 2*MAXORD+1
        ENERGY = ENERGY + EN(JORD)
        WRITE(6,'(4X,I2,8X,1E20.12,6X,1E22.14)')
     &  JORD,EN(JORD),ENERGY
      END DO
      EFINAL = ENERGY
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Final S matrix '
        CALL PRSYM(S,MAXORD+1)
        WRITE(6,*) ' Final V matrix '
        CALL PRSYM(V,MAXORD+1)
        WRITE(6,*) ' Final H0 matrix '
        CALL PRSYM(H0,MAXORD+1)
      END IF

*
      RETURN
      END
      SUBROUTINE FIFAM(FIFA)
*
*. Construct inactive + active fock matrix
*
*. On input FIFAM Should be the inactive Fock matrix, in symmetry packed form
*
* Jeppe Olsen
c      ImplICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*
CNW   DIMENSION FIFA(*)
      integer FIFA
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      CALL FIFAMS(FIFA,dbl_mb(KRHO1),IBSO,NSMOB,
     &            NTOOBS,NACOB,NTOOB,IREOST)
*
      RETURN
      END
      SUBROUTINE FIFAMS(FIFA,RHO1,IOBSM,NSMOb,lOBSM,NACOB,
     &            NORBt,ISTOb)
*
* Update inactive fock matrix with active contributions
*
*     FIFA(I,J) = FIFA(I,J) + sum(k,l) ((ij!kl)-0.5*(il!kj))*rho1(kl)
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

CNW   DIMENSION FIFA(*),RHO1(NACOB,NACOB)
      DIMENSION RHO1(NACOB,NACOB)
      integer FIFA
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
       if (ga_nodeid().eq.0) then
       WRITE(6,*) 
       WRITE(6,*) ' ======================='
       WRITE(6,*) ' Initial matrix to FIFA '
       WRITE(6,*) ' ======================='
       WRITE(6,*) 
       endif
       ISYM = 1
CNW    CALL APRBLM2(FIFA,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FIFA)
      END IF
*
*.  Assume spatial symmetric fock matrix
      ibuf = imnmx(lobsm,nsmob,2)
      ibuf = ibuf*(ibuf+1)/2
      if (.not.ma_push_get(mt_dbl,ibuf,'buffer',lbuf,kbuf)) call
     $   errquit('fifams:push_get failed', ibuf, MA_ERR)
      IJSM = 1
      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
        IF(JSM.NE.0) THEN
          istart = ij + 1     
          ii = 0
          call setvec(dbl_mb(kbuf),0.0d0,ibuf)
          DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
            DO J = IOBSM(JSM),I                          
              IP = ISTOB(I)
              JP = ISTOB(J)
               IJ= IJ + 1
               DO IA = 1, NACOB
                 DO IB = 1, NACOB
CNW  &             FIFA(IJ) = FIFA(IJ) 
CNW  &           + RHO1(IA,IB)
CNW  &           *(GTIJKL(IP,JP,IA,IB)-0.5*GTIJKL(IP,IB,IA,JP))
                   IF(RHO1(IA,IB).NE.0.0D0)        
     &             dbl_mb(kbuf+ii) = dbl_mb(kbuf+ii) + RHO1(IA,IB)
     &           *(GTIJKL(IP,JP,IA,IB)-0.5*GTIJKL(IP,IB,IA,JP))
                 END DO
               END DO
               ii = ii + 1
            END DO
          END DO
          if (ij.ge.istart)
     &    call ga_acc(FIFA,istart,ij,1,1,dbl_mb(kbuf),1,1.0d0)
        END IF
      END DO
      if (.not.ma_pop_stack(lbuf)) call
     $   errquit('fifams:pop_stack failed', lbuf, MA_ERR)
*
      IF(NTEST.NE.0) THEN
       if (ga_nodeid().eq.0) then
       WRITE(6,*) ' FI + FA in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       endif
       ISYM = 1
CNW    CALL APRBLM2(FIFA,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FIFA)
      END IF
* 
      RETURN
      END
      SUBROUTINE APRBLM_F7(A,LROW,LCOL,NBLK,ISYM)
C
C PRINT BLOCKED MATRIX
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*)
      DIMENSION LROW(NBLK),LCOL(NBLK)
C
      IBASE = 1
      WRITE(6,*)
      DO 100 IBLK = 1, NBLK
        WRITE(6,'(A,I3)') ' Block ... ',IBLK
        IF(ISYM.EQ.0) THEN
          IF(IBLK .NE. 1 ) IBASE = IBASE + LROW(IBLK-1)*LCOL(IBLK-1)
          CALL WRTMAT_F7(A(IBASE),LROW(IBLK),LCOL(IBLK),
     &         LROW(IBLK),LCOL(IBLK) )
        ELSE
          IF(IBLK .NE. 1 ) 
     &         IBASE = IBASE + LROW(IBLK-1)*(LCOL(IBLK-1)+1)/2
          CALL PRSYM(A(IBASE),LROW(IBLK))              
        END IF

  100 CONTINUE
      RETURN
      END
      SUBROUTINE APRBLM2(A,LROW,LCOL,NBLK,ISYM)
C
C PRINT BLOCKED MATRIX
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*)
      DIMENSION LROW(NBLK),LCOL(NBLK)
C
      IBASE = 1
      WRITE(6,*) ' Blocked matrix '
      WRITE(6,*) '================'
      WRITE(6,*)
      DO 100 IBLK = 1, NBLK
        WRITE(6,'(A,I3)') ' Block ... ',IBLK
        IF(ISYM.EQ.0) THEN
          IF(IBLK .NE. 1 ) IBASE = IBASE + LROW(IBLK-1)*LCOL(IBLK-1)
          CALL WRTMAT2(A(IBASE),LROW(IBLK),LCOL(IBLK),
     &         LROW(IBLK),LCOL(IBLK) )
        ELSE
          IF(IBLK .NE. 1 )
     &         IBASE = IBASE + LROW(IBLK-1)*(LCOL(IBLK-1)+1)/2
          CALL PRSYM(A(IBASE),LROW(IBLK))
        END IF

  100 CONTINUE
      RETURN
      END
      SUBROUTINE FI(FIMAT,ECC,IDOH2)
*
*. Construct inactive fockmatrix + core-core interaction energy.
*. I.e. add contributions from all orbitals
*  that belong to hole orbital spaces ( as defined by IPHGAS).
*
* Note that this is a more general definition of the
* Inactive Fockmatrix than usually used.
*
*. On input FIMAT should be the inactive Fock matrix, in symmetry packed form
*
* If I_USE_SIMTRH = 0 input and output matrices are assumed lower half packed
*                 = 1 Input and output matrices are assumed complete blocks
*
* Jeppe Olsen
*
* Revision : Dec 97 : General hole spaces
*            aug 00 : I_USE_SIMTRH switch added
c      ImplICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*
CNW   DIMENSION FIMAT(*)
      integer FIMAT
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'oper.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*
      NTEST = 00
      IF (NTEST.GE.10) THEN
        WRITE(6,*) ' FI: IPHGAS and ITPFSO '
        CALL IWRTMA(IPHGAS,1,NGAS,1,NGAS)
        CALL IWRTMA(ITPFSO,1,NTOOB,1,NTOOB)
      END IF

      IF(I_USE_SIMTRH.EQ.0) THEN
        CALL FIH(FIMAT,ECC,IBSO,NSMOB,ITPFSO,IPHGAS,NTOOBS,NTOOB,
     &           IREOST,IDOH2,NGAS)
      ELSE
        CALL FIHA(FIMAT,ECC,IBSO,NSMOB,ITPFSO,IPHGAS,NTOOBS,NTOOB,
     &           IREOST,IDOH2)
      END IF
*
      RETURN
      END
      SUBROUTINE FIHA(FI,ECC,IOBSM,NSMOB,ITPFSO,IPHGAS,LOBSM,NORBT,
     &               ISTOB,IDOH2)
*
* construct inactive fock matrix 
*
*     FI(I,J) = FI(I,J) + sum(h) (2(ij!hh)-(ih!jh))
*
* where h is summed over all hole orbitals (as declaed by IPHGAS)
* Note that this is a more general definition of the Inactive
* Fock matrix than usually used.
* (Normal realization : see FIS )
*
* Version with complete symmetry blocks
*
* Jeppe Olsen ( I admit ) 
*
* Aug 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION FI(*)
      integer FI
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
      INTEGER ITPFSO(*), IPHGAS(*)
*. To get rid of annoying and incorrect compiler warnings
      IIOFF = 0
*
      NTEST = 0
*
* Core-Core energy 
*
      ECC = 0.0D0
      IJSM = 1
*. One-electron part 
      call ga_distribution(FI,ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(FI, ilo, ihi, 1, 1, ikfi, idum)
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IIOFF = IIOFF + LOBSM(ISM-1)** 2 
        END IF
        LEN = LOBSM(ISM)
        DO I = 1, LEN                           
          II = IIOFF -1 + (I-1)*LEN + I 
CNW       IF(IPHGAS(ITPFSO(I+IOBSM(ISM)-1)).EQ.2) 
CNW  &    ECC = ECC + 2*FI(II)
          IF(IPHGAS(ITPFSO(I+IOBSM(ISM)-1)).EQ.2
     &       .and.II.ge.ilo.and.II.le.ihi) THEN
            ECC = ECC + 2*dbl_mb(ikfi+ii-ilo)
          ENDIF
        END DO
      END DO
      call ga_release(FI,ilo,ihi,1,1)
      endif
      call ga_sync()
      call ga_dgop(1, ECC, 1, '+')
      if (ga_nodeid().eq.0) WRITE(6,*) ' one-electron part to ECC ', ECC
*. Two-electron part
      IF(IDOH2.NE.0) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IOBSM(ISM), IOBSM(ISM) + LOBSM(ISM)-1
          DO J = IOBSM(JSM), IOBSM(JSM) + LOBSM(JSM)-1
              IP = ISTOB(I)
              JP = ISTOB(J)
              IF(IPHGAS(ITPFSO(I)).EQ.2.AND.IPHGAS(ITPFSO(J)).EQ.2)
     &        ECC = ECC +2*GTIJKL(IP,IP,JP,JP)-GTIJKL(IP,JP,JP,IP)
*
C?            IF(IPHGAS(ITPFSO(I)).EQ.2.AND.IPHGAS(ITPFSO(J)).EQ.2)
C?   &        THEN
C?              WRITE(6,*) ' IP, JP, (II!JJ), (IJ!JI) = ',
C?   &          IP, JP, GTIJKL(IP,IP,JP,JP),GTIJKL(IP,JP,JP,IP)
C?            END IF
*
          END DO
          END DO
        END DO
        END DO
      END IF
      call ga_sync()
c     call ga_dgop(1, ECC, 1, '+')
*
      IF(ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
      END IF
*
*.  Inactive Fock matrix
*
      IF(IDOH2.NE.0) THEN 
        ibuf = imnmx(lobsm,nsmob,2)
        ibuf = ibuf*ibuf
        if (.not.ma_push_get(mt_dbl,ibuf,'buffer',lbuf,kbuf)) call
     $     errquit('fiha:push_get failed', ibuf, MA_ERR)
        IJSM = 1
        IJ = 0
        DO ISM = 1, NSMOB
          CALL SYMCOM(2,6,ISM,JSM,IJSM)
          IF(JSM.NE.0) THEN
            istart = ij + 1
            ii = 0
            call setvec(dbl_mb(kbuf),0.0d0,ibuf)
            DO J = IOBSM(JSM),IOBSM(JSM) + LOBSM(JSM) - 1
              DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM) - 1
                IP = ISTOB(I)
                JP = ISTOB(J)
                IJ= IJ + 1
                DO KSYM = 1, NSMOB
                  DO K = IOBSM(KSYM),IOBSM(KSYM)-1+LOBSM(KSYM)
                    KP = ISTOB(K)
CNW                 IF(IPHGAS(ITPFSO(K)).EQ.2) FI(IJ) = FI(IJ) 
CNW  &            + 2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP)
                    IF(IPHGAS(ITPFSO(K)).EQ.2) THEN
                      dbl_mb(kbuf+ii) = dbl_mb(kbuf+ii) 
     &            + 2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP)
                    ENDIF
                  END DO
                END DO
                ii = ii + 1
              END DO
            END DO
            if (ij.ge.istart.and.ga_nodeid().eq.0)
     &      call ga_acc(FI,istart,IJ,1,1,dbl_mb(kbuf),1,1.0d0)
          END IF
        END DO
        if (.not.ma_pop_stack(lbuf)) call
     $     errquit('fiha:pop_stack failed', lbuf, MA_ERR)
      END IF
      call ga_sync()
*
      IF(NTEST.NE.0) THEN
*
       if (ga_nodeid().eq.0) then
       WRITE(6,*) ' FI in Symmetry blocked form '
       WRITE(6,*) ' ============================'
       WRITE(6,*) 
       endif
       ISYM = 0
CNW    CALL APRBLM2(FI,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FI)
      END IF
* 
      RETURN
      END
      SUBROUTINE FIH(FI,ECC,IOBSM,NSMOB,ITPFSO,IPHGAS,LOBSM,NORBT,ISTOB,
     &               IDOH2,NGAS)
*
* construct inactive fock matrix 
*
*     FI(I,J) = FI(I,J) + sum(h) (2(ij!hh)-(ih!jh))
*
* where h is summed over all hole orbitals (as declared by IPHGAS)
* Note that this is a more general definition of the Inactive
* Fock matrix than usually used.
* (Normal realization : see FIS )
*
* Jeppe Olsen ( I admit ) 
*
* Dec 97
* June 2010: Inactive orbitals flagged by type = 0 added.
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION FI(*)
      integer FI
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
      INTEGER ITPFSO(*), IPHGAS(*)
*. To get rid of annoying and incorrect compiler warnings
      IIOFF = 0
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN
       if (ga_nodeid().eq.0) then
       WRITE(6,*) 
       WRITE(6,*) ' Input integrals to FIH in symmetry blocked form '
       WRITE(6,*) ' ================================================'
       WRITE(6,*) 
       endif
       ISYM = 1
CNW    CALL APRBLM2(FI,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FI)
      END IF
*
* Core-Core energy 
*
      ECC = 0.0D0
      IJSM = 1
*. One-electron part 
      call ga_distribution(FI,ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(FI, ilo, ihi, 1, 1, ikfi, idum)
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IIOFF = IIOFF + LOBSM(ISM-1)*(LOBSM(ISM-1)+1)/2
        END IF
        II = IIOFF-1
        DO I = IOBSM(ISM),IOBSM(ISM)+LOBSM(ISM)-1
          II = II + (I-IOBSM(ISM)+1) 
          I_INACTIVE = 0
          IF(ITPFSO(I).GT.0.AND.ITPFSO(I).LE.NGAS) THEN
           IF (IPHGAS(ITPFSO(I)).EQ.2) I_INACTIVE = 1
          END IF
          IF(ITPFSO(I).EQ.0) I_INACTIVE = 1
C?        WRITE(6,*) ' I, I_INACTIVE = ', I,I_INACTIVE
          IF(I_INACTIVE.EQ.1) THEN
C?          WRITE(6,*) ' Contribution to ECC from I =', I
CNW         ECC = ECC + 2*FI(II)
            if (II.ge.ilo.and.II.le.ihi)
     &      ECC = ECC + 2*dbl_mb(ikfi+ii-ilo)
            
C?          WRITE(6,*) ' Updated ECC = ', ECC
          END IF
        END DO
      END DO
      call ga_release(FI,ilo,ihi,1,1)
      endif
      call ga_sync()
      call ga_dgop(1, ECC, 1, '+')
      IF(NTEST.GE.0000.and.ga_nodeid().eq.0) 
     & WRITE(6,*) ' one-electron term to Ecore ', ECC
*. Two-electron part
      IF(IDOH2.NE.0) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IOBSM(ISM), IOBSM(ISM) + LOBSM(ISM)-1
          DO J = IOBSM(JSM), IOBSM(JSM) + LOBSM(JSM)-1
*
            I_INACTIVE = 0
            IF(ITPFSO(I).GT.0.AND.ITPFSO(I).LE.NGAS) THEN
             IF (IPHGAS(ITPFSO(I)).EQ.2) I_INACTIVE = 1
            END IF
            IF(ITPFSO(I).EQ.0) I_INACTIVE = 1
*
            J_INACTIVE = 0
            IF(ITPFSO(J).GT.0.AND.ITPFSO(J).LE.NGAS) THEN
             IF (IPHGAS(ITPFSO(J)).EQ.2) J_INACTIVE = 1
            END IF
            IF(ITPFSO(J).EQ.0) J_INACTIVE = 1
C?          WRITE(6,*) ' I, J, I_INACTIVE, J_INACTIVE = ', 
C?   &                   I, J, I_INACTIVE, J_INACTIVE
*
            IP = ISTOB(I)
            JP = ISTOB(J)
*
            IF(I_INACTIVE.EQ.1.AND.J_INACTIVE.EQ.1) THEN
C?            WRITE(6,*) ' Contribution to ECC from I,J =', I,J
              ECC = ECC +2*GTIJKL(IP,IP,JP,JP)-GTIJKL(IP,JP,JP,IP)
C?            WRITE(6,*) ' Updated ECC = ', ECC
            END IF
          END DO
          END DO
        END DO
        END DO
      END IF
c     call ga_dgop(1, ECC, 1, '+')
*
      IF(NTEST.NE.0.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' FIH Core-Core interaction energy ', ECC
      END IF
*
*.  Inactive Fock matrix
*
      IF(IDOH2.NE.0) THEN 
        ibuf = imnmx(lobsm,nsmob,2)
        ibuf = ibuf*(ibuf+1)/2
        if (.not.ma_push_get(mt_dbl,ibuf,'buffer',lbuf,kbuf)) call
     $     errquit('fiha:push_get failed', ibuf, MA_ERR)
        IJSM = 1
        IJ = 0
        DO ISM = 1, NSMOB
          CALL SYMCOM(2,6,ISM,JSM,IJSM)
          IF(JSM.NE.0) THEN
            istart = ij + 1
            ii = 0
            call setvec(dbl_mb(kbuf),0.0d0,ibuf)
            DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
              DO J = IOBSM(JSM),I                          
                IP = ISTOB(I)
                JP = ISTOB(J)
                IJ= IJ + 1
                DO KSYM = 1, NSMOB
                  DO K = IOBSM(KSYM),IOBSM(KSYM)-1+LOBSM(KSYM)
                    K_INACTIVE = 0
                    IF(ITPFSO(K).GT.0.AND.ITPFSO(K).LE.NGAS) THEN
                     IF (IPHGAS(ITPFSO(K)).EQ.2) K_INACTIVE = 1
                    END IF
                    IF(ITPFSO(K).EQ.0) K_INACTIVE = 1
                    KP = ISTOB(K)
CNW                 IF(K_INACTIVE.EQ.1) FI(IJ) = FI(IJ) 
CNW  &            + 2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP)
                    IF(K_INACTIVE.EQ.1) THEN
                      dbl_mb(kbuf+ii) = dbl_mb(kbuf+ii)
     &                  +2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP)
                    ENDIF
                  END DO
                END DO
                ii = ii + 1
              END DO
            END DO
            if (ij.ge.istart.and.ga_nodeid().eq.0)
     &      call ga_acc(FI,istart,IJ,1,1,dbl_mb(kbuf),1,1.0d0)
          END IF
        END DO
        if (.not.ma_pop_stack(lbuf)) call
     $     errquit('fiha:pop_stack failed', lbuf, MA_ERR)
      END IF
      call ga_sync()
*
      IF(NTEST.NE.0) THEN
*
       if (ga_nodeid().eq.0) then
       WRITE(6,*) ' FI in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       endif
       ISYM = 1
CNW    CALL APRBLM2(FI,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FI)
      END IF
* 
      RETURN
      END
      SUBROUTINE FIHO(FI,ECC,IOBSM,NSMOB,ITPFSO,IPHGAS,LOBSM,NORBT,
     &               ISTOB,IDOH2)
*
* construct inactive fock matrix 
*
*     FI(I,J) = FI(I,J) + sum(h) (2(ij!hh)-(ih!jh))
*
* where h is summed over all hole orbitals (as declaed by IPHGAS)
* Note that this is a more general definition of the Inactive
* Fock matrix than usually used.
* (Normal realization : see FIS )
*
* Jeppe Olsen ( I admit ) 
*
* Dec 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION FI(*)
      integer FI
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
      INTEGER ITPFSO(*), IPHGAS(*)
*. To get rid of annoying and incorrect compiler warnings
      IIOFF = 0
*
      NTEST = 000
C?    WRITE(6,*) ' FIH: IPHGAS and ITPFSO '
C?    CALL IWRTMA(IPHGAS,1,2,1,2)
C?    CALL IWRTMA(ITPFSO,1,3,1,3)
*
* Core-Core energy 
*
      ECC = 0.0D0
      IJSM = 1
*. One-electron part 
      call ga_distribution(FI,ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(FI, ilo, ihi, 1, 1, ikfi, idum)
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IIOFF = IIOFF + LOBSM(ISM-1)*(LOBSM(ISM-1)+1)/2
        END IF
        II = IIOFF-1
        DO I = IOBSM(ISM),IOBSM(ISM)+LOBSM(ISM)-1
          II = II + (I-IOBSM(ISM)+1) 
CNW       IF(IPHGAS(ITPFSO(I)).EQ.2) ECC = ECC + 2*FI(II)
          IF(IPHGAS(ITPFSO(I)).EQ.2 
     &      .and.II.ge.ilo.and.II.le.ihi) THEN
            ECC = ECC + 2*dbl_mb(ikfi+ii-ilo)
          ENDIF
        END DO
      END DO
      call ga_release(FI,ilo,ihi,1,1)
      endif
      call ga_dgop(1, ECC, 1, '+')
C?    WRITE(6,*) ' one-electron part to ECC ', ECC
*. Two-electron part
      IF(IDOH2.NE.0) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IOBSM(ISM), IOBSM(ISM) + LOBSM(ISM)-1
          DO J = IOBSM(JSM), IOBSM(JSM) + LOBSM(JSM)-1
              IP = ISTOB(I)
              JP = ISTOB(J)
              IF(IPHGAS(ITPFSO(I)).EQ.2.AND.IPHGAS(ITPFSO(J)).EQ.2)
     &        ECC = ECC +2*GTIJKL(IP,IP,JP,JP)-GTIJKL(IP,JP,JP,IP)
          END DO
          END DO
        END DO
        END DO
      END IF
*
      IF(NTEST.NE.0.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
      END IF
*
*.  Inactive Fock matrix
*
      IF(IDOH2.NE.0) THEN 
        ibuf = imnmx(lobsm,nsmob,2)
        ibuf = ibuf*(ibuf+1)/2
        if (.not.ma_push_get(mt_dbl,ibuf,'buffer',lbuf,kbuf)) call
     $     errquit('fiha:push_get failed', ibuf, MA_ERR)
        IJSM = 1
        IJ = 0
        DO ISM = 1, NSMOB
          CALL SYMCOM(2,6,ISM,JSM,IJSM)
          IF(JSM.NE.0) THEN
            istart = ij + 1
            ii = 0
            call setvec(dbl_mb(kbuf),0.0d0,ibuf)
            DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
              DO J = IOBSM(JSM),I                          
                IP = ISTOB(I)
                JP = ISTOB(J)
                IJ= IJ + 1
                DO KSYM = 1, NSMOB
                  DO K = IOBSM(KSYM),IOBSM(KSYM)-1+LOBSM(KSYM)
                    KP = ISTOB(K)
CNW                 IF(IPHGAS(ITPFSO(K)).EQ.2) FI(IJ) = FI(IJ) 
CNW  &            + 2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP)
                    IF(IPHGAS(ITPFSO(K)).EQ.2) THEN
                      dbl_mb(kbuf+ii) = dbl_mb(kbuf+ii)
     &            + 2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP)
                    ENDIF
                  END DO
                END DO
                ii = ii + 1
              END DO
            END DO
            if (ij.ge.istart.and.ga_nodeid().eq.0)
     &      call ga_acc(FI,istart,IJ,1,1,dbl_mb(kbuf),1,1.0d0)
          END IF
        END DO
        if (.not.ma_pop_stack(lbuf)) call
     $     errquit('fiha:pop_stack failed', lbuf, MA_ERR)
      END IF
      call ga_sync()
*
      IF(NTEST.NE.0) THEN
*
       if (ga_nodeid().eq.0) then
       WRITE(6,*) ' FI in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       endif
       ISYM = 1
CNW    CALL APRBLM2(FI,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FI)
      END IF
* 
      RETURN
      END
      SUBROUTINE FISM_OLD(FI,ECC)
*
*. Outer routine for calculating inactive Fock-matrix and
*  core-energy. No particle-hole corrections
*
* Using Initial integrals
*
*
*. Jeppe Olsen, August 2010
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      CALL FIS_OLD(FI,ECC,IBSO,NSMOB,NINOBS,NTOOBS,IREOST)
*
      RETURN
      END
      SUBROUTINE FIS_OLD(FI,ECC,IOBSM,NSMOB,NINOBS,LOBSM,ISTOB)
*
* construct inactive fock matrix 
*
*     FI(I,J) = FI(I,J) + sum(k) (2(ij!kk)-(ik!jk))
*
* Jeppe Olsen ( I admit ) 
*
* Last revision, Nov. 2012;  Jeppe Olsen; Back to Ini-integrals
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION FI(*)
      integer FI
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*),NINOBS(*)
*. To get rid of annoying and incorrect compiler warnings
      IIOFF = 0
*
      NTEST = 000
*
* Core-Core energy 
*
      ECC = 0.0D0
      IJSM = 1
*. One-electron part 
      call ga_distribution(FI,ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(FI, ilo, ihi, 1, 1, ikfi, idum)
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IIOFF = IIOFF + LOBSM(ISM-1)*(LOBSM(ISM-1)+1)/2
        END IF
        IIL = max(ilo,iioff) - ilo 
        IIH = min((iioff+ninobs(ism)),ihi) - ilo
CNW     II = IIOFF-1
CNW     DO I = 1, NINOBS(ISM)
CNW       II = II + I
CNW       ECC = ECC + 2*FI(II)
        do i = iil, iih
          ECC = ECC + 2.0d0*dbl_mb(ikfi+i)
        END DO
      END DO
      call ga_release(FI,ilo,ihi,1,1)
      endif
      EONE = ECC
      call ga_dgop(1, EONE, 1, '+')
*. Two-electron part
C gtijkl zero if not in my range
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
          DO I = IOBSM(ISM),IOBSM(ISM) + NINOBS(ISM)-1
          DO J = IOBSM(JSM),IOBSM(JSM) + NINOBS(JSM)-1
              IP = ISTOB(I)
              JP = ISTOB(J)
              ECC = ECC +2*GTIJKL(IP,IP,JP,JP)-GTIJKL(IP,JP,JP,IP)
          END DO
          END DO
      END DO
      END DO
c     call ga_dgop(1, ECC, 1, '+')
*
      IF(NTEST.NE.0.and. ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
        WRITE(6,*) ' One- and two-electron energies ',
     &  EONE, ECC- EONE
      END IF
*
*.  Inactive Fock matrix
*
      ibuf = imnmx(lobsm,nsmob,2)
      ibuf = ibuf*(ibuf+1)/2
      if (.not.ma_push_get(mt_dbl,ibuf,'buffer',lbuf,kbuf)) call
     $   errquit('fis_old:push_get failed', ibuf, MA_ERR)
      IJSM = 1
      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
        IF(JSM.NE.0) THEN
          istart = ij + 1
          ii = 0
          call setvec(dbl_mb(kbuf),0.0d0,ibuf)
          DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
            DO J = IOBSM(JSM),I                          
              IP = ISTOB(I)
              JP = ISTOB(J)
              IJ= IJ + 1
              DO KSYM = 1, NSMOB
                DO K = IOBSM(KSYM),IOBSM(KSYM)-1+NINOBS(KSYM)
                  KP = ISTOB(K)
CNW               FI(IJ) = FI(IJ) 
CNW  &          +(2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP))
                dbl_mb(kbuf+ii) = dbl_mb(lbuf+ii)
     &          +(2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP))
                END DO
              END DO
              ii = ii + 1
            END DO
          END DO
          if (ij.ge.istart.and.ga_nodeid().eq.0)
     &    call ga_acc(FI,istart,IJ,1,1,dbl_mb(lbuf),1,1.0d0)
        END IF
      END DO
      call ga_sync()
      if (.not.ma_pop_stack(lbuf)) call
     $   errquit('fis_old:pop_stack failed', lbuf, MA_ERR)

*
      IF(NTEST.NE.0) THEN
*
       if (ga_nodeid().eq.0) then
       WRITE(6,*) ' FI in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       endif
       ISYM = 1
CNW    CALL APRBLM2(FI,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FI)
      END IF
* 
      RETURN
      END
      SUBROUTINE FISM(FI,ECC)
*
*. Outer routine for calculating inactive Fock-matrix and
*  core-energy. No particle-hole corrections
*
*
*. Jeppe Olsen, August 2010
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      CALL FIS(FI,ECC,IBSO,NSMOB,NINOBS,NTOOBS,IREOST)
*
      RETURN
      END
      SUBROUTINE FIS(FI,ECC,IOBSM,NSMOB,NINOBS,LOBSM,ISTOB)
*
* construct inactive fock matrix 
*
*     FI(I,J) = FI(I,J) + sum(k) (2(ij!kk)-(ik!jk))
*
* Jeppe Olsen ( I admit ) 
*
* Last revision, Aug. 28 2012, Jeppe Olsen, updated to modern times
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION FI(*)
      integer FI
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*),NINOBS(*)
*. To get rid of annoying and incorrect compiler warnings
      IIOFF = 0
*
      NTEST = 0000
*
* Core-Core energy 
*
      ECC = 0.0D0
      IJSM = 1
*. One-electron part 
      call ga_distribution(FI,ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(FI, ilo, ihi, 1, 1, ikfi, idum)
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IIOFF = IIOFF + LOBSM(ISM-1)*(LOBSM(ISM-1)+1)/2
        END IF
        iil = max(ilo,iioff) - ilo
        iih = min((iioff+ninobs(ism)),ihi) - ilo
CNW     II = IIOFF-1
CNW     DO I = 1, NINOBS(ISM)
CNW       II = II + I
CNW       ECC = ECC + 2*FI(II)
        do i = iil, iih
          ECC = ECC + 2.0d0*dbl_mb(ikfi+i)
        END DO
      END DO
      call ga_release(FI,ilo,ihi,1,1)
      endif
      EONE = ECC
      call ga_dgop(1, EONE, 1, '+')
*. Two-electron part
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
          DO I = IOBSM(ISM),IOBSM(ISM) + NINOBS(ISM)-1
          DO J = IOBSM(JSM),IOBSM(JSM) + NINOBS(JSM)-1
              IP = ISTOB(I)
              JP = ISTOB(J)
              ECC = ECC +2*GTIJKL_GN(IP,IP,JP,JP)-GTIJKL_GN(IP,JP,JP,IP)
          END DO
          END DO
      END DO
      END DO
c     call ga_sync()
c     call ga_dgop(1, ECC, 1, '+')
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
        WRITE(6,*) ' One- and two-electron energies ',
     &  EONE, ECC- EONE
      END IF
*
*.  Inactive Fock matrix
*
      ibuf = imnmx(lobsm,nsmob,2)
      ibuf = ibuf*(ibuf+1)/2
      if (.not.ma_push_get(mt_dbl,ibuf,'buffer',lbuf,kbuf)) call
     $   errquit('fis:push_get failed', ibuf, MA_ERR)
      IJSM = 1
      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
        IF(JSM.NE.0) THEN
          istart = ij + 1
          ii = 0
          call setvec(dbl_mb(kbuf),0.0d0,ibuf)
          DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
            DO J = IOBSM(JSM),I                          
              IP = ISTOB(I)
              JP = ISTOB(J)
              IJ= IJ + 1
              DO KSYM = 1, NSMOB
                DO K = IOBSM(KSYM),IOBSM(KSYM)-1+NINOBS(KSYM)
                  KP = ISTOB(K)
CNW               FI(IJ) = FI(IJ) 
CNW  &          +(2.0D0*GTIJKL_GN(IP,JP,KP,KP)-GTIJKL_GN(IP,KP,KP,JP))
                  dbl_mb(kbuf+ii) = dbl_mb(lbuf+ii)
     &          +(2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP))
                END DO
              END DO
              ii = ii + 1
            END DO
          END DO
          if (ij.ge.istart.and.ga_nodeid().eq.0)
     &    call ga_acc(FI,istart,IJ,1,1,dbl_mb(lbuf),1,1.0d0)
        END IF
      END DO
      call ga_sync()
      if (.not.ma_pop_stack(lbuf)) call
     $   errquit('fis:pop_stack failed', lbuf, MA_ERR)

*
      IF(NTEST.NE.0) THEN
*
       WRITE(6,*) ' FI in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       ISYM = 1
CNW    CALL APRBLM2(FI,LOBSM,LOBSM,NSMOB,ISYM)
       call ga_print(FI)
      END IF
* 
      RETURN
      END
C     SUBROUTINE APRBLM2(A,LROW,LCOL,NBLK,ISYM)
C
C PRINT BLOCKED MATRIX
C
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     DIMENSION A(*)
C     DIMENSION LROW(NBLK),LCOL(NBLK)
C
C     IBASE = 1
C     WRITE(6,*) ' Blocked matrix '
C     WRITE(6,*) '================'
C     WRITE(6,*)
C     DO 100 IBLK = 1, NBLK
C       WRITE(6,'(A,I3)') ' Block ... ',IBLK
C       IF(ISYM.EQ.0) THEN
C       IF(IBLK .NE. 1 ) IBASE = IBASE + LROW(IBLK-1)*LCOL(IBLK-1)
C       CALL WRTMAT(A(IBASE),LROW(IBLK),LCOL(IBLK),
C    &              LROW(IBLK),LCOL(IBLK) )
C      ELSE
C       IF(IBLK .NE. 1 ) 
C    &  IBASE = IBASE + LROW(IBLK-1)*(LCOL(IBLK-1)+1)/2
C       CALL PRSYM(A(IBASE),LROW(IBLK))              
C      END IF

C 100 CONTINUE
C     RETURN
C     END
      SUBROUTINE MICDV5(MV7,VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  THRES_E,CONVER,RNRM_CNV)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Special version for NROOT = 1, MAXVEC = 2 !!
*
* Input :
* =======
*        MV7: Routine for direct CI
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
*   THRES_E   : Convergence criteria for eigenvalue
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
       EXTERNAL MV7
*
       WRITE(6,*) ' MICDV5, LU3 = ', LU3
       IF(NROOT.NE.1) THEN
         WRITE(6,*) 
     &   ' MICDV5 is wrong path, since NROOT .ne. 1, NROOT =',
     &     NROOT
         STOP 20
       END IF
       IF(MAXVEC.NE.2) THEN
         WRITE(6,*) 
     &   ' MICDV5 is wrong path, since MAXVEC.ne. 2, MAXVEC=',
     &     MAXVEC
         STOP 20
       END IF
       IF(NINVEC.NE.1) THEN
         WRITE(6,*) ' Several input vectors ', NINVEC
         WRITE(6,*) ' Only the first will be used '
       END IF
*
       IPICO = 0
       IOLSTM = 1
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(6,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(6,*) ' Normal Davidson method '
*
       IF(IPRT.GE.1) THEN
         WRITE(6,*) ' Convergence threshold for eigenvalues', THRES_E
       END IF
*
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
       IROOT = 1
*
* ===================
*.Initial iteration
* ===================
*
       ITER = 1
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       IF(IPRT.GE.600) THEN
         WRITE(6,*) ' Initial C -vector'
         CALL WRTVCD(VEC1,LU1,1,LBLK)
       END IF
*
       CALL MV7(VEC1,VEC2,LU1,LU2,0,0)
*
       IF(IPRT.GE.600) THEN
         WRITE(6,*) ' Initial sigma-vector'
         CALL WRTVCD(VEC1,LU2,1,LBLK)
       END IF
C?     write(6,*) ' Enforced stop after sigma'
C?     stop ' Enforced stop after sigma'
*. Projected matrix
       APROJ(1) = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
*
       IF( IPRT .GE.3  ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,1     )
       END IF
*. Diagonalize initial projected matrix : Dimension one :simple
       EIG(1,IROOT) = APROJ(1)                          
       AVEC(1) = 1.0D0
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(6,'(5F21.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       NVEC = 1      
       ITERP = 1
*. Add shift and print out
C?     ONE = 1.0D0 
C?     CALL VECSMD(VEC1,VEC2,EIGSHF,ONE,LU1,LU2,LU3,1,LBLK)
C?     WRITE(6,*) ' Sigma vector with shift '
C?     CALL WRTVCD(VEC1,LU3,1,LBLK)
C?     write(6,*) ' Enforced stop after shifted sigma'
C?     stop ' Enforced stop after shifted sigma'
       
     
*
* ======================
*. Loop over iterations
* ======================
*
      DO 1000 ITER = 2, MAXIT+1
       ITERP = ITER -1
       IF(IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Info from iteration .... ', ITER
       END IF
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IROOT = 1
*
       EIGAPR = EIG(ITER-1,IROOT)
       FACHC = 1.0D0                      
       FACC  = -EIGAPR
       CALL VECSMD(VEC1,VEC2,FACC,FACHC,LU1,LU2,LU4,1,LBLK)
C           VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
         IF ( IPRT  .GE. 600 ) THEN
         WRITE(6,*) '  ( HX - EX ) '
         CALL WRTVCD(VEC1,LU4,1,LBLK)
       END IF
*  Strange place to put convergence but ....
       RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK) )
       RNRM(ITER-1,IROOT) = RNORM
       IF(RNORM.LT. TEST ) THEN
          CONVER = .TRUE.
          RTCNV(IROOT) = .TRUE.
       ELSE
          RTCNV(IROOT) = .FALSE.
          CONVER = .FALSE.
       END IF
       IF(ITER.GT.2.AND.
     & EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT).LT.THRES_E) CONVER = .TRUE.
       IF( ITER .EQ. MAXIT+1 .OR. CONVER ) GOTO 1001
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU3
       IF( .NOT. RTCNV(IROOT) ) THEN
         CALL REWINO(LUDIA)
         CALL REWINO(LU3)
         CALL REWINO(LU4)
*. Assuming diagonal preconditioner
         IPRECOND = 1
         CALL H0M1TD(LU3,LUDIA,LU4,LBLK,NP1+NP2+NQ,IPNTR,
     &               H0,-EIGAPR,H0SCR,XH0IX,
     &               NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)
         IF ( IPRT  .GE. 600) THEN
           WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
           CALL WRTVCD(VEC1,LU3,1,LBLK)
         END IF
*
         IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU4
           CALL REWINO(LU1)
           CALL REWINO(LU4)
           CALL REWINO(LUDIA)
*
           CALL H0M1TD(LU4,LUDIA,LU1,LBLK,Np1+Np2+NQ,
     &                 IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)

* Gamma = X(T) * (H0 - E) ** -1 * X
           GAMMA = INPRDD(VEC1,VEC2,LU1,LU4,1,LBLK)
* is X an eigen vector for (H0 - 1 ) - 1
           VNORM =
     &     SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU1,LU4,1,LBLK))
           IF(VNORM .GT. 1.0D-7 ) THEN
             IOLSAC = 1
           ELSE
             IOLSAC = 0
           END IF
           IF(IOLSAC .EQ. 1 ) THEN
             IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
             DELTA = INPRDD(VEC1,VEC2,LU1,LU3,1,LBLK)
             FACTOR = (-DELTA)/GAMMA
             IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
             IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
             CALL COPVCD(LU5,LU3,VEC1,1,LBLK)

             IF(IPRT.GE.600) THEN
               WRITE(6,*) ' Modified trial vector '
               CALL WRTVCD(VEC1,LU3,1,LBLK)
             END IF
*
            END IF
         END IF
*
*. 1.3 Orthogonalize to current vector
*
         OVLAP  = INPRDD(VEC1,VEC2,LU1,LU3,1,LBLK)
         ONE = 1.0D0
         CALL VECSMD(VEC1,VEC2,-OVLAP,ONE,LU1,LU3,
     &                LU4,1,LBLK)
*
         IF ( IPRT  .GE. 600 ) THEN
           WRITE(6,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
           CALL WRTVCD(VEC1,LU4,1,LBLK)
         END IF
*
*. 1.4 Normalize vector
*
         SCALE = INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK)
         FACTOR = 1.0D0/SQRT(SCALE)
         CALL SCLVCD(LU4,LU3,FACTOR,VEC1,1,LBLK)
         IF(IPRT.GE.600) THEN
           WRITE(6,*) '   normalized     (D-E)-1 *( HX - EX ) '
           CALL WRTVCD(VEC1,LU3,1,LBLK)
         END IF
*
       END IF
*
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
        CALL REWINO(LU3)
        CALL MV7(VEC1,VEC2,LU3,LU4,0,0)
*. Augment projected matrix
        CALL REWINO(LU1)
* <X!H! Delta>
         APROJ(2) = INPRDD(VEC1,VEC2,LU1,LU4,1,LBLK)
*<Delta!H!Delta>
         APROJ(3) = INPRDD(VEC1,VEC2,LU3,LU4,1,LBLK)
*. Diagonalize projected matrix
      CALL COPVEC(APROJ,WORK(KAPROJ),2*(2+1)/2)
C     write(6,*) ' work(aproj) '
C     call prsym(work(kaproj),2)
      CALL EIGENL(WORK(KAPROJ),AVEC,2,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR **2 + C1VAR **2 )
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
        EIG(ITER,IROOT) = WORK(KAPROJ)
*
C?     WRITE(6,*) ' APROJ(2),APROJ(3)',APROJ(2),APROJ(3)
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)')
     &   ' Eigenvalue and residual of iteration ..', ITER
         WRITE(6,'(2F21.13)') EIG(ITER,1)+EIGSHF, RNORM
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,2)
        WRITE(6,'(2X,E13.7)') EIG(ITER,1)
        CALL WRTMAT(AVEC,2,1,2,1)            
      END IF
*
*. Reset      
*
      CX = AVEC(1)
      CDELTA = AVEC(2)
*. Eigenvector
      CALL VECSMD(VEC1,VEC2,CX,CDELTA,LU1,LU3,LU5,1,LBLK)
      XNORM = INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)
      SCALE = 1.0D0/SQRT(XNORM)
      CALL SCLVCD(LU5,LU1,SCALE,VEC1,1,LBLK)
*. Sigma vector
      CXS = CX*SCALE
      CDELTAS = CDELTA*SCALE
      CALL VECSMD(VEC1,VEC2,CXS,CDELTAS,LU2,LU4,LU5,1,LBLK)
      CALL COPVCD(LU5,LU2,VEC1,1,LBLK)
*
      APROJ(1) = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
*
      IF(CONVER) GOTO 1001
 1000 CONTINUE
* ( End of loop over iterations )
 1001 CONTINUE
      ITER = ITERP 
      DO IROOT = 1, NROOT
       RNRM_CNV(IROOT) = RNRM(ITER,IROOT)
       FINEIG(IROOT) = EIG(ITER,IROOT) + EIGSHF
      END DO
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINO(LU1)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F21.10)
          IF(IPRT.GE.400) THEN
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
      END IF
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
*. Collect info for root NROOT
      
C. Test  LU3
C     WRITE(6,*) ' Test copy of LU3 in MICDV5 '
C     CALL COPVCD(LU3,LU4,VEC1,1,LBLK)
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
	  
      SUBROUTINE VC3SMD(VEC1,VEC2,FAC1,FAC2,FAC3,
     &                  LU1,LU2,LU3,LU4,IREW,LBLK)
*
*
* LU4 = FAC1*LU1 + FAC2*LU2 + FAC3*LU3
*
*
* BLocked vectors in usual format as defined by LBLK\
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
*
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU3,LBLK)
        CALL REWINE( LU4,LBLK)
      END IF
*
* LOOP OVER BLOCKS OF VECTOR
*
 1000 CONTINUE
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
          NBL3 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
          READ(LU3) NBL3
          WRITE(LU4) NBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
          CALL IFRMDS( NBL3,1,-1,LU3)
          CALL ITODS ( NBL1,1,-1,LU4)
        END IF
		
        IF( NBL1 .NE. NBL2.OR.NBL2.NE.NBL3 ) THEN
         WRITE(6,'(A,3I5)') 'DIFFERENT BLOCKSIZES IN VC3SMD ',
     &   NBL1,NBL2,NBL3
         STOP ' INCOMPATIBLE BLOCKSIZES IN VC3SMD '
        END IF

      IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
        CALL FRMDSC(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
        CALL FRMDSC(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
        IF( NBL1 .GT. 0 )
     &  CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
        CALL FRMDSC(VEC2,NBL1,KBLK,LU3,IMZERO,IAMPACK)
        ONE = 1.0D0
        IF(NBL1.GT.0)
     &  CALL VECSUM(VEC1,VEC1,VEC2,ONE,FAC3,NBL1)
        CALL TODSC(VEC1,NBL1,KBLK,LU4)
      END IF
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
      RETURN
      END
	  
      SUBROUTINE MICGCG(MV8,LU1,LU2,LU3,LU4,LU5,LUDIA,VEC1,VEC2,
     &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
     &                  LUPROJ,LUPROJ2,VFINAL,IPRT)
*
* Solve set of linear equations
*
*             AX = B
*
* with preconditioned conjugate gradient method for
* case where two complete vectors can be stored in core
*
* Initial appriximation to solution must reside on LU1
* LU2 must contain B.All files are  overwritten
*
*
* Final solution vector is stored in LU1
* A scalar w can be added to the diagonal of the preconditioner
*
* If LUPROJ .NE. 0 , the optimization subspace is restricted to be orthogonal
* to the first vector in LUPROJ.
* The vector used to orthogonalize is saved on LUPROJ2
*
* Version using blocks of vectors
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*),ERROR(MAXIT+1)
      REAL*8 INPRDD
      LOGICAL CONVER
*
      EXTERNAL MV8
*
      CALL LUCIAQENTER('MICGC')
      NTEST = 02
      NTEST = MAX(NTEST,IPRT)
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' Welcome to MICGCG '
        WRITE(6,*) ' =================='
        WRITE(6,*)
*
C?    WRITE(6,*) ' NTEST ,LU1,LU2,LU3 = ', NTEST,LU1,LU2,LU3
      END IF
      CONVER = .FALSE.
      ITER = 1
*
      LBLK = -1
*
      ONE = 1.0D0
      ONEM = -1.0D0
      ZERO = 0.0D0
*. Overlap between LUPROJ and LUPROJ2
      IF(LUPROJ.GT.0) THEN 
        X12 = INPRDD(VEC1,VEC2,LUPROJ,LUPROJ2,1,LBLK)
      ELSE
        X12 = 0.0D0
      END IF
C?    WRITE(6,*) ' MICGCG : X12 = ', X12
*
* =============
* Initial point
* =============
*
*.R = B - (A)*X on LU2
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector on LU1 '
        CALL WRTVCD(VEC1,LU1,1,LBLK)
        WRITE(6,*) ' Vector on LU2 '
        CALL WRTVCD(VEC1,LU2,1,LBLK)
      END IF
      CALL MV8(VEC1,VEC2,LU1,LU3)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector on LU3 '
        CALL WRTVCD(VEC1,LU3,1,LBLK)
      END IF
*
C          VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LU2,LU3,LU4,1,LBLK)
      CALL COPVCD(LU4,LU2,VEC1,1,LBLK)
*
*
      RNORM = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
      ERROR(1) = SQRT(RNORM)
      IF(ERROR(1).LE.TEST) THEN
*. Convergence in one shot- you are lucky or have
* supplied a vaninshing RHS
        NITER = 0
        CONVER = .TRUE.
        GOTO 1001
      END IF
*
*. Preconditioner H times initial residual, H * R on LU4
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Diagonal and input to diagonal '
       CALL WRTVCD(VEC1,LUDIA,1,LBLK)
       CALL WRTVCD(VEC1,LU2  ,1,LBLK)
       WRITE(6,*) ' SHIFT = ', W
      END IF
      CALL DMTVCD_OLD(VEC1,VEC2,LUDIA,LU2,LU4,W,1,1,LBLK)
      IF(LUPROJ.NE.0) THEN
        OVLAP = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
        FACTOR = (-OVLAP)/X12
        CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU4,LUPROJ2,LU3,1,LBLK)
        CALL COPVCD(LU3,LU4,VEC1,1,LBLK)
        OVLAP2 = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
        WRITE(6,*) ' Updated overlap of trial vector ', OVLAP2
      END IF
*. GAMMA = <R!H!R>
      GAMMA = INPRDD(VEC1,VEC2,LU2,LU4,1,LBLK)
*. P = RHO * H*R on LU3
      RHO = 1.0D0
      CALL SCLVCD(LU4,LU3,RHO,VEC1,1,LBLK)
*.S = AP on LU4
      CALL MV8(VEC1,VEC2,LU3,LU4)
*
* ====================
* Loop over iterations
* ====================
*
      NITER = 0
      DO 1000 ITER = 1, MAXIT
*
* Vectors on files :
*     X on LU1
*     R on LU2
*     P on LU3
*  S=AP on LU4
*     H on LUDIA
 
        NITER = NITER + 1
       IF ( NTEST .GE. 2 )
     & WRITE(6,*) ' INFORMATION FROM ITERATION... ',ITER
*.    D = <P!S>
        D = INPRDD(VEC1,VEC2,LU3,LU4,1,LBLK)
        C = RHO * GAMMA
        A = C/D
*.    R = R - A * S on LU2
        CALL VECSMD(VEC1,VEC2,ONE,-A,LU2,LU4,LU5,1,LBLK)
        CALL COPVCD(LU5,LU2,VEC1,1,LBLK)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Residual on LU2 '
          CALL WRTVCD(VEC1,LU2,1,LBLK)
        END IF
*.    new residual has been obtained , check for convergence
        RNORM = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
        RNORME = MAX(RNORM,0.0D0)
        ERROR(ITER+1) = SQRT(RNORME)
        IF(NTEST.GE.2) WRITE(6,*) ' Norm of residual ', ERROR(ITER+1)
*.    X = X + A * P
C?      WRITE(6,*) ' MICGCG : A = ', A
        CALL VECSMD(VEC1,VEC2,ONE,A,LU1,LU3,LU5,1,LBLK)
        CALL COPVCD(LU5,LU1,VEC1,1,LBLK)
        XNORM = INPRDD(VEC1,VEC2,LU1,LU1,1,LBLK)
        IF(NTEST.GE.5) WRITE(6,*) ' Norm of x = ', XNORM
*
        IF( SQRT(RNORME) .LT. TEST ) THEN
           CONVER = .TRUE.
           GOTO 1001
        ELSE
           CONVER = .FALSE.
*
* ============================
*. Prepare for next iteration
* ============================
*
*.H * R on LU4
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Diagonal and input to diagonal '
             CALL WRTVCD(VEC1,LUDIA,1,LBLK)
             CALL WRTVCD(VEC1,LU2  ,1,LBLK)
             WRITE(6,*) ' SHIFT = ', W
           END IF
*
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU2,LU4,W,1,1,LBLK)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Preconditioner times residual '
             CALL WRTVCD(VEC1,LU4,1,LBLK)
           END IF
           IF(LUPROJ.NE.0) THEN
             OVLAP = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
             FACTOR = (-OVLAP)/X12
             CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU4,LUPROJ2,LU5,1,LBLK)
             CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
             OVLAP2 = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
C?           WRITE(6,*) ' Updated overlap of trial vector ', OVLAP2
*. Overlap between X and LUPROJ
             OVLAP3 = INPRDD(VEC1,VEC2,LUPROJ,LU1,1,LBLK)
             WRITE(6,*) ' Overlap between LU1 and LUPROJ ', OVLAP3
           END IF
*. GAMMA = <R!H!R>
           GAMMA = INPRDD(VEC1,VEC2,LU2,LU4,1,LBLK)
           B = GAMMA/C
*. P = RHO*(H*R + B*P) on LU3
           RHO = 1.0D0
           CALL VECSMD(VEC1,VEC2,ONE,B,LU4,LU3,LU5,1,LBLK)
           CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*.S = AP on LU4
           CALL MV8(VEC1,VEC2,LU3,LU4)
*.End of prepations for next iteration
        END IF
*
 1000 CONTINUE
 1001 CONTINUE
*
      IF(CONVER) THEN
        VFINAL = ERROR(NITER+1)
      ELSE
        VFINAL = ERROR(MAXIT+1)
      END IF
      IF(NTEST .GT. 0 ) THEN
*
      IF(CONVER) THEN
       WRITE(6,1010) NITER  ,ERROR(NITER+1)
 1010  FORMAT(1H0,'  convergence was obtained in...',I3,' iterations',/,
     +        1H ,'  norm of residual..............',E13.8)
      ELSE
       WRITE(6,1020) MAXIT ,ERROR(MAXIT+1)
 1020  FORMAT(1H0,' convergence was not obtained in',I3,'iterations',/,
     +        1H ,' norm of residual...............',E13.8)
      END IF
*
      END IF
      IF(NTEST.GE. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(1H0,' solution to set of linear equations')
       CALL WRTVCD(VEC1,LU1,1,LBLK)
C?     write(6,*) ' Matrix times solutiom through another cal to MV 8'
C?     CALL MV8(VEC1,VEC2,0,0)
C?     call wrtmat(vec2,1,nvar,1,nvar)
      END IF
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(1H0,10X,'iteration point     norm of residual')
      DO 350 I=1,NITER+1
       II=I-1
       WRITE(6,1050)II,ERROR(I)
 1050  FORMAT(1H ,12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
      CALL LUCIAQEXIT('MICGC')
      RETURN
      END 
	  
      SUBROUTINE PROP_PERT(LU0,LUN,N,ISM,ISPC)
*
* Perturbation expansion of one-electron properties       
*
* It is assumed that this calculation is preceded  by 
* a call to the perturbation routine to obtain the 
* wave function corrections to the neutral state.
*
* Input       
*       LUN : File containing wave function corrections
*       LU0 : File containing reference wave funcrtion
*         N : Max order of expansion
*      ISM : Symmetry of reference state
*      ISPC : Space of referencestate
*
*
* Jeppe Olsen, April 98 ( on the train for once )
c      IMPLICIT REAL*8 (A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cintfo.inc'
      REAL *8 INPRDD   
*. Local scratch
      PARAMETER(MXNORD = 100)
*
      NTEST = 5
*
      WRITE(6,*) 
      WRITE(6,*) ' ============================ '
      WRITE(6,*) ' PROP_PERT is now in CONTROL '
      WRITE(6,*) ' ============================ '
      WRITE(6,*)
      WRITE(6,*) ' N= ', N
      IF(IRELAX.EQ.0) THEN
        WRITE(6,*) ' Property evaluated as expectation value'
      ELSE
        WRITE(6,*) ' Property evaluated as derivative '
      END IF

* a bit on files :
* LUSC36 is LUN.   
* Two additional scratch files to be used are  LUSC1 and LUSC2
* 
      LBLK = -1
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PROPPT')
*
*
*     ========================
* 1 : Local memory allocation 
*     ========================
*
*. Allocate space for two vector chunks
      CALL GET_3BLKS(KLVEC1,KLVEC2,KVEC3)
* space for one-body Density matrices through order n
      NMAT = N+1     
      LENGTH = NMAT * NTOOB ** 2
      CALL MEMMAN(KLDEN1,LENGTH,'ADDL  ',2,'DENN1 ')
*. And an extra set of density matrices
      CALL MEMMAN(KLDEN1P,LENGTH,'ADDL  ',2,'DENN1P')
*. Two-electron densities
      IF(IRELAX.EQ.0) THEN
        KLDEN2 = KLDEN1
        KLFOCK = 1
      ELSE
        LENGTH = NMAT * NTOOB ** 2 * (NTOOB**2 + 1)/2
        CALL MEMMAN(KLDEN2,LENGTH,'ADDL  ',2,'DENN2 ')
*. And relaxation terms to  one-electron density
        CALL MEMMAN(KLRELR1,NTOOB**2,'ADDL  ',2,'RELR1 ')
*. Space for Fock matrices 
        LFOCK = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
        CALL MEMMAN(KLFOCK,NMAT*LFOCK,'ADDL  ',2,'FOCK_N')
      END IF
*. A scratch matrix ( not a nice thing to say about a matrix )
      LENGTH =  2*NINT1    
      CALL MEMMAN(KLFSCR,LENGTH,'ADDL  ',2,'FSCR  ')
      CALL MEMMAN(KLSCR4,LENGTH,'ADDL  ',2,'SCR4  ')
*. S(i,j) matrix for overlap of corrections
      CALL MEMMAN(KLSIJ,(N+1)**2,'ADDL  ',2,'KLSIJ ')
      CALL MEMMAN(KLSI , N+1    ,'ADDL  ',2,'KLSI  ')
*
* =========================================================================
*.1 :  overlap of correction vectors ( intermediate normalization  assumed )
* =========================================================================
*
*. Sij(i,j) = <i!j>
      CALL REWINO(LUN)
      DO I = 1, N
*. LUN is positioned at end of vector I-1, copy vector I to LUSC1
        CALL REWINO(LUSC1)
        CALL COPVCD(LUN,LUSC1,WORK(KLVEC1),0,LBLK)
*
        CALL REWINO(LUN)
        DO J = 1, I
          IJ = I*(I-1)/2 + j
          CALL REWINO(LUSC1)
          WORK(KLSIJ-1+IJ) = 
     &    INPRDD (WORK(KLVEC1),WORK(KLVEC2),LUSC1,LUN,0,LBLK)
        END DO
      END DO
* SI(i) = sum(j=1,i-1)S(j,i-j)
      DO I = 1, N
        X = 0.0D0
        DO J = 1, I-1
          IMJ = I-J
          IJ = MAX(J,IMJ)*(MAX(J,IMJ)-1)/2+MIN(J,IMJ)
          X = X + WORK(KLSIJ-1+IJ)
        END DO
        WORK(KLSI-1+I) = X
      END DO
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' The S(i,j) Matrix '
        WRITE(6,*) ' ================= '
        CALL PRSYM(WORK(KLSIJ),N)
        WRITE(6,*)
        WRITE(6,*) ' The S(i) array '
        WRITE(6,*) ' ================= '
        CALL WRTMAT(WORK(KLSI),N,1,N,1)
      END IF
*
* ===============================================
* 2 : Construct density matrices through order N
* ===============================================
*
*
* 2a : One-body densities Rho1(N) = Sum(M=0,N) <0(M)!E!0(N-M)>
*
      ILRHO2 = 0
      LRHO2 = 0
      LRHO1 = NTOOB**2
*. No print in density matrices
      IPRDEN_SAVE = IPRDEN
      IPRDEN = 0
      I12_SAVE = I12
      I12 = 1
      DO K = 0, N
        CALL PERTDN(K,LU0,LUN,ISM,ISPC,WORK(KLVEC1),WORK(KLVEC2),
     &       WORK(KLDEN1+(K-0)*LRHO1),
     &       WORK(KLDEN2+(K-0)*LRHO2),LUSC1,LUSC2,0)
      END DO
      IPRDEN = IPRDEN_SAVE
      I12 = I12_SAVE
*
* Change the densities so the correspond to order expansion of
* normalized wf
* Rho'(n) = Rho(n) - sum(j=1,n) Si(j)Rho'(n-j)
*
      ONE = 1.0D0
      DO I = 0, N
        CALL COPVEC(WORK(KLDEN1 +(I-0)*LRHO1),
     &              WORK(KLDEN1P+(I-0)*LRHO1),LRHO1)
*
        DO J = 1, I
          FACTOR = -WORK(KLSI-1+J)
          IOFF = KLDEN1P+(I-0)*LRHO1
          JOFF = KLDEN1P+(I-J-0)*LRHO1
          CALL VECSUM(WORK(IOFF),WORK(IOFF),WORK(JOFF),
     &                ONE,FACTOR,LRHO1)
        END DO
*
C?      WRITE(6,*) ' Density correction for NORMALIZED wf '
C?      CALL WRTMAT(WORK(KLDEN1P+I*LRHO1),NTOOB,NTOOB,NTOOB,NTOOB)
      END DO
*
      IF(IRELAX.EQ.1) THEN
* Set up <0(0)!   !0(N)> densities and symmetrize
        ILRHO2 = 1
        LRHO2 = NTOOB**2*(NTOOB**2 + 1)/ 2 
        LRHO1 = NTOOB**2
*. No print in density matrices
        IPRDEN_SAVE = IPRDEN
        IPRDEN = 0
        I12_SAVE = I12
        I12 = 2
        DO K = 0, N
          CALL PERTDN(K,LU0,LUN,ISM,ISPC,WORK(KLVEC1),WORK(KLVEC2),
     &         WORK(KLDEN1+(K-0)*LRHO1),
     &         WORK(KLDEN2+(K-0)*LRHO2),LUSC1,LUSC2,1)
*. Well it was only 0.5 Times above term we wanted so 
          HALF = 0.5D0
          CALL SCALVE(WORK(KLDEN1+(K-0)*LRHO1),HALF,LRHO1)
          CALL SCALVE(WORK(KLDEN2+(K-0)*LRHO2),HALF,LRHO2)
        END DO
        IPRDEN = IPRDEN_SAVE
        I12 = I12_SAVE
* Order expansion of Fock matrix
C            GET_FN(FN,DEN1N,DEN2N,MAXN,LFOCK)
        CALL GET_FN(WORK(KLFOCK),WORK(KLDEN1),WORK(KLDEN2),N,LFOCK)
*. Restore zero order two-body densities  
        CALL COPVEC(WORK(KLDEN2),WORK(KRHO2),LRHO2)
      END IF
*
*. Properties for each order
*
      DO IORD = 0, N
        WRITE(6,*)
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' Information for order ', IORD
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        IF(IRELAX.EQ.1) THEN
*. Relaxation contribution to density
*. Restore zero order one-body densities  
          CALL COPVEC(WORK(KLDEN1),WORK(KRHO1),LRHO1)
          INOFF = KLFOCK + (IORD-0)*LRHO1
          CALL RESPDEN_FROM_F(WORK(INOFF),WORK(KLRELR1)) 
        END IF
*. 
        III = KLDEN1P + (IORD-0)*LRHO1
        CALL COPVEC(WORK(III),WORK(KRHO1),LRHO1)
*. No natural orbital analysis, so
        I_EXP_OR_TRA = 2
        CALL ONE_EL_PROP(I_EXP_OR_TRA,IRELAX,WORK(KLRELR1))
      END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PROPPT')
*
      RETURN
      END
      SUBROUTINE PERTDN
     &(N,LU0,LUN,ISM,ISPC,VEC1,VEC2,RHO1N,RHO2N,LUSC1,LUSC2,
     & I_ONLY_0N)
*
* Construct one body density matrix of order N
*
*      Jeppe + Dage, Nov. 11 1995
*                    Debugged Jan 31 '97
*
* Note : I12 added, April 98
*        I_ONLY_0N added, May 99
*
*
* If I_ONLY_0N only the <0(0)!  |0(N)> terms are included
*
c      IMPLICIT REAL*8 (A-H,O-Z)
*
*. Should not be called with ICISTR = 1
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'orbinp.inc'
C     INCLUDE 'clunit.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
*. Output 
      DIMENSION RHO1N(*),RHO2N(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'PERTDN')
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
      CALL MEMMAN(KLDEN1,LRHO1,'ADDL  ',2,'KLDEN1')
      IF(I12.EQ.2) THEN
        CALL MEMMAN(KLDEN2,LRHO2,'ADDL  ',2,'KLDEN2')
      END IF
*
      LBLK = -1
      ZERO = 0.0D0
      CALL SETVEC(RHO1N,ZERO,LRHO1)
      IF (I12.EQ.2) THEN
        CALL SETVEC(RHO2N,ZERO,LRHO2)
      END IF
*
      DO L = 0, N
C?      write(6,*) ' Will load next pair of vectors '
        NMINL = N - L
CTOBE   IF(L.LE.NMINL) THEN
*. put correction vector L and NMINL on LUSC1 and LUSC2, respectively
          IF(L.EQ.0) THEN
             CALL COPVCD(LU0,LUSC1,VEC1,1,LBLK)
          ELSE 
             CALL SKPVCD(LUN,L-1,VEC1,1,LBLK)
             CALL REWINO(LUSC1)
             CALL COPVCD(LUN,LUSC1,VEC1,0,LBLK)
          END IF
*
          IF(NMINL.EQ.0) THEN 
             CALL COPVCD(LU0,LUSC2,VEC1,1,LBLK)
          ELSE 
             CALL SKPVCD(LUN,NMINL-1,VEC1,1,LBLK)
             CALL REWINO(LUSC2)
             CALL COPVCD(LUN,LUSC2,VEC1,0,LBLK)
          END IF
C       write(6,*) ' next pair of vectors loaded '
* Do the densi
          IF(I_ONLY_0N.EQ.0.OR.L.EQ.0.OR.L.EQ.N) THEN
            LEQR = 0
            XDUM = 0.0D0
            CALL DENSI2(I12,WORK(KLDEN1),WORK(KLDEN2),VEC1,VEC2,
     &                  LUSC1,LUSC2,EXPS2R,0,XDUM,XDUM,XDUM,XDUM,1)
*
CTOBE       IF(L.NE.NMINL) THEN
*. The matrix <L! E !NMINL> was calculated, add <NMINL! E ! L> 
*. as simple transposition
CTOBE          CALL TRPAD(WORK(KLDEN),ONE,NTOOB)
CTOBE       END IF
            ONE = 1.0D0
            CALL VECSUM(RHO1N,RHO1N,WORK(KLDEN1),ONE,ONE,LRHO1)
            IF(I12.EQ.2) THEN
              CALL VECSUM(RHO2N,RHO2N,WORK(KLDEN2),ONE,ONE,LRHO2)
            END IF
CTOBE     END IF
          END IF
*         ^ End I_ONLY_0N check

      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Density matrix of order in perturbation ', N
        WRITE(6,*) ' ==========================================='
        WRITE(6,*)
        WRITE(6,*) ' One-body density '
        WRITE(6,*) ' ================ '
        CALL WRTMAT(RHO1N,NTOOB,NTOOB,NTOOB,NTOOB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-body density '
          WRITE(6,*) ' ================ '
          CALL PRSYM(RHO2N,NTOOB**2)
        END IF
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ', IDUM,'PERTDN')
*
      RETURN
      END
      SUBROUTINE MICDV4_ENLMD(MV7,VEC1,VEC2,LU1,LU2,RNRM,EIG,
     &                  FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Version using H0 + Lambda V as Sigma routine
*
* Input :
* =======
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
#include "mafdecls.fh"
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
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
       EXTERNAL MV7
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
 
       IOLSTM = 1
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
         CALL ENLMD(VEC1,VEC2,LU3,LU4)
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
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,dbl_mb(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGENL(dbl_mb(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(6,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
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
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
           CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU1,LU4,0,LBLK)
   60    CONTINUE
         IF ( IPRT  .GE. 10 ) THEN
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
           IF(IPRT.GE.3) THEN
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
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)
C               H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG,IPRECOND)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
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
             IF ( IPRT  .GE. 10 ) THEN
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
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)
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
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
                CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
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
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(6,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU3,LBLK)
           CALL SCLVCD(LU3,LU1,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
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
        CALL ENLMD(VEC1,VEC2,LU3,LU4)
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
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,dbl_mb(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGENL(dbl_mb(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = dbl_mb(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        dbl_mb(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = dbl_mb(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(6,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(6,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
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
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(6,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
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
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
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
      SUBROUTINE MICDV4_H0LVP(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Version using H0 + Lambda V as Sigma routine
*
* Input :
* =======
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
#include "mafdecls.fh"
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
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
       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(6,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
 
       IOLSTM = 1
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
         CALL H0LVP(VEC1,VEC2,LU3,LU4)
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
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,dbl_mb(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGENL(dbl_mb(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(6,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
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
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
           CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU1,LU4,0,LBLK)
   60    CONTINUE
         IF ( IPRT  .GE. 10 ) THEN
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
           IF(IPRT.GE.3) THEN
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
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)
C               H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG,IPRECOND)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
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
             IF ( IPRT  .GE. 10 ) THEN
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
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT,IPRECOND)
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
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
                CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
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
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(6,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU3,LBLK)
           CALL SCLVCD(LU3,LU1,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
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
        CALL H0LVP(VEC1,VEC2,LU3,LU4)
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
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,dbl_mb(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGENL(dbl_mb(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = dbl_mb(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        dbl_mb(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = dbl_mb(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(6,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(6,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
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
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(6,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
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
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
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
	  
      SUBROUTINE H0TVF(VEC1,VEC2,LUC,LUHC,LU0,LUSCR1,E0,LBLK)
*
* Multiply vector in LUC with H0 where H0 is defined as 
*
* H0 = (1-|0><0|) F (1-|0><0>) + E0 |0><0>
*
* Where is one-electron operator defined by WORK(KFI)
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. For communicating with sigma routine
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      include 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'chktyp.inc'
      DIMENSION VEC1(*),VEC2(*)
      REAL*8 INPRDD
*
      WRITE(6,*) ' H0TVF : LUC, LUHC, LU0, LUSCR1, E0 ',
     &                     LUC, LUHC, LU0, LUSCR1, E0
*.
      SC0 = -3006.56D0
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
      WRITE(6,*) ' MV7 will be called in a few NANOSECONDS'
      CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
      I12 = 1
      ICHKTYP = 1
      CALL MV7(VEC1,VEC2,LUSCR1,LUHC,0,0)
      ICHKTYP = 0
      I12 = 2
      CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
*. Orthogonalize LUHC to LU0
      IF(LU0.GT.0) THEN
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
      NTEST = 000
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' results from H0TVF '
        WRITE(6,*) ' ==================='
        write(6,*) ' SC0, SSIGMA0 ', SC0,SSIGMA0
      END IF
      IF(NTEST.GE.100) THEN
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
      SUBROUTINE H0LVP(VEC1,VEC2,LLUIN,LLUOUT)
*
* H0 + Lambda V times vector on LLUIN
*
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
	  
      INCLUDE 'wrkspc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      COMMON/CENOT/E0
*
* H0 + Lambda V = (1-Lambda) H0 + Lambda H
      LBLK = -1
*. H0 * LLUIN on LUSC53
       WRITE(6,*) ' E0 in H0LVP ', E0
C          H0TVF(VEC1,VEC2,LUC,LUHC,LU0,LUSCR1,E0,LBLK)
      CALL H0TVF(VEC1,VEC2,LLUIN,LUSC53,LUSC51,LUSC52,E0,LBLK)
* H * LLUIN on LUSC52
      CALL MV7(VEC1,VEC2,LLUIN,LUSC52,0,0)
* (1-Lambda)*H0 + Lambda
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      FAC1 = 1.0D0 - XLAMBDA
      FAC2 = XLAMBDA
      CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUSC53,LUSC52,LLUOUT,1,LBLK)
*
      RETURN
      END
      SUBROUTINE ENLMD(VEC1,VEC2,LLUIN,LLUOUT)
*
* H0 + Lambda V times vector on LLUIN with Epstein-Nesbet partitioning 
*
* = ((1-Lambda)H_diag + Lambda * H) * LLUIN
*
* Diagonal is assumed stored on LUDIA
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
*
* H0 + Lambda V = (1-Lambda) H0 + Lambda H
      LBLK = -1
*. H0 * LLUIN on LUSC53
C     DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
      ZERO = 0.0D0
      CALL DMTVCD(VEC1,VEC2,LLUIN,LUDIA,LUSC53,ZERO,1,0,LBLK)
* H * LLUIN on LUSC52
      CALL MV7(VEC1,VEC2,LLUIN,LUSC52,0,0)
* (1-Lambda)*H0 + Lambda
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      FAC1 = 1.0D0 - XLAMBDA
      FAC2 = XLAMBDA
      CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUSC53,LUSC52,LLUOUT,1,LBLK)
*
      RETURN
      END
      SUBROUTINE FAM(FA)
*
*. Construct active fock matrix
*
* Jeppe Olsen, July 2010
      INCLUDE 'wrkspc.inc'
*
CNW   DIMENSION FA(*)
      integer FA
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*
      CALL FAS(FA,dbl_mb(KRHO1),IBSO,NSMOB,NTOOBS,NACOB,NINOB,
     &         IREOST)
*
      RETURN
      END
*
      SUBROUTINE FAS(FA,RHO1,IBSO,NSMOB,NTOOBS,NACOB,NINOB,
     &           ISTOB)
*
* Active Fock matrix - in complete orbital space
*
*     FA(I,J) =  sum(k,l: active) ((ij!kl)-0.5*(il!kj))*rho1(kl)
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION FA(*),RHO1(NACOB,NACOB)
      DIMENSION RHO1(NACOB,NACOB)
      integer FA
      INTEGER IBSO(*),NTOOBS(*),ISTOB(*)
*
      NTEST = 00
      IF(NTEST.NE.0.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) 
       WRITE(6,*) ' ============='
       WRITE(6,*) '  FAS calling '
       WRITE(6,*) ' ============='
       WRITE(6,*) 
      END IF
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' Input one-body density matrix'
       CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
      END IF
*. Zero FA
*
*.  Assume spatial symmetric fock matrix
*
      ibuf = imnmx(lobsm,nsmob,2)
      ibuf = ibuf*(ibuf+1)/2
      if (.not.ma_push_get(mt_dbl,ibuf,'buffer',lbuf,kbuf)) call
     $   errquit('fis_old:push_get failed', ibuf, MA_ERR)
      IJSM = 1
      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
        IF(JSM.NE.0) THEN
          call setvec(dbl_mb(kbuf),0.0d0,ibuf)
          istart = ij + 1
          ii = 0
          DO I = IBSO(ISM),IBSO(ISM) + NTOOBS(ISM)-1
            DO J = IBSO(JSM),I                          
              IP = ISTOB(I)
              JP = ISTOB(J)
               IJ= IJ + 1
CNW            FA(IJ) = 0.0D0
               DO IA = NINOB+1, NACOB+NINOB
                 DO IB = NINOB+1, NACOB+NINOB
                   IF(RHO1(IA-NINOB,IB-NINOB).NE.0.0D0)        
     &             dbl_mb(kbuf+ii) = dbl_mb(kbuf+ii)
     &           + RHO1(IA-NINOB,IB-NINOB)
     &           *(GTIJKL(IP,JP,IA,IB)-0.5*GTIJKL(IP,IB,IA,JP))
                 END DO
               END DO
               ii = ii + 1
            END DO
          END DO
          if (ij.ge.istart)
     &    call ga_acc(FA,istart,IJ,1,1,dbl_mb(kbuf),1,1.0d0)
        END IF
      END DO
      if (.not.ma_pop_stack(lbuf)) call
     $   errquit('fis_old:pop_stack failed', lbuf, MA_ERR)
*
      IF(NTEST.NE.0) THEN
       if (ga_nodeid().eq.0) then
       WRITE(6,*) ' FA in Symmetry blocked form '
       WRITE(6,*) ' ============================'
       WRITE(6,*) 
       endif
       ISYM = 1
CNW    CALL APRBLM2(FA,NTOOBS,NTOOBS,NSMOB,ISYM)
       call ga_print(FA)
      END IF
* 
      RETURN
      END
      FUNCTION ECORE_TERM (IDOH1,IDOH2,IDOPH1,IDOPH2)
*
*. Construct contributions to core/inactive  energy.
*
* IDOH1 = 1:  Include terms inactive orbitals in one-electron operator
* IDOH2 = 1:  Include terms from inactive orbitals in two-electron
* operator
* IDOPH1= 1: Include terms from ph-reorganization of one-electron
* operator
* IDOPH2= 1: Include terms from ph-reorganization of two-electron
* operator
*
* If I_USE_SIMTRH = 0 input and output matrices are assumed lower half packed
*                 = 1 Input and output matrices are assumed complete blocks
*                     (not active pt)
*
* Jeppe Olsen, August 2010
*
      INCLUDE 'wrkspc.inc'
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'oper.inc'
*
*
      NTEST = 00
      IF (NTEST.GE.10) THEN
      END IF

      IF(I_USE_SIMTRH.EQ.0) THEN
        ECORE_TERM =  ECORE_TERMS(WORK(KH),
     &                IBSO,NSMOB,ITPFSO,IPHGAS,NTOOBS,NTOOB,IREOST,NGAS,
     &                NGAS,IDOH1,IDOH2,IDOPH1,IDOPH2)
      ELSE
        ECORE_TERM =  ECORE_TERMSA(WORK(KH),
     &                IBSO,NSMOB,ITPFSO,IPHGAS,NTOOBS,NTOOB,IREOST,NGAS,
     &                IDOH1,IDOH2,IDOPH1,IDOPH2)
      END IF
*
      RETURN
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' ECORE_TERM reporting '
       WRITE(6,'(A,4I4)') ' Input: IDOH1, IDOH2, IDOPH1, IDOPH2 = ', 
     &                             IDOH1, IDOH2, IDOPH1, IDOPH2
       WRITE(6,'(A,E22.15)') ' ECORE_TERM = ', ECORE_TERM
      END IF
*
      END
      FUNCTION ECORE_TERMS(H,
     &         IOBSM,NSMOB,ITPFSO,IPHGAS,LOBSM,NORBT,ISTOB,NGAS,
     &         IDOH1, IDOH2, IDOPH1, IDOPH2)
*
*. Construct contributions to core/inactive  energy.
*
* IDOH1 = 1:  Include terms inactive orbitals in one-electron operator
* IDOH2 = 1:  Include terms from inactive orbitals in two-electron
* operator
* IDOPH1= 1: Include terms from ph-reorganization of one-electron
* operator
* IDOPH2= 1: Include terms from ph-reorganization of two-electron
* operator
*
* H is the standard one-electron operator 
*
* Jeppe Olsen ( I admit ) August 2010
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(*)
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
      INTEGER ITPFSO(*), IPHGAS(*)
*. To get rid of annoying and incorrect compiler warnings
      IIOFF = 0
*
      NTEST = 00
*
      ECC = 0.0D0
      IJSM = 1
*
*. One-electron part 
*
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IIOFF = IIOFF + LOBSM(ISM-1)*(LOBSM(ISM-1)+1)/2
        END IF
        II = IIOFF-1
        DO I = IOBSM(ISM),IOBSM(ISM)+LOBSM(ISM)-1
          II = II + (I-IOBSM(ISM)+1) 
          I_INCLUDE = 0
          IF(IDOPH1.EQ.1.AND.ITPFSO(I).GT.0.AND.ITPFSO(I).LE.NGAS) THEN
           IF (IPHGAS(ITPFSO(I)).EQ.2) I_INCLUDE = 1
          END IF
          IF(IDOH1.EQ.1.AND.ITPFSO(I).EQ.0) I_INCLUDE = 1
C?        WRITE(6,*) ' I, I_INCLUDE = ', I,I_INCLUDE
          IF(I_INCLUDE.EQ.1) THEN
C?          WRITE(6,*) ' Contribution to ECC from I =', I
            ECC = ECC + 2*H(II)
C?          WRITE(6,*) ' Updated ECC = ', ECC
          END IF
        END DO
      END DO
C?    WRITE(6,*) ' one-electron part to ECC ', ECC
*
*. Two-electron part
*
      ECC2 = 0.0d0
      IF(IDOH2.EQ.1.OR.IDOPH2.EQ.1) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IOBSM(ISM), IOBSM(ISM) + LOBSM(ISM)-1
          DO J = IOBSM(JSM), IOBSM(JSM) + LOBSM(JSM)-1
*
            I_INCLUDE = 0
            IF(IDOPH2.EQ.1.AND.ITPFSO(I).GT.0.AND.ITPFSO(I).LE.NGAS)THEN
             IF (IPHGAS(ITPFSO(I)).EQ.2) I_INCLUDE = 1
            END IF
            IF(IDOH2.EQ.1.AND.ITPFSO(I).EQ.0) I_INCLUDE = 1
*
            J_INCLUDE = 0
            IF(IDOPH2.EQ.1.AND.ITPFSO(J).GT.0.AND.ITPFSO(J).LE.NGAS)THEN
             IF (IPHGAS(ITPFSO(J)).EQ.2) J_INCLUDE = 1
            END IF
            IF(IDOH2.EQ.1.AND.ITPFSO(J).EQ.0) J_INCLUDE = 1
C?          WRITE(6,*) ' I, J, I_INCLUDE, J_INCLUDE = ', 
C?   &                   I, J, I_INCLUDE, J_INCLUDE
*
            IP = ISTOB(I)
            JP = ISTOB(J)
*
            IF(I_INCLUDE.EQ.1.AND.J_INCLUDE.EQ.1) THEN
C?            WRITE(6,*) ' Contribution to ECC from I,J =', I,J
              ECC2 = ECC2 +2*GTIJKL(IP,IP,JP,JP)-GTIJKL(IP,JP,JP,IP)
C?            WRITE(6,*) ' Updated ECC = ', ECC
            END IF
          END DO
          END DO
        END DO
        END DO
      END IF
      call ga_sync()
      call ga_dgop(1, ECC2, 1, '+')
*
      ECC = ECC + ECC2
      CORE_TERMS = ECC
*
      IF(NTEST.NE.0.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
      END IF
*
      END
      FUNCTION ECORE_TERMSA(H,
     &         IOBSM,NSMOB,ITPFSO,IPHGAS,LOBSM,NORBT,ISTOB,NGAS,
     &         IDOH1, IDOH2, IDOPH1, IDOPH2)
*
      INCLUDE 'implicit.inc'
*
      WRITE(6,*) ' You have called ECORE_TERMSA to calculate '
      WRITE(6,*) ' core-energy term for sim transf H'
      WRITE(6,*) ' ECORE_TERMSA is a dummy routine so I stop '
      ECORE_TERMSA = -3006.56D0
      STOP ' ECORE_TERMSA is a dummy routine so I stop '
*
      END
      FUNCTION EXP_ONEEL_INACT(A,ISM)
*
* Calculate contributions to expectation value 
* from inactive orbitals for one-electron operator A
*
* ISM = 1 => A is packed to lower triangular symmetryblocks
* ISM = 0 => A is packed as complete symmetry blocks
*
* Jeppe Olsen, sitting in a 56 euro hotel room in Bruxelles
*              (pink doors, no table etc...)
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      NTEST = 100
*
      EXP_ONEEL_INACT = EXP_ONEEL_INACTS(A,ISM,NSMOB,NINOBS,NTOOBS)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Contribution from inactive orbitals to expectation value',
     &    EXP_ONEEL_INACT
      END IF
*
      RETURN
      END
      FUNCTION EXP_ONEEL_INACTS(A,ISM,NSMOB,NINOBS,NTOOBS)
*
* Calculate contributions to expectation value 
* from inactive orbitals for one-electron operator A
*
* ISM = 1 => A is packed to lower triangular symmetryblocks
* ISM = 0 => A is packed as complete symmetry blocks
*
* Jeppe Olsen, September 2010, sitting in a 56 euro hotel room in Bruxelles
*              (pink doors, no table etc...)
*
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*)
*. Generel input
      INTEGER NINOBS(NSMOB),NTOOBS(NSMOB)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from EXP_ONEEL_INACTS '
      END IF
*
      E = 0.0D0
*
      DO ISYM = 1,  NSMOB
        IF(ISYM.EQ.1) THEN
          IOFF = 1
        ELSE
          IF(ISM.EQ.1) THEN
            IOFF = IOFF + NTOOBS(ISYM-1)*(NTOOBS(ISYM-1)+1)/2
          ELSE
            IOFF = IOFF + NTOOBS(ISYM-1)**2
          END IF
        END IF
        WRITE(6,*) ' ISYM, IOFF =', ISYM, IOFF
        DO INAC = 1, NINOBS(ISYM)
          E = E + 2.0D0*A(IOFF -1 + INAC*(INAC+1)/2)
          WRITE(6,*) ' INAC, A(..) ', INAC, A(IOFF -1 + INAC*(INAC+1)/2)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Contribution from inactive orbitals to exp.value',
     &  E
      END IF
*
      EXP_ONEEL_INACTS = E
*
      RETURN
      END
      SUBROUTINE GET_2E_TERMS_TO_FI(F2,C,IHOLETP)
*
*. Obtain two-electron terms to inactive Fock matrix from C being
*  expansion of MO's in initial basis
*
*. Jeppe Olsen, October 2010
*
* The orbitals contributing to F2 is defined by IHOLETP:
*
* IHOLETP = 1: Only explicitly declared inactive orbitals
*         = 2: Only explicitly declared active hole-orbitals, 
*              which are not inactive orbitals
*           =3: Explicitly declared inactive + active
*               hole-orbitals
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  
*. Input
      DIMENSION C(*)
*. Output
CNW   DIMENSION F2(*)
      integer F2
*
      NTEST = 000
*
*. Assumed symmetry of density and integrals- for later 
*. generalizations
      IINTSM = 1
      IDENSM = 1
      IFSM = MULTD2H(IINTSM,IDENSM)
*
      LEN_F = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,0)
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      IDUM = -1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'2ETOFI')
*. Space for a density matrix- used also for scratch, therefore the
*  factor 2
      LENM = 2*NTOOB*NTOOB
      CALL MEMMAN(KLDINI,LENM,'ADDL  ',2,'DINI  ')
*. And a scratch matrix 
      CALL MEMMAN(KLSCR,LENM,'ADDL  ',2,'SCRM  ')
*. Obtain density matrix in initial basis
      CALL GET_D_INI_FROM_C(dbl_mb(KLDINI),C,IHOLETP)
*. Contract with two-electron integrals
      FACC = 1.0D0
      FACE = 0.5D0
C     TWO_INT_D_TERM_F(F2,DINI,FACC,FACE)
      CALL MEMCHK2('BE_TWO')
      CALL TWOINT_D_TERM_F(F2,dbl_mb(KLDINI),FACC,FACE)
      call ga_sync()
      CALL MEMCHK2('AF_TWO')
*. Transform from initial to current basis
C     TRAN_SYM_BLOC_MAT4
C    &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      CALL MEMCHK2('BE_TRA')
      CALL TRAN_SYM_BLOC_MAT4b(F2,C,C,NSMOB,NTOOBS,NTOOBS,
     &     dbl_mb(KLSCR),dbl_mb(KLDINI),0)
      call ga_sync()
      CALL MEMCHK2('AF_TRA')
CNW   CALL COPVEC(dbl_mb(KLSCR),F2,LEN_F)
      call ga_put(F2,1,LEN_F,1,1,dbl_mb(KLSCR),LEN_F)
      call ga_sync()
*
      IF(NTEST.GE.100) THEN
        if (ga_nodeid().eq.0) then
        WRITE(6,*) ' Output from GET_2E_TERMS '
        WRITE(6,*) ' ========================='
        endif
CNW     CALL APRBLM2(F2,NTOOBS,NTOOBS,NSMOB,0)
        call ga_print(F2)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'2ETOFI')
      RETURN
      END
      SUBROUTINE GET_D_INI_FROM_C(DINI,C,IHOLETP)
*
* Obtain Density matrix (in DINI) in initial basis over hole-orbitals
* from C-coefficients in C
*
* DINI(I,J) = 2*Sum_(k:hole) C(I,k) C(J,k)
*
* Hole-orbitals: IHOLETP = 1: Only explicitly declared inactive orbitals
*                        = 2: Only explicitly declared active hole-orbitals, 
*                             which are not hole-orbitals
*                         =3: Explicitly declared inactive + active
*                             hole-orbitals
*
*
      INCLUDE 'wrkspc.inc'
*. General input
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      DIMENSION C(*)
*. Specific output
      DIMENSION DINI(*)
* DINI will be calculated as complete matrix although it is a 
* symmetric matrix
*
* Jeppe Olsen, Summer of 2011
* Last Revision, Sept 4 2012, Jeppe Olsen, Correcting error for active hole-orbitals
      ZERO = 0.0D0
      ONE = 1.0D0
      TWO = 2.0D0
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_D_FROM_C' 
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' MO-INI transformation'
        CALL APRBLM2(C,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      DO  ISM = 1, NSMOB
        LORB_ISM = NTOOBS(ISM)
*
        IF(ISM.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + NTOOBS(ISM-1)**2
        END IF
        CALL SETVEC(DINI(ICOFF),ZERO,LORB_ISM**2)
*
        KOFF = -3006
        DO KTP = 0, NGAS
          IF(KTP.EQ.0) THEN
            KOFF = 1
          ELSE
            KOFF = KOFF + NOBPTS_GN(KTP-1,ISM)
          END IF
*. Should KTP be included?
          INCLUDE = 0
          IF((IHOLETP.EQ.1.OR.IHOLETP.EQ.3).AND.KTP.EQ.0) INCLUDE  = 1
          IF(IHOLETP.GE.2.AND.(1.LE.KTP.AND.KTP.LE.NGAS)) THEN
            IF(IPHGAS(KTP).EQ.2) INCLUDE = 1
          END IF
*
          IF(INCLUDE.EQ.1) THEN
*. start of orbitals relative to start of orbitals with given sym
COLD       KOFF = IOBPTS_GN(KTP,ISM) - IOBPTS_GN(0,ISM) + 1
           LORB_KTP = NOBPTS_GN(KTP,ISM)
COLD       WRITE(6,*) ' ISM, KTP, KOFF = ', ISM, KTP, KOFF
           DO KORB = KOFF, KOFF + LORB_KTP - 1
*. Update DINI(I,J) with C(I,K)C(J,K)
             DO J = 1,  LORB_ISM
               CJK = C(ICOFF-1+(KORB-1)*LORB_ISM + J)
C?             WRITE(6,*) 'ISM, J,KORB,CJK = ', ISM,J,KORB,CJK
               FAC = TWO*CJK
               CALL VECSUM(DINI(ICOFF + (J-1)*LORB_ISM),
     &                     DINI(ICOFF + (J-1)*LORB_ISM),
     &                     C(ICOFF + (KORB-1)*LORB_ISM),
     &                     ONE,FAC,LORB_ISM)
             END DO
           END DO
          END IF
*         ^ End of Ktp should be included
        END DO
*       ^ End of loop over KTP
      END DO
*     ^ End of loop over ISM
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DINI as delivered by GET_D_FROM_C' 
        WRITE(6,*) ' ================================='
        CALL APRBLM2(DINI,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE FI_FROM_INIINT(FI,CINI,H,EINAC,IHOLETP)
*
*. Calculate inactive Fock matrix from integrals in initial basis
*. and one-electron integrals in current basis
*
* IHOLETP defines which types of orbitals contribute to the 
* inactive Fock-matrix
*
* IHOLETP = 1: Explicitly declared inactive orbitals
* IHOLETP = 2: Orbitals that are declared hole-orbitals
*              through IPHGAS
* IHOLETP = 3: Combination of above two
*
*. Jeppe Olsen, October 2010
*  (Growing up in public?- Heading towards efficient treatment of 
*                inactive orbitals)
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cprnt.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  
*
*. Input: 
*. =====
* H: one-electron integrals in current MO basis, lower half of symmetry blocks
* IHOLETP: As above
* CINI: Current set of mo's expanded in initial set of orbitals
CNW   DIMENSION H(*), CINI(*)
      integer H
      dimension CINI(*)
*
*. Output:
*. ======
* FI:  Inactive Fock matrix in current basis, lower half of symmetry blocks
* EINAC: contribution to core energy from inactive/hole-orbitals
CNW   DIMENSION FI(*)
      integer FI
*
      IF(IPRINTEGRAL.GE.100) NTEST = 100
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from FI_FROM_INIINT'
        WRITE(6,*) ' ========================'
        WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN 
        WRITE(6,*) ' Input CINI '
        CALL APRBLM2(CINI,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'FIFRIN')
*. For later generalizations
      IINTSM = 1
      IDENSM = 1
      IFSM = MULTD2H(IINTSM,IDENSM)
*
      LEN_F = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,0)
      LEN_FP = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,1)
      CALL MEMMAN(KLBLF,LEN_FP,'ADDL  ',2,'LBLF  ')
*. Two-electron terms to the inactive Fock-matrix
C     GET_2E_TERMS_TO_FI(F2,C_MOINI,IHOLETP)
      CALL GET_2E_TERMS_TO_FI(FI,CINI,IHOLETP)
      call ga_sync()
*. Pack to lower half form
      CALL TRIPAK_BLKMG(FI,dbl_mb(KLBLF),1,NTOOBS,NSMOB,2)
      call ga_sync()
*. Contribution to inactive energy from two-electron interaction
       CALL GET_INA_TERM_TO_1EEXP(dbl_mb(KLBLF),E_2E,IHOLETP,1)
       E_2E = 0.5D0*E_2E
*. Contribution to inactive energy from one-electron interaction
       CALL GET_INA_TERM_TO_1EEXPb(H,E_1E,IHOLETP,1)
       EINAC = E_2E + E_1E
*. And add to one-electron Hamiltonian
       ONE = 1.0D0
CNW    CALL VECSUM(FI,H,dbl_mb(KLBLF),ONE,ONE,LEN_FP)
       call ga_put(FI,1,LEN_FP,1,1,dbl_mb(KLBLF),LEN_FP)
       call ga_sync()
       call ga_add(1.0d0,H,1.0d0,FI,FI)
*
       IF(NTEST.GE.100) THEN
         if (ga_nodeid().eq.0) then
         WRITE(6,*) ' Output from FI_FROM_INIINT: '
         WRITE(6,*) ' One- and two-electron contributions to E_INAC=',
     &               E_1E,E_2E
         WRITE(6,*) ' The inactive Fock matrix'
         endif
CNW      CALL APRBLM2(FI,NTOOBS,NTOOBS,NSMOB,1)
         call ga_print(FI)
       END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'FIFRIN')
      RETURN
      END
*
      SUBROUTINE TWOINT_D_TERM_F(F2,DINI,FACC,FACE)
*
* Obtain terms to a Fock-matrix from two-electron integrals times 
* one-body density matrix in initial basis. 
* Input density  is assumed in complete symmetryblocked form
* and complete symmetry-blocked matrix generated
*
* F2(IAL,IBE) = SUM(IGA,IDE) (FACC*(IAL,IBE!IGA, IDE)-FACE*(IAL,IDE!IGA,IBE))
*                             * D(IGA,IDE)
*
*. Jeppe Olsen, October 2010
*
*. General input
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  	  
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
*. Specific input
      DIMENSION DINI(*)
*. Local pointer to offsets in DINI
      INTEGER IDPNT(MXPOBS)
      INTEGER IDSMOS(MXPOBS)
*. Output
CNW   DIMENSION F2(*)
      integer F2
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'F2EINI')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
         WRITE(6,*) ' Entering TWOINT_D_TERM_F '
         WRITE(6,*) ' ========================= '
         WRITE(6,*) ' FACC, FACE = ',FACC, FACE
      END IF
*
*. Symmetry for integrals and densities- for later generalizations
      IINTSM = 1
      IDENSM = 1
      IFSM = MULTD2H(IINTSM,IDENSM)
*. Offsets to symmetryblocks in density matrix
C     PNT2DM(I12SM,NSMOB,NSMSX,OSXO,IPSM,JPSM,
C    &                  IJSM,ISM2,IPNTR,MXPOBS)
      CALL PNT2DM(0,NSMOB,NSMOB,MULTD2H,NTOOBS,NTOOBS,
     &            IDENSM,IDSMOS,IDPNT,MXPOBS)
*. Length of largest 2e-int symmetry-block- stored in core!
      LBLOCK = 0
      LBLKMX = 0
      DO ISM = 1, NSMOB
       LBLKMX = MAX(LBLKMX,NTOOBS(ISM))
       DO JSM = 1, NSMOB
        IJSM = MULTD2H(ISM,JSM)
        DO KSM = 1, NSMOB
         IJKSM = MULTD2H(IJSM,KSM)
         LSM = MULTD2H(IINTSM,IJKSM)
         LENGTH = NTOOBS(ISM)*NTOOBS(JSM)*NTOOBS(KSM)*NTOOBS(LSM)
         LBLOCK = MAX(LBLOCK,LENGTH)
         IF(NTEST.GE.1000) THEN
           WRITE(6,'(A,4I3,I9)')' ISM,JSM,KSM,LSM,LENGTH= ',
     &                          ISM,JSM,KSM,LSM,LENGTH
         END IF
        END DO
       END DO
      END DO
*
      CALL MEMMAN(KL2EBLK,LBLOCK,'ADDL  ',2,'E2BLK ')
      CALL MEMMAN(KFI2BLK,LBLKMX*LBLKMX,'ADDL  ',2,'FI2BLK ')
      call ga_zero(F2)
*
* FI2(IAL,IBE) = SUM(IGA,IDE) ((IAL,IBE!IGA, IDE)-0.5(IAL,IDE!IGA,IBE))
*                             * D(IGA,IDE)
      NALOB = -2810
      NBEOB = -2810
      IFOFF = -2810
      DO IALSM  = 1, NSMOB
       IBESM = MULTD2H(IFSM,IALSM)
       IABSM = IFSM
       IF(IALSM.EQ.1) THEN
         IFOFF = 1
       ELSE
         IFOFF = IFOFF + NALOB*NBEOB
       END IF
* 
       NALOB = NTOOBS(IALSM)
       NBEOB = NTOOBS(IBESM)
       ZERO = 0.0D0
CNW    CALL SETVEC(F2(IFOFF),ZERO,NALOB*NBEOB)
       CALL SETVEC(dbl_mb(KFI2BLK),ZERO,NALOB*NBEOB)
*
       DO IGASM = 1, NSMOB
        IABGSM = MULTD2H(IABSM,IGASM)
        IDESM = MULTD2H(IINTSM,IABGSM)
*
        NGAOB = NTOOBS(IGASM)
        NDEOBL = NTOOBS(IDESM)
* (NDEOBL instead of NDEOB as NDEOB is total numner of deleted orbitals
* in orbinp)
* 
        IDOFF = IDPNT(IGASM)
*. Fetch integral block
        ICOUL = 1
        IXCHNG = 1
        IKSM = 0
        JLSM = 0
*
        IDOFUSK = 0
        IFUSK = 0
        IF(IALSM.EQ.2.AND.IBESM.EQ.2.AND.IGASM.EQ.1.AND.IDESM.EQ.1) THEN
           IFUSK = 1
        END IF
        IF(IDOFUSK.EQ.1.AND.IFUSK.EQ.1) THEN
          FACC = 0.0D0
          FACE = 1.0D0
        END IF
        CALL LUCIAGETINT(dbl_mb(KL2EBLK),-1,IALSM,-1,IBESM,-1,IGASM,-1,
     &              IDESM,
     &              IXCHNG,IKSM,JLSM,ICOUL,FACC,FACE)
        IF(NTEST.GE.1000.AND.IALSM.EQ.1.AND.IBESM.EQ.1.AND.
     &    IGASM.EQ.1.AND.IDESM.EQ.1) THEN
          WRITE(6,*) ' IALSM, IBESM, IGASM, IDESM = ',
     &                 IALSM, IBESM, IGASM, IDESM 
          WRITE(6,*) ' Integrals from GETINT as (AL BE, GA DE)'
          CALL WRTMAT(dbl_mb(KL2EBLK),NALOB*NBEOB,NGAOB*NDEOBL,
     &                              NALOB*NBEOB,NGAOB*NDEOBL )
        END IF
C       GETINT_ORIG(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                  IXCHNG,IKSM,JLSM,ICOUL)
*. And multiply
        FACC = 1.0D0
        FACAB = 1.0D0
C?      WRITE(6,*) ' IFOFF, IDOFF,KL2EBLK = ', IFOFF, IDOFF,KL2EBLK
CNW     CALL MATML7(F2(IFOFF),dbl_mb(KL2EBLK),DINI(IDOFF),
        CALL MATML7(dbl_mb(KFI2BLK),dbl_mb(KL2EBLK),DINI(IDOFF),
     &              NALOB*NBEOB,1,NALOB*NBEOB,NGAOB*NDEOBL,
     &              NGAOB*NDEOBL,1,FACC,FACAB,0)
C     MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
       END DO
       if (NALOB*NBEOB.gt.0)
     &  call ga_put(F2,IFOFF,IFOFF+NALOB*NBEOB-1,1,1,dbl_mb(KFI2BLK),
     &              NALOB*NBEOB)
      END DO
      call ga_sync()
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from TWOINT_D_TERM_F'
       WRITE(6,*) ' ============================'
       WRITE(6,*) ' Contraction of 2e ints and 1e density'
CNW    CALL APRBLM2(F2,NTOOBS,NTOOBS,NSMOB,0)
       call ga_print(F2)
      END IF
*
      IF(NTEST.GE.100) WRITE(6,*) ' Leaving TWOINT_D_TERM_F '
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'F2EINI')
      RETURN
      END
      SUBROUTINE GET_INA_TERM_TO_1EEXP(F,EXPEC_INA,IHOLETP,ISYM)
* 
* Contribution from inactive terms in matrix F to one-electron
* expectation value
*
*  EXPEC = 2*sum_k F_kk, where the sum over k is over inactive/hole
*          orbitals according to IHOLETP
* F is stored in standard integral mode (symmetryblocked)
*
*. Jeppe Olsen, October 2010
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
*. Specific input
      DIMENSION F(*)
*
      NTEST = 0
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================================='
        WRITE(6,*) ' Output from  GET_INA_TERM_TO_1EEXP '
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
      END IF
*
      EXPEC_INA = 0.0D0
      IJSM = 1
*. One-electron part 
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE
          IF(ISYM.EQ.1) THEN
            IIOFF = IIOFF + NTOOBS(ISM-1)*(NTOOBS(ISM-1)+1)/2
          ELSE
            IIOFF = IIOFF + NTOOBS(ISM-1)**2
          END IF
        END IF
        II = IIOFF-1
        DO I = IBSO(ISM),IBSO(ISM)+NTOOBS(ISM)-1
          IF(ISYM.EQ.1) THEN
            II = II + (I-IBSO(ISM)+1)
           ELSE
            II = IIOFF -1 + (I-IBSO(ISM))*NTOOBS(ISM) + (I-IBSO(ISM)+1)
           END IF
          I_INACTIVE = 0
          IF(ITPFSO(I).GT.0.AND.ITPFSO(I).LE.NGAS) THEN
           IF (IHOLETP.GE.2.AND.IPHGAS(ITPFSO(I)).EQ.2) I_INACTIVE = 1
          END IF
          IF((IHOLETP.EQ.1.OR.IHOLETP.EQ.3).
     &       AND.ITPFSO(I).EQ.0) I_INACTIVE = 1
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' I, I_INACTIVE,II = ', I,I_INACTIVE,II
          IF(I_INACTIVE.EQ.1) THEN
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' Contribution to ECC from I =', I
            EXPEC_INA = EXPEC_INA + 2*F(II)
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' Updated EXPEC_INAC =', EXPEC_INA
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from GET_INA_TERM_TO_1EEXP '
        WRITE(6,*)
        WRITE(6,*) ' EXPEC_INA = ', EXPEC_INA
      END IF
*
      RETURN
      END

      SUBROUTINE GET_INA_TERM_TO_1EEXPb(F,EXPEC_INA,IHOLETP,ISYM)
* 
* Contribution from inactive terms in matrix F to one-electron
* expectation value
*
*  EXPEC = 2*sum_k F_kk, where the sum over k is over inactive/hole
*          orbitals according to IHOLETP
* F is stored in standard integral mode (symmetryblocked)
*
*. Jeppe Olsen, October 2010
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
CNW   DIMENSION F(*)
      integer F
*
      NTEST = 0
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================================='
        WRITE(6,*) ' Output from  GET_INA_TERM_TO_1EEXP '
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
      END IF
*
      EXPEC_INA = 0.0D0
      IJSM = 1
*. One-electron part 
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IF(ISYM.EQ.1) THEN
            IIOFF = IIOFF + NTOOBS(ISM-1)*(NTOOBS(ISM-1)+1)/2
          ELSE 
            IIOFF = IIOFF + NTOOBS(ISM-1)**2
          END IF
        END IF
        II = IIOFF-1
        DO I = IBSO(ISM),IBSO(ISM)+NTOOBS(ISM)-1
          IF(ISYM.EQ.1) THEN
            II = II + (I-IBSO(ISM)+1)
           ELSE
            II = IIOFF -1 + (I-IBSO(ISM))*NTOOBS(ISM) + (I-IBSO(ISM)+1)
           END IF
          I_INACTIVE = 0
          IF(ITPFSO(I).GT.0.AND.ITPFSO(I).LE.NGAS) THEN
           IF (IHOLETP.GE.2.AND.IPHGAS(ITPFSO(I)).EQ.2) I_INACTIVE = 1
          END IF
          IF((IHOLETP.EQ.1.OR.IHOLETP.EQ.3).
     &       AND.ITPFSO(I).EQ.0) I_INACTIVE = 1
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' I, I_INACTIVE,II = ', I,I_INACTIVE,II
          IF(I_INACTIVE.EQ.1) THEN
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' Contribution to ECC from I =', I
            call ga_get(F,II,II,1,1,FVAL,1)
            EXPEC_INA = EXPEC_INA + 2*FVAL
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' Updated EXPEC_INAC =', EXPEC_INA
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from GET_INA_TERM_TO_1EEXP '
        WRITE(6,*) 
        WRITE(6,*) ' EXPEC_INA = ', EXPEC_INA
      END IF
*
      RETURN
      END
      
      SUBROUTINE FA_FROM_INIINT(FA,CINI,CINIB,D,IPACK)
*
*. Calculate normal or bioorthogonal active Fock matrix from integrals in initial basis
*. and density in current basis
*
*. Jeppe Olsen, October 2010
*               Bioorthogonal extension added 
*  (Growing up in public?- Heading towards efficient treatment of 
*                inactive orbitals)
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cprnt.inc'
*
*. Input: 
*. =====
* CINI: Current set of mo's expanded in initial set of orbitals
* CINIB: Current set of bioorthogonal MO's (in standard = CINI)
* D   : Density matrix in current basis
      DIMENSION D(*), CINI(*), CINIB(*)
*
*. Output:
*. ======
* FA:  Inactive Fock matrix in current basis, lower half of symmetry blocks
CNW   DIMENSION FA(*)
      integer FA
*. Local scratch
      DIMENSION NOBPTS_L(0:6+MXPR4T,MXPOBS)
*
      NTEST = 000
      IF(IPRINTEGRAL.GE.100) NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info form FA_FROM_INI'
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' Initial density matrix '
        CALL APRBLM2(D,NACOBS,NACOBS,NSMOB,0)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'FAFRIN')
*. For later generalizations
      IINTSM = 1
      IDENSM = 1
      IFSM = MULTD2H(IINTSM,IDENSM)
*
      LEN_F = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,0)
      LEN_FP = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,1)
      CALL MEMMAN(KLCT,LEN_F,'ADDL  ',2,'CT    ')
      CALL MEMMAN(KLCBT,LEN_F,'ADDL  ',2,'CBT   ')
      CALL MEMMAN(KLDBLK,LEN_F,'ADDL  ',2,'DBLK  ')
      CALL MEMMAN(KLDINI,LEN_F,'ADDL  ',2,'DINI  ')
      CALL MEMMAN(KLSCR,2*LEN_F,'ADDL  ',2,'LSCR  ')
*. Obtain one-electron density matrix in blocked form
C     EXTR_SYMBLK_ACTMAT(AIN,AOUT,IJSM)
CERR  CALL EXTR_SYMBLK_ACTMAT(D,WORK(KLDBLK),1)
C          REORHO1(RHO1I,RHO1O,IRHO1SM,IWAY)
      CALL REORHO1(D,dbl_mb(KLDBLK),1,1)
*
*. Obtain density matrix contravariantly transformed to initial basis
*.  DINI = C_act D CB_act(T)
*
*. Extract active MO's from C and CB and save in KLSCR
*. (NACOBS is overwritten- but unchanged)
      CALL CSUB_FROM_C(CINIB,dbl_mb(KLSCR),NACOBS,NOBPTS_L,1,NGAS,0)
      CALL TRP_BLK_MAT(dbl_mb(KLSCR),dbl_mb(KLCBT),NSMOB,NTOOBS,NACOBS)
      CALL CSUB_FROM_C(CINI,dbl_mb(KLSCR),NACOBS,NOBPTS_L,1,NGAS,0)
      CALL TRP_BLK_MAT(dbl_mb(KLSCR),dbl_mb(KLCT),NSMOB,NTOOBS,NACOBS)
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,IONLY_DIM)
      CALL TRAN_SYM_BLOC_MAT4a(dbl_mb(KLDBLK),dbl_mb(KLCT),
     &     dbl_mb(KLCBT),NSMOB,NACOBS,NTOOBS,dbl_mb(KLDINI),
     &     dbl_mb(KLSCR),0)
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' Density matrix over active orbitals, initial basis'
       CALL APRBLM2(dbl_mb(KLDINI),NTOOBS,NTOOBS,NSMOB,0)
      END IF
C     TRAN_SYM_BLOC_MAT4
C    &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
*. Two-electron terms to the active Fock-matrix
      FACC = 1.0D0
      FACE = 0.5D0
      CALL TWOINT_D_TERM_F(FA,dbl_mb(KLDINI),FACC,FACE)
*. Transform to current basis- and save result in KLCT
      CALL TRAN_SYM_BLOC_MAT4b(FA,CINI,CINIB,NSMOB,NTOOBS,NTOOBS,
     &                        dbl_mb(KLCT),dbl_mb(KLSCR),0)
      call ga_sync()
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' FA in initial basis'
        CALL APRBLM2(dbl_mb(KLCT),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. Pack to lower half form
      IF(IPACK.NE.0) THEN
        CALL TRIPAK_BLKMG2(dbl_mb(KLCT),FA,1,NTOOBS,NSMOB,3)
      ELSE
CNW     CALL COPVEC(dbl_mb(KLCT),FA,LEN_F)
        call ga_put(FA,1,LEN_F,1,1,dbl_mb(KLCT),LEN_F)
      END IF
*
       IF(NTEST.GE.100) THEN
         if (ga_nodeid().eq.0) WRITE(6,*) ' The active Fock matrix'
CNW      CALL APRBLM2(FA,NTOOBS,NTOOBS,NSMOB,IPACK)
         call ga_print(FA)
       END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'FAFRIN')
      RETURN
      END
	  
      SUBROUTINE FI_FROM_INIINT_G(FI,CINI,CINIB,H,EINAC,IHOLETP,IPACK)
*
* A MOAO transformation CINI and its biorthogonal tranformation CINIB
* are given (CINIB = CINI gives standard)
*
*. Calculate standard (CINIB = CINI) or biothogonal inactive Fock 
*. matrix from integrals in initial basis
*. and one-electron integrals in current (standard or bioorthogonal)
* .basis
*
* IHOLETP defines which types of orbitals contribute to the 
* inactive Fock-matrix
*
* IHOLETP = 1: Explicitly declared inactive orbitals
* IHOLETP = 2: Orbitals that are declared hole-orbitals
*              through IPHGAS
* IHOLETP = 3: Combination of above two
*
* IPACK = 1, input and output one-electron integrals 
*             are packed
*       = 0: input and output one-electron integrals 
*            are not packed

*
*. Jeppe Olsen, July 2011

      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*
*. Input: 
*. =====
* H: one-electron integrals in current MO basis, symmetry blocks
*    packed according to IPACK
* IHOLETP: As above
* CINI: Current set of mo's expanded in initial set of orbitals
* CINIB: Current set of bioorthogonal mo's expanded in initial set of orbitals
CNW   DIMENSION H(*), CINI(*), CINIB(*)
      DIMENSION CINI(*), CINIB(*)
      integer H
*
*. Output:
*. ======
* FI:  Inactive Fock matrix in current basis, packed according to IPACK
*      
* EINAC: contribution to core energy from inactive/hole-orbitals
CNW   DIMENSION FI(*)
      integer FI
*
      NTEST = 000
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'FIFRIN')
*. For later generalizations
      IINTSM = 1
      IDENSM = 1
      IFSM = MULTD2H(IINTSM,IDENSM)
*
      LEN_F = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,0)
      LEN_FP = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,1)
      CALL MEMMAN(KLBLF,LEN_F,'ADDL  ',2,'LBLF  ')
*. Two-electron terms to the inactive Fock-matrix
C     GET_2E_TERMS_TO_FI_G(F2,C_MOINI,C_MOINIB,IHOLETP)
      CALL GET_2E_TERMS_TO_FI_G(FI,CINI,CINIB,IHOLETP)
      call ga_sync()
*. Pack to lower half form is required
      IF(IPACK.EQ.1) THEN
        CALL TRIPAK_BLKMG(FI,dbl_mb(KLBLF),1,NTOOBS,NSMOB,2)
      ELSE
CNW     CALL COPVEC(FI,dbl_mb(KLBLF),LEN_F)
        call ga_get(FI,1,LEN_F,1,1,dbl_mb(KLBLF),1)
       END IF
       call ga_sync()
*. Contribution to inactive energy from two-electron interaction
C      GET_INA_TERM_TO_1EEXP(F,EXPEC_INA,IHOLETP,ISYM)
       CALL GET_INA_TERM_TO_1EEXP(dbl_mb(KLBLF),E_2E,IHOLETP,IPACK)
       E_2E = 0.5D0*E_2E
*. Contribution to inactive energy from one-electron interaction
       CALL GET_INA_TERM_TO_1EEXPb(H,E_1E,IHOLETP,IPACK)
       EINAC = E_2E + E_1E
*. And add to one-electron Hamiltonian
       ONE = 1.0D0
       IF(IPACK.EQ.0) THEN
         LEN_MAT = LEN_F
       ELSE
         LEN_MAT = LEN_FP
       END IF
CNW    CALL VECSUM(FI,H,dbl_mb(KLBLF),ONE,ONE,LEN_MAT)
       if (IPACK.EQ.1) call ga_put(FI,1,LEN_MAT,1,1,dbl_mb(KLBLF),1)
       call ga_sync()
       call ga_add(ONE,H,ONE,FI,FI)
*
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Output from FI_FROM_INIINT_G: '
         WRITE(6,*) ' One- and two-electron contributions to E_INAC=',
     &               E_1E,E_2E
         WRITE(6,*) 
     &   ' Contribution from inactive orbitals to core-energy', EINAC
         WRITE(6,*) ' The inactive Fock matrix'
CNW      CALL APRBLM2(FI,NTOOBS,NTOOBS,NSMOB,IPACK)
CBERT    call ga_print(fi)
       END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'FIFRIN')
      RETURN
      END
	  
	  
      SUBROUTINE GET_2E_TERMS_TO_FI_G(F2,C,CB,IHOLETP)
*
*. Obtain standard or bioorthogonal
*  two-electron terms to inactive Fock matrix from C being
*  expansion of MO's in initial basis and CB being the 
*  biorthogonal transformation. Standard orthogonal expression 
*  is obtained if CB = C
*
* Output matrix is with complete symmetry blocks
*
*. Jeppe Olsen, July 2011
*
* The orbitals contributing to F2 is defined by IHOLETP:
*
* IHOLETP = 1: Only explicitly declared inactive orbitals
*         = 2: Only explicitly declared active hole-orbitals, 
*              which are not inactive orbitals
*           =3: Explicitly declared inactive + active
*               hole-orbitals
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  
*. Input
      DIMENSION C(*), CB(*)
*. Output
CNW   DIMENSION F2(*)
      integer F2
*
      NTEST = 000
*
*. Assumed symmetry of density and integrals- for later 
*. generalizations
      IINTSM = 1
      IDENSM = 1
      IFSM = MULTD2H(IINTSM,IDENSM)
*
      LEN_F = NDIM_1EL_MAT(IFSM,NTOOBS,NTOOBS,NSMOB,0)
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      IDUM = -1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'2ETOFI')
*. Space for a density matrix- used also for scratch, therefore the
*  factor 2
      LENM = 2*NTOOB*NTOOB
      CALL MEMMAN(KLDINI,LENM,'ADDL  ',2,'DINI  ')
*. And a scratch matrix 
      CALL MEMMAN(KLSCR,LENM,'ADDL  ',2,'SCRM  ')
*. Obtain density matrix in initial basis
      CALL GET_D_INI_FROM_C_G(dbl_mb(KLDINI),C,CB,IHOLETP)
*. Contract with two-electron integrals
      FACC = 1.0D0
      FACE = 0.5D0
C     TWO_INT_D_TERM_F(F2,DINI,FACC,FACE)
      CALL MEMCHK2('BE_TWO')
      CALL TWOINT_D_TERM_F(F2,dbl_mb(KLDINI),FACC,FACE)
      call ga_sync()
      CALL MEMCHK2('AF_TWO')
*. Transform from initial to current basis
C     TRAN_SYM_BLOC_MAT4
C    &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      CALL MEMCHK2('BE_TRA')
      CALL TRAN_SYM_BLOC_MAT4b(F2,C,CB,NSMOB,NTOOBS,NTOOBS,
     &     dbl_mb(KLSCR),dbl_mb(KLDINI),0)
      call ga_sync()
      CALL MEMCHK2('AF_TRA')
CNW   CALL COPVEC(dbl_mb(KLSCR),F2,LEN_F)
      call ga_put(F2,1,LEN_F,1,1,dbl_mb(KLSCR),1)
      call ga_sync()
*
      IF(NTEST.GE.100) THEN
        if (ga_nodeid().eq.0) then
        WRITE(6,*) 
     &  ' Two-electron terms to inactive Fock matrix in MO basis'
        WRITE(6,*) 
     & ' ======================================================'
        endif
CNW     CALL APRBLM2(F2,NTOOBS,NTOOBS,NSMOB,0)
        call ga_print(f2)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'2ETOFI')
      RETURN
      END
      SUBROUTINE GET_D_INI_FROM_C_G(DINI,C,CB,IHOLETP)
*
* Obtain Density matrix (in DINI) in initial basis over hole-orbitals
* from C-coefficients in C and the bioorthogonal transformation 
* matrix CB (standard is recovered for C = CB)
*
* DINI(I,J) = 2*Sum_(k:hole) C(I,k) CB(J,k)
*
* Hole-orbitals: IHOLETP = 1: Only explicitly declared inactive orbitals
*                        = 2: Only explicitly declared active hole-orbitals, 
*                             which are not hole-orbitals
*                         =3: Explicitly declared inactive + active
*                             hole-orbitals
* Jeppe Olsen, July 2011
*
*
      INCLUDE 'wrkspc.inc'
*. General input
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      DIMENSION C(*), CB(*)
*. Specific output
      DIMENSION DINI(*)
* DINI will be calculated as complete matrix although it is a 
* symmetric matrix
*
      ZERO = 0.0D0
      ONE = 1.0D0
      TWO = 2.0D0
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_D_FROM_C' 
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' MO-INI transformation'
        CALL APRBLM2(C,NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*) ' MO-INI bioorthogonal transformation'
        CALL APRBLM2(CB,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      DO  ISM = 1, NSMOB
        LORB_ISM = NTOOBS(ISM)
*
        IF(ISM.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + NTOOBS(ISM-1)**2
        END IF
        CALL SETVEC(DINI(ICOFF),ZERO,LORB_ISM**2)
*
        DO KTP = 0, NGAS
*. Should KTP be included?
          INCLUDE = 0
          IF((IHOLETP.EQ.1.OR.IHOLETP.EQ.3).AND.KTP.EQ.0) INCLUDE  = 1
          IF(IHOLETP.GE.2.AND.(1.LE.KTP.AND.KTP.LE.NGAS)) THEN
            IF(IPHGAS(KTP).EQ.2) INCLUDE = 1
          END IF
*
          IF(INCLUDE.EQ.1) THEN
*. start of orbitals relative to start of orbitals with given sym
           KOFF = IOBPTS_GN(KTP,ISM) - IOBPTS_GN(0,ISM) + 1
           LORB_KTP = NOBPTS_GN(KTP,ISM)
           DO KORB = KOFF, KOFF + LORB_KTP - 1
*. Update DINI(I,J) with C(I,K)CB(J,K)
             DO J = 1,  LORB_ISM
               CBJK = C(ICOFF-1+(KORB-1)*LORB_ISM + J)
               FAC = TWO*CBJK
               CALL VECSUM(DINI(ICOFF + (J-1)*LORB_ISM),
     &                     DINI(ICOFF + (J-1)*LORB_ISM),
     &                     C(ICOFF + (KORB-1)*LORB_ISM),
     &                     ONE,FAC,LORB_ISM)
             END DO
           END DO
          END IF
*         ^ End of Ktp should be included
        END DO
*       ^ End of loop over KTP
      END DO
*     ^ End of loop over ISM
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DINI as delivered by GET_D_FROM_C_G' 
        WRITE(6,*) ' ==================================='
        CALL APRBLM2(DINI,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE ORTH_TO_SELECT_VECTOR(X,IORTVC,S,NDIM,SCR1)
*
* Orthogonalize vectors in X to vector given by X(I,IORTVC) and normalize
* 
* Jeppe Olsen, June 2012
*
* S is input overlap matrix in full format, X is input set of vectors
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. input 
      DIMENSION S(NDIM,NDIM)
      DIMENSION X(NDIM,NDIM)
*. Scratch : vector of length NDIM
      DIMENSION SCR1(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from ORTH_TO_SELECT_VECTOR '
       WRITE(6,*) ' ============================== '
       WRITE(6,*) ' Vector orthogonalizing to ', IORTVC
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input X and S block '
        CALL WRTMAT_F7(X,NDIM,NDIM,NDIM,NDIM)
        CALL WRTMAT_F7(S,NDIM,NDIM,NDIM,NDIM)
      END IF
*
* Normalize vector IORTVC
*
      CALL MATVCB(S,X(1,IORTVC),SCR1,NDIM,NDIM,0)
      XNORM = INPROD(X(1,IORTVC),SCR1,NDIM)
      IF(XNORM.LE.0.0D0) THEN
        WRITE(6,*) ' Input vector to ORTH_TO_SELECT_VEC.. is vanishing'
        STOP       ' Input vector to ORTH_TO_SELECT_VEC.. is vanishing'
      END IF
      FACTOR = 1.0D0/SQRT(XNORM)
      CALL SCALVE(X(1,IORTVC),FACTOR,NDIM)
*. Scale also S X(1,IORTVC)
      CALL SCALVE(SCR1,FACTOR,NDIM)
*
*. Orthogonalize other vectors to IORTVC
*
      ONE = 1.0D0
      DO IVEC = 1, NDIM
       IF(IVEC.NE.IORTVC) THEN
          OVERLAP = INPROD(SCR1,X(1,IVEC),NDIM)
          CALL VECSUM(X(1,IVEC),X(1,IVEC),X(1,IORTVC),ONE,-OVERLAP,NDIM) 
       END IF
      END DO
*
*. Normalize the vectors
*
      DO IVEC = 1, NVEC
        CALL MATVCB(S,X(1,IVEC),SCR1,NDIM,NDIM,0)
        XNORM = INPROD(X(1,IVEC),SCR1,NDIM)
        IF(XNORM.LE.0.0D0) THEN
          WRITE(6,*) ' Input vector to ORTH_TO_SELECT_VEC. is vanishing'
          STOP       ' Input vector to ORTH_TO_SELECT_VEC. is vanishing'
        END IF
        FACTOR = 1.0D0/SQRT(XNORM)
        CALL SCALVE(X(1,IVEC),FACTOR,NDIM)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap matrix '
        CALL WRTMAT_F7(S,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Orthonormalized vectors '
        CALL WRTMAT_F7(X,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE ORT_MOS_TO_SELECTED_MOS(CMOAO,NORT,IORT)
*
* Orthogonalize MOs in CMO to the orbitals in CMO given by
* IORT. 
* 
* IORT is given in symmetry-type order
*
*. Jeppe Olsen, June 2012
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
*. Input
      INTEGER*4 IORT(NORT)
*. Input and output
      DIMENSION CMOAO(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ORTMOS')
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from ORT_MOS_TO_SELE... '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' MOs to be orthogonalized to '
        CALL IWRTMA(IORT,1,NORT,1,NORT)
      END IF
*
*. Obtain overlap matrix SAO in complete block form
*
      LEN = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
C                     (IHSM,NRPSM,NCPSM,NSM,IPACK)
      CALL MEMMAN(KLSAOE,LEN,'ADDL  ',2,'SAOE  ')
      CALL TRIPAK_BLKM(dbl_mb(KLSAOE),dbl_mb(KSAO),2,NTOOBS,NSMOB,3)
*
*. Scratch vector for MO of given sym
      NMAX = IMNMX(NTOOBS,NSMOB,2)
      CALL MEMMAN(KLMO,NMAX,'ADDL  ',2,'MO_VEC')
*
      DO IIIORT = 1, NORT
       IIORT = IORT(IIIORT)
       ISYM = ISMFSO(IIORT)
       IREL = IIORT-IBSO(ISYM) + 1
*
       IOFF = 1
       DO JSYM = 1, ISYM-1
        IOFF = IOFF + NTOOBS(JSYM)**2
       END DO
       N = NTOOBS(ISYM)
C      ORTH_TO_SELECT_VECTOR(X,IORTVC,S,NDIM,SCR1)
       CALL ORTH_TO_SELECT_VECTOR(CMOAO(IOFF),IREL,
     &      dbl_mb(KLSAOE+IOFF-1),N,dbl_mb(KLMO))
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output vector from ORT_MOS_TO_SELE.. '
        CALL APRBLM_F7(CMOAO,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ORTMOS')
*
      RETURN
      END
      SUBROUTINE ORT_CMO_TO_FROZEN_ORBITALS(CMOAO)
*
* Orthogonalize orbitals to the frozen orbitals
*
*. Jeppe Olsen, June 2012
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"	  
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
*
*. Specific input and output
      DIMENSION CMOAO(*)
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from ORT_CMO_TO_FROZEN_ORBITALS '
        WRITE(6,*) ' ====================================='
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ORTFRO')
*. The frozen orbitals are given in type-order in IFRZ_ORB. Reform to 
*.symmetry-order
      CALL MEMMAN(KLFRZS,NFRZ_ORB,'ADDL  ',1,'FRZ_SM')
*
      DO IFRZ_T = 1, NFRZ_ORB
        IFRZ_S = IREOTS(IFRZ_ORB(IFRZ_T))
C            ICOPVE3(IIN,IOFFIN,IOUT,IOFFOUT,NDIM)
        CALL ICOPVE3(IFRZ_S,1,int_mb(KLFRZS),IFRZ_T,1)
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Frozen orbitals in symmetry numbering '
        CALL IWRTMA(int_mb(KLFRZS),1,NFRZ_ORB,1,NFRZ_ORB)
      END IF
*. And orthonormalize
C     ORT_MOS_TO_SELECTED_MOS(CMOAO,NORT,IORT)
      CALL ORT_MOS_TO_SELECTED_MOS(CMOAO,NFRZ_ORB,int_mb(KLFRZS))
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Orbitals orthogonalized to frozen '
        CALL APRBLM_F7(CMOAO,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ORTFRO')
      RETURN
      END
      SUBROUTINE MGS4(X,S,NDIM,SCR1,THRES,NVECUT)
*
* Modified Gram-Schmidt procedure by forward orthogonalization
*
*  watch out for zero columns indicating linear dependency
*
* Jeppe Olsen, March 2013, added thres to MGS3
*
* S is input overlap matrix, X is output set of orthonormalized vectors
*
* Thres is min norm of linear independent vector- only meaningfull if all 
* initial vectors have identical norm
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. input 
      DIMENSION S(NDIM,NDIM)
*. Output
      DIMENSION X(NDIM,*)
*. Scratch : vector of length NDIM
      DIMENSION SCR1(*)
*
      NTEST = 10
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from MGS4 '
        WRITE(6,*) ' ==============='
      END IF
*
*. Initialize X to unit matrix
*
      ZERO = 0.0D0
      ONE = 1.0D0
      CALL SETVEC(X,ZERO,NDIM**2)
      CALL SETDIA(X,ONE,NDIM,0)     
C          SETDIA(MATRIX,VALUE,NDIM,IPACK)
*
      DO IVEC = 1, NDIM
*. Normalize vector IVEC
        CALL MATVCB(S,X(1,IVEC),SCR1,NDIM,NDIM,0)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
*. avoid NaN's by putting norm to at least zero
        XNORM = INPROD(X(1,IVEC),SCR1,NDIM)
        
        IF (XNORM.LE.THRES) THEN
          FACTOR = 0.0D0
        ELSE
          FACTOR = 1.0D0/SQRT(XNORM)
        END IF
        CALL SCALVE(X(1,IVEC), FACTOR, NDIM)
        CALL SCALVE(SCR1,FACTOR,NDIM)
*. Subtract X(1,IVEC) from all remaining vectors
        DO JVEC = IVEC+1,NDIM
          XSX = INPROD(SCR1,X(1,JVEC),NDIM)
          CALL VECSUM(X(1,JVEC),X(1,JVEC),X(1,IVEC),ONE,-XSX,NDIM) 
        END DO
      END DO
*
*. And remove zero vectors
*
      NVECUT = 0
      DO IVEC = 1, NDIM
        XNORM = INPROD(X(1,IVEC),X(1,IVEC),NDIM)
        IF(XNORM.GT.0.0D0) THEN
          NVECUT = NVECUT + 1
          IF(NVECUT.NE.IVEC) CALL COPVEC(X(1,IVEC),X(1,NVECUT),NDIM)
        END IF
      END DO
*
      IF(NTEST.GE.1.AND.NVECUT.NE.NDIM) THEN
        WRITE(6,*)' MGS4 reduced dim, from and to ', NDIM,NVECUT
      ELSE IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthogonalization information:'
        WRITE(6,*) ' Number of linear independent vectors ', NVECUT
        WRITE(6,*) ' Orthonormalized vectors '
        CALL WRTMAT(X,NDIM,NVECUT,NDIM,NVECUT)
      END IF
*
      RETURN
      END

c $Id$
