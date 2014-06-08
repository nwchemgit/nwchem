* Version of Oct 21
      SUBROUTINE ZERO_OFFDIAG_BLM(A,NBLOCK,LBLOCK,IPACK)
*
* Zero off-diagonal elements in a block-diagonal matrix
* if IPACK .ne. 0, the matrix is packed (inner loop over columns)
*
* Jeppe Olsen, October 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*),LBLOCK(NBLOCK)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFF = 1
       ELSE
         LPREV = LBLOCK(IBLOCK-1)
         IF(IPACK.EQ.1) THEN
           IOFF = IOFF + LPREV*(LPREV+1)/2
         ELSE
           IOFF = IOFF + LPREV*LPREV
         END IF
       END IF
       L = LBLOCK(IBLOCK)
       CALL ZERO_OFFDIAG(A(IOFF),L,IPACK)
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output matrix from ZERO_OFFDIAG_BLM'
        CALL APRBLM2(A,LBLOCK,LBLOCK,NBLOCK,IPACK)
C       CALL APRBLM2(WORK(KFOCK),NTOOBS,NTOOBS,NSMOB,ISM)
      END IF
*
      RETURN
      END 
      SUBROUTINE ZERO_OFFDIAG(A,NDIM,IPACK)
*
* Zero off-diagonal elements in matrix A
*
* Jeppe Olsen, Oct 1997
*
* IF IPACK.NE.0, the matrix is assumed packed (inner loop over columns)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*)
*
      IF(IPACK.EQ.0) THEN
        DO I = 1, NDIM
          DO J = 1, NDIM
            IF(I.NE.J) THEN
              IJ = (J-1)*NDIM+I
              A(IJ) = 0.0D0
            END IF
          END DO
        END DO
      ELSE
        DO I = 1, NDIM
          DO J = 1, I-1
            IJ = I*(I-1)/2 + J
            A(IJ) = 0.0D0
          END DO
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
*
* Invert symmetric  - hopefully nonsingular - matrix A 
* by diagonalization
*
* Jeppe Olsen, Oct 97 to check INVMAT
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output matrix
      DIMENSION A(*)       
*. Scratch matrices and vector
      DIMENSION B(*),SCR(*),VEC(*)
*
      NTEST = 00
*. Reform a to symmetric packed form
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(A,SCR,1,NDIM,NDIM)
*. Diagonalize
      CALL EIGENL(SCR,B,NDIM,0,1)
      CALL COPDIA(SCR,VEC,NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMAT(VEC,NDIM,1,NDIM,1)
      END IF
*. Invert diagonal elements - without safety at the moment
      DO I = 1, NDIM
       IF(ABS(VEC(I)).GT.1.0D-15) THEN
         VEC(I) = 1.0D0/VEC(I)
       ELSE
         VEC(I) = 0.0D0
         WRITE(6,*) ' Singular mode inactivated '
       END IF
      END DO
*. and obtain inverse matrix by transformation
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(A,B,VEC,NDIM,SCR)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverse matrix from INVERSE_BY_DIAG'
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END 
      
      SUBROUTINE PTWFNRM(LU0,LUN,N,VEC1,VEC2,LUSCR1,LUSCR2,PTNORM)
*
* Norm of wavefunction obtained as perturbation expansion
* with intermediate normalization
*
* |0> = sum (i=0,n) |0(i)>
*
* Jeppe Olsen, Oct 97
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPRDD
*. Scratch blocks
      DIMENSION VEC1(*), VEC2(*)
*. Local scratch
      PARAMETER(MXORD = 100)
      DIMENSION SCR(MXORD)       
*
      IF(N.GT.MXORD) THEN
        WRITE(6,*) ' PTWFNRM : N > MXORD : ',N,MXORD
        WRITE(6,*) ' Increase parameter MXORD '  
        STOP       ' PTWFNRM : N > MXORD  '
      END IF
*
      LBLK = -1
      ONE = 1.0D0
*
      CALL COPVCD(LU0,LUSCR1,VEC1,1,LBLK)
      CALL SETVEC(SCR,ONE,N)
      CALL MVCSMD2(LUN,SCR,ONE,LUSCR1,LUSCR2,VEC1,VEC2,N,1,LBLK)
      PTNORM = INPRDD(VEC1,VEC2,LUSCR2,LUSCR2,1,LBLK)
      PTNORM = SQRT(PTNORM)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN 
         WRITE(6,*) ' PTWFNRM: Norm of summed vector ',PTNORM
      END IF
*
      RETURN
      END 
      SUBROUTINE COPDIA(A,VEC,NDIM,IPACK)
*
* Copy diagonal of matrix A into vector VEC
*
*   IPACK = 0 : Full matrix
*   IPACK = 1 : Lower triangular matrix
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION A(*),VEC(*)
*
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          VEC(I) = A((I-1)*NDIM+I)
  100   CONTINUE
      ELSE
        DO 200 I = 1, NDIM
          VEC(I) = A(I*(I+1)/2)
  200   CONTINUE
      END IF
*
      RETURN
      END
      SUBROUTINE SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
*
* Calculate square root of positive definite symmetric matrix A
* if(ITASK .EQ. 2 ) Inverted square root matrix is also calculated
* In case of singularities in A A -1/2 is defined to have the same
* singularity
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
*
      DIMENSION A(NDIM,NDIM)
      DIMENSION ASQRT(NDIM,NDIM),AMSQRT(NDIM,NDIM)
      DIMENSION SCR(*)
*
      NTEST = 0
*
* Length of SCR should at least be 2 * NDIM ** 2 + NDIM*(NDIM+1)/2
      KLFREE = 1
*
      KLASYM = KLFREE
      KLAVAL = KLASYM
      KLFREE = KLASYM + NDIM*(NDIM+1)/2
*
      KLAVEC = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
*
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(A,SCR(KLASYM),1,NDIM,NDIM)
      CALL EIGENL(SCR(KLASYM),SCR(KLAVEC),NDIM,0,1)
      CALL COPDIA(SCR(KLASYM),SCR(KLAVAL),NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMAT(SCR(KLAVAL),NDIM,1,NDIM,1)
      END IF
*. Check for negative eigenvalues
      DO I = 1, NDIM
       IF(SCR(KLAVAL-1+I).LT.0.0D0) THEN
         WRITE(6,*) ' SQRTMT : Negative eigenvalue ', SCR(KLAVAL-1+I)
         WRITE(6,*) ' SQRTMT : I will STOP '
         STOP       ' SQRTMT : Negative eigenvalue '
       END IF
      END DO
*
      DO 100 I = 1,NDIM
        SCR(KLAVAL-1+I) = SQRT(SCR(KLAVAL-1+I))
  100 CONTINUE
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(ASQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
*
      IF(ITASK .EQ. 2 ) THEN
        DO 200 I = 1,NDIM
          IF(SCR(KLAVAL-1+I) .GT. 1.0D-13 ) then
            SCR(KLAVAL-1+I) = 1.0D0/SCR(KLAVAL-1+I)
          ELSE
            SCR(KLAVAL-1+I) = SCR(KLAVAL-1+I)
          END IF
  200   CONTINUE
        CALL XDIAXT(AMSQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
      END IF
*
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Info from SQRTMT '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' Input matrix to SQRTMT '
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Square root of matrix '
        CALL WRTMAT(ASQRT,NDIM,NDIM,NDIM,NDIM)
        IF(ITASK .EQ. 2 ) THEN
          WRITE(6,*) ' Inverse square root of matrix '
          CALL WRTMAT(AMSQRT,NDIM,NDIM,NDIM,NDIM)
        END IF
      END IF
*
      RETURN
      END

      SUBROUTINE XDIAXT(XDX,X,DIA,NDIM,SCR)
*
* Obtain XDX = X * DIA * X(Transposed)
* where DIA is an diagonal matrix stored as a vector
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XDX(NDIM,NDIM)
      DIMENSION X(NDIM,NDIM),DIA(NDIM)
      DIMENSION SCR(NDIM,NDIM)
*
* DIA * X(transposed)
      DO 100 I=1,NDIM
        CALL COPVEC(X(1,I),SCR(1,I),NDIM)
        CALL SCALVE(SCR(1,I),DIA(I),NDIM)
  100 CONTINUE
* X * DIA * X(Transposed)
      CALL MATML4(XDX,X,SCR,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,2)
*
      RETURN
      END

      SUBROUTINE EXTKOP
*
* Use Extended Koopmaans' theorem to obtain ionization energies 
*
*
* Jeppe Olsen , LUCIA version adapeted from LUCAS, Nov 1995
*
c      IMPLICIT REAL*8 (A-H,O-Z)
*
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
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"


*
      WRITE(6,*) '        *****************************************'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *  Welcome to the magic(?) world of the *'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *      Extended Koopmans Theorem        *'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *       Version of Dec 12 1996          *'
      WRITE(6,*) '        *                                       *'
      WRITE(6,*) '        *****************************************'
C     KFREEO = KFREE
      IDUM = 0
      WRITE(6,*) ' NTOOB = ', NTOOB 
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXTKOP')
      CALL MEMMAN(KDEN1SM,2*NINT1,'ADDL ',2,'DEN1SM')
      LENGTH = NTOOB ** 2
      CALL MEMMAN(KLSSQ,LENGTH,'ADDL  ',2,'SSQ   ')
      CALL MEMMAN(KLSSQI,LENGTH,'ADDL  ',2,'SSQI  ')
      CALL MEMMAN(KLFP  ,LENGTH,'ADDL  ',2,'FP    ')
      CALL MEMMAN(KLFEXP,LENGTH,'ADDL  ',2,'KLFEXP')
      CALL MEMMAN(KLSCR ,4*LENGTH,'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KFOCK ,LENGTH ,'ADDL  ',2,'KFOCK ')
      
*. Density matrices and integrals assumed in place 
*
*. Construct Fock matrix over all orbitals
      CALL FOCK_MAT(WORK(KFOCK),2)  
*. Multiply Fock matrix with a factor of 2 to improve results
*. - a deeper argument does exist
COLD  TWO = 2.0D0
*. will give problems if we reduve size of FOCL matrix
*   to symmetry blocked form
COLD  CALL SCALVE(WORK(KFOCK),TWO,NTOOB**2)
      WRITE(6,*) ' Fock matrix '
      WRITE(6,*) ' ============'
      ISM = 0
      CALL APRBLM2(WORK(KFOCK),NTOOBS,NTOOBS,NSMOB,ISM)
*. Obtain density matrix in symmetry blocks
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFFM=1
          IOFFO=1
        ELSE
          IOFFM = IOFFM + NTOOBS(ISM-1)** 2
          IOFFO = IOFFO + NTOOBS(ISM-1)
        END IF
        DO IOB = IOFFO,IOFFO + NTOOBS(ISM)-1
           LOB = NTOOBS(ISM)
           DO JOB = IOFFO,IOFFO + NTOOBS(ISM)-1
*. Corresponding type indeces
             IOBP = IREOST(IOB)
             JOBP = IREOST(JOB)
             WORK(KDEN1SM-1+IOFFM-1+(JOB-IOFFO)*LOB+IOB-IOFFO+1)
     &     = dbl_mb(KRHO1-1+(IOBP-1)*NTOOB+JOBP)
           END DO
        END DO
        WRITE(6,*) ' Symmetry packed one body density matrix '
        IMATSM = 0 
        CALL APRBLM2(WORK(KDEN1SM),NTOOBS,NTOOBS,NSMOB,IMATSM)
      END DO
*. Loop over Symmetries and  and solve F C = E S C for each symmetry
      DO 200 ISM = 1, NSMOB 
        N = NTOOBS(ISM)
        IF(ISM.EQ.1) THEN
           IBASE = 1
           IBASE2 = 1
           IBASE2P = 1
        ELSE 
           IBASE = IBASE + NTOOBS(ISM-1)
           IBASE2P = IBASE2P + NTOOBS(ISM-1)*(NTOOBS(ISM-1)+1)/2
           IBASE2 = IBASE2 + NTOOBS(ISM-1)**2                   
        END IF
*. Pack block of Fock matrix to complete form
C       CALL TRIPAK(WORK(KLFEXP),WORK(KFOCK-1+IBASE2P),2,N,N)
*. Already in complete form so
        CALL COPVEC(WORK(KFOCK-1+IBASE2),WORK(KLFEXP),N**2)
        WRITE(6,*) ' Outpacked block of FOCK '
        CALL WRTMAT(WORK(KLFEXP),n,n,n,n)

*. Solve FC = eSC with Ct S C = 1 as F'C' = e S'C'( S is one body density)
* with F' = S-(1/2)F S(-1/2), C' = S(1/2)C
*. S(-1/2),S(1/2)
        CALL SQRTMT(dbl_mb(KRHO1+IBASE2-1),N,2,WORK(KLSSQ),WORK(KLSSQI),
     &              WORK(KLSCR))
        write(6,*) ' S-1/2 matrix '
        call wrtmat(work(klssqi),n,n,n,n)
*. S(-1/2)F S(-1/2)
        CALL TRNMAT(WORK(KLFEXP),WORK(KLSSQI),WORK(KLSCR),N,N)
        write(6,*) ' Transformed matrix'
        call wrtmat(work(KLFEXP),n,n,n,n)
*. Diagonalize
        CALL TRIPAK(WORK(KLFEXP),WORK(KLSCR),1,N,N)
        CALL EIGENL(WORK(KLSCR),WORK(KLFP),N,0,1) 
        WRITE(6,*) 
     &  ' Eigenvalues of Extended Koopmans Equation for ISM = ',ISM
        DO 80 I = 1, N
          WRITE(6,'(F25.15)') WORK(KLSCR-1+I*(I+1)/2) 
   80   CONTINUE
  200 CONTINUE
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXTKOP')
      RETURN
      END
      SUBROUTINE PERTDN_OLD
     &(N,LU0,LUN,ISM,ISPC,VEC1,VEC2,RHO1N,RHO2N,LUSC1,LUSC2)
*
* Construct one body density matrix of order N
*
*      Jeppe + Dage, Nov. 11 1995
*                    Debugged Jan 31 '97
*
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
      CALL MEMMAN(KLDEN2,LRHO2,'ADDL  ',2,'KLDEN2')
*
      LBLK = -1
      ZERO = 0.0D0
      CALL SETVEC(RHO1N,ZERO,LRHO1)
      CALL SETVEC(RHO2N,ZERO,LRHO2)
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
C?      write(6,*) ' next pair of vectors loaded '
* Do the densi
          LEQR = 0
C DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR)
          I12 = 2
C?        WRITE(6,*) ' Calling DENSI2 '
          XDUM = 0.0D0
          CALL DENSI2(I12,WORK(KLDEN1),WORK(KLDEN2),VEC1,VEC2,
     &                LUSC1,LUSC2,EXPS2,0,XDUM,XDUM,XDUM,XDUM,1)
          WRITE(6,*) ' Home from DENSI2 '
C         WRITE(6,*) ' Densities fresh from DENSI2 '
C?        WRITE(6,*) ' NTOOB = ', NTOOB
C         CALL WRTMAT(WORK(KLDEN1),NTOOB,NTOOB,NTOOB,NTOOB)
C         CALL PRSYM(WORK(KLDEN2),NTOOB**2)
*
CTOBE     IF(L.NE.NMINL) THEN
*. The matrix <L! E !NMINL> was calculated, add <NMINL! E ! L> 
*. as simple transposition
CTOBE        CALL TRPAD(WORK(KLDEN),ONE,NTOOB)
CTOBE     END IF
          ONE = 1.0D0
          CALL VECSUM(RHO1N,RHO1N,WORK(KLDEN1),ONE,ONE,LRHO1)
          CALL VECSUM(RHO2N,RHO2N,WORK(KLDEN2),ONE,ONE,LRHO2)
CTOBE   END IF
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Density matrix of order in perturbation ', N
        WRITE(6,*) ' ==========================================='
        WRITE(6,*)
        WRITE(6,*) ' One-body density '
        WRITE(6,*) ' ================ '
        CALL WRTMAT(RHO1N,NTOOB,NTOOB,NTOOB,NTOOB)
        WRITE(6,*) ' Two-body density '
        WRITE(6,*) ' ================ '
        CALL PRSYM(RHO2N,NTOOB**2)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ', IDUM,'PERTDN')
*
      RETURN
      END
      SUBROUTINE EKTPERT(F,S,NDIM,NORD,EN,C,
     &                   VEC1,VEC2,VEC3,AMAT1,AMAT2,AMAT3)
*
* Perturbation expansion of generalized eigenvalue problem 
* Special version for EKT problem where there are singularities in 
* the zero order matrices.
*
* Ordering the matrices so the occupied orbitals come first, and
* then the virtual orbitals, the  zero order matrices have the form
*
*     (x *                     ) 
*     ( x*               0     )
*     (**x                     )
*     (                        )
*     (                        )
*     (                        )
*     (                        )
*     (                        )
*     (  0              0      )
*     (                        )
*     (                        )
*     (                        )
*
*
*
* The matrices in the  eigenvalue problem FC = E SC
* are separated into orders
*
* F = Sum(k=0,NORD) F(K)
* S = Sum(K=0,NORD) S(K)
*
*. Obtain corrections to energy and wawe functions 
*
*. The normalization condition used is C(K)T S(0) C(0) = 0
*
* The energy corrections become
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I) 
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
* and the wave function corrections in the occupied orbital space are
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(I=0,N-K)E(N-K-L)S(L)C(K))
*
* Whereas they read in the virtual space
*
* C(N-2) = (F(2)-E(0)S(2))-1 (-Sum(K=0,N-3)F(N-K)C(K)
*                             +Sum(K=0,N-3)Sum(I=0,N-K)E(N-K-L)S(L)C(K)
*                             + -(F(2)-E0*S(2)) C(N-2)OCC )
* Where only the 
*
* The zero order matrices F(0),S(0) are assumed diagonal
*
* Jeppe and Dage, Oct 22 1997

*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION F(NDIM**2,*),S(NDIM**2,*)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(0:NORD)
*. Scratch 
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
      DIMENSION AMAT1(NDIM,NDIM),AMAT2(NDIM,NDIM)
      DIMENSION AMAT3(NDIM,NDIM)
*
      NTEST = 1
 
*. Zero order wavefunction in virtual space
* obtained by solving - in the virtual space, the
* equations
*
* (F(2)-E(0)S(2))(virt,virt) C0(virt) = -(F(2)-E(0)S(2))C0(occ))(virt)
* 
*
*. Obtain the number of occupied orbitals by examining S(0)
      NOCC = 0
      DO I = 1, NDIM
       IF(S((I-1)*NDIM+I,1) .NE. 0) NOCC = NOCC + 1
      END DO 
      NVIRT = NDIM - NOCC
      WRITE(6,*) ' Number of occupied orbitals ', NOCC
      WRITE(6,*) ' Number of virtual  orbitals ', NVIRT
*
*. Zero order energy
*
C          MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(F,C,VEC1,NDIM,NDIM,0)
      C0FC0 = INPROD(VEC1,C,NDIM)
*
      CALL MATVCB(S,C,VEC1,NDIM,NDIM,0)
      C0SC0 = INPROD(VEC1,C,NDIM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'C(0)T F(0) C(0) = ', C0FC0
        WRITE(6,*) 'C(0)T S(0) C(0) = ', C0SC0
      END IF
      E0 = C0FC0/C0SC0
      EN(0) = E0    
*. Save diagonal of F(0) - E(0)C(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = F((I-1)*NDIM+I,1)-E0*S((I-1)*NDIM+I,1)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Zero order diagonal '
        CALL WRTMAT(VEC3,1,NDIM,1,NDIM)
      END IF
*
*. The virtual-virtual  part of the operator F(2)-E(0)S(2)
*. will be used to solve linear equations. Set up inverse
*
      ONE = 1.0D0
      CALL VECSUM(AMAT1,F(1,3),S(1,3),ONE,-E0,NDIM**2)
      DO IOCC = 1, NOCC
        DO J = 1, NDIM
          AMAT1(IOCC,J) = 0.0D0
          AMAT1(J,IOCC) = 0.0D0
        END DO
        AMAT1(IOCC,IOCC) = 0.001
      END DO
C?    write(6,*) ' Input matrix to INVMAT'
C?    CALL WRTMAT(AMAT1,NDIM,NDIM,NDIM,NDIM)
      CALL COPVEC(AMAT1,AMAT3,NDIM**2)
      CALL INVMAT(AMAT1,AMAT2,NDIM,NDIM,ISING)
*     ^ Returns inverse matrix in AMAT1
      WRITE(6,*) ' Inverse matrix obtained by INVMAT'
      CALL WRTMAT(AMAT1,NDIM,NDIM,NDIM,NDIM)
*. Inverse matrix by diagonalization
      CALL COPVEC(AMAT3,AMAT1,NDIM**2)
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
      CALL INVERT_BY_DIAG(AMAT1,AMAT2,AMAT3,VEC1,NDIM)
*. Zero occ-occ part
      DO IOCC = 1, NOCC
        DO J = 1, NDIM
          AMAT1(IOCC,J) = 0.0D0
          AMAT1(J,IOCC) = 0.0D0
        END DO
      END DO
*. -(F(2)-E(0)S(2)) Will be used in the future, reconstruct in AMAT2   
      CALL VECSUM(AMAT2,F(1,3),S(1,3),ONE,-E0,NDIM**2)
      ONEM = -1.0D0
      CALL SCALVE(AMAT2,ONEM,NDIM**2)
*. Obtain virtual part of C0
* (F(2)-E(0)S(2))(virt,virt) C0(virt) = -(F(2)-E(0)S(2))C0(occ))(virt)
      CALL MATVCB(AMAT2,C(1,1),VEC1,NDIM,NDIM,0)
      CALL MATVCB(AMAT1,VEC1,VEC2,NDIM,NDIM,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Virtual part of C0 '
        CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
      END IF
      CALL COPVEC(VEC2(1+NOCC),C(1+NOCC,1),NVIRT)
C?    WRITE(6,*) ' complete zero order correction vector '
C?    CALL WRTMAT(C(1,1),1,NDIM,1,NDIM)
*. And zero it once again ( COnstructed below )
      ZERO = 0.0D0
      CALL SETVEC(C(1+NOCC,1),ZERO,NVIRT)
*
*. And then start the iterations
      DO IORD = 1, NORD
*
* =================================================
* The (IORD-2) wf corrections in the virtual space
* =================================================
*
        IF(IORD.GE.2) THEN
          ZERO = 0.0D0
          CALL SETVEC(VEC2,ZERO,NDIM)
*. Note : Only the occupied part of C(N-2) is included in RHS of
*         expression. The virtual part was carefully zeroed !
          DO K = 0, IORD -2
            CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
            ONE = 1.0D0
            ONEM = -1.0D0
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
          END DO
*
          DO K = 0, IORD -2
            DO L = 0, IORD -K
              CALL MATVCB(S(1,L+1),C(1,K+1),VEC1,NDIM,NDIM,0)
              CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
            END DO
          END DO
*. Multiply with (E(2)-E0S(2))-1
          CALL MATVCB(AMAT1,VEC2,VEC1,NDIM,NDIM,0)
*. And save
          CALL COPVEC(VEC1(1+NOCC),C(1+NOCC,IORD-2+1),NVIRT)
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' occ+virtual part for order = ',IORD-2
            CALL WRTMAT(C(1,IORD-2+1),1,NDIM,1,NDIM)
          END IF
*
        END IF
*       ^ End of construction of virtual part of C(VIRT,IORD-2) 
*
*  =================
*. Energy correction
*  =================
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I) 
*      - Sum(I=0,N-1)Sum(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
        EN(IORD) = 0.0D0
        DO I = 1, IORD
          CALL MATVCB(F(1,I+1),C(1,IORD-I+1),VEC1,NDIM,NDIM,0)
          EN(IORD) = EN(IORD)+INPROD(C,VEC1,NDIM)
        END DO
C?      write(6,*) ' First term to En ', EN(IORD)
        DO I = 0,IORD-1
          DO J = 1, IORD-I
            CALL MATVCB(S(1,J+1),C(1,I+1),VEC1,NDIM,NDIM,0)
            EN(IORD) = EN(IORD)-EN(IORD-I-J)*INPROD(C,VEC1,NDIM)
C?          write(6,*) ' second term to EN: I J EN ',I,J,EN(IORD)
          END DO
        END DO
        EN(IORD) = EN(IORD)/C0SC0
        WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ===========================================
*. Wave function corrections, occupied of IORD
*  ===========================================
*
* The occupied part obtained from
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(L=0,N-K)E(N-K-L)S(L)C(K))
        ZERO = 0.0D0
        CALL SETVEC(VEC2,ZERO,NDIM)
        DO K = 0, IORD -1
          CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
          ONE = 1.0D0
          ONEM = -1.0D0
          CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
        END DO
C?      write(6,*) ' first term to rhs '
C?      CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
*
        DO K = 0, IORD -1
          DO L = 0, IORD -K
            CALL MATVCB(S(1,L+1),C(1,K+1),VEC1,NDIM,NDIM,0)
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
C?          write(6,*) ' second term to rhs, K,L = ', K,L
C?          CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
          END DO
        END DO
*. Check overlap with S(0) times zero order state ( should be zero )
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
        CALL MATVCB(S(1,1),C(1,1),VEC1,NDIM,NDIM,0)
        C0SSC0 = INPROD(VEC1,VEC1,NDIM)
        OVLAP = INPROD(VEC1,VEC2,NDIM)
        WRITE(6,*) ' OVLAP = ', OVLAP
        FACTOR = -OVLAP/C0SSC0
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Vector before ortho'
          CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
        END IF
        CALL VECSUM(VEC2,VEC2,VEC1  ,ONE,FACTOR,NDIM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Vector before DIAVC2'
          CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
        END IF
*. Multiply with (F(0)-E(0)S(0))-1
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
*
        JEPZAP = 1
        IF(IORD.EQ.1.AND.JEPZAP.EQ.1) THEN
          WRITE(6,*) ' First order correction zapped '
          ZERO = 0.0D0
          CALL SETVEC(VEC1,ZERO,NOCC)
        END IF
*
        CALL COPVEC(VEC1(1),C(1,IORD+1),NOCC)
*. The virtual part is still not known
        ZERO = 0.0D0
        CALL SETVEC(C(1+NOCC,IORD+1),ZERO,NVIRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Eigenfunction correction ', IORD
          CALL WRTMAT(C(1,IORD+1),1,NDIM,1,NDIM)
        END IF
      END DO
* 
      WRITE(6,*) ' Energy corrections : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) '   Order        Correction '
      WRITE(6,*) ' ==============================='
      DO IORD = 1, NORD
        WRITE(6,'(1H ,3X,I3,E25.13)')IORD,EN(IORD)
      END DO
*
      ETOT = E0
      DO IORD = 1, NORD
        ETOT = ETOT + EN(IORD)
      END DO
      WRITE(6,'(A,E25.13)') ' Zero-order energy ', E0
      WRITE(6,'(A,E25.13)') ' Sum(K=0,NORD) E(K) ', ETOT 
*
      RETURN
      END
      SUBROUTINE GENPERT(F,S,NDIM,NORD,EN,C,VEC1,VEC2,VEC3)
*
* Perturbation expansion of generalized eigenvalue problem 
*
* The matrices in the  eigenvalue problem FC = E SC
* are separated into orders
*
* F = Sum(k=0,NORD) F(K)
* S = Sum(K=0,NORD) S(K)
*
*. Obtain corrections to energy and wawe functions 
*
*. The normalization condition used is C(K)T S(0) C(0) = 0
*
* The energy corrections become
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I) 
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
* and the wave function corrections are
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(I=0,N-K)E(N-K-L)S(L)C(K))
*
* The zero order matrices F(0),S(0) are assumed diagonal
*
* Jeppe and Dage, Jan 31 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION F(NDIM**2,*),S(NDIM**2,*)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(0:NORD)
*. Scratch 
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
*
*. Zero order energy
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(F,C,VEC1,NDIM,NDIM,0)
      C0FC0 = INPROD(VEC1,C,NDIM)
*
      CALL MATVCB(S,C,VEC1,NDIM,NDIM,0)
      C0SC0 = INPROD(VEC1,C,NDIM)
*
      WRITE(6,*) 'C(0)T F(0) C(0) = ', C0FC0
      WRITE(6,*) 'C(0)T S(0) C(0) = ', C0SC0
      E0 = C0FC0/C0SC0
      EN(0) = E0    
*. Save diagonal of F(0) - E(0)C(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = F((I-1)*NDIM+I,1)-E0*S((I-1)*NDIM+I,1)
      END DO
      WRITE(6,*) ' Zero order diagonal '
      CALL WRTMAT(VEC3,1,NDIM,1,NDIM)
*. And then start the iterations
      DO IORD = 1, NORD
*
*  =================
*. Energy correction
*  =================
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I) 
*      - Sum(I=0,N-1)Sum(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
        EN(IORD) = 0.0D0
        DO I = 1, IORD
          CALL MATVCB(F(1,I+1),C(1,IORD-I+1),VEC1,NDIM,NDIM,0)
          EN(IORD) = EN(IORD)+INPROD(C,VEC1,NDIM)
        END DO
C?      write(6,*) ' First term to En ', EN(IORD)
        DO I = 0,IORD-1
          DO J = 1, IORD-I
            CALL MATVCB(S(1,J+1),C(1,I+1),VEC1,NDIM,NDIM,0)
            EN(IORD) = EN(IORD)-EN(IORD-I-J)*INPROD(C,VEC1,NDIM)
C?          write(6,*) ' second term to EN: I J EN ',I,J,EN(IORD)
          END DO
        END DO
        EN(IORD) = EN(IORD)/C0SC0
        WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ==========================
*. Wave function corrections
*  ==========================
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(L=0,N-K)E(N-K-L)S(L)C(K))
        ZERO = 0.0D0
        CALL SETVEC(VEC2,ZERO,NDIM)
        DO K = 0, IORD -1
          CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
          ONE = 1.0D0
          ONEM = -1.0D0
          CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
        END DO
C?      write(6,*) ' first term to rhs '
C?      CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
*
        DO K = 0, IORD -1
          DO L = 0, IORD -K
            CALL MATVCB(S(1,L+1),C(1,K+1),VEC1,NDIM,NDIM,0)
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
C?          write(6,*) ' second term to rhs, K,L = ', K,L
            CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
          END DO
        END DO
*. Check overlap with S(0) times zero order state ( should be zero )
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
        CALL MATVCB(S(1,1),C(1,1),VEC1,NDIM,NDIM,0)
        OVLAP = INPROD(VEC1,VEC2,NDIM)
        FACTOR = -OVLAP/C0SC0
        WRITE(6,*) ' Vector before ortho'
        CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
        CALL VECSUM(VEC2,VEC2,C(1,1),ONE,FACTOR,NDIM)
        WRITE(6,*) ' OVLAP = ', OVLAP
        WRITE(6,*) ' Vector before DIAVC2'
        CALL WRTMAT(VEC2,1,NDIM,1,NDIM)
*. Multiply with (F(0)-E(0)S(0))-1
C            DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
*
        JEPZAP = 1
        IF(IORD.EQ.1.AND.JEPZAP.EQ.1) THEN
          WRITE(6,*) ' First order correction zapped '
          ZERO = 0.0D0
          CALL SETVEC(VEC1,ZERO,NDIM)
        END IF
*
        CALL COPVEC(VEC1,C(1,IORD+1),NDIM)
*
        WRITE(6,*) ' Eigenfunction correction ', IORD
        CALL WRTMAT(C(1,IORD+1),1,NDIM,1,NDIM)
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
      ETOT = E0
      DO IORD = 1, NORD
        ETOT = ETOT + EN(IORD)
      END DO
      WRITE(6,*) ' Zero-order energy ', E0
      WRITE(6,*) ' Sum(K=0,NORD) E(K) ', ETOT 
*
      RETURN
      END


          
          





*
   
         


          

         





      
c $Id$
