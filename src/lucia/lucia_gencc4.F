      SUBROUTINE SKPVCD2(LU,NVEC,NVAR,SEGMNT,IREW,LBLK) 
*
* Skip over NVEC calls to FRMDSC
*
      INCLUDE 'implicit.inc'
*. Scratch 
      DIMENSION SEGMNT(*)
      IF(IREW.EQ.1) CALL REWINO(LU)
      DO IVEC = 1, NVEC
        CALL FRMDSC(SEGMNT,NVAR,LBLK,LU,IMZERO,IAMPACK)
      END DO
*
      RETURN
      END
*
      SUBROUTINE SYM_BLMAT(A,NBLOCK,LBLOCK)
*
* Symmetrize blocked matrix A
*
      INCLUDE 'implicit.inc'
*. Input 
      INTEGER LBLOCK(NBLOCK)
*. Input and output
      DIMENSION A(*)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Matrix to be averaged  '
        CALL APRBLM2(A,LBLOCK,LBLOCK,NBLOCK,0)
      END IF
*
      IOFF = 1
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
        IOFF = 1
       ELSE
        IOFF = IOFF + LBLOCK(IBLOCK-1)**2
       END IF
       LEN = LBLOCK(IBLOCK)
c       WRITE(6,*) ' IBLOCK, LEN, IOFF = ', IBLOCK,LEN,IOFF
       DO I = 1, LEN
        DO J = 1, LEN
          AVERAGE = 0.5D0*(A(IOFF-1+(J-1)*LEN+I) +
     &                     A(IOFF-1+(I-1)*LEN+J)   )
          A(IOFF-1+(J-1)*LEN+I) = AVERAGE
          A(IOFF-1+(I-1)*LEN+J) = AVERAGE
        END DO
       END DO
*.     ^ End of loop over I and J
      END DO
*.    ^ End of loop over blocks
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Averaged blocked matrix '
        CALL APRBLM2(A,LBLOCK,LBLOCK,NBLOCK,0)
C            APRBLM2(AOUT,LX_COL,LX_COL,NBLOCK,ISYM)
      END IF
*
      RETURN
      END
      SUBROUTINE SYM_2DENS(RHO2,NDIM,IMODE)
*
*     symmetrize the 2-electron density: 
*
*        d(ijkl) = 0.5 (d(ijkl) + d(jilk))
*
*     addressing:  ijkl = ij*(ij-1)/2 + kl  (ij.ge.kl)
*
*     if IMODE.EQ.1 ijkl goes over the full square
*
      IMPLICIT NONE

      REAL(8), INTENT(INOUT) ::
     &     RHO2(*)
      INTEGER, INTENT(IN)    ::
     &     NDIM, IMODE

      INTEGER ::
     &     II, JJ, KK, LL,
     &     IADIJ, IADJI, IADKL, IADLK, IAD1, IAD2
      REAL(8) ::
     &     AV

      IF (IMODE.EQ.0) THEN
        DO II = 1, NDIM
          DO JJ = 1, NDIM
            IADIJ = (II-1)*NDIM+JJ
            IADJI = (JJ-1)*NDIM+II

            DO KK = 1, NDIM
              DO LL = 1, NDIM
                IADKL = (KK-1)*NDIM+LL
                IADLK = (LL-1)*NDIM+KK
                
                IF (IADKL.GT.IADIJ) CYCLE
                IAD1 = IADIJ*(IADIJ-1)/2+IADKL

                IF (IADLK.LE.IADJI) THEN
                  IAD2 = IADJI*(IADJI-1)/2+IADLK
                ELSE
                  IAD2 = IADLK*(IADLK-1)/2+IADJI
                END IF
                
                AV = 0.5d0*(RHO2(IAD1)+RHO2(IAD2))
                RHO2(IAD1) = AV
                RHO2(IAD2) = AV
              END DO
            END DO
          END DO
        END DO
        
      ELSE
        DO II = 1, NDIM
          DO JJ = 1, NDIM
            IADIJ = (II-1)*NDIM+JJ
            IADJI = (JJ-1)*NDIM+II

            DO KK = 1, NDIM
              DO LL = 1, NDIM
                IADKL = (KK-1)*NDIM+LL
                IADLK = (LL-1)*NDIM+KK
                
                IAD1 = (IADIJ-1)*NDIM*NDIM+IADKL

                IAD2 = (IADJI-1)*NDIM*NDIM+IADLK
                
                AV = 0.5d0*(RHO2(IAD1)+RHO2(IAD2))
                RHO2(IAD1) = AV
                RHO2(IAD2) = AV
              END DO
            END DO
          END DO
        END DO

      END IF

      RETURN
      END
      SUBROUTINE CLEAN_2DENS(RHO2,NDIM)
*
*     clean the 2-electron correlation density from contributions
*     that are not present in the frozen core approx.
*
*     addressing:  ijkl = ij*(ij-1)/2 + kl  (ij.ge.kl)
*
      IMPLICIT NONE

      REAL(8), INTENT(INOUT) ::
     &     RHO2(*)
      INTEGER, INTENT(IN)    ::
     &     NDIM

      INTEGER ::
     &     II, JJ, KK, LL,
     &     IADIJ, IADKL, IAD, NFR

* FUSK
      NFR = 1

      DO II = 1, NDIM
        DO JJ = 1, NDIM
          IADIJ = (II-1)*NDIM+JJ

          DO KK = 1, NDIM
            DO LL = 1, NDIM
              IF (II.GT.NFR.AND.JJ.GT.NFR.AND.
     &            KK.GT.NFR.AND.LL.GT.NFR) CYCLE

              IADKL = (KK-1)*NDIM+LL
              
              IF (IADKL.GT.IADIJ) CYCLE

              IAD = IADIJ*(IADIJ-1)/2+IADKL
              
              IF (RHO2(IAD).GT.1D-12)
     &            WRITE(6,'(4i5,E15.5)') II,JJ,KK,LL,RHO2(IAD)
              RHO2(IAD) = 0d0

            END DO
          END DO
        END DO
      END DO

      RETURN
      END
      SUBROUTINE MKIR2DENS(IMODE,RHO1,RHO2,NDIM,NSPIN)
*
*     generate irredudible two-electron density
*     i.e. remove that part of RHO2 that is expressible through RHO1
*
*     addressing:  ijkl = ij*(ij-1)/2 + kl  (ij.ge.kl)
*
      IMPLICIT NONE

      INTEGER, PARAMETER ::
     &     NTEST = 000

      REAL(8), INTENT(INOUT) ::
     &     RHO1(*), RHO2(*)
      INTEGER, INTENT(IN)    ::
     &     IMODE, NDIM, NSPIN

      INTEGER ::
     &     II, JJ, KK, LL,
c     &     IADKI, IADIK,
     &     IADIJ, IADKL,
     &     IADIL, IADJK, IAD
      REAL(8) ::
     &     RHOOLD


      IF (NTEST.GT.0) THEN
        WRITE(6,*) 'MKIR2DENS AT WORK'
        WRITE(6,*) '================='
      END IF

      IF (NSPIN.NE.1) STOP 'ADAPT THIS ROUTINE'
*. remove separable part from 2-el-density:
      DO II = 1, NDIM
        DO JJ = 1, NDIM
          IADIJ = (II-1)*NDIM+JJ
            
          DO KK = 1, NDIM
c            IADKJ = (KK-1)*NDIM+JJ
            ! we assume that RHO1 is symmetric
            IADJK = (JJ-1)*NDIM+KK

            DO LL = 1, NDIM
                
              IADKL = (KK-1)*NDIM+LL
              IADIL = (II-1)*NDIM+LL
              
              IF (IADKL.GT.IADIJ) CYCLE
                
              IAD = IADIJ*(IADIJ-1)/2+IADKL
                
              RHOOLD = RHO2(IAD)

              RHO2(IAD) = RHO2(IAD) - RHO1(IADIJ)*RHO1(IADKL)
     &                              + 0.5d0*RHO1(IADIL)*RHO1(IADJK)

              IF (NTEST.GE.1000) THEN
c                IF (ABS(RHO2(IAD)).GT.1D-12.OR.ABS(RHOOLD.GT.1D-12))
                IF (ABS(RHO2(IAD)).LT.1D-12.AND.ABS(RHOOLD).GT.1D-12)
     &               WRITE(6,'(4i5,E15.5," -> ",3E15.5,2L)')
     &               II,JJ,KK,LL,RHOOLD,RHO2(IAD),
     &               RHO1(IADIJ)*RHO1(IADKL),
     &               RHO1(IADIL)*RHO1(IADJK)
     &               ,RHO2(IAD).LT.1D-12,RHOOLD.GT.1D-12
              END IF

            END DO
          END DO
        END DO
      END DO

      IF (NTEST.GE.500) THEN
        CALL SWEEPVEC(RHO2,(NDIM**2)*(NDIM**2+1)/2)
        WRITE(6,*) '2-particle density: '
        WRITE(6,*) 'index index:'
        IAD = 0
        DO II = 1, NDIM
          DO JJ = 1, NDIM
            IAD = IAD + 1
            WRITE(6,*) II,JJ,' --> ',IAD
          END DO
        END DO
        CALL PRTRLT(RHO2,NDIM**2)
      END IF

      RETURN
      END
      SUBROUTINE PRINT_2DENS(RHO2,NDIM,IMODE,THRSH)
*
*     print the 2-electron density: 
*
*     addressing:  ijkl = ij*(ij-1)/2 + kl  (ij.ge.kl)
*
      IMPLICIT NONE

      REAL(8), INTENT(INOUT) ::
     &     RHO2(*)
      REAL(8), INTENT(IN) ::
     &     THRSH
      INTEGER, INTENT(IN)    ::
     &     NDIM, IMODE

      INTEGER ::
     &     II, JJ, KK, LL,
     &     IADIJ, IADKL, IAD

      IF (IMODE.EQ.0) THEN
        DO II = 1, NDIM
          DO JJ = 1, NDIM
            IADIJ = (II-1)*NDIM+JJ
            
            DO KK = 1, NDIM
              DO LL = 1, NDIM
                IADKL = (KK-1)*NDIM+LL
                
                IF (IADKL.GT.IADIJ) CYCLE
                IAD = IADIJ*(IADIJ-1)/2+IADKL
                
                IF (RHO2(IAD).GT.THRSH)
     &               WRITE(6,'(5i5,E15.5)') II,JJ,KK,LL,IAD,RHO2(IAD)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO II = 1, NDIM
          DO JJ = 1, NDIM
            IADIJ = (II-1)*NDIM+JJ

            DO KK = 1, NDIM
              DO LL = 1, NDIM
                IADKL = (KK-1)*NDIM+LL
                
                IAD = (IADIJ-1)*NDIM*NDIM+IADKL

                IF (RHO2(IAD).GT.THRSH)
     &               WRITE(6,'(5i5,E15.5)') II,JJ,KK,LL,IAD,RHO2(IAD)

              END DO
            END DO
          END DO
        END DO

      END IF

      RETURN
      END

      SUBROUTINE GEN_FOCK_MAT(F,I12)
*
* Construct general Fock matrix using zero-order Fock matrix
*
* F(i,j) = SUM(K) F0(i,K) * RHO1(j,K)
*          + SUM(KL) RHO1(KL) A(K L i j )
*          + SUM(M,K,L) I  (i M K L ) * RHO2( j M K L )
*
*. Input
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
*
      INCLUDE 'cintfo.inc'
*. Output
      DIMENSION F(*)
*
      NTEST = 000
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'FOO   ')
*
*. In the following we will fetch FI elements through GETH1
*  so copy inactive Fock matrix from  zero order space
      CALL MEMCHK2('MARK 1')
      CALL SWAPVE(WORK(KINT1),WORK(KFIZ),NINT1)
      CALL MEMCHK2('MARK 2')

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
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT ')
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN ')
*. And a block of F
      MX2IBLK = MXTSOB ** 2
      CALL MEMMAN(KLFBLK,MX2IBLK,'ADDL  ',2,'KLFBL ')
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
            CALL SETVEC(WORK(KLFBLK),ZERO,NI*NJ)
* 1 : One-electron part 
            DO KGAS = 1, NGAS
              KSM = IJSM
              NK = NOBPTS(KGAS,KSM)
*. blocks of one-electron integrals and one-electron density
              CALL GETD1(WORK(KLDEN),JSM,JGAS,KSM,KGAS,1)
              CALL GETH1(WORK(KLINT),ISM,IGAS,KSM,KGAS)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) 
     &          ' 1-e ints for ISM IGAS KGAS ',ISM,IGAS,KGAS
                CALL WRTMAT(WORK(KLINT),NI,NK,NI,NK)
                WRITE(6,*) 
     &          ' 1-e densi for ISM JGAS KGAS ',ISM,JGAS,KGAS
                CALL WRTMAT(WORK(KLDEN),NJ,NK,NJ,NK)
              END IF
*. And then a matrix multiply( they are pretty much in fashion 
*. these days )
              CALL MATML7(WORK(KLFBLK),WORK(KLINT),WORK(KLDEN),
     &                    NI,NJ,NI,NK,NJ,NK,ONE,ONE,2)
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated block '
                 CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
               END IF
 
            END DO
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One-electron contributions'
              WRITE(6,*) ' =========================='
              CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
            END IF
            IF(I12.EQ.2) THEN
*. 2 : Two-electron part 1
            FAC = 0.5d0
            DO KLSM = 1, NSMOB
              KSM = KLSM
              LSM = KLSM
*
              DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
               
                CALL GETD1(WORK(KLDEN),KSM,KGAS,LSM,LGAS,1)

                IXCHNG = 0
                ICOUL  = 1

*. Blocks of integrals for D1(K L) (K L ! I J) 
                CALL LUCIAGETINT(WORK(KLINT),
     &               KGAS,KSM,LGAS,LSM,IGAS,ISM,JGAS,JSM,
     &               IXCHNG,0,0,ICOUL,ONE,ONE)
 
                NIJ = NI*NJ
                NKL = NK*NL
                CALL MATML7(WORK(KLFBLK),WORK(KLDEN),WORK(KLINT),
     &                        1,NIJ,1,NKL,NKL,NIJ,1D0,4D0*FAC,0)
                
c                IF(NTEST.GE.1000) THEN
c                  WRITE(6,*) ' NIJ, NKL: ',NIJ, NKL
c                  WRITE(6,*) ' INTS 1'
c                  CALL WRTMAT(WORK(KLINT),NKL,NIJ,NKL,NIJ)
c                  WRITE(6,*) ' DEN 1'
c                  CALL WRTMAT(WORK(KLDEN),NK,NL,NK,NL)
c                  WRITE(6,*) ' Updated block 1'
c                  CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
c                END IF
*. Blocks of integrals for D1(K L) (K J ! I L) 
*.  i.e. get (KJ|IL) as XINT(KL,JI)                
                ICOUL = 2

                CALL LUCIAGETINT(WORK(KLINT),
     &               KGAS,KSM,JGAS,JSM,IGAS,ISM,LGAS,LSM,
     &               IXCHNG,0,0,ICOUL,ONE,ONE)

                NIJ = NI*NJ
                NKL = NK*NL
                CALL MATML7(WORK(KLFBLK),WORK(KLDEN),WORK(KLINT),
     &                        1,NIJ,1,NKL,NKL,NIJ,1D0,-1D0*FAC,0)

c                IF(NTEST.GE.1000) THEN
c                  WRITE(6,*) ' Updated block 2'
c                  CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
c                END IF
*. Blocks of integrals for D1(K L) (K I ! L J) 
*.  i.e. get (KI|LJ) as XINT(KL,IJ)                
                ICOUL = 0

                CALL LUCIAGETINT(WORK(KLINT),
     &               KGAS,KSM,IGAS,ISM,LGAS,LSM,JGAS,JSM,
     &               IXCHNG,0,0,ICOUL,ONE,ONE)

                NIJ = NI*NJ
                NKL = NK*NL
                CALL MATML7(WORK(KLFBLK),WORK(KLDEN),WORK(KLINT),
     &                        1,NIJ,1,NKL,NKL,NIJ,1D0,-1D0*FAC,0)

c                IF(NTEST.GE.1000) THEN
c                  WRITE(6,*) ' Updated block 3'
c                  CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
c                END IF

              END DO ! LGAS
              END DO ! KGAS
            END DO ! KLSYM
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One- + two-electron contributions 1'
              WRITE(6,*) ' ==================================='
              CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
            END IF
*. 2 : Two-electron part 2
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
                CALL LUCIAGETINT(WORK(KLINT),
     &               KGAS,KSM,LGAS,LSM,IGAS,ISM,MGAS,MSM,
     &               IXCHNG,0,0,ICOUL,ONE,ONE)
                
                CALL GETD2 (WORK(KLDEN),
     &               KSM,KGAS,LSM,LGAS,JSM,JGAS,MSM,MGAS,1)
                NKL = NK*NL
                DO M = 1, NM
                  IIOFF = KLINT + (M-1)*NKL*NI
                  IDOFF = KLDEN + (M-1)*NKL*NJ
                  CALL MATML7(WORK(KLFBLK),WORK(IIOFF),WORK(IDOFF),
     &                        NI,NJ,NKL,NI,NKL,NJ,ONE,ONE,1)
                END DO
              END DO
              END DO
              END DO
            END DO
            END DO
            END IF
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One- + two-electron contributions 2'
              WRITE(6,*) ' ==================================='
              CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
            END IF
*. Block has been constructed , transfer to -complete- 
*. symmetry blocked Fock matrix
            DO J = 1, NJ
              DO I = 1, NI
C?              WRITE(6,*) 'IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1',
C?   &                      IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1
                F(IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1 ) = 
     &          WORK(KLFBLK-1+(J-1)*NI+I)
              END DO
            END DO
*
          END DO
        END DO
      END DO
*
*. And clean up 
      CALL SWAPVE(WORK(KINT1),WORK(KFIZ),NINT1)

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

      SUBROUTINE CC_DENSI(ECC,IREFDEN,
     &           CCVEC1,CCVEC2,
     &           VEC1,VEC2,
     &           LUCCAMP,LU_OMG,LU_LAMP,LUREF,LUDIA,
     &           LU_1DEN, LU_2DEN,
     &           LU_1DENR, LU_2DENR) 
*
* Obtain CC density matrices 
*
*. 1 : Obtain lagrangian multipliers 
*. 2 : Obtain CI expansion corresponding to 
*      |Lambda > = |HF> + sum_{\mu} \bar t_{\mu}\tauop_{\mu}
*      Obtain density <Lambda!a+ a|CC> 
*
* Obtain vector of Lagrange multupliers by solving 
*
* CCLAG(T) * CC_JAC = - E1
*
* Where CC_JAC is the Jacobian and E1 is the 
* vector 
*          <HF!H tau_{\mu} exp(T) |HF>
*
* Jeppe Olsen, November 2000 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      CHARACTER*6 CCTYPE
*
      INCLUDE 'crun.inc'
c      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lorr.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'opti.inc'
*. Scratch
      DIMENSION CCVEC1(*),CCVEC2(*)
      DIMENSION VEC1(*),VEC2(*)
* TEST TEST TEST
      DIMENSION OCC(100)
*
      REAL(8) INPROD
*
      lblk = -1
      NTEST = 00
*
      ONE = 1.0D0
      ZERO = 0.0D0 
      ONEM = -1.0D0
*
      CCTYPE(1:6) = 'GEN_CC'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CC_DEN')
*
*. Local memory for linear solver
*
      LEN = MXITLE**2
      CALL MEMMAN(KLAPROJ,LEN,'ADDL  ',2,'APROJ ')
      LEN = 3*MXITLE**2 + MXITLE
      CALL MEMMAN(KLSCR,LEN,'ADDL  ',2,'LSCR  ')
      LEN = MXITLE
      CALL MEMMAN(KLRNRM,LEN,'ADDL  ',2,'RNRM  ')
*
*. Obtain Right hand side of Lagrangian equations
*
      lurhs = iopen_nus('CCDEN_RHS')
      luint = 0 ! for the moment: calculate it
      CALL ZERO_ORD_RHS(CCVEC1,VEC1,VEC2,LUCCAMP,LURHS,LUINT)
*
*. Solve zero-order equations, result in CCVEC1
*
      lusc1 = iopen_nus('CCDENSCR1')
      lusc2 = iopen_nus('CCDENSCR2')
      lusc3 = iopen_nus('CCDENSCR3')
      lusc4 = iopen_nus('CCDENSCR4')
      lusc5 = iopen_nus('CCDENSCR5')
      lusc6 = iopen_nus('CCDENSCR6')

*.Lambda vector corresponds to solution for the left side of 
*.the Jacobian. Tell the program to perform left side transformations
      L_OR_R = 1
      NINVEC = 0
      CALL GENMAT_LINSOL(L_OR_R,ECC,ELAM,
     &                   CCVEC1,CCVEC2,VEC1,VEC2,
     &                   LU_LAMP,LUCCAMP,LU_OMG,LURHS,LUDIA,
     &                   LUSC2,LUSC3,LUSC4,
     &                   WORK(KLRNRM),MXITLE,N_CC_AMP,NINVEC,
     &                   MXITLE,WORK(KLAPROJ),WORK(KLSCR),
     &                   0,THRGRD)
*
      call relunit(lurhs,'delete')
*
      ICSPC = IETSPC
      ISSPC = IETSPC
*. Construct Exp(T) |HF> on LUSC6
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-20

      call vec_from_disc(ccvec1,n_cc_amp,1,lblk,luccamp)
*.  Exp(t1) !ref> on LUSC6
      CALL EXPT_REF2(LUREF,LUSC6,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &               CCVEC1,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
*. Lambda on LUHC
*. |Lambda> =   exp(-T\dagger) ( !HF> + sum_{mu} L_{mu} \tau_{mu} |HF>)
C          SIG_GCC(C,HC,LUC,LUHC,T)
*. sum(mu) L(mu) tau(mu) |HF>
      call vec_from_disc(ccvec2,n_cc_amp,1,lblk,lu_lamp)
      CALL SIG_GCC(VEC1,VEC2,LUREF,LUSC2,CCVEC2)
*. Add |HF>
C          VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      CALL VECSMD(VEC1,VEC2,ONE,ONE,LUREF,LUSC2,LUSC1,1,-1)
*. and multiply with exp(-T\dag)
      CALL SCALVE(CCVEC1,ONEM,N_CC_AMP)
*. Conjugate amplitudes (reorder)  and operators
      CALL CONJ_CCAMP(CCVEC1,1,CCVEC2)
      CALL CONJ_T
*. Exp(-T\dag), result on LUSC2
      CALL EXPT_REF2(LUSC1,LUSC2,LUSC3,LUSC4,LUSC5,XCONV,MX_TERM,
     &             CCVEC2,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
      CALL CONJ_T
*.
      IF(NTEST.GE.1000) THEN  
        WRITE(6,*) ' Lambda vector : '
        CALL WRTVCD(VEC1,LUSC2,1,-1)
      END IF
*. And calculate density in standard CI way
      LBLK = -1
      XDUM = 0.0D0
      CALL DENSI2(IDENSI,WORK(KRHO1),WORK(KRHO2),
     &     VEC1,VEC2,LUSC2,LUSC6,EXPS2,
     &     0,XDUM,XDUM,XDUM,XDUM,1)
*. Symmtrize one-body density matrix
C     SYM_BLMAT(A,NBLOCK,LBLOCK)
      CALL SYM_BLMAT(WORK(KRHO1),1,NTOOB)
*. Symmetrize two-body density matrix
      IF (IDENSI.EQ.2)
     &     CALL SYM_2DENS(WORK(KRHO2),NTOOB,0)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Symmetrized one-body CC density ' 
        CALL APRBLM2(WORK(KRHO1),NTOOB,NTOOB,1,0) 
      END IF

      IF(IDENSI.EQ.2.AND.NTEST.GE.1000) THEN
        WRITE(6,*) ' Symmetrized two-body CC density ' 
        CALL PRINT_2DENS(WORK(KRHO2),NTOOB,0,1d-1)
      END IF


      LRHO1 = NTOOB**2
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2

      CALL VEC_TO_DISC(WORK(KRHO1),LRHO1,1,LBLK,LU_1DEN)
      IF (IDENSI.EQ.2)
     &     CALL VEC_TO_DISC(WORK(KRHO2),LRHO2,1,LBLK,LU_2DEN)

      call relunit(lusc1,'delete')
      call relunit(lusc2,'delete')
      call relunit(lusc3,'delete')
      call relunit(lusc4,'delete')
      call relunit(lusc5,'delete')
      call relunit(lusc6,'delete')
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CC_DEN')
*
      RETURN
      END
      SUBROUTINE ZERO_ORD_RHS(RHSVEC,VEC1,VEC2,LUCCAMP,LU_RHS,LU_EXRF)
*
* Jeppe + Jesper, October 2000
*
* Calculate RHS for zero order Lagrange multipliers
*
* -<HF!H tau_{\mu} exp(T) |HF>
*
*. Input
*
*  CC_AMP : CC amplitudes of current CC state
*
*. Output
*
*  RHSVEC : Right hand side vector as above
* 
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc' 
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cands.inc'
*     INCLUDE 'multd2h.inc'
      INCLUDE 'cstate.inc'
*. Input on
      INTEGER LUCCAMP
*. Output
      INTEGER LU_RHS
*. Scratch
      DIMENSION RHSVEC(*)
      DIMENSION VEC1(*), VEC2(*)
*
      CHARACTER*6 CCTYPE
*
      NTEST = 00
      IDUM = 1
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ZER_RH')
*. Works only for generalized cc
      CCTYPE(1:6) = 'GEN_CC'
*
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-20
*
*
*  <HF!H tau_{\mu} exp(T) |HF>
*
*  Calculate as <1| tau_nu |2>  
*
*  with |1> = H     |HF>
*       |2> = exp T |HF>
*
*. H |HF> on LU36
      ICSM = IREFSM
      ISSM = ICSM
      ICSPC = IREFSPC
      ISSPC = IETSPC
      ICC_EXC = 0
      CALL MV7(VEC1,VEC2,LUC,LUSC36,0,0)

      IF (LU_EXRF.GT.0) THEN
        LUEXPTREF = LU_EXRF
      ELSE
        LUEXPTREF = LUHC
*. exp(T)|HF> on LUHC
        ICC_EXC = 1
        LBLK = -1
        CALL VEC_FROM_DISC(RHSVEC,N_CC_AMP,1,LBLK,LUCCAMP)
        IOPTYP = 0
        ICSPC = IETSPC
        ISSPC = IETSPC
        CALL EXPT_REF2(LUC,LUHC,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &                 RHSVEC,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,IOPTYP)
C          DEN_GCC(C,HC,LUC,LUHC,T)
      END IF
      ICSPC=IETSPC
      ISSPC=IETSPC
      CALL DEN_GCC(VEC1,VEC2,LUEXPTREF,LUSC36,RHSVEC)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '************************************'
        WRITE(6,*) '  Zero-order RHS for CC Lagrangian '
        WRITE(6,*) '************************************'
        CALL WRT_CC_VEC2(RHSVEC,LU,CCTYPE)
      END IF

      CALL SCALVE(RHSVEC,-1d0,N_CC_AMP)
      LBLK = -1
      CALL VEC_TO_DISC(RHSVEC,N_CC_AMP,1,LBLK,LU_RHS)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ZER_RH')
*
      RETURN
      END 
      SUBROUTINE GENMAT_LINSOL(L_OR_R,ECC,ELAM,
     &           VEC1,VEC2,SCRVEC1,SCRVEC2,
     &           LU1,LUCCAMP,LU_OMG,LURHS,LUDIA,
     &           LU2,LU3,LU4,
     &           RNRM,MAXIT,NVAR,NINVEC,
     &           MAXVEC,APROJ,WORK,
     &           IPRT,TEST)
*
* Iterative linear equation solver for general matrix 
*
* Modified GENMAT_EIG
*
*
* MIN version requiring two vectors in core
*
*
* Jeppe Olsen Nov 2000, from GENMAT_EIG
*
* Procedure 
*
*.   Loop over iterations 
*      If first iteration 
*        Orthogonalize vectors
*      end if
*      Find matrix in subspace :
*       Calculate matrix times new vectors  
*       augment subspace matrix 
*      Find rhs in subspace
*      Solve eqs in subspace 
*      Find solution vector and residual in full space
*      If (residual is larger than test and more its to come)
*         Obtain new direction 
*         orthogonalize direction to prev vectors
*      End if
*    End of loop over iterations
*
*
* Input :
* =======
*        LU1 : Can contain initial approximation to linear equations
*              Zero vector must not be supplied as all vectors 
*              are orthonormalized 
*        LURHS : Right hand side of eqs    
*        VEC1,VEC2 : Two vectors, each must be dimensioned to hold
*                    complete vector
*        LU2       : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NINVEC    : Number of initial vectors on LU1
*                    may be equal to zero
*        MAXVEC    : Largest allowed number of vectors
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVER
      DIMENSION VEC1(*),VEC2(*)
      REAL * 8   INPROD
*
      DIMENSION RNRM(MAXIT)
*
*. Scratch through argument list
*
*. Length of APROJ : MAXVEC**2
*  Length of WORK  : 3*MAXVEC**2 + MAXVEC
      DIMENSION APROJ(*),WORK(*)
*
      INTEGER KPIV(MAXVEC)
*
      NTEST = 00

      CALL ATIM(CPU0,WALL0)
*. Constants 
*
      ZERO = 0.0D0
      ONE  = 1.0D0
      ONEM =-1.0D0
*
      WRITE(6,*) ' Convergence threshold for residual : ', TEST
*. Division of scratch memory 
      KS = 1
      KFREE = KS + MAXVEC**2
*
      KX = KFREE
      KFREE = KX + MAXVEC**2
*
      KSCR = KFREE
      KFREE = KSCR + MAXVEC**2
*
      KRHS = KFREE 
      KFREE = KRHS + MAXVEC
*
      CONVER = .FALSE.
*. Scratch files
      LUINTM1 = -IOPEN_NUS('CCGMATSCR1')
      LUINTM2 = -IOPEN_NUS('CCGMATSCR2')
      LUINTM3 = IOPEN_NUS('CCGMATSCR3')
*.
      WRITE(6,'(">>>",A)')
     &     '  Iter.        lambda energy     res. norm     norm(L)'
      WRITE(6,'(">>>",A)')
     &     '-------------------------------------------------------'

*. Start of loop over iterations 
      NITER = 0
      DO ITER = 1, MAXIT
*
*. Orthogonalize vectors on LU1 and save vectors on LU1
*
        IF(ITER.EQ.1) THEN 
*. Orthogonalize all input vectors
          NVEC = NINVEC
          NORT = 0
        ELSE
*. Orthogonalize last vector to all previous vectors
          NORT = NVEC-1
        END IF
*
        IF(NVEC.GT.0) THEN
        
C                         ORTH_VEC(LUVEC,NVEC,NORT,VEC1,VEC2,LUSCR,LUOUT,
C    &                    S,X,SCR,ICOPY)
          IF(NVEC.GT.0) CALL ORTH_VEC(LU1,NVEC,NORT,NVAR,VEC1,VEC2,
     &                  LU3,LU4,WORK(KS),WORK(KX),WORK(KSCR),1)
*
*. Augment Projected matrix 
*
          IF(ITER.EQ.1) THEN
            NADD = NINVEC
            ISTART = 1
          ELSE
            NADD = 1
            ISTART = NVEC
          END IF
          NPREV = NVEC-NADD
*. Update Matrix * Vector products on LU2
C            SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
          CALL SKPVCD2(LU1,ISTART-1,NVAR,VEC1,1,-1)
          CALL SKPVCD2(LU2,ISTART-1,NVAR,VEC1,1,-1)
C       CALL SKPRC3(ISTART-1,LU1)
C       CALL SKPRC3(ISTART-1,LU2)
C?      WRITE(6,*) ' ITER,ISTART,NVEC=', ITER,ISTART,NVEC
          DO IVEC = ISTART,NVEC
C         CALL FRMDSC(VEC1,NVAR,-1,LUVEC,IMZERO,IAMPACK)
            IWRMOD=1
            XDUM = -1
            LUDUM = 0
            CALL JAC_T_VEC2(L_OR_R,0,IWRMOD,0,0,
     &           VEC1,VEC2,SCRVEC1,SCRVEC2,
     &           NVAR,NVAR,
     &           XDUM,XNORM,RNORM,
     &           LUCCAMP,LUDUM,LU1,LU2,LUDUM,
     &           LUINTM1,LUINTM2,LUINTM3)
            IF (LUINTM1.LT.0) LUINTM1 = ABS(LUINTM1)
            IF (LUINTM2.LT.0) LUINTM2 = ABS(LUINTM2)
c          CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
c          CALL MV10(VEC1,VEC2)
c          CALL TODSC(VEC2,NVAR,-1,LU2)

          END DO
*. Asub(I,Jnew)
          CALL REWINO(LU1)
C       CALL SKPRC3(ISTART-1,LU2)
          CALL SKPVCD2(LU2,ISTART-1,NVAR,VEC1,1,-1)
C  GET_NONSYM_SUBMAT(LUC,LUHC,NVEC,NVECP,SUBMAT,NDIM,
C    &                             VEC1,VEC2,SUBMATP)
          IF(NPREV.GT.0) CALL COPVEC(APROJ,WORK(KSCR),NPREV**2)

          CALL GET_NONSYM_SUBMAT(LU1,LU2,NVEC,NPREV,APROJ,NVAR,
     &                          VEC1,VEC2,WORK(KSCR))
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Updated subspace matrix '
            CALL WRTMAT(APROJ,NVEC,NVEC,NVEC,NVEC)
          END IF
*. project right hand vector on subspace 
c          CALL REWINO(LURHS)
C       CALL FRMDSC(VEC1,NVAR,-1,LUVEC,IMZERO,IAMPACK)
c          CALL FRMDSC(VEC2,NVAR,-1,LURHS,IMZERO,IAMPACK)
          CALL VEC_FROM_DISC(VEC2,NVAR,1,-1,LURHS)
          CALL SKPVCD2(LU1,NPREV,NVAR,VEC1,1,-1)
          DO IVEC = NPREV+1,NVEC
            CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
            WORK(KRHS-1+IVEC) = INPROD(VEC1,VEC2,NVAR)
c            WORK(KRHS-1+IVEC) = -INPROD(VEC1,VEC2,NVAR)
          END DO
          IF (NTEST.GE.100) THEN
            WRITE(6,*) ' Updated RHS in subspace'
            CALL WRTMAT(WORK(KRHS),NVEC,1,NVEC,1)
          END IF
*. Solve linear equations in subspace by inverting and multiplying
*. Save in WORK(KSCR)
          CALL COPVEC(APROJ,WORK(KX),NVEC**2)
c          CALL INVMAT(WORK(KX),WORK(KSCR),NVEC,NVEC,ISING)
c          CALL MATVCB(WORK(KX),WORK(KRHS),WORK(KSCR),NVEC,NVEC,0)

          ! factorize
          CALL DGECO(WORK(KX),NVEC,NVEC,KPIV,COND,WORK(KSCR))
          IF(NTEST.GE.150) THEN
            WRITE(6,*) ' Factorized subspace matrix '
            CALL WRTMAT(WORK(KX),NVEC,NVEC,NVEC,NVEC)
          END IF
          CALL COPVEC(WORK(KRHS),WORK(KSCR),NVEC)
          JOB = 0
          CALL DGESL(WORK(KX),NVEC,NVEC,KPIV,WORK(KSCR),JOB)
          IF (NTEST.GE.100) THEN
            WRITE(6,*) ' Solution for x in subspace (Ax = rhs)'
            CALL WRTMAT(WORK(KSCR),NVEC,1,NVEC,1)
          END IF
*
*. Matrix times solution vector and New residual 
*
          CALL SETVEC(VEC1,ZERO,NVAR)
          CALL REWINO(LU2)
          DO IVEC = 1, NVEC
            CALL FRMDSC(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK) 
            FAC = WORK(KSCR-1+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,ONE,FAC,NVAR)
          END DO
c          CALL REWINO(LURHS)
c          CALL FRMDSC(VEC2,NVAR,-1,LURHS,IMZERO,IAMPACK)
          CALL VEC_FROM_DISC(VEC2,NVAR,1,-1,LURHS)
c          ONEM = -1.0D0
          CALL VECSUM(VEC2,VEC2,VEC1,1D0,-1D0,NVAR)
c          CALL VECSUM(VEC2,VEC2,VEC1,1D0,1D0,NVAR)
          RNORM = INPROD(VEC2,VEC2,NVAR)
          NITER = NITER + 1
          RNRM(NITER) = SQRT(RNORM)
          WRITE(6,'(">>>",I5,25X,(2X,E10.4))') NITER,SQRT(RNORM)
          IF(SQRT(RNORM).GT.TEST) THEN
*. Not converged
            CONVER = .FALSE.
          ELSE
            CONVER = .TRUE.
          END IF
*
          IF(CONVER.OR.ITER.EQ.MAXIT.OR.NVEC.EQ.MAXVEC) THEN
*. Reduce subspace to single vector being current solution vector
*. First write Matrix * current solution to disc as we 
*. already have this vector in VEC1
            CALL REWINO(LU3)
            CALL TODSC(VEC1,NVAR,-1,LU3)
*. and the new solution vector
            CALL SETVEC(VEC1,ZERO,NVAR)
            CALL REWINO(LU1)
            DO IVEC = 1, NVEC
              CALL FRMDSC(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK) 
              FAC = WORK(KSCR-1+IVEC)
              CALL VECSUM(VEC1,VEC1,VEC2,ONE,FAC,NVAR)
            END DO
            IF (.NOT.CONVER) THEN
              CALL REWINO(LU1)
              CALL TODSC(VEC1,NVAR,-1,LU1)
            ELSE ! we need it in this form later:
              CALL VEC_TO_DISC(VEC1,NVAR,1,-1,LU1)
            END IF
            XLNRM = SQRT(INPROD(VEC1,VEC1,NVAR))
            IF (LU_OMG.GT.0) THEN
              CALL VEC_FROM_DISC(VEC2,NVAR,1,-1,LU_OMG)
              ELAM = ECC + INPROD(VEC1,VEC2,NVAR)
              WRITE(6,'(">>>",5X,F25.12,2(2X,E10.4))')
     &           ELAM,SQRT(RNORM),XLNRM
            ELSE
              WRITE(6,'(">>>",5X,25X,2(2X,E10.4))')
     &           SQRT(RNORM),XLNRM
            END IF
            NVEC = 1
*. Obtain Ax in vec1 for later use
            CALL REWINO(LU3)
            CALL FRMDSC(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
          END IF
*
        ELSE
          CONVER = .FALSE.
          CALL SETVEC(VEC1,ZERO,NVAR)
        END IF
*       ^ End if there is vectors in subspace
        IF(.NOT.CONVER) THEN
*. Find new trial vector  : Diag**(-1)(b-ax)
c          CALL REWINO(LURHS)
c          CALL FRMDSC(VEC2,NVAR,-1,LURHS,IMZERO,IAMPACK) 
          CALL VEC_FROM_DISC(VEC2,NVAR,1,-1,LURHS)
*. B - AX in VEC2
          CALL VECSUM(VEC2,VEC2,VEC1,1D0,-1D0,NVAR)
c          CALL VECSUM(VEC2,VEC2,VEC1,1D0,1D0,NVAR)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' New residual vector '
            CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
          END IF
*. D**(-1) (B-AX) in VEC2
c          CALL REWINO(LUDIA)
c          CALL FRMDSC(VEC1,NVAR,-1,LUDIA,IMZERO,IAMPACK) 
          CALL VEC_FROM_DISC(VEC1,NVAR,1,-1,LUDIA)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Diagonal read in '
            CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
          END IF
C              DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
          CALL DIAVC2(VEC2,VEC2,VEC1,ZERO,NVAR) 
          NVEC = NVEC + 1
*. And save on LU1
          CALL REWINO(LU1)
C         CALL SKPRC3(NVEC-1,LU1)
          CALL SKPVCD2(LU1,NVEC-1,NVAR,VEC1,1,-1)
          CALL TODSC(VEC2,NVAR,-1,LU1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' New raw direction '
            CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
          END IF
        END IF ! .NOT.CONVER
        IF(CONVER) GOTO 1001
      END DO
*
 1001 CONTINUE
c*. Obtain solution in VEC1
c      CALL REWINO(LU1)
c      CALL FRMDSC(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
*
      IF(NTEST .GT. 0 ) THEN
*
        IF(CONVER) THEN
          WRITE(6,1010) NITER  ,RNRM(NITER)
 1010  FORMAT(1H0,'  convergence was obtained in...',I3,' iterations',/,
     +        1H ,'  norm of residual..............',E13.8)
        ELSE
          WRITE(6,1020) MAXIT ,RNRM(MAXIT)
 1020  FORMAT(1H0,' convergence was not obtained in',I3,'iterations',/,
     +        1H ,' norm of residual...............',E13.8)
        END IF
*
      END IF
*
      CALL RELUNIT(ABS(LUINTM1),'DELETE')
      CALL RELUNIT(ABS(LUINTM2),'DELETE')
      CALL RELUNIT(LUINTM3,'DELETE')
*
      IF(NTEST.GE. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(1H0,' solution to set of linear equations')
       CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
      END IF
C
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(1H0,10X,'iteration point     norm of residual')
      DO 350 I=1,NITER
       WRITE(6,1050)I,RNRM(I)
 1050  FORMAT(1H ,12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
*
      CALL ATIM(CPU,WALL)
      CALL PRTIM(6,'time in linear solver',CPU-CPU0,WALL-WALL0)
*
      RETURN
      END
      SUBROUTINE OVLAP_VEC(LUVEC,NVEC,NVAR,ISTART,VEC1,VEC2,OVLAP)
*
* NVEC vectors are residing on LUVEC
* Find overlap between vectors ISTART-NVEC and all vectors
* 
* S(i,j) = <Vec(i)!Vec(j)>, i=1, NVEC, j= ISTART,NVEC
*
* Jeppe Olsen, Korshojen, November 2000
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Scratch 
      DIMENSION VEC1(*),VEC2(*)
*. Output
      DIMENSION OVLAP(NVEC,NVEC-ISTART+1)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' OVLAP_VEC speaking '
        WRITE(6,*) ' NVEC,NVAR,ISTART', NVEC,NVAR,ISTART
      END IF
*
      DO J = ISTART,NVEC
*. Read vector J in VEC1
C            SKPVCD2(LU,NVEC,NVAR,SEGMNT,IREW,LBLK) 
        CALL SKPVCD2(LUVEC,J-1,NVAR,VEC1,1,-1)
        CALL FRMDSC(VEC1,NVAR,-1,LUVEC,IMZERO,IAMPACK)
C?      WRITE(6,*) ' J vector read in '
C?      CALL WRTMAT(VEC1,1,NVAR,1,NVAR)
        CALL REWINO(LUVEC)
        DO I = 1, NVEC
          CALL FRMDSC(VEC2,NVAR,-1,LUVEC,IMZERO,IAMPACK)
C?        WRITE(6,*) ' I vector read in '
C?        CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
          OVLAP(I,J-ISTART+1)  = INPROD(VEC1,VEC2,NVAR)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap matrix from OVLAP_VEC '
        NVECACT = NVEC-ISTART+1
        CALL WRTMAT(OVLAP,NVEC,NVECACT,NVEC,NVECACT)
      END IF
*
      RETURN
      END 
      SUBROUTINE ORTH_VEC(LUVEC,NVEC,NORT,NVAR,VEC1,VEC2,LUSCR,LUOUT,
     &                    S,X,SCR,ICOPY)
*
* A set of NVEC are residing on file LUVEC.
* The first NORT vectors are orthogonal
* ortogonalize the remaining NVEC-NORT vectors to the 
* first NORT vectors and to each other 
*
* Jeppe Olsen, Korshojen, November 2000
*
      INCLUDE 'implicit.inc'
*. Scratch vectors, must each hold vector one vector from disc
      DIMENSION VEC1(*),VEC2(*)
*. Scratch matrices
      DIMENSION X(NVEC,NVEC),S(NVEC,NVEC),SCR(NVEC,NVEC)
*
      NTEST = 00
*. Set initial overlap matrix to unit matrix 
      ZERO = 0.0D0
      CALL SETVEC(S,ZERO,NVEC**2)
C SETDIA(MATRIX,VALUE,NDIM,IPACK)
      ONE = 1.0D0
      CALL SETDIA(S,ONE,NVEC,0)
* Obtain in SCR overlap matrix between vectors NORT+1-NVEC and all vectors
C     OVLAP_VEC(LUVEC,NVEC,NVAR,ISTART,VEC1,VEC2,OVLAP)
      CALL OVLAP_VEC(LUVEC,NVEC,NVAR,NORT+1,VEC1,VEC2,S(1,NORT+1))
*. Complete overlap matrix 
      DO I = NORT+1,NVEC
        DO J = 1, NORT
          S(I,J) = S(J,I)
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Overlap matrix '
        WRITE(6,*) ' ============== '
        WRITE(6,*)
        CALL WRTMAT(S,NVEC,NVEC,NVEC,NVEC)
      END IF
*
* Find orthonormalization  matrix 
*
C  MGS3(X,S,NDIM,SCR1)
      CALL MGS3(X,S,NVEC,SCR)
*. And orthonormalize 
*. The first NORT vectors are easy : just copy
      CALL REWINO(LUVEC)
      CALL REWINO(LUOUT)
      DO IVEC = 1, NORT
        CALL FRMDSC(VEC1,NVAR,-1,LUVEC,IMZERO,IAMPACK)
        CALL TODSC(VEC1,NVAR,-1,LUOUT)
      END DO
*. And vectors NORT+1 - NVEC
      DO IVEC = NORT+1,NVEC
        CALL REWINO( LUVEC)
        CALL SETVEC(VEC1,ZERO,NVAR)
        DO JVEC = 1, NVEC
          CALL FRMDSC(VEC2,NVAR,-1  ,LUVEC,IMZERO,IAMPACK)
          FACTOR =  X(JVEC,IVEC)
          CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
        END DO
        CALL TODSC(VEC1,NVAR,-1  ,LUOUT)
      END DO
*
      IF(ICOPY.EQ.1) THEN
*. Copy orthonormalized vectors over LUVEC
        CALL REWINO(LUVEC)
        CALL REWINO(LUOUT)
        DO IVEC = 1, NVEC
          CALL FRMDSC(VEC1,NVAR,-1  ,LUOUT,IMZERO,IAMPACK)
          CALL TODSC(VEC1,NVAR,-1  ,LUVEC)
        END DO
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 
     &  ' Information from ORTH_VEC : orthonormalized vectors'
        IF(ICOPY.EQ.1) THEN
           LUX = LUVEC
        ELSE 
           LUX = LUOUT
        END IF
        CALL REWINO(LUX)
        DO IVEC = 1, NVEC
          CALL FRMDSC(VEC1,NVAR,-1,LUX,IMZERO,IAMPACK)
          CALL WRTMAT(VEC1,NVAR,1,NVAR,1)
        END DO
      END IF
*
      RETURN
      END
c $Id$
