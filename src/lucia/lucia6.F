*
* lucia6.f
*
* ====================
* CI response routines lucia version
* ====================
*
* Arbitrary order perturbation theory for one or two perturbations
* Cauchy moments, effective spectra , Dispersion coefficients
*
* Jeppe Olsen
*
*
      SUBROUTINE GET_ALPHA(ALPHA,W,VEC1,VEC2)
*
* Obtain polarizability tensor for frequency W
*
* Jeppe Olsen, April 2004
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cands.inc'
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*)
*
      RETURN
      END


      SUBROUTINE SYM_FOR_OP(OPER,IXYZ,ISYM)
*
* Symmetry of perturbation operator OPER
* obtained from XYZ characters in start of OPER
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific input
      CHARACTER*8 OPER
*. General input
      INTEGER IXYZ(3)
      INCLUDE 'multd2h.inc'
*
*. Number of cartesian components
      NCART=0
      DO ICHAR = 1, 8
        IF(OPER(ICHAR:ICHAR).EQ.'X'  .OR.
     &     OPER(ICHAR:ICHAR).EQ.'Y'  .OR.
     &     OPER(ICHAR:ICHAR).EQ.'Z'      ) THEN
             NCART = NCART + 1
        ELSE
*. End of cartesian part 
          GOTO 1001
        END IF
      END DO
 1001 CONTINUE
*
      ISYM = 1
      DO ICART = 1, NCART
        IF(OPER(ICART:ICART).EQ.'X') ISYM = MULTD2H(ISYM,IXYZ(1))
        IF(OPER(ICART:ICART).EQ.'Y') ISYM = MULTD2H(ISYM,IXYZ(2))
        IF(OPER(ICART:ICART).EQ.'Z') ISYM = MULTD2H(ISYM,IXYZ(3))
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,'(A,A,A,I1)') ' Label ', OPER, ' has symmetry ', ISYM
      END IF
*
      RETURN
      END
      SUBROUTINE GET_NMLW(VEC,N,M,LW,N_MAX,M_MAX,
     &           LUCORR,LUREF,IREAD,LUOUT,IREW,LBLK,LUOUTEFF)
*
* Obtain correction vector ! N,M,LW> obtained in GNDBTFREQ
* and save on file LUOUT
* If IREAD = 0, the file is positioned 
* for reading the file, but no transfer is realized
* If IREW ne 0, LUOUT is rewinded
*
*. Jeppe in Pisa, July 17, 97
*                 Sept. 97, Changed to disc, ICISTR = 2,3
*
*                 
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)' GET_NMLW : Vector in Demand, N,M,LW :', N,M,LW
      END IF
*
      LUOUTEFF = LUCORR
      IF(IREW.NE.0) CALL REWINO(LUOUT)
*
      IF(N.EQ.0.AND.M.EQ.0) THEN
*. You asked for the reference vector and that will be
        LUOUTEFF = LUREF
        CALL REWINO(LUREF)
        IVEC = 1
        IF(IREAD.NE.0 )  THEN
          CALL COPVCDP(LUREF,LUOUT,VEC,1,LBLK)
        END IF
      ELSE
*. A real correction vector
*. Position to start of requested vector
        CALL REWINO(LUCORR)  
        IVEC = 0
        DO NP = 0, N
          IF(NP.EQ.N) THEN
           MXMP = M
          ELSE
           MXMP = M_MAX
          END IF
          DO  MP = 0, MXMP
           IF(NP.EQ.N.AND.MP.EQ.M)THEN
             MXLWP = LW-1
           ELSE
             MXLWP = NP
           END IF
           DO LWP = -NP, MXLWP,2
C?           WRITE(6,*) ' Next vector NP,MP,LWP',NP,MP,LWP
*. Skip
             IF(.NOT.(NP.EQ.0.AND.MP.EQ.0)) THEN
C?            WRITE(6,*) ' SKPVCD will be called '
              CALL SKPVCD(LUCORR,1,VEC,0,LBLK)
              IVEC = IVEC + 1
C?            write(6,*) ' Vector skipped ', IVEC
             END IF
           END DO
          END DO
        END DO
*. And then we are ready to read the vector of interest
        IF(IREAD.NE.0) THEN
          CALL COPVCDP(LUCORR,LUOUT,VEC,0,LBLK)
          IVEC = IVEC + 1
        END IF
      END IF
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Vector to be read '
        WRITE(6,*) ' N M LW =',N,M,LW
        WRITE(6,*) ' obtained as vector ',IVEC
      END IF
*
      IF(NTEST.GE.100) THEN
       IF(IREAD.NE.0) THEN 
         IF(IREW.NE.0) THEN
           WRITE(6,*) ' Vector read in '
           CALL WRTVCD(VEC,LUOUT,1,LBLK)
         ELSE
           WRITE(6,*) ' Vector is not first vector'
           WRITE(6,*) ' so it will not be printed '
         END IF
       ELSE
         WRITE(6,*) 'vector positioned at end of vector',IVEC
       END IF
      END IF
*
      RETURN
      END
      SUBROUTINE GNDBPTFREQ(V,W,IRFSM,ERF,LURF,  !LUHI0,????
     &           IVRFSM,MAXV,IWRFSM,MAXW,IDRFSM,IVSM,IWSM,
     &           LUN,VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           LUDIA,FREQ,ENM,
     &           N_AVE_OP,AVE_OP,AVINT,LAVINT,IAVE_SYM,IPRNT)
*
* Solve the double perturbation expansion arising from a 
* static perturbation(V) and a frequency dependent perturbation(W).
*
* Susceptibilities corresponding to expectation values of 
* the N_AVE_OP operators defined by AVINT are also constructed
* 
* The perturbations or the squared perturbations are assumed 
* to be symmetric
*
* Jeppe Olsen, adapted from LUCAS codes
* The wave function corrections reads
*
* d(n,m,lw) = -(e[2] - lw I)-1 *
*      ( Vw[2] d(n-1,m,(l-1)w)
*       +Vo[2] d(n,m-1,kw)
*       +V-w[2]d(n-1,m,(l+1)w)
*       -sum(ld=1,n;ls=1,m-1;ll)Vo[1]t  d(ld,ls,lw)d(n-ld,m-ls-1,(l-ll)w)
*       -sum(ld=1,n-1;ls=1,m;ll)Vw[1]t  d(ld,ls,lw)d(n-ld-1,m-ls,(l-ll-1)w)
*       -sum(ld=1,n.1;ls=1,m;ll)V-w[1]t d(ld,ls,lw)d(n-ld-1,m-ls,(l-ll+1)w)
*
* Where 
* n : is the order in the frequency dependent expansion
* m : is the order in tha static expansion
* w : is the frequency (FREQ)
*
* 
*
* Input
* =====
* V : one-electron integrals defining  V
* W : one-electron integrals defining  W
* IRFSM : symmetry of reference vector
* LURF : file containing reference vector
* IVRFSM : symmetry of V times reference vector
* IWRFSM : symmetry of W times reference vector
* IDRFSM : symmetry of VW times reference vector
* LUN : file number for file to contain perturbation vectors
* ERF   : Reference energy
* MAXV : Order in V through which the equations should be solved
* MAXW : Order in W through which the equations should be solved
* VEC1,VEC2 : Scratch vectors , should be able to hold a complete vector
* LU1, LU2, LU3, LU4, LU5, LU6, LU7  : scratch files
*
* N_AVE_OP : Number of operators for which expectation values are expanded
*            (A- operators)
* AVE_OP   : Labels of above operators
* AVINT    : Integrals of A operators AVINT(*,IAV) contains integrals
*            for A operator number IAV
* LAVINT   : Row length of AVINT
* IAVE_SYM : Symmetry of A operators
*
* FREQ
*
* Output
* ======
* LUN : contains the MAXORD perturbation vectors
*       Correction nm is stored as record (n-1)*MAXV + m
* ENM : Contains the energy corrections
*
* Internal links
* ===============
* The actual solutions of linear eqs and the multiplicartion
* of V times a vector is realized through calls to BVEC and HMWITV
* /CANDS/ and /OPER( is used to communicate to these routines
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL * 8 INPRDD
      DIMENSION V(*),VEC1(*),VEC2(*)
      DIMENSION ENM(MAXV+1,MAXW+1)
*
      DIMENSION ISMVW(2,2)
      PARAMETER(MAXORD = 10)
      DIMENSION V01D(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
      DIMENSION VW1D(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. Inner products <d!d>
      DIMENSION DD(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. Inner products <d!W!d>
      DIMENSION DWD(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. And susceptibilities
      DIMENSION CHI(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. A-operators
      CHARACTER*8 AVE_OP(N_AVE_OP) 
      DIMENSION  AVINT(LAVINT,N_AVE_OP),IAVE_SYM(N_AVE_OP)
C    &           N_AVE_OP,AVE_OP,AVINT,LAVINT,IAVE_SYM)
*  V ** n W ** m !ref> : Symmetry is  ISMVW(MOD(n,2)+1,MOD(m,2)+1)
      INCLUDE 'cands.inc'
*. Local scratch 
      DIMENSION SCR(10*MAXORD)
      INCLUDE 'oper.inc'
      INCLUDE 'multd2h.inc'
*. Transfer of shift 
      INCLUDE 'cshift.inc'
*. Transfer of zero order energy
      COMMON/CENOT/E0
*. Still testing after all these years
      NTESTL =   1  
      WRITE(6,*) ' Input Print flag ', IPRNT
      NTEST = MAX(NTESTL,IPRNT)
*. redefine LUN
      LUN = 41
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' LU1, LU2, LU3, LU4, LU5, LU6, LU7 '
        WRITE(6,*)   LU1, LU2, LU3, LU4, LU5, LU6, LU7  
        WRITE(6,*) ' LUN LURF ', LUN,LURF
      END IF
*
      X = 99999999
*. I think there is a problem with V01D, and VW1D, so initialize with outrageous values 
      CALL SETVEC(V01D,X,(MAXORD+1)*(MAXORD+1)*(2*MAXORD+1))
      CALL SETVEC(VW1D,X,(MAXORD+1)*(MAXORD+1)*(2*MAXORD+1))
C     WRITE(6,*) ' FREQ = ', FREQ
*
      E0 = 0.0D0
      LBLK = -1
*
      NNNVEC = 0
      ISMVW(1,1) = IRFSM
      ISMVW(2,1) = IVRFSM
      ISMVW(1,2) = IWRFSM
      ISMVW(2,2) = IDRFSM
*
      IF(NTEST.GE.5) THEN 
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*************************************************'
      WRITE(6,*) '*                                               *'
      WRITE(6,*) '* Welcome to GNDBPTFREQ                         *'
      WRITE(6,*) '* General double perturbation calculations      *'
      WRITE(6,*) '* Involving a static and a dynamic perturbation *'
      WRITE(6,*) '*                                               *'
      WRITE(6,*) '*                      Best wishes from         *'
      WRITE(6,*) '*                      Antonio and Jeppe        *'
      WRITE(6,*) '*                                               *'
      WRITE(6,*) '*************************************************'
      WRITE(6,*)
      WRITE(6,*) 
      END IF
C     WRITE(6,*)  ' LUN and LURF ',LUN,LURF
 
*.  Real perturbation
      IRC = 1
*. For a imaginary pert set IRC to 2
*
*
*. The vectors P V0!0> and P Vw!0> ( P = 1 - !0><0!) 
*  will be used a number of times in the
*. following. Construct these guys and store on LU5
*
*
      CALL REWINO(LU5)
*
* V !0> as first vector on LU5
*
      ICSM = IRFSM
      ISSM = IVRFSM
      IPERTOP = 0
      IAPR = 0
C          BVEC(B,IBSM,LUC,LUB,VEC1,VEC2)
      CALL BVEC(V,IVSM,LURF,LU5,VEC1,VEC2)
      IF(IVSM.NE.1) THEN
       V00 = 0.0D0
      ELSE
       V00 = INPRDD(VEC1,VEC2,LU5,LURF,1,LBLK)
       V1NORM = INPRDD(VEC1,VEC2,LU5,LU5,1,LBLK)
       IF(NTEST.GE.10) THEN
         WRITE(6,*) ' V1NORM = ', V1NORM
         WRITE(6,*) ' V00 = ', V00
       END IF
       ONE = 1.0D0
       FACTOR = -V00
       CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU5,LURF,LU4,1,LBLK)
       CALL COPVCDP(LU4,LU5,VEC1,1,LBLK)
      END IF
*
* W !0> as second vector on LU5
*
C?    WRITE(6,*) ' IWSM = ', IWSM
      ISSM  = IWRFSM
      CALL BVEC(W,IWSM,LURF,LU4,VEC1,VEC2)
      IF(IWSM.NE.1) THEN
       W00 = 0.0D0
       CALL REWINO(LU4)
       CALL COPVCDP(LU4,LU5,VEC1,0,LBLK)
      ELSE
       W00 = INPRDD(VEC1,VEC2,LU4,LURF,1,LBLK)
C?     WRITE(6,*) ' W00 = ', W00
       ONE = 1.0D0
       FACTOR = -W00
       CALL REWINO(LU4)
       CALL REWINO(LURF)
       CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU4,LURF,LU5,0,LBLK)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The two vectors on LU5 '
        CALL REWINO(LU5)
        DO IVEC = 1, 2
           WRITE(6,*) ' Vector ', IVEC
           CALL WRTVCD(VEC1,LU5,0,LBLK) 
        END DO
      END IF
 
*
*. A note on the organization of the correction vectors on disc :
*  The vectors are stored as Loop over freq dep op, Loop over static op,
*  Loop over freq. components
 
      DO N = 0, MAXV
*
       IF(N.EQ.0) THEN
        MMIN = 1
       ELSE
        MMIN = 0
       END IF
*
       DO M = MMIN, MAXW
*. Symmetry of this correction  vector
         MNSM  = ISMVW(MOD(N,2)+1,MOD(M,2)+1)
*. Frequency components
         DO LW = -N,N,2
*
          WRITE(6,*)
          WRITE(6,*) ' Correction vector will be solved for '
          WRITE(6,*)
          WRITE(6,*) '    Order in dynamic perturbation ', N 
          WRITE(6,*) '    Order in static  perturbation ', M 
          WRITE(6,*) '    Frequency component           ', LW
          WRITE(6,*) '    Symmetry of correction vector ', MNSM
          WRITE(6,*)  
*
          ISSM = ISMVW(MOD(N,2)+1,MOD(M,2)+1)

* =========================================================================
* 1 : construct (Vw-Vw00)!N-1,M,LW-1>,(Vw-Vw00)!N-1,M,LW+1>,(W-W00)!N,M-1,LW> ( and save on LU4)
* =========================================================================
          CALL REWINO(LU4)
*. Components !N-1,M,LW-1>, !N-1,M,LW+1> corresponds to allowed freq.
          IF(N.GE.1.AND.(ABS(LW-1).LE.N-1)) THEN
            IMIN = 1
          ELSE
            IMIN = 0
          END IF
*
          IF(N.GE.1.AND.(ABS(LW+1).LE.N-1)) THEN
            IPLUS= 1
          ELSE
            IPLUS= 0
          END IF
C?        WRITE(6,*) ' N LW IPLUS',N,LW,IPLUS
          ICSM = ISMVW(MOD(N-1,2)+1,MOD(M,2)+1)
*
* (Vw - <0!Vw!0>) !N-1,M,LW-1>( on LU4)
*
          IF(IMIN.EQ.1) THEN
            IF(NTEST.GE.100) WRITE(6,*) ' IMIN active '
* : Place !0(N-1,M,LW-1)> on LU6 
            CALL GET_NMLW(VEC1,N-1,M,LW-1,MAXV,MAXW,LUN,LURF,1,
     &                    LU6,1,LBLK,LUOUTEFF)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Correction vector !0(N-1,M,LW-1)> read '
            END IF
            IF(NTEST.GE.100) THEN
              CALL WRTVCD(VEC1,LU6,1,LBLK) 
            END IF
            I12 = 1
            CALL BVEC(V,IVSM,LU6,LU7,VEC1,VEC2)
*. - <0!Vw!0> !N-1,M,LW-1>
C           IF(V00.NE.0.0D0 .AND. N+M.NE.1) THEN
            IF(V00.NE.0.0D0 ) THEN
              ONE = 1.0D0
              FACTOR = -V00
              CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU7,LU6,LU4,1,LBLK)
            ELSE
              CALL COPVCDP(LU7,LU4,VEC1,1,LBLK)
            END IF
*
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' V times Correction vector on LU4 '
              CALL WRTVCD(VEC1,LU4,1,LBLK) 
            END IF
          END IF
*
* Vw!N-1,M,LW+1> on LU4
*
          IF(IPLUS.EQ.1) THEN
* : Place !0(N-1,M,LW+1)> in VEC1
            IF(NTEST.GE.100) WRITE(6,*) ' IPLUS active '
            CALL GET_NMLW(VEC1,N-1,M,LW+1,MAXV,MAXW,LUN,LURF,1,
     &                    LU6,1,LBLK,LUOUTEFF)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Correction vector !0(N-1,M,LW+1)> read '
            END IF
            IF(NTEST.GE.100) THEN
              CALL WRTVCD(VEC1,LU6,1,LBLK)    
            END IF
            I12 = 1
            CALL BVEC(V,IVSM,LU6,LU7,VEC1,VEC2)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' V times Correction vector'
              CALL WRTVCD(VEC2,LU7,1,LBLK)
            END IF
*. - <0!Vw!0> !N-1,M,LW+1>
C           IF(V00.NE.0.0D0 .AND. N+M .NE. 1) THEN
            IF(V00.NE.0.0D0 ) THEN
              ONE = 1.0D0
              FACTOR = -V00
              CALL REWINO(LU7)
              CALL REWINO(LU6)
              CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU7,LU6,LU4,0,LBLK)
            ELSE
              CALL REWINO(LU7)
              CALL COPVCDP(LU7,LU4,VEC1,0,LBLK)
            END IF
          END IF
*
*  W!N,M-1,LW> on LU4
*
          IF(M.GT.0.AND.ABS(LW).LE.N) THEN
            INOT = 1
          ELSE
            INOT = 0
          END IF
          IF(INOT.EQ.1) THEN
            IF(NTEST.GE.100) WRITE(6,*) ' INOT active '
            ICSM = ISMVW(MOD(N,2)+1,MOD(M-1,2)+1)
*. !N,M-1,LW> on LU6
            CALL GET_NMLW(VEC1,N,M-1,LW,MAXV,MAXW,LUN,LURF,1,
     &                    LU6,1,LBLK,LUOUTEFF)
            IF(NTEST.GE.100) THEN
             WRITE(6,*) ' !n,m-1> read in'
            END IF
            IF(NTEST.GE.100) THEN
              CALL WRTVCD(VEC1,LU6,1,LBLK)    
            END IF
C?          RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?          write(6,*) ' norm of input vector to BVEC ',RHSNORM
            I12 = 1    
            CALL REWINO(LU7)
            CALL BVEC(W,IWSM,LU6,LU7,VEC1,VEC2)
C?          RHSNORM = INPRDD(VEC1,VEC2,LU7,LU7,1,LBLK)
C?          WRITE(6,*) ' Norm of RHS after BVEC ', RHSNORM
           
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' W times Correction vector on LU7 '
              CALL WRTVCD(VEC2,LU7,1,LBLK)  
            END IF
*. - <0!W!0> !N,M-1,LW>
C           IF(W00.NE.0.0D0 .AND. N+M .NE. 1) THEN
            IF(W00.NE.0.0D0 ) THEN
              ONE = 1.0D0
              FACTOR = -W00
              CALL REWINO(LU7)
              CALL REWINO(LU6)
              CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU7,LU6,LU4,0,LBLK)
            ELSE
C?            write(6,*) ' route 1 '
              CALL REWINO(LU7)
              CALL COPVCDP(LU7,LU4,VEC1,0,LBLK)
            END IF
            RHSNORM = INPRDD(VEC1,VEC2,LU7,LU7,1,LBLK)
C?          WRITE(6,*) ' Norm of vector on LU7  ', RHSNORM
          END IF
          IF(NTEST.GE.100) THEN
            write(6,*) ' v,w times vectors calculated '
          END IF
*. Add the above vectors and save on LU6
          NVEC = IMIN + IPLUS + INOT
C?        write(6,*) ' imin iplus, inot', imin,iplus,inot
C?          RHSNORM = INPRDD(VEC1,VEC2,LU4,LU4,1,LBLK)
C?          WRITE(6,*) ' Norm of first vector on lu4  ', RHSNORM
          DO IVEC = 1, NVEC
            SCR(IVEC) = 1.0D0
          END DO
          ZERO = 0.0D0
          CALL MVCSMD2(LU4,SCR,ZERO,LU6,LU7,VEC1,VEC2,
     &                 NVEC,1,LBLK)
          IF(NTEST.GE.100) THEN
           WRITE(6,*) ' rhs before sum V01D* D '
           CALL WRTVCD(VEC1,LU6,1,LBLK)       
          END IF
C?        RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?        WRITE(6,*) ' Norm of RHS before sum V01D*D ', RHSNORM
* 
*       -sum(ld=1,n;ls=1,m-1;llw)Vo[1]t d(ld,ls,llw)d(n-ld,m-ls-1,lw-llw)
*
*. d(ld,ls,lw) must be of the same symmety as Vo=W
*. d(n-ld,m-ls-1,(k-l)w) will then have the same symmetry as d(n,m,*)
*
* Add up on LU6, using LU7 as scratch
          ID1SM  = IWRFSM
          ID2SM = MNSM
*
C?        write(6,*) ' TERM 1'
*. Avoid occurance of zero order vector
          DO LD = 0, N
            IF(LD.EQ.0) THEN
              LSMIN = 1
            ELSE
              LSMIN = 0
            END IF
            IF(LD.EQ.N) THEN
              LSMAX = M-2
            ELSE
              LSMAX = M-1
            END IF
            DO LS = LSMIN, LSMAX
              DO LLW = -LD,LD, 2
*. within frequency bounds ?
              IF(ABS(LW-LLW).LE.N-LD) THEN
                IF(NTEST.GE.5) THEN
                  WRITE(6,*) ' Within frequency bounds for LD,LS,LLW=',
     &            LD,LS,LLW
                END IF
                ID1PSM  = ISMVW(MOD(LD,2)+1,MOD(LS,2)+1)
                ID2PSM  = ISMVW(MOD(N-LD,2)+1,MOD(M-LS-1,2)+1)
                IF(ID1PSM.EQ.ID1SM.AND.ID2PSM.EQ.ID2SM) THEN
                  IF(NTEST.GE.5) THEN
                    WRITE(6,*) 
     &              ' term (Vo[1]t d(ld,ls,llw))d(n-ld,m-ls-1,lw-llw)',
     &              ' with right symmetry'
                    WRITE(6,*)
     &              'LD LS',LD,LS
                  END IF
*. Position LUN for reading D(n-ld,m-ls-1,(k-l)w)
                  CALL GET_NMLW(VEC2,N-LD,M-LS-1,LW-LLW,
     &                 MAXV,MAXW,LUN,LURF,0,LU7,1,LBLK,LUOUTEFF)
CW                IF(LD.EQ.N-LD.AND.LS.EQ.M-LS-1.AND.LLW.EQ.LW-LLW)
CW   &            THEN
CW                  COEF = -2.0D0*V01D(LD,LS,LLW)
CW                ELSE
                    COEF = -V01D(LD,LS,LLW)
CW                END IF
C?                write(6,*) ' COEF = ',coef
                  ONE = 1.0D0
                  CALL REWINO(LU6)
                  CALL REWINO(LU7)
                  CALL VECSMDP(VEC1,VEC2,ONE,COEF,LU6,LUOUTEFF,
     &            LU7,0,LBLK)
                  CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
                END IF
              END IF
              END DO
            END DO
          END DO
*
*       -sum(ld=1,n-1;ls=1,m;llw)Vw[1]t d(ld,ls,lw)d(n-ld-1,m-ls,lw-llw-1)
*
C?        write(6,*) ' TERM 2 '
          ID1SM  = IVRFSM
          DO LD = 0, N-1
            IF(LD.EQ.0) THEN
              LSMIN = 1
            ELSE
              LSMIN = 0
            END IF
            IF(LD.EQ.N-1) THEN
              LSMAX = M-1
            ELSE
              LSMAX = M
            END IF
            DO LS = LSMIN, LSMAX
              DO LLW = -LD,LD,2
*. within frequency bounds ?
              IF(ABS(LW-LLW-1).LE.N-LD-1) THEN
                IF(NTEST.GE.5) THEN
                  WRITE(6,*) ' Within frequency bounds for LD,LS,LLW=',
     &            LD,LS,LLW
                END IF
                ID1PSM  = ISMVW(MOD(LD,2)+1,MOD(LS,2)+1)
                ID2PSM  = ISMVW(MOD(N-LD-1,2)+1,MOD(M-LS,2)+1)
                IF(ID1PSM.EQ.ID1SM.AND.ID2PSM.EQ.ID2SM) THEN
                  IF(NTEST.GE.5) THEN
                    WRITE(6,*) 
     &              ' term (Vw[1]t d(ld,ls,lw))d(n-ld-1,m-ls,lw-llw-1)',
     &              ' with right symmetry'
                    WRITE(6,*)
     &              'LD LS',LD,LS
                  END IF
*. Position LUN at start of D(n-ld-1,m-ls,(k-l)w)
                  CALL GET_NMLW(VEC2,N-LD-1,M-LS,LW-LLW-1,
     &                 MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
CW                IF(LD.EQ.N-LD-1.AND.LS.EQ.M-LS.AND.LLW.EQ.LW-LLW-1)
CW   &            THEN
CW                  COEF = -2.0D0*VW1D(LD,LS,LLW)
CW                ELSE
                    COEF = -VW1D(LD,LS,LLW)
CW                END IF
C?                write(6,*) ' COEF = ',coef
                  ONE = 1.0D0
                  CALL REWINO(LU6)
                  CALL REWINO(LU7)
                  CALL VECSMDP(VEC1,VEC2,ONE,COEF,LU6,
     &                 LUOUTEFF,LU7,0,LBLK)
                  CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
                END IF
              END IF
              END DO
            END DO
          END DO
*
*       -sum(ld=1,n-1;ls=1,m;ll)(V-w[1]t d(ld,ls,lw))d(n-ld-1,m-ls,(l-ll+1)w)
*
C?        write(6,*) ' TERM 3 '
          DO LD = 0, N-1
            IF(LD.EQ.0) THEN
              LSMIN = 1
            ELSE
              LSMIN = 0
            END IF
            IF(LD.EQ.N-1) THEN
              LSMAX = M-1
            ELSE
              LSMAX = M
            END IF
            DO LS = LSMIN, LSMAX
              DO LLW = -LD,LD,2 
*. within frequency bounds ?
              IF(ABS(LW-LLW+1).LE.N-LD-1) THEN
                IF(NTEST.GE.5) THEN
                  WRITE(6,*) ' Within frequency bounds for LD,LS,LLW=',
     &            LD,LS,LLW
                END IF
                ID1PSM  = ISMVW(MOD(LD,2)+1,MOD(LS,2)+1)
                ID2PSM  = ISMVW(MOD(N-LD-1,2)+1,MOD(M-LS,2)+1)
                IF(ID1PSM.EQ.ID1SM.AND.ID2PSM.EQ.ID2SM) THEN
                  IF(NTEST.GE.5) THEN
                    WRITE(6,*) 
     &              'term (V-w[1]t d(ld,ls,lw))d(n-ld-1,m-ls,lw-llw+1)'
     &              ,' with right symmetry'
                    WRITE(6,*)
     &              'LD LS',LD,LS
                  END IF
*. Position LUN for reading D(n-ld-1,m-ls,(k-l)w)
                  CALL GET_NMLW(VEC2,N-LD-1,M-LS,LW-LLW+1,
     &                 MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
*. We assume that the perturbation is real, 
                  IF(IRC.EQ.1) THEN
                    COEF = -VW1D(LD,LS,LLW)
                  ELSE
                    COEF =  VW1D(LD,LS,LLW)
                  END IF
CW                IF(LD.EQ.N-LD-1.AND.LS.EQ.M-LS.AND.LLW.EQ.LW-LLW+1)
CW   &            THEN
CW                  COEF = 2.0D0*COEF             
CW                END IF
C?                write(6,*) ' COEF = ',coef
                  ONE = 1.0D0
                  CALL REWINO(LU6)
                  CALL REWINO(LU7)
                  CALL VECSMDP(VEC1,VEC2,ONE,COEF,LU6,
     &                 LUOUTEFF,LU7,0,LBLK)
                  CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
                END IF
              END IF
              END DO
            END DO
          END DO
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' After the 3 terms '
            CALL WRTVCD(VEC1,LU6,1,LBLK)
          END IF
*. We have now assambled the right hand side vector
*. Well there was a minus in front of (e[2]-wI)**-1 
          ONEM = -1.0D0
          CALL SCLVCD(LU6,LU7,ONEM,VEC1,1,LBLK)
CSEPT29   CALL COPVCDP(LU7,LU6,VEC ,1,LBLK)
          CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
          IF(NTEST.GE.5) THEN
             WRITE(6,*) ' Assembling right hand site vector finished'
          END IF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Right hand side vector before projection '
            CALL WRTVCD(VEC1,LU6,1,LBLK)
          END IF
*. Symmetry of right hand side is also ISSM
          ICSM = ISSM
C?        RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?        WRITE(6,*) ' Norm of RHS 1 ', RHSNORM
*. project !0> component out, save result on LU6
          IF(ICSM.EQ.IRFSM) THEN
            OVLAP = INPRDD(VEC1,VEC2,LURF,LU6,1,LBLK)
C?          write(6,*) ' ovlap ', OVLAP
            FACTOR = -OVLAP
            ONE = 1.0D0
            CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU6,LURF,LU7,1,LBLK)
            CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
          END IF
C?        RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?        WRITE(6,*) ' Norm of RHS ', RHSNORM
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Right hand side vector after projection '
            CALL WRTVCD(VEC1,LU6,1,LBLK)
          END IF
*.save in LU1 and solve linear eq .
          CALL REWINO(LU1)
          CALL COPVCDP(LU6,LU1,VEC1,1,LBLK)
*. Should !0> be printed out in each micro
          IF(ICSM.EQ.IRFSM) THEN
            LUPROJ = LURF
            IPROJ = 1
          ELSE
            LUPROJ = 0
            IPROJ = 0
          END IF
*. Change core energy to include -lw *FREQ
          SHIFT = -ERF + LW*FREQ
C?        WRITE(6,*) ' Shift in GNDB.. ', SHIFT
          I12 = 2
C?        WRITE(6,*) ' LUPROJ and SHIFT', LUPROJ,SHIFT
          CALL HINTV(LU1,LU2,SHIFT,SHIFT,VEC1,VEC2,LBLK,LUPROJ,LUPROJ)
*. Solution to equations are hiding on LU2, transfer to LUN
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' NEW Correction vector'
            CALL WRTVCD(VEC1,LU2,1,LBLK)
          END IF
*. Add to our collection of correction vectors
          CALL GET_NMLW(VEC1,N,M,LW,
     &         MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
          CALL REWINO(LU2)
          CALL COPVCDP(LU2,LUN,VEC1,0,LBLK)
*. Test overlap with reference
          IF(ICSM.EQ.IRFSM) THEN
            OVLAP = INPRDD(VEC1,VEC2,LURF,LU2,1,LBLK)
C?          write(6,*) ' overlap between reference and correction ',
C?   &      IORD,OVLAP
          END IF
*. Obtain the V[1]t D(n,m,lw) arrays
*
          CALL REWINO(LU5)
*
C?        write(6,*) ' MNSM IVRFSM IWRFSM ', MNSM,IVRFSM,IWRFSM
          IF(MNSM.EQ.IVRFSM) THEN
            CALL REWINO(LU2)
            VW1D(N,M,LW) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
          ELSE
            CALL SKPVCD(LU5,1,VEC1,0,LBLK)
            VW1D(N,M,LW) = 0.0D0
          END IF
          IF(NTEST.GE.5) 
     &    write(6,*) '  VW1D(N,M,LW) =  ', VW1D(N,M,LW)
*
          IF(MNSM.EQ.IWRFSM) THEN
            CALL REWINO(LU2)
            V01D(N,M,LW) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
          ELSE
            CALL SKPVCD(LU5,1,VEC1,0,LBLK)
            V01D(N,M,LW) = 0.0D0
          END IF
          IF(NTEST.GE.5) 
     &    write(6,*) ' V01D(N,M,LW) =  ', V01D(N,M,LW)
*
         END DO
*        ^ End of loop over frequency components
        END DO
*       ^ End of loop over W
       END DO
*       ^ End of loop over V
*
       IF(NTEST.GE.5) THEN
*
       WRITE(6,*)
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' V0[1]t D(ld,ls,lw) : '
       WRITE(6,*) ' ==================== '
       WRITE(6,*)
       DO LD = 0, MAXV
         IF(LD.EQ.0) THEN
           LSMIN = 1
         ELSE
           LSMIN = 0
         END IF
         DO LS = LSMIN, MAXW
           DO LW = -LD,LD,2
             WRITE(6,'(A,3I3,A,E18.10)') 
     &       ' ld, ls, lw = ',LD,LS,LW,'   ,   ',V01D(LD,LS,LW)
           END DO
         END DO
       END DO
*
       WRITE(6,*)
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' Vw[1]t D(ld,ls,lw) : '
       WRITE(6,*) ' ==================== '
       WRITE(6,*)
       DO LD = 0, MAXV
         IF(LD.EQ.0) THEN
           LSMIN = 1
         ELSE
           LSMIN = 0
         END IF
         DO LS = LSMIN, MAXW
           DO LW = -LD,LD,2
             WRITE(6,'(A,3I3,A,E18.10)') 
     &       ' ld, ls, lw = ',LD,LS,LW,'   ,   ',Vw1D(LD,LS,LW)
           END DO
         END DO
       END DO
*
       END IF
*      ^ End of print section
*
*. And then the susceptibilities
*
* The n'th ,m'th order susceptibility with w the operator to be studied
*  is <<w;v,v...,w,w..>>f1,f2,f3,...fn,0,0,0, 
* where f1,f2,f3 ... fn equals +/- freq. It depends only upon the total 
* frequency  ftot= sum(i=1,n)fi and can be written as
*
* Khi(ftot,f,f,f,f,f...,0,0,0,...)w,v,v,v,....,w,w,...
*
*. The susceptibility is the term on the Taylor expansion of 
* <0!w!0> with the right frequencies and perturbations
* and by insisting on symmetry on all indeces after the 
*
*. Inner products <d!d>
*
*. Currently we only asemble up to MAXV,MAXW
      DO LD = 0,MAXV
        DO LS = 0,MAXW
          DO LW = -LD,LD,2     
            IF(NTEST.GE.3) THEN
              WRITE(6,*)
              write(6,*) ' Terms to DD for LD LS LW =',LD,LS,LW
              WRITE(6,*)
            END IF
            DD(LD,LS,LW) = 0.0D0
*. Obtain <d|d>(ld,ls,lw) = 
*  sum(ldl,lsl,lwl)dt(ldl,lsl,-lwl)d(ld-ldl,ls-lsl,lw-lwl)
            DO LDL = 0, LD
              DO LSL = 0, LS
                LDR = LD - LDL
                LSR = LS - LSL
*. Identical symmetries
                ILSM = ISMVW(MOD(LDL,2)+1,MOD(LSL,2)+1)
                IRSM = ISMVW(MOD(LDR,2)+1,MOD(LSR,2)+1)
                IF(LDR.LE.MAXV.AND.LSR.LE.MAXW.AND.ILSM.EQ.IRSM) THEN
                 DO LWL = -LDL,LDL,2
                  LWR  = LW -LWL
*. Other frequency component allowed ?
                  IF(ABS(LWR).LE.LDR) THEN
*. Fetch D(LDL,LSL,-LWL) and save on LU7
                    CALL GET_NMLW(VEC2,LDL,LSL,-LWL,         
     &                   MAXV,MAXW,LUN,LURF,1,LU7,1,LBLK,LUOUTEFF)
*. Position at start of D(LDR,LSR,LWR)
                    CALL GET_NMLW(VEC1,LDR,LSR, LWR,         
     &                   MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
                    CALL REWINO(LU7)
                    CONT = INPRDD(VEC1,VEC2,LUOUTEFF,LU7,0,LBLK)
*
                    IF(NTEST.GE.3) THEN
                      WRITE(6,*) ' D(ldl,lsl,-lwl)D(ldr,lsr,lwr) for '
                      WRITE(6,*) '  ldl lsl lwl ldr lsr lwr :',
     &                LDL,LSL,LWL,LDR,LSR,LWR
                      WRITE(6,*) '    is ',CONT
                    END IF
                    DD(LD,LS,LW) = DD(LD,LS,LW)+CONT
                  END IF
                 END DO
                END IF
*               ^ End of check of correct symmetries
              END DO
*             ^ End of loop over LSL
            END DO
*           ^ End of loop over LSR
          IF(NTEST.GE.3) WRITE(6,*) ' DD(LD,LS,LW)= ',DD(LS,LS,LW)
          END DO
*         ^ End of loop over LW
        END DO
*       ^ End of loop over LS
      END DO
*     ^ End of loop over LD
      WRITE(6,*)
      WRITE(6,*) ' ==================== '
      WRITE(6,*) '    DD(ld,ls,lw) : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      DO LD = 0, MAXV
        IF(LD.EQ.0) THEN
          LSMIN = 1
        ELSE
          LSMIN = 0
        END IF
        DO LS = LSMIN, MAXW
          DO LW = -LD,LD,2
            WRITE(6,'(A,3I3,A,E18.10)') 
     &      ' ld, ls, lw = ',LD,LS,LW,'   ,   ',DD(LD,LS,LW)
          END DO
        END DO
      END DO
*. We will now calculate elements <d!A!d> and the 
*  corresponding response functions.
* This involves a loop over A- operators so :
      DO I_AVE_OP = 1, N_AVE_OP
        WRITE(6,*)
        WRITE(6,'(A)') 
     &  ' *****************************************'
        WRITE(6,'(A,A)') 
     &  ' Info for A - operator ', AVE_OP(I_AVE_OP)
        WRITE(6,'(A)') 
     &  ' *****************************************'
        WRITE(6,*)
*
        IASM = IAVE_SYM(I_AVE_OP)
C?      WRITE(6,*) ' IASM = ', IASM
*
*. Inner products <d!A!d>
*
* Notation A => W
      DO LD = 0,MAXV
        DO LS = 0,MAXW
             
          DO LW = -LD,LD,2     
            DWD(LD,LS,LW) = 0.0D0
*. Obtain <d|w!d>(ld,ls,lw) = 
*  sum(ldl,lsl,lwl)dt(ldl,lsl,-lwl)d(ld-ldl,ls-lsl,lw-lwl)
            IF(NTEST.GE.5) THEN
              WRITE(6,*) 
     &        ' DAD under construction for LD,LS,LW=',LD,LS,LW
            END IF
*
            DO LDL = 0, LD
              DO LSL = 0, LS
                LDR = LD - LDL
                LSR = LS - LSL
*. Identical symmetries
                ILSM = ISMVW(MOD(LDL,2)+1,MOD(LSL,2)+1)
                IRSM = ISMVW(MOD(LDR,2)+1,MOD(LSR,2)+1)
*W*right
C               IWRSM = ISMVW(MOD(LDR,2)+1,MOD(LSR+1,2)+1)
                IWRSM = MULTD2H(IRSM,IASM)
C?              WRITE(6,*) 'IWRSM = ', IWRSM
                IF(LDR.LE.MAXV.AND.LSR.LE.MAXW.AND.ILSM.EQ.IWRSM) THEN
                 DO LWL = -LDL,LDL,2
                  LWR  = LW -LWL
*. Other frequency component allowed ?
                  IF(ABS(LWR).LE.LDR) THEN
*. Fetch D(LDR,LSR,LWR) on LU7
                    CALL GET_NMLW(VEC1,LDR,LSR, LWR,         
     &                   MAXV,MAXW,LUN,LURF,1,LU7,1,LBLK,LUOUTEFF)
*. A time D(LDR,LSR,LWR) on LU6
                    ICSM = ISMVW(MOD(LDR,2)+1,MOD(LSR,2)+1)
C                   ISSM= ISMVW(MOD(LDR,2)+1,MOD(LSR+1,2)+1)
                    ISSM = MULTD2H(ICSM,IASM)
C?                  WRITE(6,*) ' ISSM = ', ISSM
*
                    I12 = 1
                    CALL BVEC(AVINT(1,I_AVE_OP),IASM,LU7,LU6,
     &                        VEC1,VEC2)
C                        BVEC(B,IBSM,LUC,LUB,VEC1,VEC2)
*. Fetch D(LDL,LSL,-LWL) on LU7
                    CALL GET_NMLW(VEC1,LDL,LSL,-LWL,         
     &                   MAXV,MAXW,LUN,LURF,1,LU7,1,LBLK,LUOUTEFF)
                    CONT = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                    IF(NTEST.GE.3) THEN
                     WRITE(6,*) ' Contribution to DAD : '
                     WRITE(6,*) '    LDL LSL LWL',LDL,LSL,LWL
                     WRITE(6,*) '    LDR LSR LWR',LDR,LSR,LWR
                     WRITE(6,*) '    CONT = ',CONT
                    END IF
                    DWD(LD,LS,LW) = DWD(LD,LS,LW)+CONT
                  END IF
                 END DO
                END IF
*               ^ End of check of correct symmetries
              END DO
*             ^ End of loop over LSL
            END DO
*           ^ End of loop over LSR
          END DO
*         ^ Enf of loop over LW
        END DO
*       ^ End of loop over LS
      END DO
*     ^ End of loop over LD
      WRITE(6,*)
      WRITE(6,*) ' ==================== '
      WRITE(6,*) '    DAD(ld,ls,lw) : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      DO LD = 0, MAXV
        IF(LD.EQ.0) THEN
          LSMIN = 1
        ELSE
          LSMIN = 0
        END IF
        DO LS = LSMIN, MAXW
          DO LW = -LD,LD,2
            WRITE(6,'(A,3I3,A,E18.10)') 
     &      ' ld, ls, lw = ',LD,LS,LW,'   ,   ',DWD(LD,LS,LW)
          END DO
        END DO
      END DO
*
* **********************************
* . And then : The susceptibilities
* **********************************
*
* The khi's are to obtained as expansion of (<0!A!0>/<0!0>) 
*
* It is slightly inconvenient to expand 1/(<0!0>) to arbitary order.
*. Instead we multiply both sides by <0!0>
*  Khi*<0!0> = <0!A!0> and expand this
*
      DO LD = 0, MAXV
       DO LS = 0, MAXW
         DO LW = -LD,LD,2
         IF(NTEST.GE.3) THEN
           WRITE(6,*) ' Susceptibility for LD,LS,LW =',LD,LS,LW
           WRITE(6,*) ' ========================================='
           WRITE(6,*)
         END IF
*. Obtain chi(ld,ls,-lw) = 
*<0!A!0>(ld,ls,lw) -
*sum(ldp,lsp,lwp,ldp+lsp.ge.0)chi(ldp,lsp,-lwp)<0!0>(ld-ldp,ls-lsp,lw-lwp)
         CHI(LD,LS,-LW) = DWD(LD,LS,LW)
         IF(NTEST.GE.3)THEN
            WRITE(6,*) 
     &      ' Initial term <0!A!0>(ld,ls,lw) ',CHI(LD,LS,-LW)
            WRITE(6,'(A)')  
     &      '  Terms -chi(ldp,lsp,-lwp) <0!0>(ld-ldp,ls-lsp,lw-lwp) : '
         END IF
         DO LDP = 0, LD
*
          IF(LDP.EQ.LD) THEN
           LSPMAX = LS-1
          ELSE
           LSPMAX = LS
          END IF
*
          DO LSP = 0,LSPMAX
           DO LWP = -LDP,LDP
*. Allowed frequency splitting
            IF(ABS(LW-LWP).LE.(LD-LDP)) THEN
              TERM =- CHI(LDP,LSP,-LWP)*DD(LD-LDP,LS-LSP,LW-LWP)
              CHI(LD,LS,-LW) = CHI(LD,LS,-LW) + TERM
              IF(NTEST.GE.3) THEN
                WRITE(6,*)
     &          '   ldp lsp -lwp = ', LDP,LSP,-LWP ,' is ', term
              END IF
            END IF
           END DO
          END DO
*         ^ End of LSP,LSW loops
         END DO
*.       ^ End of LDP loop, end thereby end of terms for given susceptibility
         IF(NTEST.GE.3) THEN
         WRITE(6,*) ' Unscaled Susceptibility for LD,LS,-LW = ',
     &              LD,LS,-LW,  ' is ',   CHI(LD,LS,-LW)
         END IF
*
         END DO
*        ^ End of loop over LW
       END DO
*      ^ End of loop over LS
      END DO
*     ^ End of loop over LD
*
*
*. What we have calculated now is the term in a powerseries expansion
*  of given orders and frequency. Three scalings must be performed
* 1 : susceptibilities are taylor coefficients, multiply with ( ld+ls)!
* 2 : Susceptibilities like <<A;Vd,Vs>> and <<A;Vs,Vd>> are identical
*     and only their sum was obtained. Divide by number of different
*     components to obtain individual susceptibilities
*     but separate.
* 3 : There can be several frequence components underlying a givrn total
*     frequence, f.ex, second order zero freq corresponds to +w,-w
*     and -w +w. These two should be considered separate susceptibilities.
       WRITE(6,*)
       WRITE(6,*) ' ================================= '
       WRITE(6,*) ' Response functions Chi(ld,ls,lw) : '
       WRITE(6,*) ' ================================= '
       WRITE(6,*)
       WRITE(6,*)
       DO LD = 0, MAXV
         IF(LD.EQ.0) THEN
           LSMIN = 1
         ELSE
           LSMIN = 0
         END IF
         DO LS = LSMIN, MAXW
           DO LW = -LD,LD,2
*
             FAC1 = IFAC(LD+LS)
*. Number of ld, ls combinations
             FAC2 = IBION(LD+LS,LD)
*. Number of components with frequency up
             NUP = (LD+LW)/2
             NCOMB = IBION(LD,NUP)
             FAC = FAC1/FLOAT(NCOMB)/FAC2
C?           WRITE(6,*) ' FAC1, FAC2, NCOMB ',FAC1,FAC2, NCOMB
             CHI(LD,LS,LW) = FAC * CHI(LD,LS,LW)
*
             WRITE(6,'(A,3I3,A,E18.10)') 
     &       ' ld, ls, lw = ',LD,LS,LW,'   ,   ',CHI(LD,LS,LW)
           END DO
         END DO
       END DO
*
       END DO
*      ^ End of loop over A-operators 
*
      WRITE(6,*) ' Returning from GNDBPTFREQ '
*
      RETURN
      END
      SUBROUTINE BVEC(B,IBSM,LUC,LUB,VEC1,VEC2)
*
* Construct B defined as
*
*       B = B !0> (no projection) 
*
* Where B is a one-electron operator with symmetry IBSM
* and !0> is assumed stored on LUC
*
* Jeppe Olsen , Jan. 1  1990
*               Sept97 : Modified to LUCIA
*
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      DIMENSION B(*)
      integer kvec1,kvec2
      INCLUDE 'glbbas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cands.inc'
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from BVEC '
        WRITE(6,*) ' ==============='
        WRITE(6,*)
        WRITE(6,*) ' ICSM, IBSM = ', ICSM, IBSM
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input vector in BVEC '
        LBLK = -1
        CALL WRTVCD(VEC1,LUC,1,LBLK)
      END IF
*. Ensure that symmetries are consistent
      ISSMP = MULTD2H(IBSM,ICSM)
      IF(ISSMP .NE. ISSM ) THEN
        WRITE(6,*) ' Inconsistent symmetries in BVEC'
        WRITE(6,*) ' ICSM, ISSM, IBSM ', ICSM,ISSM,IBSM
        STOP       ' Inconsistent symmetries in BVEC'
      END IF
*. Get integrals in place
      CALL SWAPVE(B,KINT1,NTOOB**2)
*. Just one-electron operator
      I12 = 1
*. No additional approximations
      IAPR = 0
      IPERTOP = 0
*. And : DO IT
      CALL MV7(VEC1,VEC2,LUC,LUB,0,0)
*. Restore order
      CALL SWAPVE(B,KINT1,NTOOB**2)
*
      RETURN
      END
      SUBROUTINE CI_RESPONS(LU1,LU2,LU3,LU4,LU5,LU6,LU7,LUC,LUDIA,
     &                      VEC1,VEC2,ENOT,CCALC)
*
* Master routine for LUCIA calculation of properties
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cands.inc'
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*)
*. Local scratch 
       INTEGER IAVE_SYM(20)
*
      NTEST = 000
*
      IDUM = -1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'RESPON')
      IRFSM = ICSM
*
C?    write(6,*) ' ENOT IN CIRESP ', ENOT
*
      WRITE(6,*) ' Number of CI response calculations ', NRESP
*
      DO ICALC= 1, NRESP
       WRITE(6,*)
       WRITE(6,'(A)')
     & '  ********************************************** '
       WRITE(6,'(A,I3)')
     & '  Information about response calculation ..',ICALC
       WRITE(6,'(A)')
     & '  ********************************************** '
       WRITE(6,*)
*
*
**  General procedure for for double perturbation theory
*                 through arbitrary order
*
         
          WRITE(6,*)
          WRITE(6,*) ' ================================================'
          WRITE(6,*) ' General expansion in two external perturbations'
          WRITE(6,*) ' ================================================'
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*) 
     &    ' Operators used for expectation values ( A-operators) '
          WRITE(6,*)
     &    ' ====================================================='
          WRITE(6,*)
          DO IAVE = 1, N_AVE_OP
            WRITE(6,'(1H ,5X,A)') AVE_OP(IAVE)
          END DO
          WRITE(6,*)
          WRITE(6,*) ' Perturbations : '
          WRITE(6,*) ' =============== '
          WRITE(6,*)
          WRITE(6,*) '     Operator 1 : ', RESP_OP(1,ICALC)
          WRITE(6,*) '     Operator 2 : ', RESP_OP(2,ICALC)
          WRITE(6,*) 
     &    '     Max order in operator 1 : ', MAXORD_OP(1,ICALC)
          WRITE(6,*) 
     &    '     Max order in operator 2 : ', MAXORD_OP(2,ICALC)
          IF(RESP_W(ICALC).EQ.0.0D0) THEN
            WRITE(6,*) '     Static perturbations'
          ELSE
            WRITE(6,*) 
     &      '    frequency of second operator  ',RESP_W(ICALC)
          END IF
          WRITE(6,*)
*. Space for one_electron integrals
          CALL MEMMAN(KINTOP1,  NTOOB**2,'ADDL  ',2,'INTOP1') !done
          CALL MEMMAN(KINTOP2,  NTOOB**2,'ADDL  ',2,'INTOP2') !done
          CALL MEMMAN(KLSCR  ,5*NTOOB**2,'ADDL  ',2,'RSPSCR') !done
          CALL MEMMAN(KINTAV,N_AVE_OP*NTOOB **2,'ADDL  ',2,'INTAVE') !d
*. Energy corrections
          NCORR = (MAXORD_OP(1,ICALC)+1)*(MAXORD_OP(2,ICALC)+1)
          CALL MEMMAN(KLENM,NCORR,'ADDL  ',2,'ENM   ') !done
*. Symmetry of perturbations
C     SYM_FOR_OP(OPER,IXYX,ISYM)
*. W : Op1, V : Op 2
          CALL SYM_FOR_OP(RESP_OP(1,ICALC),IXYZSYM,IWSM)
          CALL SYM_FOR_OP(RESP_OP(2,ICALC),IXYZSYM,IVSM)
*. Symmetry of average operators
          DO I_AVE_OP = 1, N_AVE_OP
            CALL SYM_FOR_OP(AVE_OP(I_AVE_OP),IXYZSYM,
     &                      IAVE_SYM(I_AVE_OP))
          END DO
*. Perturbations times reference
* W !0>
          IWRFSM = MULTD2H(IWSM,IRFSM)
* V !0>
          IVRFSM = MULTD2H(IVSM,IRFSM)
* VW !0>
          IWVRFSM = MULTD2H(IWSM,IVRFSM)
*. Obtain property integrals
          CALL GET_PROP_PERMSM(RESP_OP(1,ICALC),IPERMSM)
          CALL GET_PROPINT(dbl_mb(KINTOP1),IWSM,RESP_OP(1,ICALC),
     &                     dbl_mb(KLSCR),NTOOBS,NTOOBS,NSMOB,1,IPERMSM )
          CALL GET_PROP_PERMSM(RESP_OP(2,ICALC),IPERMSM)
          CALL GET_PROPINT(dbl_mb(KINTOP2),IVSM,RESP_OP(2,ICALC),
     &                     dbl_mb(KLSCR),NTOOBS,NTOOBS,NSMOB,1,IPERMSM )
          DO I_AVE_OP = 1, N_AVE_OP
            KINTAV2 = KINTAV + (I_AVE_OP-1)*NTOOB**2
            CALL GET_PROP_PERMSM(AVE_OP(I_AVE_OP),IPERMSM)
            CALL GET_PROPINT(dbl_mb(KINTAV2),IAVE_SYM(I_AVE_OP),
     &           AVE_OP(I_AVE_OP),dbl_mb(KLSCR),NTOOBS,NTOOBS,NSMOB,1,
     &           IPERMSM ) 
          END DO
*
          MAXW = MAXORD_OP(1,ICALC)
          MAXV = MAXORD_OP(2,ICALC)
C?        WRITE(6,*) ' MAXW, MAXV ', MAXW,MAXV
*
          FREQ = RESP_W(ICALC)
          NTOOB2 = NTOOB ** 2
          IF(FREQ.EQ.0.0D0) THEN
*. Static perturbation
            ZERO = 0.0D0
            CALL GNDBPTFREQ(dbl_mb(KINTOP2),dbl_mb(KINTOP1),IRFSM,
     &           ENOT,LUC, ! LUHI0, ????
     &           IVRFSM,MAXV,IWRFSM,MAXW,IWVRFSM,IVSM,IWSM,
     &           LUN,VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           LUDIA,ZERO,dbl_mb(KLENM),
     &           N_AVE_OP,AVE_OP,dbl_mb(KINTAV),NTOOB2,IAVE_SYM,
     &           IPRRSP)
          ELSE
*. Frequency dependent perturbation
            CALL GNDBPTFREQ(dbl_mb(KINTOP2),dbl_mb(KINTOP1),IRFSM,
     &           ENOT,LUC, ! LUHI0, ????
     &           IVRFSM,MAXV,IWRFSM,MAXW,IWVRFSM,IVSM,IWSM,
     &           LUN,VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           LUDIA,FREQ,dbl_mb(KLENM),
     &           N_AVE_OP,AVE_OP,dbl_mb(KINTAV),NTOOB2,IAVE_SYM,
     &           IPRRSP)
        END IF
*       ^ End of switch between static/dynamic perturbation
      END DO
*.    ^ End of loop over calculations
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'RESPON')
      RETURN
      END
c $Id$
