*
* Fusk routines for new spin densities by MCLR
*
* A set of non-redundant 'singlet single excitations with spinflip'
*  are defined as
*
*  S_ij = aT(1,-1)_ij S+ + bT(1,0) S0 + c T(1,1)_ij S-
*
*  S_ij = a+ia ajb S- + (a+ia aja - a+ib ajb ) Sz + a+ib aja S+
*
*. Accordingly we can write S_ij as the singlet operator
*
* S_ij = sum(k) (-1/2 e_ij,kk - e_ik,kj )+ 3/2 E_ij 
*
* The 2-e part of the above operator is not symmetric in particle 1,2 
* so we rewrite
*
* S_ij = sum(k) (-1/4 e_ij,kk -1/4 e_kk,ij - 1/2 e_ik,kj - 1/2 e_kj,ik) 
*      + 3/2 E_ij
*
* 
*
* and the antisymmetric combination is
*
*

*
* The S operators have the adjoint property
*
* (S_ij) \dag = S_ji
*
* 1 : Set 
*
* Jeppe + Anders, Feb. 99
*
      SUBROUTINE SP_MCLR(ISM,ISPC)
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc' 
      INCLUDE 'cands.inc'
*
      COMMON/VEC12/KVEC1,KVEC2
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SP_MCL')
      WRITE(6,*)  
      WRITE(6,*) ' Welcome to MCLR_SP '
      WRITE(6,*) 
      WRITE(6,*) ' Symmetry and space ', ISM,ISPC
*
*.Transfer to CANDS
      ICSM = ISM
      ISSM = ISM
      ICSPC = ISPC
      ISSPC = ISPC
*. Scratch space for CI - behind the curtain 
       CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)
*
*.1 : Number of non-redundant operators of this symmetry and
*.    indeces of non-redundant operators
* ( it may be that we must subtract some, but ... )
*
*. Divide orbital types into inactive, active, secondary
      CALL GASSPC2(I_IAD,IGSOCCX(1,1,1)) 
*. Space for excitation index matrix and Type-type excitation vector
      CALL MEMMAN(KLTTACT,NGAS*NGAS,'ADDL  ',1,'TTACT ')
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
*. Matrix giving nonredundant type-type excitations
      CALL NONRED_TT_EXC(WORK(KLTTACT),IDUM,0)
*. Orbital excitation indeces and number of orbital excitations
C          NONRED_OO_EXC(NOOEXC,IOOEXC,IPHSM,ITTACT)
      CALL NONRED_OO_EXC(NOOEXC,WORK(KLOOEXC),1,WORK(KLTTACT))
*. Space for E1 vector and E2 matrix, Kappa parameters
      WRITE(6,*) ' Number of nonredundant orbital excitations ',
     &             NOOEXC
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1VEC ')
      CALL MEMMAN(KLE2,(NOOEXC+1)*NOOEXC/2,'ADDL  ',2,'E2MAT ')
*. Complete Hessian, two copies
      CALL MEMMAN(KLE2E,NOOEXC*NOOEXC,'ADDL  ',2,'E2MATE')
      CALL MEMMAN(KLE2E2,NOOEXC*NOOEXC,'ADDL  ',2,'E2MAT2')
*. Kappa parameters
      CALL MEMMAN(KLKAPS,NOOEXC,'ADDL  ',2,'KAPS  ')


*
*. 2: Construct : <0![H,T+mu - T mu, T+ vu - T vu ]]0>
*     (through a sequence of SD calculations (sic))
      CALL E2S_FUSK(WORK(KLE2),NOOEXC,WORK(KLOOEXC))
*. 3: Construct : <0![H,T+mu - T mu, ]]0>
      CALL E1S_FUSK(WORK(KLE2),NOOEXC,WORK(KLOOEXC))
*. 4 : Construct metric  <0![T+mu - T mu, T+ vu - T vu ]]0>
*. 5 : Test for singularities and transform to nonsingular basis
*. 6 : Invert E2    
*. Expand to complete form
C       TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(WORK(KLE2E),WORK(KLE2),2,NOOEXC,NOOEXC)
*          INVMAT(A,B,MATDIM,NDIM,ISING)
      CALL INVMAT(WORK(KLE2E),WORK(KLE2E2),NOOEXC,NOOEXC,ISING)
*. 7 : Solve Kappa S = -E[2] ** E[1]
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)|
      CALL MATVCB(WORK(KLE2E),WORK(KLE1),WORK(KLKAPS),NOOEXC,
     &            NOOEXC,0)
      ONEM = -1.0D0
      CALL SCALVE(WORK(KLKAPS),ONEM,NOOEXC)
      WRITE(6,*) ' Kappa S rotation parameters ' 
      CALL WRTMAT(WORK(KLKAPS),1,NOOEXC,1,NOOEXC)

*. 8 : set up modifications to spin density : <0![Kappa S,T]0>
*
      RETURN
      END
      SUBROUTINE E2S_FUSK(E2S,NOOEXC,IOOEXC)
*
* Set up E2S for singlet spinflip operators 
*
* Outrageous sloppy routine !!!!           
*
* The author would like to remain anonymous
*
* Jan/Feb 99
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER IOOEXC(NTOOB,NTOOB)
*. Output : Lower half packed column
      DIMENSION E2S(*)
*
*. Loop over orbital pairs in 1 and 2
*
      DO I = 1, NTOOB
       DO J = 1, NTOOB
        IJNUM = IOOEXC(I,J)
        IF(IJNUM.GT.0) THEN
         DO K = 1, NTOOB
          DO L = 1, NTOOB
           KLNUM = IOOEXC(K,L)
           IF(KLNUM.GT.0.AND.IJNUM.GE.KLNUM) THEN
            WRITE(6,*) ' Hessian element will be calculated '
            WRITE(6,*) ' (I,J), (K,L) ', I,J,K,L
            CALL E2S_ELEMENT(E2SE,I,J,K,L)
            IJKL = IJNUM*(IJNUM-1)/2 + KLNUM
            E2S(IJKL) = E2SE
           END IF
          END DO
         END DO
*        ^ End of loop over KL
        END IF
       END DO
      END DO
*     ^ End of loop over IJ
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' E2fusk matrix '
        CALL PRSYM(E2S,NOOEXC)
      END IF
*
      RETURN
      END
*
      SUBROUTINE E1S_FUSK(E1S,NOOEXC,IOOEXC)
*
* Set up E1S for singlet spinflip operators 
*
* Outrageous sloppy routine !!!!           
*
* The author would like to remain anonymous
*
* Jan/Feb 99
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER IOOEXC(NTOOB,NTOOB)
*. Output  
      DIMENSION E1S(*)
*
*. Loop over non-redundant orbital excitations 
*
      DO I = 1, NTOOB
       DO J = 1, NTOOB
        IJNUM = IOOEXC(I,J)
        IF(IJNUM.GT.0) THEN
          CALL Z_E1S_ELEMENT(E1SE,I,J)
          E1S(IJNUM) = E1SE
        END IF
       END DO
      END DO
*     ^ End of loop over IJ
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' E1S Vector'
        CALL WRTMAT(E1S,1,NOOEXC1,1,NOOEXC)
      END IF
*
      RETURN
      END
      SUBROUTINE  Z_E1S_ELEMENT(E1SE,I,J)
*
* Construct element <0![H,S_ij-S_ji]!0> = 2 <0!H(S_ij-S_ji)!0>
*
* J.O. Feb. 99
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'glbbas.inc'
      COMMON/VEC12/KVEC1,KVEC2
      REAL * 8 INPRDD
*
      LBLK = -1
*
* Fake integral list for excitation S_ij - S_ji
*
      CALL SET_SP_MCLR_INTLIST(I,J)
*. Tell integral program to work with MCLR fusk integrals
      I_SOP = 1 
      ICC_EXC = 0
*. |0> is on LUC, generate (S_ij-S_ji) on LUSC1
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUC,LUSC1,0,0)
*. H( (S_ij-S_ji) !0> on LUSC2
      I_SOP = 0
      ICC_EXC = 0
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC2,0,0)
*. Complete by taking inner product with <0!
      E1SE = 2.0D0*
     &       INPRDD(WORK(KVEC1),WORK(KVEC2),LUC,LUSC2,1,LBLK)
*
      NTEST = 100
      IF( NTEST.GE.100) THEN
        WRITE(6,*) ' S gradient for i,j = ', I,J, E1SE
      END IF
*
      RETURN
      END
*
      SUBROUTINE E2S_ELEMENT(E2SE,I,J,K,L)
*
* Set up Hessian element <0![[H,S_ij-S_ji,S_kl-S_lk]|0>
*                       = 2 <0![[H,S_ij-S_ji] S_kl-S_lk|0>
*                       + 2 <0![[H,S_kl-S_lk] S_ij-S_ji|0>
*
* J.O. Feb. 1, 1999
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'clunit.inc' 
      COMMON/VEC12/KVEC1,KVEC2
      REAL * 8 INPRDD
*
      LBLK = -1
*
* ==================================
* a : <0![[H,S_ij-S_ji] S_kl-S_lk|0>
* ==================================
*
*( S_kl - S_lk ) |0> on file LUSC1
      CALL SET_SP_MCLR_INTLIST(K,L)
      I_SOP = 1 
      ICC_EXC = 0
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUC,LUSC1,0,0)
* (S_ij - S_ji) ( S_kl - S_lk ) |0>  on file LUSC2
      CALL SET_SP_MCLR_INTLIST(I,J)
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC2,0,0)
*H (S_ij-S_ji) ( S_kl - S_lk ) |0> on file LUSC3 
      I_SOP = 0
      ICC_EXC = 0
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC2,LUSC3,0,0)
*. H ( S_kl - S_lk ) |0>  on LUSC2
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUC,LUSC2,0,0)
* (S_ij-S_ji)  H ( S_kl - S_lk ) |0> on LUH
      I_SOP = 1
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC2,LUH,0,0)
*[H, (S_ij-S_ji)] ( S_kl - S_lk ) |0> on LUSC1
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC3,LUSC1,0,0)
      ONE = 1.0D0
      ONEM = -1.0D0
C          VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      CALL VECSMD(WORK(KVEC1),WORK(KVEC2),ONE,ONEM,LUSC3,LUH,LUSC1,
     &             1,LBLK)
*. Inner product with <0!  
      E2A = 2.0D0*INPRDD(WORK(KVEC1),WORK(KVEC2),LUC,LUSC1,1,LBLK)
*
*
* ==================================
* b : <0![[H,S_kl-S_lk] S_ij-S_ji|0>
* ==================================
*
*( S_ij - S_ji ) |0> on file LUSC1
      CALL SET_SP_MCLR_INTLIST(I,J)
      I_SOP = 1 
      ICC_EXC = 0
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUC,LUSC1,0,0)
* (S_kl - S_lk) ( S_ij - S_ji ) |0>  on file LUSC2
      CALL SET_SP_MCLR_INTLIST(K,L)
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC2,0,0)
*H (S_kl-S_lk) ( S_ij - S_ji ) |0> on file LUSC3 
      I_SOP = 0
      ICC_EXC = 0
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC2,LUSC3,0,0)
*. H ( S_ij - S_ji ) |0>  on LUSC2
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUC,LUSC2,0,0)
* (S_kl-S_lk)  H ( S_ij - S_ji ) |0> on LUH
      I_SOP = 1
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC2,LUH,0,0)
*[H, (S_kl-S_lk)] ( S_ij - S_ji ) |0> on LUSC1
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC3,LUSC1,0,0)
      ONE = 1.0D0
      ONEM = -1.0D0
C          VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      CALL VECSMD(WORK(KVEC1),WORK(KVEC2),ONE,ONEM,LUSC3,LUH,LUSC1,
     &             1,LBLK)
*. Inner product with <0!  
      E2B = 2.0D0*INPRDD(WORK(KVEC1),WORK(KVEC2),LUC,LUSC1,1,LBLK)
*
      E2SE = E2A + E2B
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Hessian element for (i,j), (k,l) ', I,J,K,K 
        WRITE(6,*) ' Sum and the two tems ', E2SE,E2A,E2B
      END IF
*
      RETURN 
      END
      SUBROUTINE MCLR_SPINDEN(KAPPA,S)
*
* Obtain first order correction to spin density matrix
*
* Input Kappa, S are first order orbital and state  respons 
*
* Currently only orbital response is calculated
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
*
      RETURN
      END
      SUBROUTINE SET_SP_MCLR_INTLIST(K,L)
      RETURN 
      END

      
      
      
      
   

    
  


      



      
 

c $Id$
