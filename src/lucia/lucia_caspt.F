      SUBROUTINE PTFOCK(LU0,LUN,N,ISM,ISPC)
*
* Perturbation expansion of general fock matrix
* through arbitrary order.
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
* Jeppe, June 98
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
*
      NTEST = 500
*
      WRITE(6,*) 
      WRITE(6,*) ' ======================== '
      WRITE(6,*) ' PTFOCK is now in CONTROL '
      WRITE(6,*) ' ======================== '
      WRITE(6,*)

* a bit on files :
* LUSC36 is LUN.   
* Two additional scratch files to be used are  LUSC1 and LUSC2
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PTFOCK')
*
*
*     ========================
* 1 : Local memory allocation 
*     ========================
*
*. Allocate space for two vector chunks
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(IPRCIX.GE.1)
     &WRITE(6,*) ' PTEKT  : LCSBLK = ', LCSBLK
      LSCR1 = MAX(LCSBLK,LSCR1)
      CALL MEMMAN(KLVEC1,LSCR1,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,LSCR1,'ADDL  ',2,'VEC2  ')
* one-body Density matrices through order n
      NMAT = N+1     
      LENGTH = NMAT * NTOOB ** 2
      CALL MEMMAN(KLDEN1,LENGTH,'ADDL  ',2,'DENN1 ')
*. in symmetry packed form
      CALL MEMMAN(KLDEN1P,LENGTH,'ADDL  ',2,'DENN1P')
*. and Two-body density matrices
       LENGTH = NMAT * NTOOB ** 2*(NTOOB**2+1)/2
      CALL MEMMAN(KLDEN2,LENGTH,'ADDL  ',2,'DENN2 ')
*. Fock matrices through order n
      LENGTH = NMAT * 2*NINT1    
      CALL MEMMAN(KLFN,LENGTH,'ADDL  ',2,'F(N)  ')
*. A scratch matrix ( not a nice thing to say about a matrix )
      LENGTH =  2*NINT1    
      CALL MEMMAN(KLFSCR,LENGTH,'ADDL  ',2,'FSCR  ')
*. Space for Fock matrices, density matrix belonging to a 
*  given symmetry
      LENGTH = NMAT * 2*NINT1    
      CALL MEMMAN(KLFNS,LENGTH,'ADDL  ',2,'F(N)S ')
      LENGTH = NMAT * 2*NINT1    
      CALL MEMMAN(KLDNS,LENGTH,'ADDL  ',2,'D(N)S ')
*. 
      CALL MEMMAN(KLSCR4,2*NTOOB**2,'ADDL  ',2,'SCR4  ')
      CALL MEMMAN(KLSCR5,2*NTOOB**2,'ADDL  ',2,'SCR5  ')
      CALL MEMMAN(KLSCR6,2*NTOOB**2,'ADDL  ',2,'SCR6  ')
*
* ===============================================
* 2 : Construct density matrices through order N
* ===============================================
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
*. No print in density matrices
      IPRDEN_SAVE = IPRDEN
      IPRDEN = 0
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLSCR4),ZERO,LRHO1)
      I12_SAVE = I12
      I12 = 2
      DO K = 0, N
        CALL PERTDN(K,LU0,LUN,ISM,ISPC,WORK(KLVEC1),WORK(KLVEC2),
     &       WORK(KLDEN1+(K-0)*LRHO1),
     &       WORK(KLDEN2+(K-0)*LRHO2),LUSC1,LUSC2,0)
      END DO
      IPRDEN = IPRDEN_SAVE
      I12 = I12_SAVE
*
      WRITE(6,*) ' Memtest 1 : '
      CALL MEMCHK
*
* =====================================================
* 3 : Construct for all orders Fock matrix with unpartioned Hamiltonian 
* =====================================================
*
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB **2*(NTOOB**2+1)/2
      LFOCK = 2*NINT1
*
      ONE = 1.0D0
      ONEM = -1.0D0
      DO K = 0, N
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLFN+(K-0)*LFOCK),ZERO,LFOCK)
*. Full Hamiltonian with K order density
        CALL COPVEC(WORK(KLDEN1+(K-0)*LRHO1),WORK(KRHO1),LRHO1)
        CALL COPVEC(WORK(KLDEN2+(K-0)*LRHO2),WORK(KRHO2),LRHO2)
        CALL FOCK_MAT(WORK(KLFSCR),2)
        CALL COPVEC(WORK(KLFSCR),WORK(KLFN+(K-0)*LFOCK),LFOCK)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 'Correction to Fock matrix of order =',K
          CALL APRBLM2(WORK(KLFN+(K-0)*LFOCK),
     &                 NOCOBS,NOCOBS,NSMOB,0)
        END IF
*
      END DO
*
      IF(NTEST.GE.100) THEN
* Accumulate corrections to Fock matrix
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLSCR4),ZERO,LFOCK)
        ONE = 1.0D0
        DO K = 0, N
          CALL VECSUM(WORK(KLSCR4),WORK(KLSCR4),
     &         WORK(KLFN+(K-0)*LFOCK),ONE,ONE,LFOCK)
        END DO
*
        WRITE(6,*) '  sum(k) Fock(k) '
        WRITE(6,*) ' =============== '
        CALL APRBLM2(WORK(KLSCR4),NOCOBS,NOCOBS,NSMOB,0)
      END IF
      
*
*.Finito
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',1,'PTFOCK')
*
      RETURN
      END
c $Id$
