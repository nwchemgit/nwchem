      SUBROUTINE vib_HMASS(HESST,ATMASS)
* $Id: vib_hmass.F,v 1.1 1995-07-29 23:19:32 d3e129 Exp $
C
C  This routine mass weights and scales the Hessian matrix
C       The scaling is done to avoid numerical problems in the
C       diagonalization routine
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HESST(NHESST) ! lower triangular Hessian
      DIMENSION ATMASS(NATOM) ! mass of the atoms
      COMMON /HESS/ NATOM,NAT3,NHESS,NHESST   ! Hessian information
C
C      set up function for locating i,j elements packed canonically as ij
C
      ISYM2(I,J)=MAX(I,J)*((MAX(I,J))-1)/2 + MIN(I,J)
      IATOM(I)  = (I+2)/3   ! function for coordinate I is on atom IATOM
C
      DO 00100 II = 1,NAT3  ! loop over coordinates
         JJEND = II
         IATII = IATOM(II)  ! coordinate II is for atom IATII
         DO 00100 JJ = 1,JJEND   ! loop over coordinates
            IDUM = ISYM2(II,JJ)  ! get canonical index
            IATJJ = IATOM(JJ)   ! coordinate JJ is for atom IATJJ
            FACT = SQRT(ATMASS(IATII))*SQRT(ATMASS(IATJJ)) ! mass weight
            HESST(IDUM) = HESST(IDUM)/FACT  ! weight Hessian
00100 CONTINUE
      SCALE = 1.0D03                      ! Hessian scaling factor
*dscal
      call dscal(nhesst,scale,hesst,1) ! Scale Hessian for diagonaization
      WRITE(6,9000)
      CALL vib_PRNTPD(HESST,NHESST,NAT3,6)

 9000 FORMAT(/,10X,'-----------------------------------------------',
     &                                                       '-----',
     &  /,10X,'MASS-WEIGHTED NUCLEAR HESSIAN (Hartree/Bohr/Bohr/Kg)',
     &  /,10X,'----------------------------------------------------',
     &       //)
      RETURN
      END
