      SUBROUTINE vib_HMASS(HESST,ATMASS,first_pass)
* $Id: vib_hmass.F,v 1.3 1997-04-11 09:26:40 d3e129 Exp $
C
C  This routine mass weights and scales the Hessian matrix
C       The scaling is done to avoid numerical problems in the
C       diagonalization routine
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
      INTEGER NAT, NAT3, NHESS, NHESST
      COMMON /CVIB_HESS/ NAT,NAT3,NHESS,NHESST   ! HESSIAN INFORMATION
c
      double precision HESST(NHESST) ! lower triangular Hessian
      double precision ATMASS(NAT) ! mass of the atoms
      logical first_pass
c
      double precision fact, scale
      integer ii, jj, jjend, iatii, iatjj, idum
C
C      set up function for locating i,j elements packed canonically as ij
C
      integer i, j, isym2, iatom
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
      if (first_pass) then
        WRITE(6,9000)
        CALL vib_PRNTPD(HESST,NHESST,NAT3,6)
      endif

 9000 FORMAT(/,10X,'-----------------------------------------------',
     &                                                       '-----',
     &  /,10X,'MASS-WEIGHTED NUCLEAR HESSIAN (Hartree/Bohr/Bohr/Kg)',
     &  /,10X,'----------------------------------------------------',
     &       //)
      RETURN
      END
