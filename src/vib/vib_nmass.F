      SUBROUTINE vib_NMASS(VECTOR,ATMASS)
* $Id: vib_nmass.F,v 1.2 1995-10-21 09:16:13 d3e129 Exp $
C
C  This routine unmass weights the normal mode vectors.
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
      INTEGER NATOM, NAT3, NHESS, NHESST
      COMMON /cvib_HESS/ NATOM,NAT3,NHESS,NHESST   ! Hessian information
      
      double precision VECTOR(NAT3,NAT3) ! matrix of normal mode eigenvectors
      double precision ATMASS(NATOM) ! mass of the atoms
c
      DOUBLE PRECISION FACT
      INTEGER II, JJ, IATII
C
C      set up function for locating i,j elements packed canonically as ij
C
      integer iatom, i
      IATOM(I)  = (I+2)/3   ! function for coordinate I is on atom IATOM
C
      DO 00100 II = 1,NAT3  ! loop over coordinates
         IATII = IATOM(II)  ! coordinate II is for atom IATII
         FACT = SQRT(ATMASS(IATII))   !mass weight
         DO 00100 JJ = 1,NAT3   ! loop over coordinates
            VECTOR(II,JJ) = VECTOR(II,JJ)/FACT  ! "unweight normal mode"
00100 CONTINUE
      RETURN
      END
