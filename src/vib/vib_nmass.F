      SUBROUTINE vib_NMASS(VECTOR,ATMASS)
* $Id$
C
C  This routine unmass weights the normal mode vectors.
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
      INTEGER NATOM, NAT3, NHESS, NHESST
      COMMON /cvib_HESS/ NATOM,NAT3,NHESS,NHESST   ! Hessian information
      
      double precision VECTOR(NAT3,NAT3) ! matrix of normal mode eigenvectors
      double precision ATMASS(NATOM) ! mass of the atoms
c
      double precision mass_ii
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
         mass_ii = ATMASS(IATII)
         if (abs(mass_ii).lt.1.0d-01) mass_ii = 1.0d05
*        FACT = SQRT(ATMASS(IATII))   !mass weight
         FACT = SQRT(mass_ii)   !mass weight
         DO 00100 JJ = 1,NAT3   ! loop over coordinates
            VECTOR(II,JJ) = VECTOR(II,JJ)/FACT  ! "unweight normal mode"
00100 CONTINUE
      RETURN
      END
      subroutine vib_vecnormal(vectors,nvec)
      implicit none
#include "util.fh"
*.. normalize the eigen-vectors 
      integer nvec
      double precision vectors(nvec,nvec)
c
      double precision vnorm
      integer col
      do col = 1,nvec
        vnorm = ddot(nvec,vectors(1,col),1,vectors(1,col),1)
        vnorm = sqrt(1.0d00/vnorm)
        call dscal(nvec,vnorm,vectors(1,col),1)
      enddo
      end
      subroutine vib_vecphase(coord,vectors,nvec)
      implicit none
#include "util.fh"
c
c...  Fix the phase of the normal modes. To do this we compute the
c...  dot-product of every normal mode with the nuclear coordinates
c...  (which provide a fixed reference point) and we insist that this
c...  shall be non-negative.
c
      integer nvec
      double precision vectors(nvec,nvec)
      double precision coord(nvec)
c
      double precision vnorm
      integer col
      do col = 1,nvec
        vnorm = ddot(nvec,coord,1,vectors(1,col),1)
        if (vnorm.lt.0.0d0) then
          vnorm = -1.0d0
          call dscal(nvec,vnorm,vectors(1,col),1)
        endif
      enddo
      end
