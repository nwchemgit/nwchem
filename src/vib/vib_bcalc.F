       SUBROUTINE vib_BCALC(IZMAT,BMTRX,HINT,ZMAT,SCR,
     &                 HESSIAN,EIGVEC,EIGVAL)
C*MODULE VIB  *DECK BCALC
* $Id$
C
      IMPLICIT none
#include "nwc_const.fh"
C
      INTEGER MXATM
      PARAMETER (MXATM=nw_max_atom)
      double precision zan,c
      COMMON /cvib_INFOA /  ZAN(MXATM),C(3,MXATM)
c
      LOGICAL LINEAR
      integer nzmat, nzvar, nvar
      COMMON /cvib_ZMAT  / NZMAT,NZVAR,NVAR,LINEAR
c
      integer nrow, ndeg
      COMMON /cvib_VBDCMP/ NROW,NDEG
c
      integer nat, nat3, nhess, nhesst
      COMMON /cvib_HESS/   NAT,NAT3,NHESS,NHESST
c
      integer izmat(nzmat)
      double precision bmtrx(nat3, nat3)
      double precision hint(nat3, nat3)
      double precision zmat(8*nat3)
      double precision scr(*)
      double precision hessian(nat3, nat3)
      double precision eigvec(nat3, nat3)
      double precision eigval(nat3)
      WRITE(6,9010)
      WRITE(6,9020)
 9010 FORMAT(///,10X,'------------------------------------------',/
     &           10X,'INTERNAL COORDINATE DESCRIPTION OF HESSIAN',/
     &           10X,'------------------------------------------')
 9020 FORMAT(///,10X,'------------------------------------------------',
     &         /,10X,'NORMAL MODE EIGENVECTORS IN INTERNAL COORDINATES',
     &         /,10X,'------------------------------------------------',
     &         /,10X,'      (Frequencies expressed in cm-1) ')
      END
