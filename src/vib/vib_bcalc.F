       SUBROUTINE vib_BCALC(IZMAT,BMTRX,HINT,ZMAT,SCR,
     &                 HESSIAN,EIGVEC,EIGVAL)
C*MODULE VIB  *DECK BCALC
* $Id: vib_bcalc.F,v 1.5 1999-07-14 00:29:13 d3e129 Exp $
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
c:old      DIMENSION IZMAT(NZMAT),BMTRX(NAT3,NAT3),HINT(NAT3,NAT3),
c:old     &          ZMAT(8*NAT3),SCR(1)
c:old      DIMENSION HESSIAN(NAT3,NAT3),EIGVEC(NAT3,NAT3),EIGVAL(NAT3)
C
C
*      CALL GETARRAY(10,ZAN,INTOWP(NAT),301,IDUM)
*      CALL GETARRAY(10,C,INTOWP(3*NAT),342,IDUM)
C
C     ----- CALCULATE THE B MATRIX -----
C
*      CALL vib_BMAT(BMTRX,NZVAR,NVAR,NAT,NAT3,NHESS,IZMAT,ZMAT,LINEAR)
C
*      CALL vib_TRANSP(BMTRX,NAT3,NAT3)
*      CALL vib_CEBDAB(HINT,NAT3,HESSIAN,NAT3,BMTRX)
      WRITE(6,9010)
*      NR = NAT3 - NDEG
*      NC = NAT3 - NDEG
*      CALL vib_PRDATA(HINT,DUM,NAT3,NAT3,NR,NC,.FALSE.,NPRI)
C
C     ----- WRITE OUT NORMAL MODES AND FREQUENCIES IN TERMS
C     ----- OF INTERNAL COORDINATES.
C
*      CALL vib_TRANSP(BMTRX,NAT3,NAT3)
*      CALL MATMUL(BMTRX,EIGVEC,HINT,
*     &            NAT3,NAT3,NAT3,NAT3,NAT3,NAT3)
      WRITE(6,9020)
C
*      CALL vib_PRDATA(HINT,EIGVAL,NAT3,NAT3,NR,NAT3,.TRUE.,NPRI)
C
 9010 FORMAT(///,10X,'------------------------------------------',/
     &           10X,'INTERNAL COORDINATE DESCRIPTION OF HESSIAN',/
     &           10X,'------------------------------------------')
 9020 FORMAT(///,10X,'------------------------------------------------',
     &         /,10X,'NORMAL MODE EIGENVECTORS IN INTERNAL COORDINATES',
     &         /,10X,'------------------------------------------------',
     &         /,10X,'      (Frequencies expressed in cm-1) ')
      END
