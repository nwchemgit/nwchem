c  $Id$
c
c  These are a bunch of utility routines that were in hessian/rhf_hessian.F
c  that were moved here since the hessian code no longer depends on them.
c  It may be that this needs to go to the util directory in the future, but
c  a lot of the functionality is not just utility routines.
c
      SUBROUTINE HND_REWFIL(NFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REWIND NFT
      RETURN
      END
c
      subroutine hnd_nwhnd_tran(dnwc,dhnd,ndim)
      implicit double precision (a-h,o-z)
#include "global.fh"
      parameter (mxnbf =2048)
      common/hnd_nwtohnd/inw_to_hnd(mxnbf)
      common/hnd_iofile/ir,iw,ip
      common/hnd_facntoh/fac_nwthnd(mxnbf)
      dimension dnwc(ndim,ndim),dhnd(ndim,ndim)
c
      logical out
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... nwhnd_tran ... '
      endif
c
c     ----- matrices tran. from nwchem to hondo shell order-----
c
      do j=1,ndim
        do i=1,ndim

           dhnd(inw_to_hnd(i),inw_to_hnd(j))=
     &     dnwc(i,j)/(fac_nwthnd(i)*fac_nwthnd(j))

        enddo
      enddo
c
      return
      end
c
      SUBROUTINE HND_DAREAD(IDAF,IODA,IX,NX,IDAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_DAFNAV/IXDDAF(2048),NAV10,NAV20
      DIMENSION IODA(1),IX(1)
      DIMENSION IXSDAF(1024)
      EQUIVALENCE (IXSDAF(1),IXDDAF(1))
      DATA IXMAX /1024/
C
      IF(IDAF.EQ.10) NAV10=IODA(IDAR+0)
      IF(IDAF.EQ.20) NAV20=IODA(IDAR+0)
      MAXIX=IXMAX*ISINGL
      LDAR =   NX*ISINGL
C
      MAX=0
   10 MIN=MAX+1
      MAX=MAX+MAXIX
      IF(MAX.GT.LDAR) MAX=LDAR
      IF(IDAF.EQ.10.AND.ISINGL.EQ.1) READ(IDAF,REC=NAV10) IXSDAF
      IF(IDAF.EQ.10.AND.ISINGL.EQ.2) READ(IDAF,REC=NAV10) IXDDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.1) READ(IDAF,REC=NAV20) IXSDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.2) READ(IDAF,REC=NAV20) IXDDAF
      DO 20 I=MIN,MAX
   20 IX(I)=IXDDAF(I-MIN+1)
      IF(IDAF.EQ.10) NAV10=NAV10+1
      IF(IDAF.EQ.20) NAV20=NAV20+1
      IF(MAX.LT.LDAR) GO TO 10
      RETURN
      END
C
      SUBROUTINE HND_DAWRIT(IDAF,IODA,IX,NX,IDAR,NAV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      PARAMETER (MXIODA=255)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MACHIN/ISINGL,NBITS
      COMMON/HND_DAFNAV/IXDDAF(2048),NAV10,NAV20
      COMMON/HND_DAFREC/LENDAR(MXIODA)
      DIMENSION ERRMSG(3)
      DIMENSION IODA(1),IX(1)
      DIMENSION IXSDAF(1024)
      EQUIVALENCE (IXSDAF(1),IXDDAF(1))
      DATA ERRMSG /'PROGRAM ','STOP IN ','-DAWRIT-'/
      DATA IXMAX  /1024/
C
      OUT=.FALSE.
C
      MAXIX=IXMAX*ISINGL
      LDAR =   NX*ISINGL
      IF(IODA(IDAR+0).NE.0) GO TO 100
C
C     ----- FIRST WRITE -----
C
      IODA(IDAR+0)=NAV
      IF(IDAF.EQ.10) NAV10=NAV
      IF(IDAF.EQ.20) NAV20=NAV
C
      IF(IDAF.EQ.10) THEN
         LENDAR(IDAR+0)=NX
         IF(OUT) THEN
            WRITE(IW,9998) IDAR,NX
         ENDIF
      ENDIF
C
      MAX=0
   10 MIN=MAX+1
      MAX=MAX+MAXIX
      IF(MAX.GT.LDAR) MAX=LDAR
      DO 20 I=MIN,MAX
   20 IXDDAF(I-MIN+1)=IX(I)
      IF(IDAF.EQ.10.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV10) IXSDAF
      IF(IDAF.EQ.10.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV10) IXDDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV20) IXSDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV20) IXDDAF
      IF(IDAF.EQ.10) NAV10=NAV10+1
      IF(IDAF.EQ.20) NAV20=NAV20+1
      IF(MAX.LT.LDAR) GO TO 10
C
      IF(IDAF.EQ.10) NAV=NAV10
      IF(IDAF.EQ.20) NAV=NAV20
      RETURN
C
C     ----- REWRITE -----
C
  100 CONTINUE
      IF(IDAF.EQ.10) NAV10=IODA(IDAR+0)
      IF(IDAF.EQ.20) NAV20=IODA(IDAR+0)
C
      IF(IDAF.EQ.10) THEN
         IF(OUT) THEN
            WRITE(IW,9997) IDAR,NX
         ENDIF
         NX0=LENDAR(IDAR+0)
         IF(NX.GT.NX0) THEN
            IF(OUT) THEN
               WRITE(IW,9999) IDAR,NX,NX0
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ELSEIF(NX.LT.NX0) THEN
            IF(OUT) THEN
               WRITE(IW,9999) IDAR,NX,NX0
            ENDIF
         ENDIF
      ENDIF
C
      MAX=0
  110 MIN=MAX+1
      MAX=MAX+MAXIX
      IF(MAX.GT.LDAR) MAX=LDAR
      DO 120 I=MIN,MAX
  120 IXDDAF(I-MIN+1)=IX(I)
      IF(IDAF.EQ.10.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV10) IXSDAF
      IF(IDAF.EQ.10.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV10) IXDDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.1) WRITE(IDAF,REC=NAV20) IXSDAF
      IF(IDAF.EQ.20.AND.ISINGL.EQ.2) WRITE(IDAF,REC=NAV20) IXDDAF
      IF(IDAF.EQ.10) NAV10=NAV10+1
      IF(IDAF.EQ.20) NAV20=NAV20+1
      IF(MAX.LT.LDAR) GO TO 110
C
      RETURN
 9999 FORMAT(' INCONSISTENT RECORD LENGTH FOR -IDAR- = ',I4,/,
     1       ' NX,NX0 = ',2I10)
 9998 FORMAT(' FIRST-WRITE, IDAR = ',I5,' LENGTH = ',I10)
 9997 FORMAT('    RE-WRITE, IDAR = ',I5,' LENGTH = ',I10)
      END
C
      subroutine hnd_hndnw_tran(dhnd,dnwc,ndim)
      implicit double precision (a-h,o-z)
#include "global.fh"
      parameter (mxnbf =2048)
      common/hnd_hndtonw/ihnd_to_nw(mxnbf)
      common/hnd_iofile/ir,iw,ip
      dimension dnwc(ndim,ndim),dhnd(ndim,ndim)
      common/hnd_facntoh/fac_nwthnd(mxnbf)
c
      logical out
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... hndnw_tran ... '
      endif
c
c     ----- matrices tran. from hondo to nwchem shell order-----
c
      do j=1,ndim
        do i=1,ndim
           dnwc(ihnd_to_nw(i),ihnd_to_nw(j))=
     &     dhnd(i,j)*(fac_nwthnd(ihnd_to_nw(i))*
     &                fac_nwthnd(ihnd_to_nw(j)) )
        enddo
      enddo
c
      return
      end
c
      SUBROUTINE HND_DSXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
#include "hnd_whermt.fh"
c
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
     1                               ,NI,NJ
     1                              ,CX,CY,CZ
      DIMENSION MIN(7),MAX(7)
      DATA MIN  /1,2,4, 7,11,16,22/
      DATA MAX  /1,3,6,10,15,21,28/
      DATA ZERO /0.0D+00/
C
      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      NPTS=(NI+NJ-2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 16 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)*T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (7,6,5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    6 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    7 GO TO (15,14,13,12,11,10,9,8),NJ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   12 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   13 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   14 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   15 CONTINUE
      XINT=XINT+PX
      YINT=YINT+PY
      ZINT=ZINT+PZ
   16 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_PREAD(IJK,XX,IX,NXX,MAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PACK2E
      COMMON/HND_PCKLAB/LABSIZ
      COMMON/HND_PCKOPT/NHEX,NTUPL,PACK2E
      COMMON/HND_INTOPT/NOPK,NOK,NOSQUR
      COMMON/HND_INTFIL/NINTMX
      DIMENSION XX(MAX),IX(LABSIZ*MAX)
      IF(PACK2E) GO TO 10
      READ(IJK) XX,IX,NX
      NXX=NX
      RETURN
   10 CALL HND_READPK(IJK,XX,XX,       NXX,NH,IERR,IEND)
      CALL HND_READPN(IJK,IX,IX,LABSIZ*NXX,NT,IERR,IEND)
      RETURN
      END
C
      SUBROUTINE HND_PRTRLAB(D,LAB,N)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "nwc_const.fh"
C
C     ----- PRINT OUT A TRIANGULAR MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION D(1),DD(10)
      character*8 LAB
      DIMENSION LAB(1)
C
      LIST=1
      IF(LIST.EQ.1) MAX=7
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (IW,9008)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 160 J = 1,N
      K = 0
      DO 140 I = IMIN,IMAX
      K = K+1
      II = MAX0( I, J)
      JJ = MIN0( I, J)
      IJ = (II*(II-1))/2 + JJ
  140 DD(K) = D(IJ)
      IF(LIST.EQ.1) WRITE(IW,9148) J,LAB(J),(DD(I),I = 1,K)
  160 CONTINUE
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(/)
 9128 FORMAT(14X,7(6X,I3,6X))
 9148 FORMAT(I3,1X,A8,2X,7F15.10)
      END
C
      SUBROUTINE HND_PRTRL(D,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
C
C     ----- PRINT OUT A TRIANGULAR MATRIX WITH LABELS
C
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXBFN=3072)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
      DIMENSION D(1),DD(10)
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (IW,9008)
c%%
C     IF(LIST.EQ.0) WRITE (IW,9028) (BFLAB(I),I = IMIN,IMAX)
C     IF(LIST.EQ.1) WRITE (IW,9128) (BFLAB(I),I = IMIN,IMAX)
C     IF(LIST.EQ.2) WRITE (IW,9228) (BFLAB(I),I = IMIN,IMAX)
c%%
      WRITE (IW,9008)
      DO 160 J = 1,N
      K = 0
      DO 140 I = IMIN,IMAX
      K = K+1
      II = MAX0( I, J)
      JJ = MIN0( I, J)
      IJ = (II*(II-1))/2 + JJ
  140 DD(K) = D(IJ)
c%%
C     IF(LIST.EQ.0) WRITE (IW,9048) J,BFLAB(J),(DD(I),I = 1,K)
C     IF(LIST.EQ.1) WRITE (IW,9148) J,BFLAB(J),(DD(I),I = 1,K)
C     IF(LIST.EQ.2) WRITE (IW,9248) J,BFLAB(J),(DD(I),I = 1,K)
c%%
  160 CONTINUE
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(15X,10(2X,A8,1X))
 9048 FORMAT(I5,2X,A8,10F11.5)
 9128 FORMAT(15X,7(4X,A8,3X))
 9148 FORMAT(I5,2X,A8,7F15.10)
 9228 FORMAT(15X,7(4X,A8,3X))
 9248 FORMAT(I5,2X,A8,7E15.8)
      END
c
      SUBROUTINE HND_SPRTR(D,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
      DIMENSION D(1),DD(10)
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 160 J = 1,N
      K = 0
      DO 140 I = IMIN,IMAX
      K = K+1
C     II = MAX0( I, J)
C     JJ = MIN0( I, J)
C     IJ = (II*(II-1))/2 + JJ
      IJ = (J-1)*N + I
  140 DD(K) = D(IJ)
      IF(LIST.EQ.0) WRITE (IW,9048) J,(DD(I),I = 1,K)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(DD(I),I = 1,K)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(DD(I),I = 1,K)
  160 CONTINUE
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
      END
c
      SUBROUTINE HND_SECOND(SEC,WSEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_CLOCK0/SEC0,WSEC0
      COMMON/HND_CLOCKS/DAT,TIM,CPUSEC,TOTSEC,DELCPU,DELTOT
      DIMENSION DIAL(6)
      EQUIVALENCE (DIAL(1),DAT)
C
C     ----- RETURN ELAPSED - CPU- TIME IN SECONDS IN - SEC- -----
C     ----- RETURN ELAPSED -WALL- TIME IN SECONDS IN -WSEC- -----
C
      CALL HND_SYSCLK(DIAL)
       SEC= CPUSEC
       SEC= SEC- SEC0
      WSEC= TOTSEC
      WSEC=WSEC-WSEC0
      RETURN
      END
C
      SUBROUTINE HND_JKWRYS(RWV,ABV,NUMG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_rys.fh"
      DIMENSION RWV(2,NUMG,1),ABV(5,1)
C
      IF(MROOTS.GT.5) GO TO 100
C
C     ----- MROOTS .LE. 5 -----
C
      DO 20 NG=1,NUMG
      XX=ABV(5,NG)
      IF(MROOTS.LE.3) CALL HND_RT123
      IF(MROOTS.EQ.4) CALL HND_ROOT4
      IF(MROOTS.EQ.5) CALL HND_ROOT5
      DO 10 NR=1,MROOTS
      RWV(1,NG,NR)=U(NR)
      RWV(2,NG,NR)=W(NR)
   10 CONTINUE
   20 CONTINUE
      RETURN
C
  100 CONTINUE
C
C     ----- MROOTS .GT. 5 -----
C
      DO 120 NG=1,NUMG
      YY=ABV(5,NG)
      CALL HND_DROOT
      DO 110 NR=1,MROOTS
      RWV(1,NG,NR)=U9(NR)
      RWV(2,NG,NR)=W9(NR)
  110 CONTINUE
  120 CONTINUE
      RETURN
      END
C
      SUBROUTINE HND_JKBCDF(B00,B01,B10,C00,D00,F00,DIJ,DKL,
     1                  ABV,CV,RWV,NUMG,NROOTS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/HND_SHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      DIMENSION B00(NUMG,NROOTS,3),B01(NUMG,NROOTS,3),B10(NUMG,NROOTS,3)
      DIMENSION C00(NUMG,NROOTS,3)
      DIMENSION D00(NUMG,NROOTS,3)
      DIMENSION F00(NUMG,NROOTS,3),DIJ(NUMG,NROOTS,3),DKL(NUMG,NROOTS,3)
      DIMENSION ABV(5,1),CV(18,1)
      DIMENSION RWV(2,NUMG,NROOTS)
      DATA PT5,ONE /0.5D+00,1.0D+00/
C
      DO 40 NR=1,NROOTS
      DO 30 NG=1,NUMG
      AA =ABV(1,NG)
      BB =ABV(2,NG)
      RHO=ABV(3,NG)
      QAB=ABV(4,NG)
      UU =RHO*RWV(1,NG,NR)
      WW =    RWV(2,NG,NR)
      AAUU=AA+UU
      BBUU=BB+UU
      F00(NG,NR,1)=WW*QAB
      F00(NG,NR,2)=ONE
      F00(NG,NR,3)=ONE
      DUM2=PT5/(AA*BB+UU*(AA+BB))
      AUDUM=AAUU*DUM2
      BUDUM=BBUU*DUM2
       UDUM=  UU*DUM2
      B00(NG,NR,1)= UDUM
      B00(NG,NR,2)= UDUM
      B00(NG,NR,3)= UDUM
      B01(NG,NR,1)=AUDUM
      B01(NG,NR,2)=AUDUM
      B01(NG,NR,3)=AUDUM
      B10(NG,NR,1)=BUDUM
      B10(NG,NR,2)=BUDUM
      B10(NG,NR,3)=BUDUM
       UDUM= UDUM+ UDUM
      IF(MMAXS) GO TO 10
      AUDUM=AUDUM+AUDUM
      D00(NG,NR,1)= UDUM*CV( 1,NG) + AUDUM*CV( 2,NG)
      D00(NG,NR,2)= UDUM*CV( 3,NG) + AUDUM*CV( 4,NG)
      D00(NG,NR,3)= UDUM*CV( 5,NG) + AUDUM*CV( 6,NG)
   10 IF(NMAXS) GO TO 20
      BUDUM=BUDUM+BUDUM
      C00(NG,NR,1)= UDUM*CV( 8,NG) + BUDUM*CV( 7,NG)
      C00(NG,NR,2)= UDUM*CV(10,NG) + BUDUM*CV( 9,NG)
      C00(NG,NR,3)= UDUM*CV(12,NG) + BUDUM*CV(11,NG)
   20 CONTINUE
C
   30 CONTINUE
   40 CONTINUE
C
      DO 60 NR=1,NROOTS
      DO 50 NG=1,NUMG
      DIJ(NG,NR,1)=CV(13,NG)
      DIJ(NG,NR,2)=CV(14,NG)
      DIJ(NG,NR,3)=CV(15,NG)
      DKL(NG,NR,1)=CV(16,NG)
      DKL(NG,NR,2)=CV(17,NG)
      DKL(NG,NR,3)=CV(18,NG)
   50 CONTINUE
   60 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKGNMS(GNM,NG,NMAX,MMAX,
     1 B00,B01,B10,C00,D00,F00)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/HND_SHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      DIMENSION GNM(NG,NMAX,MMAX)
      DIMENSION C00(NG),D00(NG),F00(NG)
      DIMENSION B00(NG,1),B01(NG,1),B10(NG,1)
C
C     ----- G(0,0) -----
C
      DO 10 IG=1,NG
      GNM(IG,1,1)=F00(IG)
   10 CONTINUE
      IF(NMAXS.AND.MMAXS) RETURN
      IF(NMAXS) GO TO 30
C
C     ----- G(1,0) = C00 * G(0,0) -----
C
      DO 20 IG=1,NG
      GNM(IG,2,1)=C00(IG)*GNM(IG,1,1)
   20 CONTINUE
C
   30 CONTINUE
      IF(MMAXS) GO TO 60
C
C     ----- G(0,1) = D00 * G(0,0) -----
C
      DO 40 IG=1,NG
      GNM(IG,1,2)=D00(IG)*GNM(IG,1,1)
   40 CONTINUE
      IF(NMAXS) GO TO 60
C
C     ----- G(1,1) = B00 * G(0,0) + D00 * G(1,0) -----
C
      DO 50 IG=1,NG
      GNM(IG,2,2)=B00(IG,1)*GNM(IG,1,1)+D00(IG)*GNM(IG,2,1)
   50 CONTINUE
C
   60 CONTINUE
      MAX=MAX0(NMAX-1,MMAX-1)
      DO 70 M=2,MAX
      DO 70 IG=1,NG
   70 B00(IG,M)=B00(IG,M-1)+B00(IG,1)
C
      IF(NMAXP) GO TO 120
C
C     ----- G(N+1,0) = N * B10 * G(N-1,0) + C00 * G(N,0) -----
C
      DO 80 N=2,NMAX-1
      DO 80 IG=1,NG
      B10(IG,N)=B10(IG,N-1)+B10(IG,1)
   80 CONTINUE
      DO 90 N=2,NMAX-1
      DO 90 IG=1,NG
      GNM(IG,N+1,1)=B10(IG,N-1)*GNM(IG,N-1,1)+C00(IG)*GNM(IG,N,1)
   90 CONTINUE
      IF(MMAXS) GO TO 110
C
C     ----- G(N,1) = N * B00 * G(N-1,0) + D00 * G(N,0) -----
C
      DO 100 N=2,NMAX-1
      DO 100 IG=1,NG
      GNM(IG,N+1,2)=B00(IG,N)*GNM(IG,N,1)+D00(IG)*GNM(IG,N+1,1)
  100 CONTINUE
C
  110 CONTINUE
C
  120 CONTINUE
      IF(MMAXP) GO TO 170
C
C     ----- G(0,M+1) = M * B01 * G(0,M-1) + D00 * G(O,M) -----
C
      DO 130 M=2,MMAX-1
      DO 130 IG=1,NG
      B01(IG,M)=B01(IG,M-1)+B01(IG,1)
  130 CONTINUE
      DO 140 M=2,MMAX-1
      DO 140 IG=1,NG
      GNM(IG,1,M+1)=B01(IG,M-1)*GNM(IG,1,M-1)+D00(IG)*GNM(IG,1,M)
  140 CONTINUE
      IF(NMAXS) GO TO 160
C
C     ----- G(1,M) = M * B00 * G(0,M-1) + C00 * G(0,M) -----
C
      DO 150 M=2,MMAX-1
      DO 150 IG=1,NG
      GNM(IG,2,M+1)=B00(IG,M)*GNM(IG,1,M)+C00(IG)*GNM(IG,1,M+1)
  150 CONTINUE
C
  160 CONTINUE
  170 IF(NMAXP.OR.MMAXP) RETURN
C
C     ----- G(N+1,M) = N * B10 * G(N-1,M  ) -----
C                    +     C00 * G(N  ,M  )
C                    + M * B00 * G(N  ,M-1)
C
      DO 180 M=2,MMAX-1
      DO 180 N=2,NMAX-1
      DO 180 IG=1,NG
      GNM(IG,N+1,M+1)=B10(IG,N-1)*GNM(IG,N-1,M+1)+
     1                C00(IG    )*GNM(IG,N  ,M+1)+
     2                B00(IG,M  )*GNM(IG,N  ,M  )
  180 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKGNMV(GNM,NG,NMAX,MMAX,
     1 B00,B01,B10,C00,D00,F00)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/HND_SHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      DIMENSION GNM(NG,NMAX,MMAX)
      DIMENSION C00(NG),D00(NG),F00(NG)
      DIMENSION B00(NG,1),B01(NG,1),B10(NG,1)
C
C     ----- G(0,0) -----
C
      DO 10 IG=1,NG
      GNM(IG,1,1)=F00(IG)
   10 CONTINUE
      IF(NMAXS.AND.MMAXS) RETURN
      IF(NMAXS) GO TO 30
C
C     ----- G(1,0) = C00 * G(0,0) -----
C
      DO 20 IG=1,NG
      GNM(IG,2,1)=C00(IG)*GNM(IG,1,1)
   20 CONTINUE
C
   30 CONTINUE
      IF(MMAXS) GO TO 60
C
C     ----- G(0,1) = D00 * G(0,0) -----
C
      DO 40 IG=1,NG
      GNM(IG,1,2)=D00(IG)*GNM(IG,1,1)
   40 CONTINUE
      IF(NMAXS) GO TO 60
C
C     ----- G(1,1) = B00 * G(0,0) + D00 * G(1,0) -----
C
      DO 50 IG=1,NG
      GNM(IG,2,2)=B00(IG,1)*GNM(IG,1,1)+D00(IG)*GNM(IG,2,1)
   50 CONTINUE
C
   60 CONTINUE
      MAX=MAX0(NMAX-1,MMAX-1)
      DO 70 IG=1,NG
      DO 70 M=2,MAX
   70 B00(IG,M)=B00(IG,M-1)+B00(IG,1)
C
      IF(NMAXP) GO TO 120
C
C     ----- G(N+1,0) = N * B10 * G(N-1,0) + C00 * G(N,0) -----
C
      DO 80 IG=1,NG
      DO 80 N=2,NMAX-1
      B10(IG,N)=B10(IG,N-1)+B10(IG,1)
   80 CONTINUE
      DO 90 IG=1,NG
      DO 90 N=2,NMAX-1
      GNM(IG,N+1,1)=B10(IG,N-1)*GNM(IG,N-1,1)+C00(IG)*GNM(IG,N,1)
   90 CONTINUE
      IF(MMAXS) GO TO 110
C
C     ----- G(N,1) = N * B00 * G(N-1,0) + D00 * G(N,0) -----
C
      DO 100 IG=1,NG
      DO 100 N=2,NMAX-1
      GNM(IG,N+1,2)=B00(IG,N)*GNM(IG,N,1)+D00(IG)*GNM(IG,N+1,1)
  100 CONTINUE
C
  110 CONTINUE
C
  120 CONTINUE
      IF(MMAXP) GO TO 170
C
C     ----- G(0,M+1) = M * B01 * G(0,M-1) + D00 * G(O,M) -----
C
      DO 130 IG=1,NG
      DO 130 M=2,MMAX-1
      B01(IG,M)=B01(IG,M-1)+B01(IG,1)
  130 CONTINUE
      DO 140 IG=1,NG
      DO 140 M=2,MMAX-1
      GNM(IG,1,M+1)=B01(IG,M-1)*GNM(IG,1,M-1)+D00(IG)*GNM(IG,1,M)
  140 CONTINUE
      IF(NMAXS) GO TO 160
C
C     ----- G(1,M) = M * B00 * G(0,M-1) + C00 * G(0,M) -----
C
      DO 150 IG=1,NG
      DO 150 M=2,MMAX-1
      GNM(IG,2,M+1)=B00(IG,M)*GNM(IG,1,M)+C00(IG)*GNM(IG,1,M+1)
  150 CONTINUE
C
  160 CONTINUE
  170 IF(NMAXP.OR.MMAXP) RETURN
C
C     ----- G(N+1,M) = N * B10 * G(N-1,M  ) -----
C                    +     C00 * G(N  ,M  )
C                    + M * B00 * G(N  ,M-1)
C
      DO 180 IG=1,NG
      DO 180 N=2,NMAX-1
      DO 180 M=2,MMAX-1
      GNM(IG,N+1,M+1)=B10(IG,N-1)*GNM(IG,N-1,M+1)+
     1                C00(IG    )*GNM(IG,N  ,M+1)+
     2                B00(IG,M  )*GNM(IG,N  ,M  )
  180 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKXYZS(GIJKL,HIJKL,GNKL,HNKL,FNKL,GNM,HNM,
     1 NG,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     2 EXPNDI,EXPNDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EXPNDI,EXPNDK
      DIMENSION GIJKL(NG*NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION HIJKL(NG*NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION  GNKL(NG,NLMAX,NKMAX,NMAX)
      DIMENSION  HNKL(NG*NLMAX*NKMAX,NMAX)
      DIMENSION  FNKL(NG*NLMAX*NKMAX*NMAX)
      DIMENSION   GNM(NG,NMAX,MMAX)
      DIMENSION   DIJ(NG)
      DIMENSION   DKL(NG)
C
C     ----- G(N,K,L) -----
C
      IF(EXPNDK) GO TO 40
C
      DO 30 NK=1,NKMAX
      DO 10 NL=1,NLMAX
      DO 10  N=1,NMAX
      DO 10 IG=1,NG
   10 GNKL(IG,NL,NK,N)=GNM(IG,N,NL)
      IF(NK.EQ.NKMAX) GO TO 30
      MAX=MMAX-NK
      DO 20  M=1,MAX
      DO 20  N=1,NMAX
      DO 20 IG=1,NG
   20 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   30 CONTINUE
C
      GO TO 100
   40 CONTINUE
C
      DO 70 NL=1,NLMAX
      DO 50 NK=1,NKMAX
      DO 50  N=1,NMAX
      DO 50 IG=1,NG
   50 GNKL(IG,NL,NK,N)=GNM(IG,N,NK)
      IF(NL.EQ.NLMAX) GO TO 70
      MAX=MMAX-NL
      DO 60  M=1,MAX
      DO 60  N=1,NMAX
      DO 60 IG=1,NG
   60 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   70 CONTINUE
C
  100 CONTINUE
C
C     ----- G(I,J,K,L) -----
C
      IF(EXPNDI) GO TO 140
C
      DO 130 NI=1,NIMAX
      DO 110 IGLKJ=1,NG*NLMAX*NKMAX*NJMAX
  110 HIJKL(IGLKJ,NI)=FNKL(IGLKJ)
      IF(NI.EQ.NIMAX) GO TO 130
      MAX=NMAX-NI
      DO 120  N=1,MAX
      DO 120 NK=1,NKMAX
      DO 120 NL=1,NLMAX
      DO 120 IG=1,NG
  120 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  130 CONTINUE
C
      RETURN
  140 CONTINUE
C
      DO 170 NJ=1,NJMAX
      DO 150 NI=1,NIMAX
      DO 150 IGLK=1,NG*NLMAX*NKMAX
  150 GIJKL(IGLK,NJ,NI)=HNKL(IGLK,NI)
      IF(NJ.EQ.NJMAX) GO TO 170
      MAX=NMAX-NJ
      DO 160  N=1,MAX
      DO 160 NK=1,NKMAX
      DO 160 NL=1,NLMAX
      DO 160 IG=1,NG
  160 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  170 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_JKXYZV(GIJKL,HIJKL,GNKL,HNKL,FNKL,GNM,HNM,
     1 NG,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     2 EXPNDI,EXPNDK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EXPNDI,EXPNDK
      DIMENSION GIJKL(NG*NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION HIJKL(NG*NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION  GNKL(NG,NLMAX,NKMAX,NMAX)
      DIMENSION  HNKL(NG*NLMAX*NKMAX,NMAX)
      DIMENSION  FNKL(NG*NLMAX*NKMAX*NMAX)
      DIMENSION   GNM(NG,NMAX,MMAX)
      DIMENSION   DIJ(NG)
      DIMENSION   DKL(NG)
C
C     ----- G(N,K,L) -----
C
      IF(EXPNDK) GO TO 40
C
      DO 30 NK=1,NKMAX
      DO 10 IG=1,NG
      DO 10 NL=1,NLMAX
      DO 10  N=1,NMAX
   10 GNKL(IG,NL,NK,N)=GNM(IG,N,NL)
      IF(NK.EQ.NKMAX) GO TO 30
      MAX=MMAX-NK
      DO 20 IG=1,NG
      DO 20  M=1,MAX
      DO 20  N=1,NMAX
   20 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   30 CONTINUE
C
      GO TO 100
   40 CONTINUE
C
      DO 70 NL=1,NLMAX
      DO 50 IG=1,NG
      DO 50 NK=1,NKMAX
      DO 50  N=1,NMAX
   50 GNKL(IG,NL,NK,N)=GNM(IG,N,NK)
      IF(NL.EQ.NLMAX) GO TO 70
      MAX=MMAX-NL
      DO 60 IG=1,NG
      DO 60  N=1,NMAX
      DO 60  M=1,MAX
   60 GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   70 CONTINUE
C
  100 CONTINUE
C
C     ----- G(I,J,K,L) -----
C
      IF(EXPNDI) GO TO 140
C
      DO 130 NI=1,NIMAX
      DO 110 IGLKJ=1,NG*NLMAX*NKMAX*NJMAX
  110 HIJKL(IGLKJ,NI)=FNKL(IGLKJ)
      IF(NI.EQ.NIMAX) GO TO 130
      MAX=NMAX-NI
      DO 120 IG=1,NG
      DO 120 NL=1,NLMAX
      DO 120 NK=1,NKMAX
      DO 120  N=1,MAX
  120 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  130 CONTINUE
C
      RETURN
  140 CONTINUE
C
      DO 170 NJ=1,NJMAX
      DO 150 IGLK=1,NG*NLMAX*NKMAX
      DO 150 NI=1,NIMAX
  150 GIJKL(IGLK,NJ,NI)=HNKL(IGLK,NI)
      IF(NJ.EQ.NJMAX) GO TO 170
      MAX=NMAX-NJ
      DO 160 IG=1,NG
      DO 160 NL=1,NLMAX
      DO 160 NK=1,NKMAX
      DO 160  N=1,MAX
  160 GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+GNKL(IG,NL,NK,N+1)
  170 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE HND_READPK(IS,XX,YY,NXX,NH,IERR,IEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-READPK-'/
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
      END
c
      SUBROUTINE HND_READPN(IS,IX,IY,NNX,NT,IERR,IEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-READPN-'/
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
      END
c
      SUBROUTINE HND_SYSCLK(DIAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TEN3=1.0D+03)
C
C     ----- -CPU- AND -WALL- TIME -----
C
C     DIAL(1) = DATE
C     DIAL(2) = TIME
C     DIAL(3) = CPUSEC
C     DIAL(4) = TOTSEC
C     DIAL(5) = DELCPU
C     DIAL(6) = DELTOT
C
      DIMENSION DIAL(6)
*     REAL*4 CPUTIM,ETIME_
*     TYPE TB_TYPE
*        SEQUENCE
*        REAL*4 USRTIM
*        REAL*4 SYSTIM
*     END TYPE
c     TYPE (TB_TYPE) ETIME_STRUCT
c     CPUTIM =ETIME_(ETIME_STRUCT)
c     DIAL(3)=CPUTIM
*     TOTTIM =TIMEF()/TEN3
*     DIAL(4)=TOTTIM
      RETURN
      END
