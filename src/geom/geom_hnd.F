      SUBROUTINE HND_ICOORD(GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "util.fh"
c
      logical status
      integer geom
      logical  geom_zmt_opn_fil
C
C     ---                 CONSTRUCT -Z- MATRIX                     ---
C     ---             CALCULATE ATOM-ATOM DISTANCES                ---
C     --- -NZMAT-  =  NUMBER OF ELEMENTS OF -IZMAT- GIVEN IN $ZMAT ---
C     --- -NZVAR-  =  NUMBER OF INTERNAL COORDINATES FROM -IZMAT-  ---
C     --- -NVAR-   =  (3N-6) OR (3N-5) OR -NZMOD-  .LE. -NZVAR-    ---
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL ZDONE
      LOGICAL DBUG
      LOGICAL SOME
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM)
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION R(MXATOM)
      DATA ONE   /1.0D+00/
      DATA UNIT  /0.52917715D+00/
C
      DBUG=.FALSE.
      SOME=.FALSE.
      SOME=SOME.OR.DBUG
      IF(DBUG) THEN
         WRITE(IW,*) 'IN HND_ICOORD'
         WRITE(IW,*) 'IFZMAT',IFZMAT
      ENDIF
C
C     ----- in -nwchem- , some MUST BE true -----
C
      SOME=.TRUE.
C
      IF(IFZMAT.NE.0) THEN
C
C     ----- LOOK FOR -Z- MATRIX DATA AND READ IT IN -----
C
         IF(SOME) THEN
            CALL HND_RDZMAT(NAT)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,*) 'IFZMAT',IFZMAT
         WRITE(IW,*) 'NZMAT ',NZMAT
      ENDIF
      IF(IFZMAT.NE.0.AND.NZMAT.GT.0) THEN
C
C     -----    THERE IS A -Z- MATRIX    -----
C     ----- CALCULATE -B- AND -B**(-1)- -----
C
         STATUS=GEOM_ZMT_OPN_FIL(.TRUE.,.TRUE.,.FALSE.,.FALSE.)
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('GEOM_ZMT_OPN_FIL?',911)
         ENDIF
         CALL HND_BANDBI(GEOM,NAT,C)
****         CALL HND_ZMTPRT(GEOM)  ! RJH - no need to print every time
      ELSEIF(IFZMAT.NE.0.AND.NZMAT.EQ.0) THEN
C
C     -----       THERE IS NO -Z- MATRIX     -----
C     ----- TRY TO CREAT ONE ; IF SUCCESSFUL -----
C     -----    CALCULATE -B- AND -B**(-1)-   -----
C
         CALL HND_ZMTMAK(ZDONE,geom,
     $        util_print('zmat','print_default'))
         IF(ZDONE) THEN
            STATUS=GEOM_ZMT_OPN_FIL(.TRUE.,.TRUE.,.FALSE.,.FALSE.)
            IF(.NOT.STATUS) THEN
               CALL ERRQUIT('GEOM_ZMT_OPN_FIL?',911)
            ENDIF
            CALL HND_BANDBI(GEOM,NAT,C)
****  CALL HND_ZMTPRT(GEOM)
         ELSE
C
C     -----  DEFINE ROTATIONS+TRANSLATIONS  -----
C           BY MEANS OF HND_ECKARTs CONDITIONS
C
            NCART=3*NAT
            NZVAR=3*NAT
            NVAR =3*NAT
            NZMAT=0
            NCOOR=MAX(NCART,NZVAR)
            CALL HND_TRNROT(NAT,C,NCART,NCOOR)
C
C     ----- PRINT INTERNUCLEAR DISTANCES -----
C
            IF(SOME) THEN
               IPASS=1
               FAC=UNIT
               WRITE(IW,9995)
   10          CONTINUE
               JMAX=0
   20          JMIN=JMAX+1
               JMAX=JMAX+10
               IF(JMAX.GT.NAT) JMAX=NAT
               WRITE(IW,9999)
               WRITE(IW,9998) (ANAM(J),BNAM(J),J=JMIN,JMAX)
               WRITE(IW,9999)
               DO I=1,NAT
                  DO J=JMIN,JMAX
                     RR=(C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2
     1                                    +(C(3,I)-C(3,J))**2
                     R(J)= SQRT(RR)*FAC
                  ENDDO
                  WRITE(IW,9997) I,ANAM(I),BNAM(I),(R(J),J=JMIN,JMAX)
               ENDDO
               IF(JMAX.LT.NAT) GO TO 20
               IPASS=2
               IF(IPASS.EQ.2) RETURN
               FAC=ONE
               WRITE(IW,9996)
               GO TO 10
            ENDIF
C
         ENDIF
      ELSEIF(IFZMAT.EQ.0) THEN
C
C     -----  DEFINE ROTATIONS+TRANSLATIONS  -----
C           BY MEANS OF HND_ECKARTs CONDITIONS
C
         NCART=3*NAT
         NZVAR=3*NAT
         NVAR =3*NAT
         NCOOR=MAX(NCART,NZVAR)
         CALL HND_TRNROT(NAT,C,NCART,NCOOR)
C
C     ----- PRINT INTERNUCLEAR DISTANCES -----
C
         IF(SOME) THEN
            IPASS=1
            FAC=UNIT
            WRITE(IW,9995)
  110       CONTINUE
            JMAX=0
  120       JMIN=JMAX+1
            JMAX=JMAX+10
            IF(JMAX.GT.NAT) JMAX=NAT
            WRITE(IW,9999)
            WRITE(IW,9998) (ANAM(J),BNAM(J),J=JMIN,JMAX)
            WRITE(IW,9999)
            DO I=1,NAT
               DO J=JMIN,JMAX
                  RR=(C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2
     1                                 +(C(3,I)-C(3,J))**2
                  R(J)= SQRT(RR)*FAC
               ENDDO
               WRITE(IW,9997) I,ANAM(I),BNAM(I),(R(J),J=JMIN,JMAX)
            ENDDO
            IF(JMAX.LT.NAT) GO TO 120
            IPASS=2
            IF(IPASS.EQ.2) RETURN
            FAC=ONE
            WRITE(IW,9996)
            GO TO 110
         ENDIF
C
      ENDIF
      RETURN
 9999 FORMAT(/)
 9998 FORMAT(15X,10(1X,A8,A2))
 9997 FORMAT(1X,I3,1X,A8,A2,10(1X,F8.4,2X))
 9995 FORMAT(/,10X,'INTERNUCLEAR DISTANCES (ANGS.)',/,10X,30(1H-))
 9996 FORMAT(/,10X,'INTERNUCLEAR DISTANCES ( A.U.)',/,10X,30(1H-))
      END
      SUBROUTINE HND_TRNROT(NAT,C,NCART,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
C
      LOGICAL OUT 
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
C     ----- MEMORY POINTERS -----
C
C     I10=1
C     I20=I10+NCART*NCART
C     I30=I20+NCART*3     
C     I40=I30+NCART*3     
C     I50=I40+NDIM*NDIM  
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NCART*NCART,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NCART*3    ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NCART*3    ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I30  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NDIM*NDIM  ,'MEM I40 ',I_I40,I40))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I40  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I50 ',I_I50,I50))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I50  FAILED',911)
C
C     ----- CREATE HND_ECKARTs PROJECTOR AND WRITE ON -94- -----
C
      CALL HND_ECKART(NAT,C,DBL_MB(I10),DBL_MB(I20),
     &                      DBL_MB(I30),DBL_MB(I20),NCART,DBUG)
C
      DO J=1,NCART
         DO I=1,NCART
            DBL_MB(I+(J-1)*NDIM+I40-1)=DBL_MB(I+(J-1)*NCART+I10-1)
         ENDDO
      ENDDO
C
      NDUM=NDIM*NDIM    
      IF(OUT) THEN
         WRITE(IW,9999)
         CALL HND_PRSQ(DBL_MB(I40),NCART,NCART,NDIM)
      ENDIF
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I50))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I50  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I40))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I40  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -HND_TRNROT- PROJECTOR = ')
      END
      SUBROUTINE HND_ECKART(NAT,C,P,T,R,TR,NCART,OUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,1)
      DIMENSION G(3)
      DIMENSION P(NCART,*),T(3,NAT,*),R(3,NAT,*),TR(NCART,*)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
      DATA TOL  /1.0D-09/ 
C
      DO JTR=1,6
         DO ICART=1,NCART
            TR(ICART,JTR)=ZERO
         ENDDO
      ENDDO
C
C     ----- DEFINE CENTER OF MASS ( UNIT MASS ) -----
C
      DO IXYZ=1,3
         DUM=ZERO
         DO IAT=1,NAT
            DUM=DUM+C(IXYZ,IAT)
         ENDDO
         DUM=DUM/DBLE(NAT)
         G(IXYZ)=DUM
      ENDDO
C
C     ----- DEFINE TRANSLATIONS -----
C
      DO IXYZ=1,3
         DO IAT=1,NAT
            T(IXYZ,IAT,IXYZ)=ONE
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9997)
         CALL HND_PRSQ(TR,3,NCART,NCART)
      ENDIF
C
C     ----- DEFINE ROTATIONS -----
C
      DO IXYZ=1,3
         IF(IXYZ.EQ.1) THEN
            DO IAT=1,NAT
               R(2,IAT,IXYZ)=-C(3,IAT)+G(3)
               R(3,IAT,IXYZ)= C(2,IAT)-G(2)
            ENDDO
         ELSEIF(IXYZ.EQ.2) THEN
            DO IAT=1,NAT
               R(1,IAT,IXYZ)= C(3,IAT)-G(3)
               R(3,IAT,IXYZ)=-C(1,IAT)+G(1)
            ENDDO
         ELSEIF(IXYZ.EQ.3) THEN
            DO IAT=1,NAT
               R(1,IAT,IXYZ)=-C(2,IAT)+G(2)
               R(2,IAT,IXYZ)= C(1,IAT)-G(1)
            ENDDO
         ENDIF
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9996)
         CALL HND_PRSQ(TR,6,NCART,NCART)
      ENDIF
C
C     ----- ORTHONORMALIZE -----
C
      NTR=0
      DO ITR=1,6
         DUM=dnrm2(NCART,TR(1,ITR),1)
         IF(ABS(DUM).GT.TOL) THEN
            NTR=NTR+1
            DUM=ONE/DUM
            CALL dscal(NCART,DUM,TR(1,ITR),1)
            IF(ITR.LT.6) THEN
               DO JTR=ITR+1,6
                  DUM= -ddot(NCART,    TR(1,ITR),1,TR(1,JTR),1)
                  CALL daxpy(NCART,DUM,TR(1,ITR),1,TR(1,JTR),1)
               ENDDO
            ENDIF
         ELSE
            DO ICART=1,NCART
               TR(ICART,ITR)=ZERO
            ENDDO
         ENDIF
      ENDDO
      IF(OUT.AND.NTR.NE.6) THEN
         WRITE(IW,9999) NTR
      ENDIF
      IF(OUT) THEN
         WRITE(IW,9995)
         CALL HND_PRSQ(TR,6,NCART,NCART)
      ENDIF
C
C     ----- CREATE PROJECTOR ONTO INTERNAL SPACE -----
C
      DO JCART=1,NCART
         DO ICART=1,NCART
            P(ICART,JCART)=ZERO
         ENDDO
         P(JCART,JCART)=ONE
         DO ICART=1,NCART
            DUM=ZERO
            DO ITR=1,6
               DUM=DUM+TR(ICART,ITR)*TR(JCART,ITR)
            ENDDO
            P(ICART,JCART)=P(ICART,JCART)-DUM
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9998)
         CALL HND_PRSQ(P,NCART,NCART,NCART)
      ENDIF
C
      RETURN
 9999 FORMAT(/,' IN -HND_ECKART- , ONLY ',I1,' TRAN+ROT WERE FOUND .')
 9998 FORMAT(/,' HND_ECKART PROJECTOR',/,' ----------------')
 9997 FORMAT(/,' TRANSLATIONS',/,' ------------')
 9996 FORMAT(/,' TRANS. + ROT.',/,' -------------')
 9995 FORMAT(/,' TRANS. + ROT. (ORTHONORM)',/,' ',25(1h-))
      END
      SUBROUTINE HND_RDZMAT(NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- READ IN Z-MATRIX INPUT -----
C
      PARAMETER (    NA=10 )
      PARAMETER (MXNIJS=500)
      PARAMETER (MXIZMT=6000)
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL OUT
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      COMMON/HND_ZMTDAT/ZMAT(MXCOOR),IZMAT(MXIZMT),NIZMAT
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION SIJ(MXNIJS),IJS(2,MXNIJS)
      DIMENSION A(3,NA)
      DIMENSION ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-RDZMAT-'/
      DATA LINEAR  /0/
      DATA NZMOD   /0/
      DATA NATURAL /2/
      DATA NPFLG   /0/
      DATA ZERO    /0.0D+00/
C
C     NAMELIST /ZMAT/ IZMAT,IJS,SIJ,A,LINEAR,NZMOD,NATURAL,NPFLG
C     NAMELIST /ZMT/  IZMAT,IJS,SIJ,A,LINEAR,NZMOD,NATURAL,NPFLG
C
      NZMAT=0
      NZVAR=3*NAT
      NVAR =3*NAT
      IF(NIZMAT.LE.0) THEN
         DO I=1,MXIZMT
            IZMAT(I)=0
         ENDDO
      ENDIF
      DO I=1,MXNIJS
         IJS(1,I)=0
         IJS(2,I)=0
         SIJ(I)=ZERO
      ENDDO
      DO J=1,NA
         DO I=1,3
            A(I,J)=ZERO
         ENDDO
      ENDDO
      NATZMT=-1
C
C     ----- READ NAMELIST -$ZMAT- OR NAMELIST -$ZMT- -----
C           BOTH PROVIDE THE SAME INPUT DATA. ONLY
C           $ZMAT HAS PRIORITY OVER $ZMT IF BOTH ARE
C           GIVEN IN INPUT DATA.
C
C     REWIND IR
C     READ(IR,ZMAT,END=100,ERR=100)
C     GO TO 110
C 100 CONTINUE
C     REWIND IR
C     READ(IR, ZMT,END=140,ERR=140)
C 110 CONTINUE
C
      IF(NPFLG.NE.1.AND.NPFLG.NE.2) THEN
         NPFLG=0
      ENDIF
      IPRZMT=NPFLG
      IF(NATURAL.NE.1.AND.NATURAL.NE.2) THEN
         NATURAL=2
      ENDIF
      IF(NATURAL.GT.0) THEN
         NATZMT=NATURAL
      ENDIF
C
      OUT=.FALSE.
      OUT=OUT.OR.NPRINT.EQ.1
      OUT=OUT.OR.IPRZMT.GE.1
      IF(OUT) THEN
         WRITE(IW,*) 'IN HND_RDZMAT'
      ENDIF
C
C     ----- INCLUDE SYMMETRICALLY EQUIVALENT -ZMAT- IF NEEDED -----
C
      CALL HND_ZMTSYM
C
C     ----- SCAN AND ANALYZE IZMAT -----
C
      NZMAT = 0
      DO I = 1,MXIZMT
         IF (IZMAT(I) .NE. 0) THEN
            NZMAT = NZMAT + 1
         ENDIF
      ENDDO
      IF(NZMAT.EQ.0) RETURN
C
      IF (OUT) WRITE (IW,9048)
      ICOUNT = 0
      NCOUNT = 0
  120 CONTINUE
      ICOUNT = ICOUNT + 1
      IF(ICOUNT.GT.NZMAT) GO TO 130
C
C     ----- CHECK FOR A VALID TYPE -----
C
C      1 ..... STRETCH
C      2 ..... BEND
C      3 ..... TORSION (DIHEDRAL ANGLE BETWEEN PLANES HAVING
C                       TWO ATOMS IN COMMON )
C      4 ..... OUT OF PLANE BEND
C      5 ..... LINEAR BEND
C      6 ..... DIHEDRAL ANGLE BETWEEN PLANES HAVING ONE ATOM
C              IN COMMON
C
      IF (IZMAT(ICOUNT) .LE. 0 .OR. IZMAT(ICOUNT) .GE. 8) THEN
         WRITE (IW,9008)
         WRITE (IW,9068) IZMAT(ICOUNT)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF (IZMAT(ICOUNT) .EQ. 1) INC = 2
      IF (IZMAT(ICOUNT) .EQ. 2) INC = 3
      IF (IZMAT(ICOUNT) .EQ. 3) INC = 4
      IF (IZMAT(ICOUNT) .EQ. 4) INC = 4
      IF (IZMAT(ICOUNT) .EQ. 5) INC = 4
      IF (IZMAT(ICOUNT) .EQ. 6) INC = 5
      IF (IZMAT(ICOUNT) .EQ. 6) INC = 5
      IF (IZMAT(ICOUNT) .EQ. 7) INC = 4 ! Bond angle + orientation
C
C     ----- MAKE SURE THE ATOM NUMBERS ARE VALID -----
C
      JDUM = INC
      IF (IZMAT(ICOUNT) .EQ. 5) JDUM = INC - 1
      DO J = 1,JDUM
         IF (IZMAT(ICOUNT+J) .GT. NAT) THEN
            WRITE (IW,9088) IZMAT(ICOUNT+J)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      IF (IZMAT(ICOUNT).EQ.5.AND.IZMAT(ICOUNT+INC).GT.NA) THEN
         WRITE (IW,9088) IZMAT(ICOUNT+INC)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF (IZMAT(ICOUNT).NE.5) NCOUNT=NCOUNT+1
      IF (IZMAT(ICOUNT).EQ.5) NCOUNT=NCOUNT+2
C
      IF (OUT) THEN
         WRITE (IW,9108) IZMAT(ICOUNT),(IZMAT(ICOUNT+J),J = 1,INC)
      ENDIF
      ICOUNT = ICOUNT + INC
      GO TO 120
  130 CONTINUE
C
C     ----- SAVE IZMAT ON DIRECT ACCESS -----
C     -----   -IZMAT-  ON -96-          -----
C     -----   -A-      ON -90-          -----
C
C
C     ----- SET -NVAR- AND -NZVAR- -----
C
      IF(LINEAR.EQ.0) NVAR=3*NAT-6
      IF(LINEAR.NE.0) NVAR=3*NAT-5
      IF( NZMOD.NE.0) NVAR=NZMOD
      NZVAR=NCOUNT
      IF(NZVAR.LT.NVAR) THEN
         WRITE(IW,9028)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(NZVAR.GT.MXZMAT) THEN
         WRITE(IW,9038)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- READ MATRIX FOR SYMMETRIC INTERNAL COORDINATES -----
C                 ( ALSO FOR LINEAR DEPENDENCIES )
C
      CALL HND_ZMTLIN(NZVAR,MXNIJS,IJS,SIJ)
C
C     ----- FINAL ADJUSTMENT OF -NATZMT- IF NEEDED -----
C
      IF(NATZMT.LT.0) THEN
         NATZMT=2
      ENDIF
      RETURN
C
  140 CONTINUE
C
C     ----- FINAL ADJUSTMENT OF -NATZMT- IF NEEDED -----
C
      IF(OUT) THEN
         WRITE(IW,9018)
      ENDIF
      IF(NATZMT.LT.0) THEN
         NATZMT=2
      ENDIF
      RETURN
 9108 FORMAT(5X,6I10)
 9088 FORMAT(' INVALID ATOM NUMBER ',I5,' . STOP')
 9068 FORMAT(' INVALID ZMAT TYPE ',I5,' . STOP')
 9048 FORMAT(/,10X,'INPUT Z-MATRIX ELEMENTS',
     1       /,5X,6X,'TYPE',9X,'I',9X,'J',9X,'K',9X,'L',9X,'M')
 9038 FORMAT(' NUMBER OF INTERNAL COORDINATES IS GREATER THAN ',
     1       ' -MXZMAT- . STOP ' )
 9028 FORMAT(' NUMBER OF INTERNAL COORDINATES IS LESS THAN (3N-6)',
     1       ' OR -NZMOD- .STOP')
 9018 FORMAT(' NO -$ZMAT- OR -$ZMT- NAMELISTS FOUND; PROCEED ...')
 9008 FORMAT(' ERROR IN $ZMAT INPUT . STOP')
      END
      SUBROUTINE HND_ZMTSYM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE HND_ZMTLIN(NZVAR,MXNIJS,IJS,SIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "rtdb.fh"
c
      logical status
      logical geom_zmt_put_zsym
      logical geom_zmt_opn_fil
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (  LBUF=9  )
      LOGICAL SOME
      LOGICAL OUT
      LOGICAL DBUG
      LOGICAL SYM
      CHARACTER*8 ERRMSG
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION IJS(2,1),SIJ(1)
      DIMENSION IBUF(LBUF),BUF(LBUF)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTLIN-'/
      DATA TENM1  /1.0D-01/
      DATA TOL    /1.0D-08/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE. 
      OUT =OUT.OR.DBUG
      OUT =OUT.OR.NPRINT.EQ.1 
      OUT =OUT.OR.IPRZMT.GE.1 
      SOME=       NPRINT.NE.-5
      SOME=SOME.OR.OUT
      IF(OUT) THEN
         WRITE(IW,*) 'IN HND_ZMTLIN'
      ENDIF
C
      NDIM=MAX(3*NAT,NZVAR)
C
C     ----- MEMORY POINTERS -----
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NDIM *NDIM ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_ZMTLIN, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_ZMTLIN, MALLOC OF I20  FAILED',911)
C
      DO IJ=1,NDIM *NDIM 
         DBL_MB(IJ+I10-1)=ZERO
      ENDDO
      DO I=1,NZVAR        
         II=I+NDIM*(I-1)
         DBL_MB(II+I10-1)=ONE
      ENDDO
      IF(DBUG) THEN
         CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
      ENDIF
C
C     ----- CHECK TO SEE IF WE HAVE A -LIN- MATRIX -----
C
      SYM=.FALSE.
      NIJS=0
      DO N=1,MXNIJS
         I=IJS(1,N)
         J=IJS(2,N)
         S=SIJ(  N)
         IF(      I) 100,50,10
   10    IF(NZVAR-I) 100,20,20
   20    IF(      J) 100,50,30
   30    IF(NZVAR-J) 100,40,40
   40    CONTINUE
         NIJS=NIJS+1
         DO K=1,NZVAR
            KJ=K+NDIM*(J-1)
            DBL_MB(KJ+I10-1)=ZERO
         ENDDO
         SYM=.TRUE.
   50    CONTINUE
      ENDDO
C
C     ----- SET AND PRINT -LIN- MATRIX -----
C
      IF(SYM) THEN
         MAXJ=0
         DO N=1,NIJS
            I=IJS(1,N)
            J=IJS(2,N)
            S=SIJ(  N)
            IJ=I+NDIM*(J-1)
            DBL_MB(IJ+I10-1)=S
            IF(J.GT.MAXJ) THEN
               MAXJ=J
            ENDIF
         ENDDO   
      ENDIF
      IF(DBUG) THEN
         WRITE(IW,9993)
         CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
      ENDIF
C
      IF(.NOT.SYM) THEN
C
C     ----- NO SYMMETRY COORDINATES GIVEN -----
C
         DO I=1,NZVAR
            II=I+NDIM*(I-1)
            DBL_MB(II+I10-1)=ONE
         ENDDO
      ELSE
C
C     ----- SYMMETRY COORDINATES ARE GIVEN -----
C
         DO J=1,MAXJ 
C
C     -----           NORMALIZE            -----
C
            DUM=ZERO
            DO K=1,NZVAR
               KJ=K+NDIM*(J-1)
               DUM=DUM+DBL_MB(KJ+I10-1)**2
            ENDDO
            IF(DUM.LT.TOL) THEN
               WRITE(IW,9998) J
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            DUM=ONE/ SQRT(DUM)
            DO K=1,NZVAR
               KJ=K+NDIM*(J-1)
               DBL_MB(KJ+I10-1)=DBL_MB(KJ+I10-1)*DUM
            ENDDO
C
C     -----         ORTHONORMALIZE         -----
C
            IF(J.LT.MAXJ) THEN
               DO I=J+1,NZVAR
                  DUM=ZERO
                  DO K=1,NZVAR
                     KI=K+NDIM*(I-1)
                     KJ=K+NDIM*(J-1)
                     DUM=DUM+DBL_MB(KJ+I10-1)*DBL_MB(KI+I10-1)
                  ENDDO
                  DO K=1,NZVAR
                     KI=K+NDIM*(I-1)
                     KJ=K+NDIM*(J-1)
                     DBL_MB(KI+I10-1)=DBL_MB(KI+I10-1)
     &                               -DBL_MB(KJ+I10-1)*DUM
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         IF(DBUG) THEN
            WRITE(IW,9993)
            CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
         ENDIF
C
C     -----    COMPLETE THE -LIN- MATRIX   -----
C      
         IF(MAXJ.LT.NZVAR) THEN
            I=0 
            DO J=MAXJ+1,NZVAR
   60          I=I+1
               IF(I.GT.NZVAR) THEN
                  WRITE(IW,9992) J
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               DO K=1,NZVAR
                  KJ=K+(J-1)*NDIM 
                  DBL_MB(KJ+I10-1)=ZERO
               ENDDO
               IJ=I+(J-1)*NDIM
               DBL_MB(IJ+I10-1)=ONE 
C
C     -----         ORTHOGONALIZE        -----
C
               DO L=1,J-1
                  DUM=ZERO
                  DO K=1,NZVAR
                     KJ=K+(J-1)*NDIM 
                     KL=K+(L-1)*NDIM 
                     DUM=DUM+DBL_MB(KL+I10-1)*DBL_MB(KJ+I10-1)
                  ENDDO
                  DO K=1,NZVAR
                     KJ=K+(J-1)*NDIM 
                     KL=K+(L-1)*NDIM 
                     DBL_MB(KJ+I10-1)=DBL_MB(KJ+I10-1)
     &                               -DBL_MB(KL+I10-1)*DUM
                  ENDDO
               ENDDO
C
C     -----            NORMALIZE         -----
C
               DUM=ZERO
               DO K=1,NZVAR
                  KJ=K+(J-1)*NDIM 
                  DUM=DUM+DBL_MB(KJ+I10-1)*DBL_MB(KJ+I10-1)
               ENDDO
               IF(DUM.GT.TOL) THEN
                  DUM=ONE/SQRT(DUM)
                  DO K=1,NZVAR
                     KJ=K+(J-1)*NDIM 
                     DBL_MB(KJ+I10-1)=DBL_MB(KJ+I10-1)*DUM            
                  ENDDO
               ELSE
                  GO TO 60
               ENDIF
            ENDDO
         ENDIF
         IF(OUT) THEN
            WRITE(IW,9993)
            IF(DBUG) THEN
               CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
            ENDIF
         ENDIF
C
C     ----- PRINT -LIN- MATRIX -----
C
         IF(SOME) THEN
            WRITE(IW,9997)
            DO J=1,NZVAR
               NBUF=1
               N=0
               DO I=1,NZVAR
                  IJ=I+NDIM*(J-1)
                  DUM=DBL_MB(IJ+I10-1)
                  IF( ABS(DUM).GE.TENM1) THEN
                     N=N+1
                     IBUF(N)=I
                      BUF(N)=DUM
                      IF(N.GE.LBUF) THEN
                   IF(NBUF.EQ.1) WRITE(IW,9995) J,(BUF(M),IBUF(M),M=1,N)
                   IF(NBUF.GT.1) WRITE(IW,9994)   (BUF(M),IBUF(M),M=1,N)
                         NBUF=NBUF+1
                         N=0
                      ENDIF
                  ENDIF
               ENDDO
               IF(N.GT.0) THEN
                  IF(NBUF.EQ.1) WRITE(IW,9995) J,(BUF(M),IBUF(M),M=1,N)
                  IF(NBUF.GT.1) WRITE(IW,9994)   (BUF(M),IBUF(M),M=1,N)
               ENDIF
            ENDDO
C
         ENDIF
C
      ENDIF
      GO TO 1000
C
C     ----- ERRORS -----
C
  100 CONTINUE
      WRITE(IW,9999) I,J,S
      CALL HND_HNDERR(3,ERRMSG)
C
 1000 CONTINUE
C
C     ----- SAVE SYMMETRY MATRIX -----
C           -S-       AT  -92-
C
      STATUS=GEOM_ZMT_OPN_FIL(.FALSE.,.FALSE.,.TRUE.,.TRUE.)
      IF(.NOT.STATUS) THEN
         CALL ERRQUIT('GEOM_ZMT_OPN_FIL?',911)
      ENDIF
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF(.NOT.GEOM_ZMT_PUT_ZSYM(DBL_MB(I10),NZVAR,NDIM))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZIND FAILED',0)
      CALL GA_SYNC
      STATUS=RTDB_PARALLEL(.TRUE.)
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_ZMTLIN, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_ZMTLIN, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' INPUT ERROR IN -ZMTLIN- . I,J,S = ',2I5,F10.5)
 9998 FORMAT(' NORMALIZATION ERROR IN -ZMTLIN- FOR VECTOR J = ',I5)
 9997 FORMAT(/,10X,'SYMMETRIC INTERNAL COORDINATES',/,10X,30(1H-))
 9995 FORMAT(' COORD. = ',I3,9(1X,F5.1,2H (,I3,1H)))
 9994 FORMAT(13X,9(1X,F5.1,2H (,I3,1H)))
 9993 FORMAT(/,10X,'INT. COORD. SYM. MATRIX',/,10X,23(1H-))
 9992 FORMAT(' IN -ZMTLIN- THE ',I4,'-TH INDEPENDENT VECTOR',
     1       ' COULD NOT BE FOUND. STOP')
      END
      SUBROUTINE HND_BANDBI(GEOM,NAT,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "geom.fh"
      INTEGER GEOM
C
C     ----- CALCULATE -B- AND -BI- MATRICES -----
C
      LOGICAL DBUG,OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION C(3,1)
C
C     ----- RETURN IF NO INTERNAL COORDINATES -----
C
      IF(NZMAT.EQ.0) RETURN
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      OUT =OUT.OR.NPRINT.EQ.1
      OUT =OUT.OR.IPRZMT.GE.1
      IF(OUT) THEN
         WRITE(IW,9999) NZMAT,NZVAR,NVAR,NAT
      ENDIF
C
      NCART=3*NAT
      NX1  =NCART
      NQ1  =NZVAR
      ND1  =MAX(NX1,NQ1)
      NZ1  =NZMAT
C
C     ----- MEMORY POINTERS -----
C
C     I10=1
C     I20=I10+NX1*NQ1
C     I30=I20+    NQ1
C     I40=I30+ABS(NZ1)  
C     LAST1=I40-1
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NQ1        ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZ1)   ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I30  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I40 ',I_I40,I40))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I40  FAILED',911)
C
C     I11=I10
C     I21=I11+NX1*NQ1
C     I31=I21+NX1*NQ1
C     I41=I31+ND1*ND1
C     I51=I41+ND1*ND1
C     I61=I51+ND1
C     I71=I61+ND1*ND1
C     I81=I71+ND1*ND1
C     I91=I81+ND1
C     I01=I91+ND1
C     LAST2=I01-1
C     LAST=MAX0(LAST1,LAST2)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I11 ',I_I11,I11))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I11  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I21 ',I_I21,I21))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I21  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I31 ',I_I31,I31))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I31  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I41 ',I_I41,I41))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I41  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1        ,'MEM I51 ',I_I51,I51))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I51  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I61 ',I_I61,I61))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I61  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I71 ',I_I71,I71))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I71  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1        ,'MEM I81 ',I_I81,I81))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I81  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1        ,'MEM I91 ',I_I91,I91))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I91  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I01 ',I_I01,I01))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I01  FAILED',911)
C
C     ----- CALCULATE -B- MATRIX -----
C
      CALL HND_BMAT(GEOM,NZVAR,NCART,
     1              DBL_MB(I10),DBL_MB(I20),DBL_MB(I30),
     2              NZMAT,C)
      IF(OUT) THEN
         WRITE(IW,9998)
      ENDIF
C
C     ----- CALCULATE THE B INVERSE MATRIX -----
C
      CALL HND_BINVR(NVAR,NZVAR,NCART,DBL_MB(I11),DBL_MB(I21),
     1        DBL_MB(I31),DBL_MB(I41),DBL_MB(I51),DBL_MB(I61),
     2        DBL_MB(I71),DBL_MB(I81),DBL_MB(I91),ND1)
      IF(OUT) THEN
         WRITE(IW,9997)
      ENDIF
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I01))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I01  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I91))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I91  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I81))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I81  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I71))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I71  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I61))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I61  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I51))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I51  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I41))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I41  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I31))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I31  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I21))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I21  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I11))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I11  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I40))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I40  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -BANDBI- NZMAT,NZVAR,NVAR,NAT = ',4I5)
 9998 FORMAT(' IN -BANDBI- AFTER -BMAT- , BEFORE -BINVR- ')
 9997 FORMAT(' IN -BANDBI- AFTER -BINVR- ')
      END
      SUBROUTINE HND_BMAT(GEOM,NZVAR,NCART,B,ZMAT,IZMAT,NZMAT,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
      logical status
      integer geom
C
C     ----- CONSTRUCT THE B MATRIX -----
C
      PARAMETER (    NA=10 )
      LOGICAL     DBUG
      LOGICAL     OUT
      LOGICAL     GEOM_ZMT_GET_IZMAT
      LOGICAL     GEOM_ZMT_GET_NIZMAT
      LOGICAL     GEOM_ZMT_PUT_ZMAT 
      LOGICAL     GEOM_ZMT_PUT_BMAT 
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION B(NCART,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION A(3,NA)
      DIMENSION C(3,*)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- BMAT -'/
      DATA PI2    /6.28318530717958D+00/
      DATA DEGREE /360.0D+00/
      DATA BOHR   /0.52917715D+00/
      DATA ZERO   /0.0D+00/
      DATA FOUR   /4.0D+00/
C
      DBUG  =.FALSE.
      OUT   =.FALSE. 
      OUT   =OUT.OR.DBUG
      OUT   =OUT.OR.NPRINT.EQ.1
      OUT   =OUT.OR.IPRZMT.GE.1
C
      PIO2  =PI2/FOUR
      RTOD  =PI2/DEGREE
C
C     ----- READ -IZMAT- -----
C          -A-     AT -90-
C          -IZMAT- AT -96-
C
      NDUM=3*NA
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF( GA_NODEID().EQ.0) THEN
         IF(.NOT.GEOM_ZMT_GET_NIZMAT(GEOM,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_NIZMAT FAILED',0)
         IF(.NOT.GEOM_ZMT_GET_IZMAT(GEOM,IZMAT,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_IZMAT FAILED',0)
         CALL GA_BRDCST(20,NZMAT,MITOB(1     ),0)
         CALL GA_BRDCST(21,IZMAT,MITOB(NZMAT ),0)
      ELSE 
         CALL GA_BRDCST(20,NZMAT,MITOB(1     ),0)
         CALL GA_BRDCST(21,IZMAT,MITOB(NZMAT ),0)
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
      IF(OUT) THEN
         WRITE(IW,9997) 
         WRITE(IW,9996) (IZMAT(I),I=1,NZMAT)
      ENDIF
C
C     ----- ZERO OUT THE B MATRIX -----
C
      DO J = 1,NZVAR
         DO I = 1,NCART
            B(I,J) = ZERO
         ENDDO
      ENDDO
C
      IADD=1
      I   =0
      DO IZVAR=1,NZVAR
      I    =I+1            
      ITYPE=IZMAT(IADD)
      IF(OUT) THEN
         WRITE(IW,9999) IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART
      ENDIF
      IF(ITYPE.EQ.1) THEN
C
C     ----- BOND STRETCH -----
C
         CALL HND_BSTR(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             C,B,NCART,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 3
      ELSEIF(ITYPE.EQ.2 .or. itype.eq.7) THEN
C
C     ----- ANGLE BEND -----
C
         CALL HND_BEND(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 4
         if (itype .eq. 7) iadd = iadd + 1 ! Additional orientation
      ELSEIF(ITYPE.EQ.3) THEN
C
C     ----- TORSION -----
C
         CALL HND_TORS(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.4) THEN
C
C     ----- OUT OF PLANE ANGLE BEND -----
C
         CALL HND_OPLA(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,PIO2,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.5) THEN
C
C     ----- LINEAR ANGLE BEND -----
C
         CALL HND_LIBE(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,A,NA,BOHR)
         ZMAT(I  ) = EQVAL
         ZMAT(I+1) = EQVAL
         I=I+1
         IADD=IADD+5
      ELSEIF(ITYPE.EQ.6) THEN
C
C     ----- DIHEDRAL ANGLE BETWEEN TWO PLANES SHARING ONE ATOM -----
C
         CALL HND_DIHPLA(EQVAL,I,IZMAT(IADD+1),C,B,NCART,RTOD)
         ZMAT(I)=EQVAL
         IADD=IADD+6
C
      ELSE
         WRITE(IW,9995) IZVAR,IADD,ITYPE
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF    
C
      ENDDO
C
      IF(OUT) THEN
         CALL HND_PRTZMT(NZVAR,IZMAT,ZMAT,DUM,0)
         WRITE(IW,9998)
         IF(DBUG) THEN
            CALL HND_PRSQ(B,NZVAR,NCART,NCART)
         ENDIF
      ENDIF
C
C     ----- SAVE -B- AND -Z- ON -IDAF- -----
C          -A-                                 AT -90-
C          -B-  (pristine)                     AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C                                              AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C
      NDUM=NCART*NZVAR
C
C     ----- STORE FOR -NWCHEM- -----
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF(.NOT.GEOM_ZMT_PUT_ZMAT(GEOM,ZMAT,NZVAR))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZMAT FAILED',0)
      IF(.NOT.GEOM_ZMT_PUT_BMAT(B   ,NCART,NZVAR,NCART))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_BMAT FAILED',0)
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
      RETURN
 9999 FORMAT(' IN -BMAT- IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART = ',
     1       6I5)
 9998 FORMAT(/,' -B- MATRIX (ROWS = CARTESIAN COORDINATE.',
     1         ' COLUMNS = INTERNAL COORDINATES.)')
 9997 FORMAT(' IN -BMAT- IZMAT = ')
 9996 FORMAT(12I5)
 9995 FORMAT(' ERROR IN -BMAT- , IZVAR,IADD,ITYPE = ',3I5)
      END
      SUBROUTINE HND_BINVR(NVAR,NZVAR,NCART,B,BINV,
     1                 G,V,E,H,S,T,IA,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- CONSTRUCT THE INVERSE OF THE -B- MATRIX -----
C
C     B            NCART *NZVAR
C     BINV         NCART *NZVAR
C     G            NZVAR *NZVAR
C     V            NZVAR *NZVAR
C     E            NZVAR
C     H       MAX( NZVAR *NZVAR , NCART*NCART )
C     S            NZVAR *NZVAR
C     T            NZVAR
C     IA      MAX( NZVAR        , NCART       )
C
      LOGICAL     OUT
      LOGICAL     DBUG
      LOGICAL     NATURL
      LOGICAL     SYMTRC
      LOGICAL     INDEP
      LOGICAL     GEOM_ZMT_PUT_BINV 
      LOGICAL     GEOM_ZMT_GET_BMAT
      LOGICAL     GEOM_ZMT_PUT_ZIND
      LOGICAL     GEOM_ZMT_PUT_ZSYM
      LOGICAL     GEOM_ZMT_GET_ZSYM
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      DIMENSION B(NCART,*),BINV(NCART,*)
      DIMENSION G(*),V(NDIM,*),E(*),IA(*)
      DIMENSION      S(NDIM,*),T(*)
      DIMENSION      H(NDIM,*)     
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- BINVR-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TM08   /1.0D-08/
C
      DBUG=.FALSE.
      DBUG=DBUG.OR.IPRZMT.GE.2
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      OUT =OUT.OR.NPRINT.EQ.1
      OUT =OUT.OR.IPRZMT.EQ.1
C
      DO I=1,NDIM           
         IA(I)=(I*(I-1))/2
      ENDDO
C
C     --- NATURL .. FINAL BASIS = NATURAL   INTERNAL        COORD. ---
C     --- SYMTRC .. FINAL BASIS = SYMMETRIC INTERNAL        COORD. ---
C     --- INDEP  .. FINAL BASIS = (NON-REDUNDANT+REDUNDANT) COORD. ---
C
      NATURL=NATZMT.EQ.0
      SYMTRC=NATZMT.EQ.1
      INDEP =NATZMT.EQ.2              
C
      IF(OUT) THEN
         WRITE(IW,9984) NZVAR,NVAR,NCART,NDIM
         WRITE(IW,9999) NATURL,SYMTRC,INDEP,NATZMT
         IF(NATURL) WRITE(IW,9998)
         IF(SYMTRC) WRITE(IW,9997)
         IF(INDEP ) WRITE(IW,9996)
      ENDIF
C
C     ----- READ -B- MATRIX FOR NATURAL INTERNAL   COORDINATES -----
C     -----  AND -S- MATRIX FOR                                -----
C     -----  TRANSFORMATION  TO SYMMETRIC INTERNAL COORDINATES -----
C
      IF(.NOT.GEOM_ZMT_GET_BMAT(B   ,NCART,NZVAR,NCART))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_BMAT FAILED',0)
      IF(OUT) THEN
         WRITE(IW,9995)
         IF(DBUG) THEN
            CALL HND_PRSQ(B,NZVAR,NCART,NCART)
         ENDIF
      ENDIF
C
C     -----   TRANSFORM -B- TO BASIS OF    -----
C     ----- SYMMETRIC INTERNAL COORDINATES -----
C
      IF(SYMTRC.OR.INDEP) THEN
         IF(.NOT.GEOM_ZMT_GET_ZSYM(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_ZSYM FAILED',0)
         DO I=1,NCART
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NZVAR
                  DUM=DUM+B(I,K)*S(K,J)
               ENDDO
               T(J)=DUM
            ENDDO
            DO J=1,NZVAR
               B(I,J)=T(J)
            ENDDO
         ENDDO
         IF(OUT) THEN
            WRITE(IW,9994)
            WRITE(IW,9993)
            IF(DBUG) THEN
               CALL HND_PRSQ(B,NZVAR,NCART,NCART)
            ENDIF
         ENDIF
      ENDIF
C
C          -A-                                 AT -90-
C          -B-  (maybe ... B tilde = B * S )   AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C                                              AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C         NOTE THAT AT THIS POINT WE MAY BE WORKING WITH 
C           B(tilde) WHICH IS THE -S- TRANSFORMED OF B
C                      B(tilde) = B * S
C
C     -----      CONSTRUCT THE -G- MATRIX       -----
C                      G = B(T) * B
C             ACTUALLY WORK WITH ( - B(T) * B )
C           ( EASY WAY TO RE-ORDER THE VECTORS )
C
      DO I=1,NZVAR
         DO J=1,I
            IJ=IA(I)+J
            DUM=ZERO
            DO K=1,NCART
               DUM=DUM+B(K,I)*B(K,J)
            ENDDO
            G(IJ)=-DUM
         ENDDO
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,9992)
         CALL HND_PRTR(G,NZVAR)
      ENDIF
C
C     ----- GET THE ( NON-REDUNDANT + REDUNDANT ) COORDINATES -----
C
      CALL HND_DIAGIV(G,V,E,IA,NZVAR,NZVAR,NDIM) 
C
C     -----  WE WORKED WITH ( - B(T) * B ) TO GET THE CORRECT ORDER  -----
C           FOR THE VECTORS. WE NEED TO CHANGE THE SIGN OF THE E's
C
      DO I=1,NZVAR
         E(I)=-E(I)
      ENDDO
C
C     -----    IF WE ARE TO WORK IN THE SPACE OF     -----
C     -----  (NON-REDUNDANT+REDUNDANT) COORDINATES   -----
C     -----   WE MUST MAINTAIN A PHASE CONSISTENCY   -----
C     -----   ... IF NOT ... WE ARE IN TROUBLES ...  -----
C     -----  SOME (PATHOLOGICAL) CASES WILL NOT WORK -----
C
      DO J=1,NZVAR 
         IMAX=0
         VMAX=ZERO
         DO I=1,NZVAR
            IF(ABS(V(I,J)).GT.VMAX) THEN
               VMAX=ABS(V(I,J))
               IMAX=I
            ENDIF
         ENDDO
         IF(V(IMAX,J).LT.ZERO) THEN
            DO I=1,NZVAR
               V(I,J)=-V(I,J)
            ENDDO
         ENDIF
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,9991)
         WRITE(IW,9990)
         CALL HND_PREV(V,E,NZVAR,NZVAR,NDIM) 
      ENDIF
C
C     ----- CHECK FOR SINGULARITY -----
C
      NINDEP=0
      NREDUN=0
      DO I=1,NZVAR
         IF( ABS(E(I)).LE.TM08) THEN
            E(I)  =ZERO
            NREDUN=NREDUN+1
         ELSE
            E(I)  =ONE/E(I)
            NINDEP=NINDEP+1
         ENDIF
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9989) NINDEP,NREDUN,NZVAR,NVAR
      ENDIF
C
C     IF NREDUN.EQ.0 , WE CAN WORK WITH NATURAL COORDINATES -----
C
      IF(NREDUN.EQ.0) THEN
         NATZMT=0
         NATURL=NATZMT.EQ.0
         SYMTRC=NATZMT.EQ.1
         INDEP =NATZMT.EQ.2              
         IF(OUT) THEN
            WRITE(IW,9979)
         ENDIF
      ENDIF
C
      IF(NATURL.OR.SYMTRC) THEN
         DO J=1,NZVAR
            DO I=1,NZVAR
               S(I,J)=ZERO
            ENDDO
            S(J,J)=ONE
         ENDDO
         IF(.NOT.GEOM_ZMT_PUT_ZIND(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZIND FAILED',0)
      ENDIF
C
C     ----- IF WORKING BASIS IS (NON-REDUNDANT+REDUNDANT) -----
C     -----         TRANSFORM -B- TO THAT BASIS           -----
C
      IF(INDEP) THEN
         IF(.NOT.GEOM_ZMT_GET_ZSYM(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_ZSYM FAILED',0)
         DO I=1,NZVAR
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NZVAR
                  DUM=DUM+S(I,K)*V(K,J)
               ENDDO
               T(J)=DUM
            ENDDO
            DO J=1,NZVAR
               S(I,J)=T(J)
            ENDDO
         ENDDO
         IF(.NOT.GEOM_ZMT_PUT_ZSYM(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZSYM FAILED',0)
         IF(OUT) THEN
            WRITE(IW,9980)
            IF(DBUG) THEN
               CALL HND_PRSQ(S,NZVAR,NZVAR,NDIM)
            ENDIF
         ENDIF
C
         DO J=NINDEP+1,NZVAR
            DO I=1,NZVAR
               V(I,J)=ZERO
            ENDDO
         ENDDO
C
         DO I=1,NCART
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NZVAR
                  DUM=DUM+B(I,K)*V(K,J)
               ENDDO
               T(J)=DUM
            ENDDO
            DO J=1,NZVAR
               B(I,J)=T(J)
            ENDDO
         ENDDO
         IF(OUT) THEN
            WRITE(IW,9985)
            IF(DBUG) THEN
               CALL HND_PRSQ(B,NZVAR,NCART,NCART)
            ENDIF
         ENDIF
C
         DO J=1,NINDEP
            DO I=1,NZVAR
               V(I,J)=ZERO
            ENDDO
            V(J,J)=ONE 
         ENDDO
      ENDIF
C
C     ----- CALCULATE  K * K(T) PROJECTOR -----
C
      DO I=1,NZVAR
         DO J=1,NZVAR
            DUM=ZERO
            DO K=1,NINDEP
               DUM=DUM+V(I,K)*V(J,K)
            ENDDO
            H(I,J)=DUM
         ENDDO
      ENDDO
      NDUM=NDIM*NDIM  
      IF(OUT) THEN
         IF(DBUG) THEN
            WRITE(IW,9981)
            CALL HND_PRSQ(H,NZVAR,NZVAR,NDIM) 
         ENDIF
      ENDIF
C
C     ----- SET  ... NVAR = NZVAR ... IF WE WORK -----
C           WITH THE SYMMETRIC INTERNAL SPACE OR
C           WITH THE NATURAL   INTERNAL SPACE
C
      IF(NATURL.OR.SYMTRC) THEN
         NVAR=NZVAR
      ELSE
         IF(OUT) THEN
            IF(NINDEP.NE.NVAR) THEN
               WRITE(IW,9988) NZVAR,NINDEP,NREDUN,NVAR
            ENDIF
            IF(NREDUN.NE.NZVAR-NVAR) THEN
               WRITE(IW,9988) NZVAR,NINDEP,NREDUN,NVAR
            ENDIF
         ENDIF
      ENDIF
C
C     ----- GENERALIZED INVERSE OF -G- IN THE BASIS OF -----
C     ----- NATURAL/SYMMETRIC INTERNAL COORDINATES  OR -----
C     ----- (NON-REDUNDANT+REDUNDANT)  COORDINATES     -----
C
      DO I=1,NZVAR
         DO J=1,NZVAR
            DUM=ZERO
            DO K=1,NZVAR
               DUM=DUM+V(I,K)*E(K)*V(J,K)
            ENDDO
            H(I,J)=DUM
         ENDDO
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,9987)
         CALL HND_PRSQ(H,NZVAR,NZVAR,NDIM) 
      ENDIF
C
C     -----   BINV(T) = B * ( ( B(T)*B ) **(-1) )(T)  -----
C
      DO I=1,NCART
         DO J=1,NZVAR
            DUM=ZERO
            DO K=1,NZVAR
               DUM=DUM+B(I,K)*H(K,J)
            ENDDO
            BINV(I,J)=DUM
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9982)
         WRITE(IW,9986)
         IF(DBUG) THEN
            CALL HND_PRSQ(BINV,NZVAR,NCART,NCART)
         ENDIF
      ENDIF
      NDUM=NCART*NVAR
C
C     ----- STORE FOR -NWCHEM- -----
C
      IF(.NOT.GEOM_ZMT_PUT_BINV(BINV,NCART,NZVAR,NCART))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_BINV FAILED',0)
C
C     -----        DEBUG CHECK OF INVERSE             -----
C     -----         H = Bi * B = K * K(T)             -----
C     -----         WHERE K IS THE MATRIX             -----
C     -----    OF (NON-REDUNDANT) COORDINATES         -----
C     -----          H IS NOT EQUAL TO   1            -----
C     ----- UNLESS THERE ARE NO REDUNDANT COORDINATES -----
C
      IF(DBUG) THEN
         NDUM= NCART*NZVAR
         DO I=1,NZVAR
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NCART
                  DUM=DUM+BINV(K,I)*B(K,J)
               ENDDO
               H(I,J)=DUM
            ENDDO
         ENDDO
         WRITE(IW,9983)
         CALL HND_PRSQ(H,NZVAR,NZVAR,NDIM) 
      ENDIF
C
C     ----- NEXT TIME AROUND WE WILL WORK WITH -----
C           WHAT ARE NOW SYM.INTERN.COORD.
C
      NATZMT=1
C
C          -A-                                 AT -90-
C          -B-  (B tilde)                      AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C          -P = K * K(T) PROJECTOR             AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C     -NZMAT-   = N.OF. -IZMAT- ELEMENTS GIVEN IN -$ZMAT-
C     -NZVAR-   = N.OF. INTERNAL COORDINATES DEFINED IN -IZMAT-
C     -NVAR-    = (3N-6), OR (3N-5), OR -NZMOD- , .LE. NZVAR
C
      RETURN
 9999 FORMAT(/,' IN -BINVR- -NATURL-, -SYMTRC- , -INDEP- , -NATZMT- = ',
     1       3L5,I4)
 9998 FORMAT(/,' FINAL BASIS WILL BE =    NATURAL   INTERNAL    ')
 9997 FORMAT(/,' FINAL BASIS WILL BE =    SYMMETRIC INTERNAL    ')
 9996 FORMAT(/,' FINAL BASIS WILL BE = (NON-REDUNDANT+REDUNDANT)')
 9995 FORMAT(/,' -B - IN NATURAL   INTERNAL        COORDINATES ')
 9994 FORMAT(/,' -B - IN NATURAL   INTERNAL        COORDINATES ',
     1                                             'TRANSFORMED',
     2          '     TO SYMMETRIC INTERNAL        COORDINATES ')
 9993 FORMAT(/,' -B - IN SYMMETRIC INTERNAL        COORDINATES ')
 9992 FORMAT(/,' -G - IN SYMMETRIC INTERNAL        COORDINATES ')
 9991 FORMAT(/,'         (NON-REDUNDANT+REDUNDANT) COORDINATES = ',
     1            'EIGENVECTORS OF -G- ')
 9990 FORMAT(/,' -V - OF (NON-REDUNDANT+REDUNDANT) COORDINATES ')
 9989 FORMAT(/,' NUM. OF LIN.INDEP. COORD. (NINDEP) = ',I5,
     1       /,' NUM. OF REDUNDANT  COORD. (NREDUN) = ',I5,
     2       /,'                            -NZVAR- = ',i5,
     3       /,'                            -NVAR-  = ',I5)
 9988 FORMAT(/,' -NINDEP- SHOULD BE EQUAL TO -NVAR-, BUT IT IS NOT.',
     1         ' STOP. NZVAR,NINDEP,NREDUN,NVAR = ',4I5)
 9987 FORMAT(/,' - G**(-1) - ') 
 9986 FORMAT(/,' -Bi- IN ..... CHOSEN .....        COORDINATES ')
 9985 FORMAT(/,' -B - IN (NON-REDUNDANT+REDUNDANT) COORDINATES ')
 9984 FORMAT(/,' IN -BINVR- NZVAR,NVAR,NCART,NDIM = ',4I5)
 9983 FORMAT(/,' CHECK OF ( Bi * B = K * K(T) )')
 9982 FORMAT(/,' -Bi(T) = B * ( G**(-1) )(T)')
 9981 FORMAT(/,' K * K(T) PROJECTOR ')         
 9980 FORMAT(/,' -ZIND- FROM NATURAL TO (INDEP+DEP) COORDINATES')
 9979 FORMAT(/,' COORDINATE BASIS CHANGED BACK TO -NATURAL-')
      END
      SUBROUTINE HND_BANDZ(GEOM,NAT,C,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "geom.fh"
C
C     ----- CALCULATE -B- AND -BI- MATRICES -----
C
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION C(3,*),Z(*)
C
C     ----- RETURN IF NO INTERNAL COORDINATES -----
C
      IF(NZMAT.EQ.0) RETURN
C
      OUT=.FALSE.
      OUT=OUT.OR.NPRINT.EQ.1
      OUT=OUT.OR.IPRZMT.GE.1
      IF(OUT) THEN
         WRITE(IW,9999) NZMAT,NZVAR,NVAR,NAT
      ENDIF
C
      NCART=3*NAT
      NX1  =NCART
      NQ1  =NZVAR
      ND1  =MAX(NX1,NQ1)
      NZ1  =NZMAT
C
C     ----- MEMORY POINTERS -----
C
C     I10=1
C     I20=I10+NX1*NQ1
C     I30=I20+    NQ1
C     I40=I30+ABS(NZ1)  
C     LAST1=I40-1
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NQ1        ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZ1)   ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I30  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I40 ',I_I40,I40))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I40  FAILED',911)
C
C     ----- CALCULATE -B- MATRIX -----
C
      CALL HND_BZ(GEOM,NZVAR,NCART,
     1            DBL_MB(I10),DBL_MB(I20),DBL_MB(I30),
     2            NZMAT,C)
      DO I=1,NZVAR
         Z(I)=DBL_MB(I+I20-1)
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,9998)
      ENDIF
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I40))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I40  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -BANDZ- NZMAT,NZVAR,NVAR,NAT = ',4I5)
 9998 FORMAT(' IN -BANDZ- AFTER -BZ- ')
      END
      SUBROUTINE HND_BZ(GEOM,NZVAR,NCART,B,ZMAT,IZMAT,NZMAT,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer geom
C
C     ----- CONSTRUCT THE B MATRIX -----
C
      PARAMETER (    NA=10 )
      LOGICAL     OUT
      LOGICAL     GEOM_ZMT_GET_IZMAT
      LOGICAL     GEOM_ZMT_GET_NIZMAT
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION B(NCART,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION A(3,NA)
      DIMENSION C(3,*)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- BMAT -'/
      DATA PI2    /6.28318530717958D+00/
      DATA DEGREE /360.0D+00/
      DATA BOHR   /0.52917715D+00/
      DATA ZERO   /0.0D+00/
      DATA FOUR   /4.0D+00/
C
      OUT   =.FALSE. 
      OUT   =OUT.OR.NPRINT.EQ.1
      OUT   =OUT.OR.IPRZMT.GE.1
C
      PIO2  =PI2/FOUR
      RTOD  =PI2/DEGREE
C
C     ----- READ -IZMAT- -----
C          -A-     AT -90-
C          -IZMAT- AT -96-
C
      NDUM=3*NA
      IF(.NOT.GEOM_ZMT_GET_NIZMAT(GEOM,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_NIZMAT FAILED',0)
      IF(.NOT.GEOM_ZMT_GET_IZMAT(GEOM,IZMAT,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_IZMAT FAILED',0)
      IF(OUT) THEN
         WRITE(IW,9997) 
         WRITE(IW,9996) (IZMAT(I),I=1,NZMAT)
      ENDIF
C
C     ----- ZERO OUT THE B MATRIX -----
C
      DO J = 1,NZVAR
         DO I = 1,NCART
            B(I,J) = ZERO
         ENDDO
      ENDDO
C
      IADD=1
      I   =0
      DO IZVAR=1,NZVAR
      I    =I+1            
      ITYPE=IZMAT(IADD)
      IF(OUT) THEN
         WRITE(IW,9999) IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART
      ENDIF
      IF(ITYPE.EQ.1) THEN
C
C     ----- BOND STRETCH -----
C
         CALL HND_BSTR(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             C,B,NCART,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 3
      ELSEIF(ITYPE.EQ.2 .or. itype.eq.7) THEN
C
C     ----- ANGLE BEND -----
C
         CALL HND_BEND(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 4
         if (itype .eq. 7) iadd = iadd + 1 ! Additional orientation
      ELSEIF(ITYPE.EQ.3) THEN
C
C     ----- TORSION -----
C
         CALL HND_TORS(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.4) THEN
C
C     ----- OUT OF PLANE ANGLE BEND -----
C
         CALL HND_OPLA(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,PIO2,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.5) THEN
C
C     ----- LINEAR ANGLE BEND -----
C
         CALL HND_LIBE(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,A,NA,BOHR)
         ZMAT(I  ) = EQVAL
         ZMAT(I+1) = EQVAL
         I=I+1
         IADD=IADD+5
      ELSEIF(ITYPE.EQ.6) THEN
C
C     ----- DIHEDRAL ANGLE BETWEEN TWO PLANES SHARING ONE ATOM -----
C
         CALL HND_DIHPLA(EQVAL,I,IZMAT(IADD+1),C,B,NCART,RTOD)
         ZMAT(I)=EQVAL
         IADD=IADD+6
C
      ELSE
         WRITE(IW,9995) IZVAR,IADD,ITYPE
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF    
C
      ENDDO
C
      IF(OUT) THEN
         CALL HND_PRTZMT(NZVAR,IZMAT,ZMAT,DUM,0)
         WRITE(IW,9998)
         CALL HND_PRSQ(B,NZVAR,NCART,NCART)
      ENDIF
C
C     ----- SAVE -B- AND -Z- ON -IDAF- -----
C          -A-                                 AT -90-
C          -B-  (pristine)                     AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C                                              AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C
      NDUM=NCART*NZVAR
C
      RETURN
 9999 FORMAT(' IN -BZ- IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART = ',
     1       6I5)
 9998 FORMAT(/,' -B- MATRIX (ROWS = CARTESIAN COORDINATE.',
     1         ' COLUMNS = INTERNAL COORDINATES.)')
 9997 FORMAT(' IN -BZ- IZMAT = ')
 9996 FORMAT(12I5)
 9995 FORMAT(' ERROR IN -BZ- , IZVAR,IADD,ITYPE = ',3I5)
      END
      SUBROUTINE HND_ZMTPRT(GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
C
C     ----- PRINT THE CURRENT Z MATRIX -----
C
      LOGICAL STATUS
      INTEGER GEOM
      LOGICAL GEOM_ZMT_GET_IZMAT
      LOGICAL GEOM_ZMT_GET_NIZMAT
      LOGICAL GEOM_ZMT_GET_ZMAT
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
c
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999) NZMAT,NZVAR,NVAR
      ENDIF
C
      IF(NZMAT.EQ.0) RETURN
C
C     ----- MEMORY POINTERS -----
C
C     I10 = 1
C     I20 = I10 + NZVAR
C     I30 = I20 + ABS(NZMAT)
C     LAST = I30 - 1
      IF(.NOT.MA_PUSH_GET(MT_DBL,NZVAR      ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_ZMTPRT, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZMAT) ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_ZMTPRT, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_ZMTPRT, MALLOC OF I30  FAILED',911)
C
C     ----- GET IZMAT AND ZMAT FROM DIRECT ACCESS -----
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF( GA_NODEID().EQ.0) THEN 
         IF(.NOT.GEOM_ZMT_GET_NIZMAT(GEOM,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_NIZMAT FAILED',0)
         IF(.NOT.GEOM_ZMT_GET_IZMAT(GEOM,DBL_MB(I20),NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_IZMAT FAILED',0)
         IF(.NOT.GEOM_ZMT_GET_ZMAT(GEOM,DBL_MB(I10),NZVAR))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_ZMAT FAILED',0)
         CALL GA_BRDCST(22,DBL_MB(I20),MDTOB(NZMAT),0)
         CALL GA_BRDCST(23,DBL_MB(I10),MDTOB(NZVAR),0)
      ELSE
         CALL GA_BRDCST(22,DBL_MB(I20),MDTOB(NZMAT),0)
         CALL GA_BRDCST(23,DBL_MB(I10),MDTOB(NZVAR),0)
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
C     ----- PRINT THE Z MATRIX -----
C
      IF( GA_NODEID().EQ.0.OR.DBUG) THEN 
         CALL HND_PRTZMT(NZVAR,DBL_MB(I20),DBL_MB(I10),DBL_MB(I30),0)
      ENDIF
      CALL GA_SYNC()
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_ZMTPRT, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_ZMTPRT, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_ZMTPRT, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -ZMTPRT- NZMAT, NZVAR, NVAR = ',3i5)
      END
      SUBROUTINE HND_ZMTPRG(GRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
C
C     ----- PRINT THE CURRENT Z MATRIX AND ITs GRADIENT -----
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION GRAD(1)
      DATA ZERO /0.0D+00/
C
      IF(NZMAT.EQ.0) RETURN
C
C     ----- MEMORY POINTERS -----
C
C     I10 = 1
C     I20 = I10 + NZVAR
C     I30 = I20 + ABS(NZMAT)
C     LAST = I30 - 1
      IF(.NOT.MA_PUSH_GET(MT_DBL,NZVAR      ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_ZMTPRG, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZMAT) ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_ZMTPRG, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_ZMTPRG, MALLOC OF I30  FAILED',911)
C
C     ----- GET IZMAT AND ZMAT FROM DIRECT ACCESS -----
C
C
C     ----- PRINT THE Z MATRIX -----
C
      IF(NZVAR.GT.NVAR) THEN
         DO I=NVAR+1,NZVAR
            GRAD(I)=ZERO
         ENDDO
      ENDIF
      CALL HND_PRTZMT(NZVAR,DBL_MB(I20),DBL_MB(I10),GRAD,1)
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_ZMTPRG, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_ZMTPRG, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_ZMTPRG, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
      END
      SUBROUTINE HND_PRTZMT(NZVAR,IZMAT,ZMAT,ZGRAD,IPRGRD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT THE Z-MATRIX AND ITs GRADIENT -----
C
C     ----- IPRGRD = 0  PRINT JUST THE Z MATRIX
C                  = 1  PRINT THE Z MATRIX AND ITs GRADIENT
C
      CHARACTER*8 RTYPE
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      DIMENSION IZMAT(1),ZMAT(1),ZGRAD(1)
      DIMENSION RTYPE(7)
      DATA RTYPE/'STRETCH ','BEND    ','TORSION ','PLA.BEND',
     1           'LIN.BEND','DIH.PLA ','BEND(2) '/
      IF (IPRGRD.EQ.0) WRITE (IW,9008)
      IF (IPRGRD.EQ.1) WRITE (IW,9028)
      IADD = 1
      I = 0
   80 I = I + 1
      IF(I.GT.NZVAR) GO TO 220
      ITYPE = IZMAT(IADD)
      GO TO (100,120,140,160,180,190,121),ITYPE
C
C     ----- BOND STRETCH -----
C
  100 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9048) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9048) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),ZMAT(I),ZGRAD(I)
      IADD = IADD + 3
      GO TO 200
C
C     ----- ANGLE BEND -----
C
  120 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9068) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9068) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),ZMAT(I),ZGRAD(I)
      IADD = IADD + 4
      GO TO 200
C
C     ----- ANGLE BEND ----- + orientation (type 7)
C
 121  CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9069) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9069) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I),ZGRAD(I)
      IADD = IADD + 5
      GO TO 200
C
C     ----- TORSION -----
C
  140 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I),ZGRAD(I)
      IADD = IADD + 5
      GO TO 200
C
C     ----- OUT OF PLANE ANGLE BEND -----
C
  160 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I),ZGRAD(I)
      IADD = IADD + 5
      GO TO 200
C
C     ----- LINEAR ANGLE BEND -----
C
  180 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I),ZGRAD(I)
      I = I + 1
      IADD = IADD + 5
      GO TO 200
C
C     ----- DIHEDRAL ANGLE BETWEEN PLANES WITH ONE COMMON ATOM -----
C
  190 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9108) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),IZMAT(IADD+5),
     2    ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9108) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),IZMAT(IADD+5),
     2    ZMAT(I),ZGRAD(I)
      IADD = IADD + 6
C
  200 CONTINUE
      GO TO 80
  220 CONTINUE
C
      IF(IPRGRD.EQ.1) THEN
         IF(NATZMT.EQ.0) WRITE(IW,9128)
         IF(NATZMT.EQ.1) WRITE(IW,9148)
         IF(NATZMT.EQ.2) WRITE(IW,9168)
      ENDIF
      RETURN
 9008 FORMAT(/,22X,20(1H-),/,22X,20Hinternal coordinates,/,22X,20(1H-),
     1       /,2X,3Hno.,2X,8H  type  ,
     2         2X,6H -i-  ,6H -j-  ,6H -k-  ,6H -l-  ,6H -m-  ,
     3         5X,7Hvalue  ,
     3       /,2X,57(1H-))
 9028 FORMAT(/,22X,20(1H-),/,22X,20Hinternal coordinates,/,22X,20(1H-),
     1       /,2X,3Hno.,2X,8H  type  ,
     2         2X,6H -i-  ,6H -j-  ,6H -k-  ,6H -l-  ,6H -m-  ,
     3         5X,7Hvalue  ,4X,8Hgradient,
     4       /,2X,69(1H-))
 9048 FORMAT(1X,I4,2X,A8,2X,2(I3,3X),18X,F12.5,F12.7)
 9068 FORMAT(1X,I4,2X,A8,2X,3(I3,3X),12X,F12.5,F12.7)
 9069 FORMAT(1X,I4,2X,A8,2X,3(I3,3X),6X,i3,3x,F12.5,F12.7)
 9088 FORMAT(1X,I4,2X,A8,2X,4(I3,3X), 6X,F12.5,F12.7)
 9108 FORMAT(1X,I4,2X,A8,2X,5(I3,3X),    F12.5,F12.7)
 9128 FORMAT(7X,'                  gradients are for natural ',
     1          'internal coordinates')
 9148 FORMAT(7X,'                gradients are for symmetric ',                 
     1          'internal coordinates')
 9168 FORMAT(7X,'gradients are for (non-redundant+redundant) ',
     1          'internal coordinates')
      END
      SUBROUTINE HND_BSTR(EQVAL,NOINT,I,J,C,B,NCART,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A
C          BOND STRETCH AS DEFINED BY WILSON (SEE WDC P.55) -----
C
      DIMENSION C(3,1),B(NCART,1)
      DIMENSION RIJ(3)
      DATA ZERO /0.0D+00/
      DIJSQ = ZERO
      DO 100 M = 1,3
      RIJ(M) = C(M,J)-C(M,I)
  100 DIJSQ = DIJSQ+RIJ(M)*RIJ(M)
      DO 120 M = 1,3
      NOCOL1 = 3*(I-1)+M
      NOCOL2 = 3*(J-1)+M
      B(NOCOL1,NOINT) = -RIJ(M)/ SQRT(DIJSQ)
  120 B(NOCOL2,NOINT) = -B(NOCOL1,NOINT)
      EQVAL = SQRT(DIJSQ)*BOHR
      RETURN
      END
      SUBROUTINE HND_BEND(EQVAL,NOINT,I,J,K,C,B,NCART,RTOD,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS OF A
C          VALENCE ANGLE BENDING COORDINATE AS DEFINED BY WILSON.
C          SEE WDC P. 56 -----
C
C     -----I AND K ARE THE NUMBERS OF THE END ATOMS.  J IS THE
C          NUMBER OF THE CENTRAL ATOM -----
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,1),B(NCART,1)
      DIMENSION RJI(3),RJK(3),EJI(3),EJK(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-BEND  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DO 100 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      RJK(M) = C(M,K)-C(M,J)
      DJISQ = DJISQ+RJI(M)*RJI(M)
  100 DJKSQ = DJKSQ+RJK(M)*RJK(M)
      DJI =  SQRT(DJISQ)
      DJK =  SQRT(DJKSQ)
      DOTJ = ZERO
      DO 120 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
  120 DOTJ = DOTJ+EJI(M)*EJK(M)
      IF (ONE- ABS(DOTJ)) 200,200,140
  140 continue
      SINJ = SQRT(ONE-DOTJ*DOTJ)
      DJI=DJI*BOHR
      DJK=DJK*BOHR
      DO 160 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NOINT) = ((DOTJ*EJI(M)-EJK(M)))/(DJI*SINJ)
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NOINT) = ((DOTJ*EJK(M)-EJI(M)))/(DJK*SINJ)
      NOCOL3 = 3*(J-1)+M
  160 B(NOCOL3,NOINT) = -B(NOCOL1,NOINT)-B(NOCOL2,NOINT)
      DOT = ZERO
      DO 180 M = 1,3
  180 DOT = DOT+EJI(M)*EJK(M)
      EQVAL = ACOS(DOT)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
      GO TO 220
  200 CONTINUE
      WRITE(IW,9008)
      CALL HND_HNDERR(3,ERRMSG)
  220 CONTINUE
      RETURN
 9008 FORMAT(' ROUNDOFF ERROR IN -BEND- STOP ')
      END
      SUBROUTINE HND_TORS(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR THE
C           TORSION AS DEFINED BY WILSON.  SEE WDC P60.
C           NOTE THAT THE SIGN OF THE SECOND TERM ON THE RIGHT HAND
C           SIDE IN EQUATION (22) IS INCORRECT AS PRINTED ON P60 OF
C           WDC.
C
C     ----- THE ORDER OF ATOMS IS I--J--K--L. IF THE OBSERVER LOOKS
C           ALONG THE VECTOR J-->K WITH J NEARER THE OBSERVER, THEN
C           THE CLOCKWISE ROTATION OF J-->I WHICH SUPERPOSES J--> I
C           WITH K-->L IS A NEGATIVE TORSION ANGLE.
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,1),B(NCART,1)
      DIMENSION RIJ(3),RJK(3),RKL(3)
      DIMENSION EIJ(3),EJK(3),EKL(3)
      DIMENSION CR1(3),CR2(3),D(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-TORS  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOLRD     /1.0001D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
C     ----- CALCULATING UNIT VECTORS ALONG I-->J, J-->K AND K-->L.
C
      DIJSQ = ZERO
      DJKSQ = ZERO
      DKLSQ = ZERO
      DO 120 M = 1,3
         RIJ(M) = C(M,J)-C(M,I)
         DIJSQ  = DIJSQ+RIJ(M)*RIJ(M)
         RJK(M) = C(M,K)-C(M,J)
         DJKSQ  = DJKSQ+RJK(M)*RJK(M)
         RKL(M) = C(M,L)-C(M,K)
         DKLSQ  = DKLSQ+RKL(M)*RKL(M)
  120 CONTINUE
      DIJ =  SQRT(DIJSQ)
      DJK =  SQRT(DJKSQ)
      DKL =  SQRT(DKLSQ)
      DO 180 M = 1,3
         EIJ(M) = RIJ(M)/DIJ
         EJK(M) = RJK(M)/DJK
         EKL(M) = RKL(M)/DKL
  180 CONTINUE
C
C                   -->   --->   --->     -->   --->   --->
C     ----- FORMING CR1 = I--J X J--K AND CR2 = J--K X K--L
C
      CR1(1) = EIJ(2)*EJK(3)-EIJ(3)*EJK(2)
      CR1(2) = EIJ(3)*EJK(1)-EIJ(1)*EJK(3)
      CR1(3) = EIJ(1)*EJK(2)-EIJ(2)*EJK(1)
      CR2(1) = EJK(2)*EKL(3)-EJK(3)*EKL(2)
      CR2(2) = EJK(3)*EKL(1)-EJK(1)*EKL(3)
      CR2(3) = EJK(1)*EKL(2)-EJK(2)*EKL(1)
C
C     ----- CALCULATING COS AND SIN OF ANGLES I-J-K AND J-K-L.
C
      DOTPJ = ZERO
      DOTPK = ZERO
      DO 220 M = 1,3
         DOTPJ = DOTPJ-EIJ(M)*EJK(M)
         DOTPK = DOTPK-EJK(M)*EKL(M)
  220 CONTINUE
      IF (ONE- ABS(DOTPJ)) 320,340,240
  240 IF (ONE- ABS(DOTPK)) 360,380,260
  260 continue
      SINPJ =  SQRT(ONE-DOTPJ*DOTPJ)
      SINPK =  SQRT(ONE-DOTPK*DOTPK)
C
C     ----- CALCULATING THE B MATRIX ELEMENTS.
C
      DIJ=DIJ*BOHR
      DJK=DJK*BOHR
      DKL=DKL*BOHR
      DO 280 M = 1,3
         SMI = -CR1(M)/(DIJ*SINPJ*SINPJ)
         NOCOL1 = 3*(I-1)+M
         B(NOCOL1,NOINT) = SMI
         F1 = (CR1(M)*(DJK-DIJ*DOTPJ))/(DJK*DIJ*SINPJ*SINPJ)
         F2 = (DOTPK*CR2(M))/(DJK*SINPK*SINPK)
         SMJ = F1-F2
         NOCOL2 = 3*(J-1)+M
         B(NOCOL2,NOINT) = SMJ
         SML = CR2(M)/(DKL*SINPK*SINPK)
         NOCOL3 = 3*(L-1)+M
         B(NOCOL3,NOINT) = SML
         NOCOL4 = 3*(K-1)+M
         B(NOCOL4,NOINT) = -(SMI+SMJ+SML)
  280 CONTINUE
C
C     ----- CALCULATING THE MAGNITUDE OF THE TORSION ANGLE, WHICH IS
C           BETWEEN 0.0 AND 180.0 DEGREES.
C
      DOT = ZERO
      DO 300 M = 1,3
         DOT = DOT+(CR1(M)*CR2(M))/(SINPJ*SINPK)
  300 CONTINUE
      IF( ABS(DOT) .GT. TOLRD) GO TO 400
      IF( ABS(DOT) .GT. ONE  ) DOT = SIGN(ONE,DOT)
      EQVAL = ACOS(DOT)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
C
C     ----- CALCULATING THE TORSION ANGLE WITH THE CONVENTIONAL SIGN.
C
C                       -PI .LE. TORSION ANGLE .LE. +PI
C
C               -->  -->   -->
C           HERE D = CR1 X CR2
C
      D(1) = CR1(2)*CR2(3)-CR1(3)*CR2(2)
      D(2) = CR1(3)*CR2(1)-CR1(1)*CR2(3)
      D(3) = CR1(1)*CR2(2)-CR1(2)*CR2(1)
      EQVAL =   SIGN(EQVAL,(D(1)*EJK(1)+D(2)*EJK(2)+D(3)*EJK(3)))
      RETURN
C
C     ----- ERROR SITUATIONS.
C
  320 WRITE (IW,9999) I,J,K,DOTPJ
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  340 WRITE (IW,9998) I,J,K
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  360 WRITE (IW,9999) J,K,L,DOTPK
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  380 WRITE (IW,9998) J,K,L
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  400 WRITE (IW,9997) I,J,K,L
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(' COS ',I3,' - ',I3,' - ',I3,' = ',1PE12.5,'?',
     1       ' ROUNDOFF ERROR. RUN STOPPED.')
 9998 FORMAT(' ATOMS ',I3,1X,I3,1X,I3,' ARE COLLINEAR. RUN STOPPED.')
 9997 FORMAT(' ROUNDOFF ERROR IN CALCULATING TORSION ANGLE FOR',
     1       I3,' - ',I3,' - ',I3,' - ',I3,'. RUN STOPPED.')
      END
      SUBROUTINE HND_OPLA(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,PIO2,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR AN
C          OUT OF PLANE WAGGING COORDINATE AS DEFINED BY WILSON.
C          SEE WDC P58. -----
C     -----I IS THE END ATOM, J IS THE APEX ATOM, AND
C          K AND L ARE THE ANCHOR ATOMS -----
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION B(NCART,1),C(3,1)
      DIMENSION RJI(3),RJK(3),RJL(3)
      DIMENSION EJI(3),EJK(3),EJL(3)
      DIMENSION C1(3),C2(3),C3(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-OPLA  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DJLSQ = ZERO
      DO 100 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      DJISQ  = DJISQ+RJI(M)*RJI(M)
      RJK(M) = C(M,K)-C(M,J)
      DJKSQ  = DJKSQ+RJK(M)*RJK(M)
      RJL(M) = C(M,L)-C(M,J)
  100 DJLSQ  = DJLSQ+RJL(M)*RJL(M)
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DJL = SQRT(DJLSQ)
      DO 120 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
  120 EJL(M) = RJL(M)/DJL
      C1(1) = EJK(2)*EJL(3)-EJK(3)*EJL(2)
      C1(2) = EJK(3)*EJL(1)-EJK(1)*EJL(3)
      C1(3) = EJK(1)*EJL(2)-EJK(2)*EJL(1)
      C2(1) = EJL(2)*EJI(3)-EJL(3)*EJI(2)
      C2(2) = EJL(3)*EJI(1)-EJL(1)*EJI(3)
      C2(3) = EJL(1)*EJI(2)-EJL(2)*EJI(1)
      C3(1) = EJI(2)*EJK(3)-EJI(3)*EJK(2)
      C3(2) = EJI(3)*EJK(1)-EJI(1)*EJK(3)
      C3(3) = EJI(1)*EJK(2)-EJI(2)*EJK(1)
      DOTI = ZERO
      DO 140 M = 1,3
  140 DOTI = DOTI+EJK(M)*EJL(M)
      IF (ONE- ABS(DOTI)) 240,240,160
  160 SINI = SQRT(ONE-DOTI*DOTI)
      DET = EJI(1)*C1(1)+EJI(2)*C1(2)+EJI(3)*C1(3)
      SINT = DET/SINI
      IF (ONE- ABS(SINT)) 260,260,180
  180 COST = SQRT(ONE-SINT*SINT)
      TANT = SINT/COST
      DJI=DJI*BOHR
      DJK=DJK*BOHR
      DJL=DJL*BOHR
      DO 200 M = 1,3
      NOCOL1 = 3*(I-1)+M
      SMI = ((C1(M)/(COST*SINI))-(TANT*EJI(M)))/DJI
      B(NOCOL1,NOINT) = SMI
      NOCOL2 = 3*(K-1)+M
      SMK = ((C2(M)/(COST*SINI))-((TANT*(EJK(M)-DOTI*EJL(M)))/(SINI*
     1     SINI)))/DJK
      B(NOCOL2,NOINT) = SMK
      NOCOL3 = 3*(L-1)+M
      SML = ((C3(M)/(COST*SINI))-((TANT*(EJL(M)-DOTI*EJK(M)))/(SINI*
     1     SINI)))/DJL
      B(NOCOL3,NOINT) = SML
      NOCOL4 = 3*(J-1)+M
      B(NOCOL4,NOINT) = -(SMI+SMK+SML)
  200 CONTINUE
      DOT = ZERO
      DO 220 M = 1,3
  220 DOT = DOT+C1(M)*EJI(M)/SINI
      THETA = ACOS(DOT)
      THETA = PIO2-THETA
      EQVAL = THETA/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
      GO TO 260
  240 CONTINUE
      WRITE (IW,9008)
      CALL HND_HNDERR(3,ERRMSG)
  260 CONTINUE
      RETURN
 9008 FORMAT(' ROUNDOFF ERROR IN -OPLA- STOP ')
      END
      SUBROUTINE HND_LIBE(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,AA,NA,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A
C           PAIR OF PERPENDICULAR LINEAR BENDING COORDINATES.  SEE
C           FERIGLE AND MEISTER, J CHEM PHYS 19,982(1951) -----
C
C     ----- I AND K ARE THE END ATOMS AND J IS THE CENTRAL ATOM.
C           NO2 IS THE NUMBER OF THE SECOND INTERNAL COORDINATE,
C           I.E. THE ONE IN A PLANE PERPENDICULAR TO THE FIRST. A(I)
C           GIVES THE CARTESIAN COORDINATES OF A POINT IN SPACE SUCH
C           THAT THE VECTOR FROM ATOM J TO POINT A DEFINES
C           THE ORIENTATION OF THE COORDINATES IN SPACE.
C           THE FIRST OF THE TWO INTERNAL COORDINATES IS IN THE I-A-K
C           PLANE AND THE SECOND IN A PLANE PERPENDICULAR TO THE FIRST
C           THROUGH POINTS I,J, AND K -----
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION B(NCART,1),C(3,1)
      DIMENSION AA(3,1),A(3)
      DIMENSION RJI(3),RJK(3)
      DIMENSION UNIT(3),UP(3),UN(3)
      DIMENSION EJI(3),EJK(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-LIBE  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
      DATA PT0001    /0.0001D+00/
C
      IF(L.LE.NA) GO TO 120
      WRITE(IW,9048) L
      CALL HND_HNDERR(3,ERRMSG)
  120 CONTINUE
      A(1)=AA(1,L)
      A(2)=AA(2,L)
      A(3)=AA(3,L)
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DAJSQ = ZERO
      DO 140 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      DJISQ  = DJISQ+RJI(M)*RJI(M)
      RJK(M) = C(M,K)-C(M,J)
      DJKSQ  = DJKSQ+RJK(M)*RJK(M)
      UN(M)  = A(M)-C(M,J)
  140 DAJSQ  = DAJSQ+UN(M)*UN(M)
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DAJ = SQRT(DAJSQ)
      DOTJ = ZERO
      DOTP = ZERO
      DO 160 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
      UNIT(M) = UN(M)/DAJ
      DOTJ = DOTJ+EJI(M)*EJK(M)
  160 DOTP = DOTP+EJI(M)*UNIT(M)
      TEST = ( ABS(DOTJ)-ONE)
      IF (PT0001- ABS(TEST)) 280,180,180
  180 CONTINUE
      TEST = ( ABS(DOTP)-ONE)
      IF (PT0001- ABS(TEST)) 200,300,300
  200 CONTINUE
      UP(1) = EJK(2)*UNIT(3)-EJK(3)*UNIT(2)
      UP(2) = EJK(3)*UNIT(1)-EJK(1)*UNIT(3)
      UP(3) = EJK(1)*UNIT(2)-EJK(2)*UNIT(1)
      DUM=ONE/ SQRT(UP(1)**2+UP(2)**2+UP(3)**2)
      UP(1)=UP(1)*DUM
      UP(2)=UP(2)*DUM
      UP(3)=UP(3)*DUM
      UNIT(1) = UP(2)*EJK(3)-UP(3)*EJK(2)
      UNIT(2) = UP(3)*EJK(1)-UP(1)*EJK(3)
      UNIT(3) = UP(1)*EJK(2)-UP(2)*EJK(1)
      DJI=DJI*BOHR
      DJK=DJK*BOHR
C
C     ----- FIRST COMPONENT -----
C
      DO 220 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NOINT) = -UNIT(M)/DJI
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NOINT) = -UNIT(M)/DJK
      NOCOL3 = 3*(J-1)+M
  220 B(NOCOL3,NOINT) = (ONE/DJI+ONE/DJK)*UNIT(M)
      EQVAL = ACOS(DOTJ)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
C
C     ----- SECOND COMPONENT -----
C
      NO2=NOINT+1
      DO 260 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NO2) = -UP(M)/DJI
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NO2) = -UP(M)/DJK
      NOCOL3 = 3*(J-1)+M
  260 B(NOCOL3,NO2) = (ONE/DJI+ONE/DJK)*UP(M)
      GO TO 320
  280 CONTINUE
      WRITE (IW,9008) I,J,K,DOTJ
      CALL HND_HNDERR(3,ERRMSG)
  300 CONTINUE
      WRITE (IW,9028) I,J,K,L
      CALL HND_HNDERR(3,ERRMSG)
  320 CONTINUE
      RETURN
 9008 FORMAT(' IN -LIBE- I,J,K ARE NOT COLLINEAR. I,J,K,DOTJ = ',
     1 3I5,F15.9)
 9028 FORMAT(' IN -LIBE- DUMMY ATOM IS COLLINEAR WITH I,J,K.',
     1 ' I,J,K,L = ',4I5)
 9048 FORMAT(' IN -LIBE- DUMMY ATOM IS OUT OF RANGE. L = ',I5)
      END
      SUBROUTINE HND_DIHPLA(DIHANG,NOINT,IZMAT,CARTC,BMAT,NCART,
     1                      DTORAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION IZMAT(5), CARTC(3,1), BMAT(NCART,1),
     1          A(3), B(3), C(3), D(3), E1(3), E2(3), E3(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG   /'PROGRAM ','STOP IN ','-DIHPLA-'/
      DATA ZERO,ONE /0.0D+00,1.0D+00/
      DATA TOL      /1.0D-06/
C
C
C     COMPUTE THE B MATRIX AND DIHEDRAL ANGLE BETWEEN 5 ATOMS
C     FORMING TWO PLANES HAVING ONLY ONE ATOM IN COMMON.
C
C     WRITTEN BY STEPHEN T. ELBERT, AMES LABORATORY, ISU, 2/23/83
C
C     ON ENTRY,
C     IW     - PRINT MESSAGE FILE
C     NOINT  - NUMBER OF THE INTERNAL COORDINATE BEING EVALUATED
C     NCART  - NUMBER OF COORDINATES (3 * NATOMS)
C     IZMAT  - LIST OF THE FIVE ATOMS INVOLVED
C     CARTC  - CARTESIAN COORDINATES OF THE ATOMS
C     DTORAD - CONVERSION FACTOR (RADIANS TO DEGREES)
C
C     ON EXIT,
C     DIHANG - DIHEDRAL ANGLE, IN DEGREES
C     BMAT   - THE NOINTth COLUMN IS COMPUTED BY THIS ROUTINE
C
C     THE FIRST THREE ATOMS DEFINE THE FIRST PLANE,
C     THE LAST THREE ATOMS DEFINE THE SECOND PLANE,
C     THE MIDDLE ATOM IS IN BOTH PLANES.
C
      IATOM = IZMAT(1)
      JATOM = IZMAT(2)
      KATOM = IZMAT(3)
      LATOM = IZMAT(4)
      MATOM = IZMAT(5)
C
C     DEFINE VECTORS FROM THE PIVOT ATOM TO THE OTHER ATOMS
C     A= I-K, B=J-K, C=L-K, D=M-K
C
      DO 110 IXYZ = 1,3
         A(IXYZ) = CARTC(IXYZ,IATOM) - CARTC(IXYZ,KATOM)
         B(IXYZ) = CARTC(IXYZ,JATOM) - CARTC(IXYZ,KATOM)
         C(IXYZ) = CARTC(IXYZ,LATOM) - CARTC(IXYZ,KATOM)
         D(IXYZ) = CARTC(IXYZ,MATOM) - CARTC(IXYZ,KATOM)
  110 CONTINUE
C
C     E1 = A X B   IS PERPENDICULAR TO FIRST PLANE
C
      E1(1) = A(2)*B(3) - A(3)*B(2)
      E1(2) = A(3)*B(1) - A(1)*B(3)
      E1(3) = A(1)*B(2) - A(2)*B(1)
      E1MAG =  SQRT( E1(1)**2 + E1(2)**2 + E1(3)**2)
      IF(E1MAG.LT.TOL) GO TO 810
C
C     E2 = C X D   IS PERPENDICULAR TO SECOND PLANE
C
      E2(1) = C(2)*D(3) - C(3)*D(2)
      E2(2) = C(3)*D(1) - C(1)*D(3)
      E2(3) = C(1)*D(2) - C(2)*D(1)
      E2MAG =  SQRT( E2(1)**2 + E2(2)**2 + E2(3)**2)
      IF(E2MAG.LT.TOL) GO TO 820
C
C     NORMALIZE E1 AND E2
C
      E1MAGI = ONE/E1MAG
      E2MAGI = ONE/E2MAG
      DO 120 IXYZ=1,3
         E1(IXYZ) = E1(IXYZ) * E1MAGI
         E2(IXYZ) = E2(IXYZ) * E2MAGI
  120 CONTINUE
C
C     THE DIRECTION OF E2 IS DEFINED SUCH THAT WHEN THE PLANES ARE
C     ROTATED ALONG THE INTERSECTION E3 = E1 X E2 SO THEY COINCIDE
C     (DIHEDRAL ANGLE=0), THE ATOMS 1 AND 4 LIE ON THE SAME SIDE OF E3.
C     (IF ATOMS 1 AND/OR 4 LIE ON E3, THEN ATOMS 2 AND/OR 5 ARE USED)
C
      ADOTE2 = A(1)*E2(1) + A(2)*E2(2) + A(3)*E2(3)
      IF(ADOTE2.LT.TOL) ADOTE2 = B(1)*E2(1) + B(2)*E2(2) + B(3)*E2(3)
      CDOTE1 = C(1)*E1(1) + C(2)*E1(2) + C(3)*E1(3)
      IF(CDOTE1.LT.TOL) CDOTE1 = D(1)*E1(1) + D(2)*E1(2) + D(3)*E1(3)
C
C     POSSIBLY CHANGE SENSE OF E2 SO A (OR B) WILL BE ON THE SAME
C     SIDE OF E1 X E3 AS C (OR D)
C
      IF(ADOTE2*CDOTE1 .GT.ZERO) GO TO 140
         ADOTE2 = -ADOTE2
         DO 130 IXYZ=1,3
            E2(IXYZ) = -E2(IXYZ)
  130    CONTINUE
  140 CONTINUE
C
C      EVALUATE DIHEDRAL ANGLE    0 < ABS(DIHANG) < 180
C      ANGLE IS POSITIVE FOR SECOND PLANE 'ABOVE' ATOM I
C      ANGLE IS NEGATIVE FOR SECOND PLANE 'BELOW' ATOM I
C
      E1DE2 = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
      IF((      ABS(E1DE2)-ONE)  .GT.ZERO) GO TO 840
      IF(( ABS( ABS(E1DE2)-ONE)) .LT.TOL ) GO TO 830
      DIHANG = ACOS(E1DE2)
      IF(ADOTE2.GT.ZERO) DIHANG = -DIHANG
      SINDI = ONE / SIN(DIHANG)
      DIHANG = DIHANG/DTORAD
C
C     E3 = E1 X E2 IS ALONG INTERSECTION LINE OF THE PLANES
C
      E3(1) = SINDI * (E1(2)*E2(3) - E1(3)*E2(2))
      E3(2) = SINDI * (E1(3)*E2(1) - E1(1)*E2(3))
      E3(3) = SINDI * (E1(1)*E2(2) - E1(2)*E2(1))
C
      F1 =  E1MAGI * (B(1)*E3(1) + B(2)*E3(2) + B(3)*E3(3))
      F2 = -E1MAGI * (A(1)*E3(1) + A(2)*E3(2) + A(3)*E3(3))
      F4 = -E2MAGI * (D(1)*E3(1) + D(2)*E3(2) + D(3)*E3(3))
      F5 =  E2MAGI * (C(1)*E3(1) + C(2)*E3(2) + C(3)*E3(3))
C
      I = 3*IATOM - 3
      J = 3*JATOM - 3
      K = 3*KATOM - 3
      L = 3*LATOM - 3
      M = 3*MATOM - 3
C
      DO 160 N=1,3
         B1 = E1(N) * F1
         B2 = E1(N) * F2
         B4 = E2(N) * F4
         B5 = E2(N) * F5
         B3 = -(B1+B2+B4+B5)
         BMAT(NOINT,I+N) = B1
         BMAT(NOINT,J+N) = B2
         BMAT(NOINT,K+N) = B3
         BMAT(NOINT,L+N) = B4
         BMAT(NOINT,M+N) = B5
  160 CONTINUE
      RETURN
C                                        ERROR MESSAGES
  810 WRITE(IW,910)
      GO TO 880
  820 WRITE(IW,920)
      GO TO 880
  830 WRITE(IW,930)
      GO TO 880
  840 WRITE(IW,940) E1DE2
  880 CONTINUE
      WRITE(IW,950) NOINT,IATOM,JATOM,KATOM,LATOM,MATOM
      CALL HND_HNDERR(3,ERRMSG)
C
  910 FORMAT(1X,'---- ERROR, FIRST THREE ATOMS ARE COLLINEAR')
  920 FORMAT(1X,'---- ERROR, LAST THREE ATOMS ARE COLLINEAR ')
  930 FORMAT(1X,'---- ERROR, ALL FIVE ATOMS ARE COPLANAR ',
     1          '     THE -B- MATRIX COLUMN CAN NOT BE DEFINED. STOP')
  940 FORMAT(1X,'---- COSIN OF DIHEDRAL ANGLE=',F16.10,
     1          ' IS GREATER THAN ONE')
  950 FORMAT(1X,'THIS MESSAGE COMES TO YOU FROM -DIHPLA-, WHERE ',
     1          'INTERNAL COORDINATE',I4,' USES ATOMS',5I4)
      END
      SUBROUTINE HND_ZMTMAK(ZDONE,geom,oprint)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "util.fh"
      integer  geom
      logical  oprint
      parameter (max_zcoord=16)
c
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)    
      PARAMETER (MXCOOR=1500)    
      PARAMETER (MXIZMT=6000)
      PARAMETER (MXBOND= 64)
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXANGS= 8*MXATOM)
      PARAMETER (MXTORS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      PARAMETER (MXLINB= 8*MXATOM)
      PARAMETER (MXLNBA=10       )
      PARAMETER (MXSEG=64)
      LOGICAL SOME,OUT,DBUG
      LOGICAL DBUG1
      LOGICAL BADZ
      LOGICAL REMOVE
      LOGICAL NUBOND
      LOGICAL MODBND
      LOGICAL MODANG
      LOGICAL MODTOR
      LOGICAL MODOOP
      LOGICAL MODLNB
      LOGICAL ENDATM
      LOGICAL ENDMOD
      LOGICAL  DONE
      LOGICAL IDONE
      LOGICAL ZDONE
      logical geom_autoz_info_get
      external geom_autoz_info_get
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)  
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTDAT/ZMAT(MXCOOR),IZMAT(MXIZMT),NIZMAT
      DIMENSION XYZLNB(3,MXLNBA)
      DIMENSION IJBOND(2,MXBOND)
      DIMENSION IJKANG(3,MXANGS)
      DIMENSION IJKLTO(4,MXTORS)
      DIMENSION IJKLOP(4,MXOOPA)
      DIMENSION IJKLNB(4,MXLINB)
      DIMENSION IJBNDS(2,MXBNDS)
      DIMENSION IJMODS(2,MXBNDS)
      DIMENSION NIBOND(MXATOM),IIBOND(MXATOM)
      DIMENSION NIMODS(MXATOM),IIMODS(MXATOM)
      DIMENSION ENDATM(MXATOM)
      DIMENSION ENDMOD(MXATOM)
      DIMENSION  IDONE(MXATOM)
      DIMENSION IATSEG(MXATOM,MXSEG)
      DIMENSION NAMSEG(       MXSEG)
      DIMENSION LENSEG(       MXSEG)
      DIMENSION RCOV(103)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTMAK-'/
      DATA ZERO   /0.0D+00/
C
C      COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
C      EXCEPT VAN DER WAALS RADII FOR HE,NE,AR,KR (SAME SOURCE),
C      AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT,RN
C
      DATA (RCOV(NUCZ),NUCZ=1,2)
     1  /0.40D+00,1.22D+00/     ! RJH Changed 0.3 for H to 0.4 for H2
      DATA (RCOV(NUCZ),NUCZ=3,10)
     1  /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     2   0.70D+00,0.66D+00,0.71D+00,1.60D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     1  /1.66D+00,1.36D+00,1.25D+00,1.17D+00,
     2   1.10D+00,1.04D+00,0.99D+00,1.91D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     1  /2.03D+00,1.74D+00,
     2   1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     3  1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     4   1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.98D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     1  /2.22D+00,1.92D+00,
     2   1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     3   1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     4   1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,2.09D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     1  /2.35D+00,1.98D+00,
     2   1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D+00,
     3   1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D+00,
     4   1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     5   1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     6   1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,2.20D+00/
      DATA (RCOV(I),i=87,103) 
     1  /3.24d00,2.68d00,2.25d00,2.16d00,1.93d00,
     2   1.66d00,1.57d00,1.81d00,2.21d00,1.43d00,
     3   1.42d00,1.40d00,1.39d00,1.38d00,1.37d00,
     4   1.36d00,1.34d00/
C
      DIST(IAT,JAT)=SQRT((C(1,IAT)-C(1,JAT))**2+
     1                   (C(2,IAT)-C(2,JAT))**2+
     2                   (C(3,IAT)-C(3,JAT))**2) * TOANGS
      BADZ(NUCZ) = NUCZ.LE.0.OR.NUCZ.GT.103
c
      DBUG1=.FALSE.
      DBUG =.FALSE.
      OUT  =.FALSE.
      OUT  =OUT.OR.DBUG
      SOME =.TRUE.
      SOME =SOME.OR.OUT
      IF(SOME .and. oprint) THEN
         WRITE(IW,8888)
      ENDIF
c
      ncvrpass = 0
10101 ncvrpass = ncvrpass + 1
C
      IF(DBUG) THEN
         DO IAT=1,NAT
            WRITE(IW,9976) IAT,NUC(IAT),C(1,IAT),C(2,IAT),C(3,IAT)
         ENDDO
      ENDIF
C
      NZVAR =0
      NIZMAT=0
      DO I=1,MXIZMT
         IZMAT(I)=0
      ENDDO
      DO I=1,MXCOOR
         ZMAT(I)=ZERO
      ENDDO
      DO MBNDS=1,MXBNDS
         IJBNDS(1,MBNDS)=0
         IJBNDS(2,MBNDS)=0
      ENDDO
      DO MBOND=1,MXBOND
         IJBOND(1,MBOND)=0
         IJBOND(2,MBOND)=0
      ENDDO
      DO MANGS=1,MXANGS
         DO I=1,3
            IJKANG(I,MANGS)=0
         ENDDO
      ENDDO
      DO MTORS=1,MXTORS
         DO I=1,4
            IJKLTO(I,MTORS)=0
         ENDDO
      ENDDO
      DO MOOPA=1,MXOOPA
         DO I=1,4
            IJKLOP(I,MOOPA)=0
         ENDDO
      ENDDO
      DO MLINB=1,MXLINB
         DO I=1,4
            IJKLNB(I,MLINB)=0
         ENDDO
      ENDDO
      DO MLNBA=1,MXLNBA
         DO I=1,3
            XYZLNB(I,MLNBA)=ZERO
         ENDDO
      ENDDO
      CVR_SCALING=1.20D+00
C
C     ----- READ IMPOSED BONDS, ANGLES, TORSIONS -----    
C     ----- O-O-PLANE, LINEAR BENDS, ... IF ANY  -----
C
c           read in zcoord data from -rtdb-
c    CVR_SCALING,IJBOND,IJKANG,IJKLTO,IJKLOP,IJKLNB
c
c
      if(.not. geom_autoz_info_get(geom,
     $     cvr_scaling, ijbond, ijkang, ijklto, ijklop, 
     $     ijklnb)) then
         continue               ! No zcoord data input
      endif
c
      if (ncvrpass .eq. 1) CVFAC=CVR_SCALING ! Use default on first pass
c
      IF(DBUG) THEN
         WRITE(IW,*) 'CVFAC = ',CVFAC
      ENDIF
C
      NUMBND=0
      DO MBOND=1,MXBOND
         ICON=IJBOND(1,MBOND)
         JCON=IJBOND(2,MBOND)
         IF(ICON.GT.0.AND.JCON.GT.0) THEN
            NUMBND=NUMBND+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODBND=NUMBND.NE.0
      IF(DBUG) THEN
         IF(MODBND) THEN
            WRITE(IW,8887)
         ELSE
            WRITE(IW,8886)
         ENDIF
      ENDIF
C
      NUMANG=0
      DO MANGS=1,MXANGS
         ICON=IJKANG(1,MANGS)
         JCON=IJKANG(2,MANGS)
         KCON=IJKANG(3,MANGS)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.KCON.GT.0) THEN
            NUMANG=NUMANG+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND.KCON.EQ.0) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODANG=NUMANG.NE.0
      IF(DBUG) THEN
         IF(MODANG) THEN
            WRITE(IW,8885)
         ELSE
            WRITE(IW,8884)
         ENDIF
      ENDIF
C
      NUMTOR=0
      DO MTORS=1,MXTORS
         ICON=IJKLTO(1,MTORS)
         JCON=IJKLTO(2,MTORS)
         KCON=IJKLTO(3,MTORS)
         LCON=IJKLTO(4,MTORS)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.
     1      KCON.GT.0.AND.LCON.GT.0     ) THEN
            NUMTOR=NUMTOR+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND. 
     1          KCON.EQ.0.AND.LCON.EQ.0     ) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODTOR=NUMTOR.NE.0
      IF(DBUG) THEN
         IF(MODTOR) THEN
            WRITE(IW,8883)
         ELSE
            WRITE(IW,8882)
         ENDIF
      ENDIF
C
      NUMOOP=0
      DO MOOPA=1,MXOOPA
         ICON=IJKLOP(1,MOOPA)
         JCON=IJKLOP(2,MOOPA)
         KCON=IJKLOP(3,MOOPA)
         LCON=IJKLOP(4,MOOPA)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.
     1      KCON.GT.0.AND.LCON.GT.0     ) THEN
            NUMOOP=NUMOOP+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND.
     1          KCON.EQ.0.AND.LCON.EQ.0     ) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODOOP=NUMOOP.NE.0
      IF(DBUG) THEN
         IF(MODOOP) THEN
            WRITE(IW,8881)
         ELSE
            WRITE(IW,8880)
         ENDIF
      ENDIF
C
      NUMLNB=0
      DO MLINB=1,MXLINB
         ICON=IJKLNB(1,MLINB)
         JCON=IJKLNB(2,MLINB)
         KCON=IJKLNB(3,MLINB)
         LCON=IJKLNB(4,MLINB)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.
     1      KCON.GT.0.AND.LCON.GT.0     ) THEN
            NUMLNB=NUMLNB+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND.
     1          KCON.EQ.0.AND.LCON.EQ.0     ) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODLNB=NUMLNB.NE.0
      IF(DBUG) THEN
         IF(MODLNB) THEN
            WRITE(IW,8879)
         ELSE
            WRITE(IW,8878)
         ENDIF
      ENDIF
C
C     ----- SET UP CONNECTIVITY TABLE, INCLUDING FORCED BONDS -----
C
      IF(DBUG) THEN
         DO IAT=1,NAT
            WRITE(IW,9976) IAT,NUC(IAT),C(1,IAT),C(2,IAT),C(3,IAT)
         ENDDO
      ENDIF
C
      MBNDS=0
      DO IAT=1,NAT
      IIBOND(IAT)=MBNDS
C
C     ----- FIRST CHECK FORCED BONDS -----
C
         IF(MODBND) THEN
            DO MBOND=1,MXBOND
               NUBOND=.FALSE.
               IF(IJBOND(1,MBOND).EQ.IAT) THEN
                  MBNDS=MBNDS+1
                  IJBNDS(1,MBNDS)=IAT
                  IJBNDS(2,MBNDS)=IJBOND(2,MBOND)
                  NUBOND=.TRUE.
               ELSEIF(IJBOND(2,MBOND).EQ.IAT) THEN
                  MBNDS=MBNDS+1
                  IJBNDS(1,MBNDS)=IAT
                  IJBNDS(2,MBNDS)=IJBOND(1,MBOND)
                  NUBOND=.TRUE.
               ENDIF
               IF(NUBOND.AND.DBUG) THEN
                  WRITE(IW,9975) MBNDS,IJBNDS(1,MBNDS),IJBNDS(2,MBNDS)
               ENDIF
C
C     ----- CHECK FOR DUPLICATES AND ELIMINATE -----
C
               IF(NUBOND.AND.MBNDS.GT.1) THEN
                  M=MBNDS-1
                  DO I=1,M
                     IF(IJBNDS(1,I).EQ.IJBNDS(1,MBNDS).AND.
     1                  IJBNDS(2,I).EQ.IJBNDS(2,MBNDS)     ) THEN
                        MBNDS=MBNDS-1
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
C
C     ----- DISCOVER CONNECTIVITY TO OTHER ATOMS -----
C
         DO JAT=1,NAT
            IF(JAT.NE.IAT) THEN
               IZ=NUC(IAT)
               JZ=NUC(JAT)
               IF(BADZ(IZ).OR.BADZ(JZ)) THEN
                  WRITE(IW,9993) IAT,JAT
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               NUBOND=.FALSE.
               RIJ=DIST(IAT,JAT)
               RCV=(RCOV(IZ)+RCOV(JZ))
               IF(RIJ.LE.RCV*CVFAC) THEN
                  IF(MBNDS.EQ.MXBNDS) THEN
                     WRITE(IW,*) 'TOO MANY BONDS FOUND',
     1                           ' MBNDS, MXBNDS = ',MBNDS, MXBNDS
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  MBNDS=MBNDS+1
                  IJBNDS(1,MBNDS)=IAT
                  IJBNDS(2,MBNDS)=JAT
                  NUBOND=.TRUE.
                  IF(DBUG1) THEN
                     WRITE(IW,9992) IAT,JAT,IZ,JZ,RIJ,RCV,MBNDS,
     1                              IJBNDS(1,MBNDS),IJBNDS(2,MBNDS),
     2                              NUBOND
                  ENDIF
C
C     ----- CHECK FOR DUPLICATES AND ELIMINATE -----
C
                  IF(NUBOND.AND.MBNDS.GT.1) THEN
                     REMOVE=.FALSE.
                     M=MBNDS-1
                     DO I=1,M
                        IF(IJBNDS(1,I).EQ.IJBNDS(1,MBNDS).AND.
     1                     IJBNDS(2,I).EQ.IJBNDS(2,MBNDS)     ) THEN
                           IF(DBUG) THEN
                              WRITE(IW,9974) MBNDS,I
                           ENDIF
                           REMOVE=.TRUE.
                        ENDIF
                     ENDDO
                     IF(REMOVE) THEN
                        IJBNDS(1,MBNDS)=0
                        IJBNDS(2,MBNDS)=0
                        MBNDS=MBNDS-1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
C
         NIBOND(IAT)=MBNDS-IIBOND(IAT)
         IF(NIBOND(IAT).GT.1) THEN
            ENDATM(IAT)=.FALSE.
         ELSE
            ENDATM(IAT)=.TRUE.
         ENDIF
      ENDDO
      NBNDS=MBNDS
C
      IF(NBNDS.LE.0) THEN
         WRITE(IW,9999)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- DEBUG PRINTING OF CONNECTIVITY -----
C
      IF(DBUG1) THEN
         WRITE(IW,9998) NBNDS
         WRITE(IW,9997)
         DO IAT=1,NAT
            WRITE(IW,9996) IAT,NIBOND(IAT),IIBOND(IAT),ENDATM(IAT)
         ENDDO
         DO IAT=1,NAT
            WRITE(IW,9995) IAT,NIBOND(IAT)
            IF(NIBOND(IAT).GT.0) THEN
               IIBND1=IIBOND(IAT)+1
               IIBND2=IIBOND(IAT)+NIBOND(IAT)
               WRITE(IW,9994) (IJBNDS(1,I),IJBNDS(2,I),I=IIBND1,IIBND2)
            ENDIF
         ENDDO
      ENDIF
C
C
C     -----    TERMINAL ATOMS HAVE BEEN DETECTED   -----
C           TRY FOR SEGMENTS ( RINGS AND BRIDGES )
C
      DO I=1,NBNDS
         IJMODS(1,I)=IJBNDS(1,I)
         IJMODS(2,I)=IJBNDS(2,I)
      ENDDO
      MXCONN=0
      DO IAT=1,NAT
         ENDMOD(IAT)=ENDATM(IAT)
         IIMODS(IAT)=IIBOND(IAT)
         NIMODS(IAT)=NIBOND(IAT)
         IF(NIMODS(IAT).GT.MXCONN) THEN
            MXCONN=NIMODS(IAT)
         ENDIF
         DO MSEG=1,MXSEG
            IATSEG(IAT,MSEG)=0
         ENDDO
      ENDDO
      DO MSEG=1,MXSEG
         NAMSEG(MSEG)=0
         LENSEG(MSEG)=0
      ENDDO
C
C     -----     TAKE AWAY -END ATOM- CONNECTIVITY       -----
C           RING AND BRIDGE CONNECTIVITY SHOULD BE LEFT
C
      NPASS=MXCONN
      DO IPASS=1,NPASS
         IF(DBUG) THEN
            WRITE(IW,9991) IPASS,NPASS
         ENDIF
         DO IAT=1,NAT
            IF(NIMODS(IAT).EQ.1) THEN
               DO I=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
                  KAT=IJMODS(1,I)
                  JAT=IJMODS(2,I)
                  IJMODS(1,I)=0
                  IJMODS(2,I)=0
                  IF(KAT.NE.0) THEN
                     IF(NIMODS(JAT).GE.1) THEN
                        DO J=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                           IF(IJMODS(2,J).EQ.IAT) THEN
                              IJMODS(1,J)=0
                              IJMODS(2,J)=0
                              NIMODS(JAT)=NIMODS(JAT)-1
                           ENDIF
                        ENDDO
                     ELSE
                        WRITE(IW,9978) JAT,IAT
                        CALL HND_HNDERR(3,ERRMSG)
                     ENDIF
                     IF(NIMODS(JAT).EQ.1) THEN
                        ENDMOD(JAT)=.TRUE.
                     ENDIF
                  ENDIF
               ENDDO
               NIMODS(IAT)=NIMODS(IAT)-1
            ENDIF
         ENDDO
C
         IF(DBUG) THEN
            DO IAT=1,NAT
               WRITE(IW,9996) IAT,NIMODS(IAT),IIMODS(IAT),ENDMOD(IAT)
            ENDDO
            DO IAT=1,NAT
               WRITE(IW,9995) IAT,NIMODS(IAT)
               IF(NIBOND(IAT).GT.0) THEN
                  IIMOD1=IIBOND(IAT)+1
                  IIMOD2=IIBOND(IAT)+NIBOND(IAT)
                  WRITE(IW,9994) (IJMODS(1,I),IJMODS(2,I),
     1                                                 I=IIMOD1,IIMOD2)
               ENDIF
            ENDDO
         ENDIF
C
      ENDDO
C
C     ----- CHECK IF THERE WILL BE MORE TO DO ... -----
C
      DONE=.TRUE.
      DO IAT=1,NAT
         IF(NIMODS(IAT).GT.0) THEN
            IDONE(IAT)=.FALSE.
         ELSE
            IDONE(IAT)=.TRUE. 
         ENDIF
         DONE=DONE.AND.IDONE(IAT)
      ENDDO
      DONE=.TRUE.
      IF(DONE) THEN
         IF(DBUG) THEN
            WRITE(IW,9980)
         ENDIF
C
C     ----- CREATE -IZMAT- -----
C
         NZVAR =0
         NIZMAT=0
         DO I=1,MXIZMT
            IZMAT(I)=0
         ENDDO
         DO I=1,NBNDS
            IJMODS(1,I)=IJBNDS(1,I)
            IJMODS(2,I)=IJBNDS(2,I)
         ENDDO
         DO IAT=1,NAT
            ENDMOD(IAT)=ENDATM(IAT)
            IIMODS(IAT)=IIBOND(IAT)
            NIMODS(IAT)=NIBOND(IAT)
         ENDDO
         CALL HND_ZMTBLD(ZDONE,NZVAR,NIZMAT,IZMAT,C,
     1                   NAT,NBNDS,IJMODS,IIMODS,NIMODS,NUMBND,
     2                   IJKANG,NUMANG,IJKLTO,NUMTOR,
     3                   IJKLOP,NUMOOP,IJKLNB,NUMLNB,XYZLNB)
c
c     Attempt to recover from a failure by increasing cvfac.
c
         if (.not. zdone) then
            if (cvfac .lt. 2.5d0) then
               write(6,9876) ncvrpass, cvfac, cvfac+0.1d0
 9876          format(//i5,' AUTOZ failed with cvr_scaling =',f4.1,
     $              ' increasing to',f4.1)
               cvfac = cvfac + 0.1d0
               goto 10101
            endif
         endif
C
C     ----- SAVE -IZMAT(1:NZMAT)-   -----
C           SAVE -XYZLNB(1:MXLNBA)-
C
         NZMAT=NIZMAT
         NLNBA=3*MXLNBA
C
         IF(SOME .and. oprint) THEN
            WRITE(IW,8877) ZDONE
         ENDIF
         RETURN
      ENDIF
C
C     ----- NOW COMPRESS MODIFIED CONNECTIVITY TABLE -----
C
      DO IAT=1,NAT
         IF(NIMODS(IAT).GT.0) THEN
            NUMODS=0
            DO I=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
               IF(IJMODS(1,I).NE.0.AND.IJMODS(2,I).NE.0) THEN
                  NUMODS=NUMODS+1
                  IJMODS(1,NUMODS)=IJBNDS(1,I)
                  IJMODS(2,NUMODS)=IJBNDS(2,I)
               ENDIF
            ENDDO
            IF(NIMODS(IAT).NE.NUMODS) THEN
               WRITE(IW,9985) IAT,NIMODS(IAT),NUMODS
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            IIBOND(IAT)=IIMODS(IAT)
            NIBOND(IAT)=IIMODS(IAT)
            ENDATM(IAT)=ENDMOD(IAT)
         ENDIF
      ENDDO
C
      IF(DBUG) THEN
         DO IAT=1,NAT
            WRITE(IW,9996) IAT,NIBOND(IAT),IIBOND(IAT),ENDATM(IAT)
         ENDDO
         DO IAT=1,NAT
            WRITE(IW,9995) IAT,NIMODS(IAT)
            IF(NIMODS(IAT).GT.0) THEN
               IIMOD1=IIMODS(IAT)+1
               IIMOD2=IIMODS(IAT)+NIMODS(IAT)
               WRITE(IW,9994) (IJMODS(1,I),IJMODS(2,I),I=IIMOD1,IIMOD2)
            ENDIF
         ENDDO
      ENDIF
C
C     ----- NOW WE HAVE ONLY ATOMS WITH REMAINING CONNECTIVITY -----
C           OF 2 OR MORE. THEY BELONG TO RINGS AND BRIDGES.
C
      NUCONN=0
      DO IAT=1,NAT
         IF(NIMODS(IAT).GT.NUCONN) THEN
            NUCONN=NIMODS(IAT)
         ENDIF
      ENDDO
C
      IF(NUCONN.LT.2) THEN
C
C     ----- WE SHOULD HAVE ONLY ( > 2 ) CONNECTIVITIES -----
C
         WRITE(IW,9989) NUCONN
         CALL HND_HNDERR(3,ERRMSG)
      ELSEIF(NUCONN.GT.3) THEN
C
C     ----- -NUCONN- > 3 .... NOT TREATED YET -----
C
         WRITE(IW,9989) NUCONN
         CALL HND_HNDERR(3,ERRMSG)
C
      ELSEIF(NUCONN.EQ.3) THEN
C
C     ----- -NUCONN- = 3 .... SEVERAL RINGS , FUSED OR NOT -----
C                             BRIDGES AS WELL.
C           FIRST , LIST OUT ALL THE CONNECTIVITY SEGMENTS
C
         IF(DBUG) THEN
            WRITE(IW,9984)
         ENDIF
C
         NSEG  =0
         IANCHR=0
  100    IANCHR=IANCHR+1
         IF(IANCHR.GT.NAT) THEN
            GO TO 120
         ENDIF
         IF(NIMODS(IANCHR).NE.3) THEN
            GO TO 100
         ENDIF
            DO I=1,3
               NSEG=NSEG+1
               IF(NSEG.GT.MXSEG) THEN
                  WRITE(IW,9983) NSEG
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
C
               JATOM=1
               IATSEG(JATOM,NSEG)=IANCHR
C
               IAT=IANCHR
               JAT=IJMODS(2,IIMODS(IAT)+I)
  110          JATOM=JATOM+1
               IATSEG(JATOM,NSEG)=JAT
               IF(NIMODS(JAT).EQ.2) THEN
                  IF(IJMODS(2,IIMODS(JAT)+1).NE.IAT) THEN
                     KAT=IJMODS(2,IIMODS(JAT)+1)
                  ELSE
                     KAT=IJMODS(2,IIMODS(JAT)+2)
                  ENDIF
                  IAT=JAT
                  JAT=KAT
                  GO TO 110
               ELSEIF(NIMODS(JAT).EQ.3) THEN
                  LENSEG(NSEG)=JATOM
                  IF(JAT.EQ.IANCHR) THEN
                     NAMSEG(NSEG)=1
                  ELSE
                     NAMSEG(NSEG)=2
                  ENDIF
               ELSE
                  WRITE(IW,9982) JAT
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ENDDO
C
  120    CONTINUE
         IF(DBUG) THEN
            DO MSEG=1,NSEG
               WRITE(IW,9981) MSEG,LENSEG(MSEG)
               WRITE(IW,9979) (IATSEG(IATOM,MSEG),IATOM=1,LENSEG(MSEG))
            ENDDO
         ENDIF
C
      ELSEIF(NUCONN.EQ.2) THEN
C
C     ----- -NUCONN- = 2 .... DISCONNECTED RINGS -----
C
         IF(DBUG) THEN
            WRITE(IW,9990)
         ENDIF
C
         NSEG  =0
         IANCHR=0
  200    IANCHR=IANCHR+1
         IF(IANCHR.GT.NAT) THEN
            GO TO 220
         ENDIF
         IF(NIMODS(IANCHR).NE.2) THEN
            GO TO 200
         ENDIF
            NSEG =NSEG+1
            IATOM=1
            IATSEG(IATOM,NSEG)=IANCHR
C
            IAT=IANCHR
            JAT=IJMODS(2,IIMODS(IAT)+1)
  210       IATOM=IATOM+1
            IATSEG(IATOM,NSEG)=JAT
            IF(NIMODS(JAT).EQ.2) THEN
               IF(IJMODS(2,IIMODS(JAT)+1).NE.IAT) THEN
                  KAT=IJMODS(2,IIMODS(JAT)+1)
               ELSE
                  KAT=IJMODS(2,IIMODS(JAT)+2)
               ENDIF
               IF(KAT.NE.IANCHR) THEN
                  IAT=JAT
                  JAT=KAT
                  GO TO 210
               ENDIF
            ELSE
               WRITE(IW,9986) JAT
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            LENSEG(NSEG)=IATOM
C
         IF(DBUG) THEN
            DO MSEG=1,NSEG
               WRITE(IW,9987) LENSEG(MSEG)
               WRITE(IW,9979) (IATSEG(IATOM,MSEG),IATOM=1,LENSEG(MSEG))
            ENDDO
         ENDIF
  220    CONTINUE
C
      ENDIF
C
      RETURN
 9999 FORMAT(' SOMETHING IS WRONG ... NO BONDS COULD BE DETECTED. STOP')
 9998 FORMAT(' TOTAL NUMBER OF BONDS = ',I5)
 9997 FORMAT(1X,'ATOM   NUM. OF BONDS   FIRST BOND  ENDATM',/,
     1       1X,41(1H-))
 9996 FORMAT(I5,6X,I5,5X,6X,I5,2X,4X,L1,3X)
 9995 FORMAT(' CONNECTIVITY FOR ATOM -IAT- = ',I4,' NUM.BONDS = ',I5)
 9994 FORMAT(6(1H(,I4,1H,,I4,1H)))
 9993 FORMAT(' NO COVALENT RADIUS SPECIFIED FOR ONE OF THESE ATOMS = ',
     1       2I5)
 9992 FORMAT(' IAT,JAT,IZ,JZ,RIJ,RCV,MBNDS= ',2I4,2I4,2F10.3,3I5,L4)
 9991 FORMAT(/,' PASS NUM. = ',I3,' OUT OF ',I3,' PASS(ES) ')
 9990 FORMAT(/,' SEARCHING FOR A SINGLE RING ...')
 9989 FORMAT(/,' ILLEGAL -MXCONN- HERE, -MXCONN- = ',I3,' . STOP.')
 9988 FORMAT(/,' THIS ROUTINE DOES NOT HANDLE -MXCONN- .GT. 3. STOP.')
 9987 FORMAT(/,' A ',I2,'-MEMBER RING WAS FOUND. ATOMS ARE = ',/,20I4)
 9986 FORMAT(/,' ATOM ',I3,' SHOULD HAVE TWO CONNECTIONS,',
     1         ' AND ONLY TWO. STOP.')
 9985 FORMAT(/,' ERROR DURING COMPRESSION OF CONNECTIVITY TABLE',
     1         ' FOR ATOM = ',I4,/,' NIMODS(IAT), NUMODS = ',2I5)
 9984 FORMAT(/,' SEARCHING FOR MULTIPLE RINGS AND BRIDGES ...')
 9983 FORMAT(/,' TOO MANY SEGMENTS FROM 3-CONNECTED ATOMS.',
     1         ' NSEG, MXSEG = ',2I4)
 9982 FORMAT(/,' ATOM ',I3,' SHOULD HAVE AT LEAST TWO CONNECTIONS.',
     1         ' STOP.')
 9981 FORMAT(/,' SEARCHING SEGMENT NO. = ',I4)
 9980 FORMAT(/,' ALL 1-CONNECTIVITY COMPLETED.')
 9979 FORMAT(16I4)
 9978 FORMAT(/,' ATOM ',I4,' HAS NO CONNECTIVITY LEFT,',/,
     1         ' BUT THE CODE IS TRYING TO REMOVE CONNECTION',
     2         ' TO ATOM ',I4,' . STOP.')
 9977 FORMAT(' INCORRECT DATA IN NAMELIST -COORDI- . STOP .')
 9976 FORMAT(2I5,3F10.6)
 9975 FORMAT(' FORCED BOND INSERTED = ',3I5)
 9974 FORMAT(' BOND REMOVED = ',I5,' EQUIVALENT TO = ',I5)
 8888 FORMAT(/,10X,6(1H-),/,10X,'auto-z',/,10X,6(1H-))
 8887 FORMAT(' SOME IMPOSED BONDS           IN -ZMTMAK- ')
 8886 FORMAT('   NO IMPOSED BONDS           IN -ZMTMAK- ')
 8885 FORMAT(' SOME IMPOSED ANGLES          IN -ZMTMAK- ')
 8884 FORMAT('   NO IMPOSED ANGLES          IN -ZMTMAK- ')
 8883 FORMAT(' SOME IMPOSED TORSIONS        IN -ZMTMAK- ')
 8882 FORMAT('   NO IMPOSED TORSIONS        IN -ZMTMAK- ')
 8881 FORMAT(' SOME IMPOSED O-O-PLANE-BENDS IN -ZMTMAK- ')
 8880 FORMAT('   NO IMPOSED O-O-PLANE-BENDS IN -ZMTMAK- ')
 8879 FORMAT(' SOME IMPOSED LINEAR-BENDS    IN -ZMTMAK- ')
 8878 FORMAT('   NO IMPOSED LINEAR-BENDS    IN -ZMTMAK- ')
 8877 FORMAT(' zdone = ',L4)
      END
      SUBROUTINE HND_ZMTBLD(ZDONE,NZVAR,NIZMAT,IZMAT,C,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,NUMBND,
     2                      IJKANG,NUMANG,IJKLTO,NUMTOR,
     3                      IJKLOP,NUMOOP,
     4                      IJKLNB,NUMLNB,XYZLNB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXANGS= 8*MXATOM)
      PARAMETER (MXTORS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      PARAMETER (MXLINB= 8*MXATOM)
      LOGICAL     DBUG
      LOGICAL     OUT
      LOGICAL     ZDONE
      LOGICAL     DOOOPA
      LOGICAL     DOLINB
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION IJKANG(3,*),IJKLTO(4,*),IJKLOP(4,*)
      DIMENSION IJKLNB(4,*),XYZLNB(3,*)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTBLD-'/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
      IF(NAT.LE.1) THEN
         ZDONE=.TRUE.
         RETURN
      ENDIF
C
C     ----- BONDS -----
C
      CALL HND_ZMTYP1(NZVAR,NIZMAT,IZMAT,
     1                NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,NUMBND)
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,9995) NZVAR
      ENDIF
C
C     ----- SOME CHECKS FOR DIATOMICS AND TRIATOMICS -----
C
      IF(NAT.EQ.2) THEN
         IF(NZVAR.EQ.1) THEN
            ZDONE=.TRUE. 
            RETURN
         ELSE
            NZVAR=0
            ZDONE=.FALSE.
            RETURN
         ENDIF
      ENDIF
      IF(NAT.EQ.3) THEN
         IF(NZVAR.NE.2.AND.NZVAR.NE.3) THEN
            NZVAR=0
            ZDONE=.FALSE.
            RETURN
         ENDIF
      ENDIF
C
C     ----- ANGLES -----
C
      DOLINB=.FALSE.
      CALL HND_ZMTYP2(NZVAR,NIZMAT,IZMAT,C,
     1                NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKANG,NUMANG,
     2                DOLINB)
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,9994) NZVAR
      ENDIF
C
C     ----- MORE CHECKS FOR TRIATOMICS -----
C
      IF(NAT.EQ.3) THEN
         IF(DOLINB) THEN
            NZVAR=0
            ZDONE=.FALSE.
            RETURN
         ENDIF
      ENDIF
C
C     ----- TORSIONS -----
C
      CALL HND_ZMTYP3(NZVAR,NIZMAT,IZMAT,C,
     1                NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKLTO,NUMTOR)
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,9993) NZVAR
      ENDIF
C
c     RJH.  Michel says OOP were working, but generated far
c     .     too many degrees of freedom to be useful.
c
c     .     The linear bends are not working but a fix 
c     .     could be to include a point(s) at some distance
c     .     and just use the regular bends.
c
c     .     In the meantime just punt after the torsions
c     .     and then analyze the generated coordinates to
c     .     make sure that we have at least (3n-6) coordinates.
c
c$$$      DOOOPA=.FALSE.
c$$$      IF(DOOOPA) THEN
c$$$C
c$$$C     ----- OUT-OF-PLANE BENDS -----
c$$$C
c$$$         CALL HND_ZMTYP4(NZVAR,NIZMAT,IZMAT,
c$$$     1                   NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKLOP,NUMOOP)
c$$$         IF(DBUG) THEN
c$$$            WRITE(IW,9998)
c$$$            WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
c$$$            WRITE(IW,9996)
c$$$            WRITE(IW,9992) NZVAR
c$$$         ENDIF
c$$$      ENDIF
c$$$C
c$$$C     ----- LINEAR BENDS -----
c$$$C
c$$$      IF(DOLINB) THEN
c$$$         CALL HND_ZMTYP5(NZVAR,NIZMAT,IZMAT,C,
c$$$     1                   NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
c$$$     2                   IJKLNB,NUMLNB,XYZLNB,IJKANG,NUMANG)
c$$$         IF(DBUG) THEN
c$$$            WRITE(IW,9998)
c$$$            WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
c$$$            WRITE(IW,9996)
c$$$            WRITE(IW,9991) NZVAR
c$$$         ENDIF
c$$$         ZDONE=.FALSE.
c$$$         RETURN
c$$$      ENDIF
C
C     ----- DONE -----
C
      IF(NZVAR.LT.3*NAT-6) THEN
         ZDONE=.FALSE.
      ELSE
         ZDONE=.TRUE.
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTBLD- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' AFTER -ZMTYP1- , NZVAR = ',I5)
 9994 FORMAT(' AFTER -ZMTYP2- , NZVAR = ',I5)
 9993 FORMAT(' AFTER -ZMTYP3- , NZVAR = ',I5)
 9992 FORMAT(' AFTER -ZMTYP4- , NZVAR = ',I5)
 9991 FORMAT(' AFTER -ZMTYP5- , NZVAR = ',I5)
      END
      SUBROUTINE HND_ZMTYP1(NZVAR,NIZMAT,IZMAT,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,NUMBND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IONE=1)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      LOGICAL     DBUG
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION IZMAT(*)   
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP1-'/
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.1) RETURN
C
C     ----- BONDS -----
C
      NBND=NBOND
      DO IBND=1,NBND
         II=IJBOND(1,IBND)
         IJ=IJBOND(2,IBND)
         IF(II.NE.0) THEN
            IF(IBND.LT.NBND) THEN
               DO JBND=IBND+1,NBND
                  JI=IJBOND(1,JBND)
                  JJ=IJBOND(2,JBND)
                  IF(JI.NE.0) THEN
                     IF((JI.EQ.II.AND.JJ.EQ.IJ).OR.        
     1                  (JI.EQ.IJ.AND.JJ.EQ.II)    ) THEN
                        IJBOND(1,JBND)=0
                        IJBOND(2,JBND)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
             IZMAT(NIZMAT+1)=IONE
             IZMAT(NIZMAT+2)=II
             IZMAT(NIZMAT+3)=IJ   
            NIZMAT=NIZMAT+3
            NZVAR =NZVAR +1
            NUMBND=NUMBND+1
         ENDIF
      ENDDO
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP1- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
      END
      SUBROUTINE HND_ZMTYP2(NZVAR,NIZMAT,IZMAT,C,
     1           NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKANG,NUMANG,
     2           DOLINB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ITWO=2)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXANGS= 8*MXATOM)
      LOGICAL     DBUG
      LOGICAL     DOLINB
      LOGICAL     LINBND
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJKANG(3,*)     
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP2-'/
      DATA TOL    /1.0D-08/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      DISTSQ(IAT,JAT)=(C(1,IAT)-C(1,JAT))**2+
     1                (C(2,IAT)-C(2,JAT))**2+
     2                (C(3,IAT)-C(3,JAT))**2
      DOTPRD(IAT,JAT,KAT)=(C(1,IAT)-C(1,JAT))*(C(1,KAT)-C(1,JAT))+
     1                    (C(2,IAT)-C(2,JAT))*(C(2,KAT)-C(2,JAT))+
     2                    (C(3,IAT)-C(3,JAT))*(C(3,KAT)-C(3,JAT))
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
C     ----- ANGLES -----
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.2) RETURN
C
      NANG=0
      DO IAT=1,NAT
         DO IJBND=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
            ICON=IJBOND(1,IJBND)
            JCON=IJBOND(2,IJBND)
            IF(ICON.NE.IAT) THEN
               WRITE(IW,9995) IJBND,IAT,ICON,JCON
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            JAT=JCON
            IF(NIBOND(JAT).GT.1) THEN
               DO JKBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                  JCON=IJBOND(1,JKBND)
                  KCON=IJBOND(2,JKBND)
                  IF(JCON.NE.JAT) THEN
                     WRITE(IW,9995) JKBND,JAT,JCON,KCON
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  KAT=KCON
                  IF(KAT.NE.IAT) THEN
                     IF(NANG.EQ.MXANGS) THEN
                        WRITE(IW,*) 'TOO MANY BOND ANGLES',
     1                              ' NANG,MXANG = ',NANG,MXANGS
                        CALL HND_HNDERR(3,ERRMSG)
                     ENDIF
                     NANG=NANG+1
                     IJKANG(1,NANG)=IAT
                     IJKANG(2,NANG)=JAT
                     IJKANG(3,NANG)=KAT
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF(DBUG) THEN
         IF(NANG.GT.1) THEN
            WRITE(IW,9994)
            DO MANG=1,NANG
               WRITE(IW,9993) MANG,IJKANG(1,MANG),
     1                             IJKANG(2,MANG),IJKANG(3,MANG)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
C     ----- ELIMINATE DUPLICATES AND CLEAN UP -----
C
      IF(NANG.GE.2) THEN
         DO IANG=1,NANG-1
            II=IJKANG(1,IANG)
            IJ=IJKANG(2,IANG)
            IK=IJKANG(3,IANG)
            IF(IJKANG(1,IANG).NE.0) THEN
               DO JANG=IANG+1,NANG
                  JI=IJKANG(1,JANG)
                  JJ=IJKANG(2,JANG)
                  JK=IJKANG(3,JANG)
                  IF(JJ.EQ.IJ) THEN
                     IF((II.EQ.JI.AND.IK.EQ.JK).OR.
     1                  (II.EQ.JK.AND.IK.EQ.JI)    ) THEN
                        IJKANG(1,JANG)=0
                        IJKANG(2,JANG)=0
                        IJKANG(3,JANG)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      IANG=0
  100 IANG=IANG+1
      IF(IANG.GT.NANG) GO TO 110
         IF(IJKANG(1,IANG).EQ.0) THEN
            IF(IANG.LT.NANG) THEN
               DO JANG=IANG+1,NANG
                  IJKANG(1,JANG-1)=IJKANG(1,JANG)
                  IJKANG(2,JANG-1)=IJKANG(2,JANG)
                  IJKANG(3,JANG-1)=IJKANG(3,JANG)
               ENDDO
               IANG=IANG-1
            ENDIF
            NANG=NANG-1
         ELSE
            RIJSQ=DISTSQ(IJKANG(1,IANG),IJKANG(2,IANG))
            RIJ  =SQRT(RIJSQ)                          
            RJKSQ=DISTSQ(IJKANG(2,IANG),IJKANG(3,IANG))
            RJK  =SQRT(RJKSQ)                           
            RIKSQ=DISTSQ(IJKANG(1,IANG),IJKANG(3,IANG))
            COSB =(RIJSQ+RJKSQ-RIKSQ)/(TWO*RIJ*RJK)
            DUM  =ABS(COSB)-ONE
            IF(DUM.GT.ZERO) THEN
               IF(DUM.LE.TOL) THEN
                  LINBND=.TRUE.
                  DOLINB=DOLINB.OR.LINBND
               ELSE
                  WRITE(IW,*) 'SOMETHING IS WRONG IN -ZMTYP2-'
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ELSE
               DUM=-DUM
               IF(DUM.LE.TOL) THEN
                   LINBND=.TRUE.
                   DOLINB=DOLINB.OR.LINBND
               ELSE
                  LINBND=.FALSE.
               ENDIF
            ENDIF
            IF(.NOT.LINBND) THEN
                IZMAT(NIZMAT+1)=ITWO
                IZMAT(NIZMAT+2)=IJKANG(1,IANG)
                IZMAT(NIZMAT+3)=IJKANG(2,IANG)
                IZMAT(NIZMAT+4)=IJKANG(3,IANG)
               NIZMAT=NIZMAT+4
               NZVAR =NZVAR +1
               NUMANG=NUMANG+1
            ENDIF
         ENDIF
         GO TO 100
  110 CONTINUE
C
      IF(DBUG) THEN
         IF(NANG.GT.1) THEN
            WRITE(IW,9994)
            DO MANG=1,NANG
               WRITE(IW,9993) MANG,IJKANG(1,MANG),
     1                             IJKANG(2,MANG),IJKANG(3,MANG)
            ENDDO
            IF(DOLINB) THEN
               WRITE(IW,9991)
            ENDIF
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP2- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' SOMETHING WRONG WITH BOND ',I2,' OF ATOM ',I4,
     1       ' WITH -ICON-, -JCON- = ',2I5)
 9994 FORMAT(' TABLE OF ANGLES = ')
 9993 FORMAT(4I5)
 9992 FORMAT(' NO BOND ANGLES FOUND ')
 9991 FORMAT(' COLLINEAR ATOMS FOUND. DO -LINEAR BENDS- ... ')
      END
      SUBROUTINE HND_ZMTYP3(NZVAR,NIZMAT,IZMAT,C,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
     2                                       IJKLTO,NUMTOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ITHREE=3)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXTORS= 8*MXATOM)
      LOGICAL     DBUG
      LOGICAL     LINIJK,LINJKL
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJKLTO(4,*)      
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP3-'/
      DATA TOL    /1.0D-08/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      DISTSQ(IAT,JAT)=(C(1,IAT)-C(1,JAT))**2+
     1                (C(2,IAT)-C(2,JAT))**2+
     2                (C(3,IAT)-C(3,JAT))**2
      DOTPRD(IAT,JAT,KAT)=(C(1,IAT)-C(1,JAT))*(C(1,KAT)-C(1,JAT))+
     1                    (C(2,IAT)-C(2,JAT))*(C(2,KAT)-C(2,JAT))+
     2                    (C(3,IAT)-C(3,JAT))*(C(3,KAT)-C(3,JAT))
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
C     ----- TORSIONS -----
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.3) RETURN
C
      NTOR=0
      DO IAT=1,NAT
         DO IJBND=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
            ICON=IJBOND(1,IJBND)
            JCON=IJBOND(2,IJBND)
            IF(ICON.NE.IAT) THEN
               WRITE(IW,9995) IJBND,IAT,ICON,JCON
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            JAT=JCON
            IF(NIBOND(JAT).GT.1) THEN
               DO JKBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                  JCON=IJBOND(1,JKBND)
                  KCON=IJBOND(2,JKBND)
                  IF(JCON.NE.JAT) THEN
                     WRITE(IW,9995) JKBND,JAT,JCON,KCON
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  KAT=KCON
                  IF(KAT.NE.IAT) THEN
                     IF(NIBOND(KAT).GT.1) THEN
                        DO KLBND=IIBOND(KAT)+1,IIBOND(KAT)+NIBOND(KAT)
                           KCON=IJBOND(1,KLBND)
                           LCON=IJBOND(2,KLBND)
                           IF(KCON.NE.KAT) THEN
                              WRITE(IW,9995) KLBND,KAT,KCON,LCON
                              CALL HND_HNDERR(3,ERRMSG)
                           ENDIF
                           LAT=LCON
                           IF(LAT.NE.JAT) THEN
                              IF(NTOR.EQ.MXTORS) THEN
                                 WRITE(IW,*) 'TOO MANY TORSIONS',
     1                                       ' NTOR, MXTORS = ',NTOR,
     2                                       MXTORS
                                 CALL HND_HNDERR(3,ERRMSG)
                              ENDIF
                              if ( iat.ne.jat .and. ! RJH from (h2o)5 test
     $                             iat.ne.kat .and. 
     $                             iat.ne.lat .and.
     $                             jat.ne.kat .and.
     $                             jat.ne.lat .and.
     $                             kat.ne.lat) then
                                 NTOR=NTOR+1
                                 IJKLTO(1,NTOR)=IAT
                                 IJKLTO(2,NTOR)=JAT
                                 IJKLTO(3,NTOR)=KAT
                                 IJKLTO(4,NTOR)=LAT
                              endif
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF(DBUG) THEN
         IF(NTOR.GT.1) THEN
            WRITE(IW,9994)
            DO MTOR=1,NTOR
               WRITE(IW,9993) MTOR,IJKLTO(1,MTOR),IJKLTO(2,MTOR),
     1                             IJKLTO(3,MTOR),IJKLTO(4,MTOR)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
C     ----- ELIMINATE DUPLICATES AND CLEAN UP -----
C
      IF(NTOR.GE.2) THEN
         DO ITOR=1,NTOR-1
            II=IJKLTO(1,ITOR)
            IJ=IJKLTO(2,ITOR)
            IK=IJKLTO(3,ITOR)
            IL=IJKLTO(4,ITOR) 
            IF(IJKLTO(1,ITOR).NE.0) THEN
               DO JTOR=ITOR+1,NTOR
                  JI=IJKLTO(1,JTOR)
                  JJ=IJKLTO(2,JTOR)
                  JK=IJKLTO(3,JTOR)
                  JL=IJKLTO(4,JTOR)
                  IF(JJ.EQ.IJ.AND.JK.EQ.IK) THEN
                     IF(II.EQ.JI.AND.IL.EQ.JL) THEN
                        IJKLTO(1,JTOR)=0
                        IJKLTO(2,JTOR)=0
                        IJKLTO(3,JTOR)=0
                        IJKLTO(4,JTOR)=0
                     ENDIF
                  ELSEIF(JJ.EQ.IK.AND.JK.EQ.IJ) THEN
                     IF(II.EQ.JL.AND.IL.EQ.JI) THEN
                        IJKLTO(1,JTOR)=0
                        IJKLTO(2,JTOR)=0
                        IJKLTO(3,JTOR)=0
                        IJKLTO(4,JTOR)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      ITOR=0
  100 ITOR=ITOR+1
      IF(ITOR.GT.NTOR) GO TO 110
         IF(IJKLTO(1,ITOR).EQ.0) THEN
            IF(ITOR.LT.NTOR) THEN
               DO JTOR=ITOR+1,NTOR
                  IJKLTO(1,JTOR-1)=IJKLTO(1,JTOR)
                  IJKLTO(2,JTOR-1)=IJKLTO(2,JTOR)
                  IJKLTO(3,JTOR-1)=IJKLTO(3,JTOR)
                  IJKLTO(4,JTOR-1)=IJKLTO(4,JTOR)
               ENDDO
               ITOR=ITOR-1
            ENDIF
            NTOR=NTOR-1
         ELSE
            RIJSQ=DISTSQ(IJKLTO(1,ITOR),IJKLTO(2,ITOR))
            RIJ  =SQRT(RIJSQ)
            RJKSQ=DISTSQ(IJKLTO(2,ITOR),IJKLTO(3,ITOR))
            RJK  =SQRT(RJKSQ)
            RKLSQ=DISTSQ(IJKLTO(3,ITOR),IJKLTO(4,ITOR))
            RKL  =SQRT(RKLSQ)
            RIKSQ=DISTSQ(IJKLTO(1,ITOR),IJKLTO(3,ITOR))
            RJLSQ=DISTSQ(IJKLTO(2,ITOR),IJKLTO(4,ITOR))
            COSBJ=(RIJSQ+RJKSQ-RIKSQ)/(TWO*RIJ*RJK)
            COSBK=(RJKSQ+RKLSQ-RJLSQ)/(TWO*RJK*RKL)
            DUM  =ABS(COSBJ)-ONE
            IF(DUM.GT.ZERO) THEN
               IF(DUM.LE.TOL) THEN
                  LINIJK=.TRUE.
               ELSE
                  WRITE(IW,*) 'SOMETHING IS WRONG IN -ZMTYP2-'
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ELSE
               DUM=-DUM
               IF(DUM.LE.TOL) THEN
                   LINIJK=.TRUE.
               ELSE
                  LINIJK=.FALSE.
               ENDIF
            ENDIF
            DUM  =ABS(COSBK)-ONE
            IF(DUM.GT.ZERO) THEN
               IF(DUM.LE.TOL) THEN
                  LINJKL=.TRUE.
               ELSE
                  WRITE(IW,*) 'SOMETHING IS WRONG IN -ZMTYP2-'
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ELSE
               DUM=-DUM
               IF(DUM.LE.TOL) THEN
                   LINJKL=.TRUE.
               ELSE
                  LINJKL=.FALSE.
               ENDIF
            ENDIF
            IF(.NOT.LINIJK.AND..NOT.LINJKL) THEN
                IZMAT(NIZMAT+1)=ITHREE
                IZMAT(NIZMAT+2)=IJKLTO(1,ITOR)
                IZMAT(NIZMAT+3)=IJKLTO(2,ITOR)
                IZMAT(NIZMAT+4)=IJKLTO(3,ITOR)
                IZMAT(NIZMAT+5)=IJKLTO(4,ITOR)
               NIZMAT=NIZMAT+5
               NZVAR =NZVAR +1
            ENDIF
         ENDIF
         GO TO 100
  110 CONTINUE
C
      IF(DBUG) THEN
         IF(NTOR.GT.1) THEN
            WRITE(IW,9994)
            DO MTOR=1,NTOR
               WRITE(IW,9993) MTOR,IJKLTO(1,MTOR),IJKLTO(2,MTOR),
     1                             IJKLTO(3,MTOR),IJKLTO(4,MTOR)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP3- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' SOMETHING WRONG WITH BOND ',I2,' OF ATOM ',I4,
     1       ' WITH -ICON-, -JCON- = ',2I5)
 9994 FORMAT(' TABLE OF TORSIONS = ')
 9993 FORMAT(5I5)
 9992 FORMAT(' NO TORSIONS FOUND ')
      END
      SUBROUTINE HND_ZMTYP4(NZVAR,NIZMAT,IZMAT,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
     2                                       IJKLOP,NUMOOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IFOUR=4)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      LOGICAL     DBUG
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION IZMAT(*)   
      DIMENSION IJKLOP(4,*)      
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP4-'/
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
C     ----- OUT-OF-PLANE-BENDS -----
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.3) RETURN
C
      NOOP=0
      DO IAT=1,NAT
         DO IJBND=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
            ICON=IJBOND(1,IJBND)
            JCON=IJBOND(2,IJBND)
            IF(ICON.NE.IAT) THEN
               WRITE(IW,9995) IJBND,IAT,ICON,JCON
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            JAT=JCON
            IF(NIBOND(JAT).GT.2) THEN
               DO JKBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                  JCON=IJBOND(1,JKBND)
                  KCON=IJBOND(2,JKBND)
                  IF(JCON.NE.JAT) THEN
                     WRITE(IW,9995) JKBND,JAT,JCON,KCON
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  KAT=KCON
                  IF(KAT.NE.IAT) THEN
                     DO JLBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                        JCON=IJBOND(1,JLBND)
                        LCON=IJBOND(2,JLBND)
                        IF(JCON.NE.JAT) THEN
                           WRITE(IW,9995) JLBND,JAT,JCON,LCON
                           CALL HND_HNDERR(3,ERRMSG)
                        ENDIF
                        LAT=LCON
                        IF(LAT.NE.IAT.AND.LAT.NE.KAT) THEN
                           NKCON=NIBOND(KAT)
                           NLCON=NIBOND(LAT)
                           IF((NKCON.EQ.1.AND.NLCON.EQ.1).OR.
     1                        (NKCON.GT.1.AND.NLCON.GT.1)    ) THEN
                              IF(NOOP.EQ.MXOOPA) THEN
                              WRITE(IW,*) 'TOO MANY OUT-OF-PLANE BENDS',
     1                                    ' NOOP, MXOOPA = ',NOOP,MXOOPA
                                 CALL HND_HNDERR(3,ERRMSG)
                              ENDIF
                              NOOP=NOOP+1
                              IJKLOP(1,NOOP)=IAT
                              IJKLOP(2,NOOP)=JAT
                              IJKLOP(3,NOOP)=KAT
                              IJKLOP(4,NOOP)=LAT
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF(DBUG) THEN
         IF(NOOP.GT.1) THEN
            WRITE(IW,9994)
            DO MOOP=1,NOOP
               WRITE(IW,9993) MOOP,IJKLOP(1,MOOP),IJKLOP(2,MOOP),
     1                             IJKLOP(3,MOOP),IJKLOP(4,MOOP)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
C     ----- ELIMINATE DUPLICATES AND CLEAN UP -----
C
      IF(NOOP.GE.2) THEN
         DO IOOP=1,NOOP-1
            II=IJKLOP(1,IOOP)
            IJ=IJKLOP(2,IOOP)
            IK=IJKLOP(3,IOOP)
            IL=IJKLOP(4,IOOP) 
            IF(IJKLOP(1,IOOP).NE.0) THEN
               DO JOOP=IOOP+1,NOOP
                  JI=IJKLOP(1,JOOP)
                  JJ=IJKLOP(2,JOOP)
                  JK=IJKLOP(3,JOOP)
                  JL=IJKLOP(4,JOOP)
                  IF(JI.EQ.II.AND.JJ.EQ.IJ) THEN
                     IF((JK.EQ.IK.AND.JL.EQ.IL).OR.
     1                  (JK.EQ.IL.AND.JL.EQ.IK)    ) THEN
                        IJKLOP(1,JOOP)=0
                        IJKLOP(2,JOOP)=0
                        IJKLOP(3,JOOP)=0
                        IJKLOP(4,JOOP)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      IOOP=0
  100 IOOP=IOOP+1
      IF(IOOP.GT.NOOP) GO TO 110
         IF(IJKLOP(1,IOOP).EQ.0) THEN
            IF(IOOP.LT.NOOP) THEN
               DO JOOP=IOOP+1,NOOP
                  IJKLOP(1,JOOP-1)=IJKLOP(1,JOOP)
                  IJKLOP(2,JOOP-1)=IJKLOP(2,JOOP)
                  IJKLOP(3,JOOP-1)=IJKLOP(3,JOOP)
                  IJKLOP(4,JOOP-1)=IJKLOP(4,JOOP)
               ENDDO
               IOOP=IOOP-1
            ENDIF
            NOOP=NOOP-1
         ELSE
             IZMAT(NIZMAT+1)=IFOUR 
             IZMAT(NIZMAT+2)=IJKLOP(1,IOOP)
             IZMAT(NIZMAT+3)=IJKLOP(2,IOOP)
             IZMAT(NIZMAT+4)=IJKLOP(3,IOOP)
             IZMAT(NIZMAT+5)=IJKLOP(4,IOOP)
            NIZMAT=NIZMAT+5
            NZVAR =NZVAR +1
         ENDIF
         GO TO 100
  110 CONTINUE
C
      IF(DBUG) THEN
         IF(NOOP.GT.1) THEN
            WRITE(IW,9994)
            DO MOOP=1,NOOP
               WRITE(IW,9993) MOOP,IJKLOP(1,MOOP),IJKLOP(2,MOOP),
     1                             IJKLOP(3,MOOP),IJKLOP(4,MOOP)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP4- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' SOMETHING WRONG WITH BOND ',I2,' OF ATOM ',I4,
     1       ' WITH -ICON-, -JCON- = ',2I5)
 9994 FORMAT(' TABLE OF OUT-OF-PLANE-BENDS = ')
 9993 FORMAT(5I5)
 9992 FORMAT(' NO OUT-OF-PLANE-ANGLE FOUND ')
      END
      SUBROUTINE HND_ZMTYP5(NZVAR,NIZMAT,IZMAT,C,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
     2                      IJKLNB,NUMLNB,XYZLNB,IJKANG,NUMANG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IFOUR=4)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      LOGICAL     DBUG
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJKLNB(4,*)      
      DIMENSION XYZLNB(3,*)
      DIMENSION IJKANG(3,*)           
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP5-'/
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
         WRITE(IW,9991)
      ENDIF
C
C     ----- LINEAR-BENDS -----
C
      NLBND=0
C
      IF(NAT.LE.3) RETURN
C
      IF(DBUG) THEN
         IF(NLBND.GT.1) THEN
            WRITE(IW,9994)
            DO MLBND=1,NLBND
               WRITE(IW,9993) MLBND,IJKLNB(1,MLBND),IJKLNB(2,MLBND),
     1                              IJKLNB(3,MLBND),IJKLNB(4,MLBND)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP5- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9994 FORMAT(' TABLE OF LINEAR-BENDS = ')
 9993 FORMAT(5I5)
 9992 FORMAT(' NO LINEAR-BENDS FOUND ')
 9991 FORMAT(' LINEAR-BENDS NOT SET UP CURRENTLY . STOP . ')
      END
      subroutine hnd_dparsc(a,la,c,lc)
      implicit double precision (a-h,o-z)
      character*(*) a
      character*(*) c
      character*1   blk
      data  blk /' '/
      do ic=1,lc
         c(ic:ic)=blk
      enddo
      if(lc.le.(la-2)) then
         c(1:lc)=a(2:lc+1)
      else
         c(1:lc)=a(2:la-1)
      endif
      return
      end
      subroutine hnd_dparsi(a,la,n)
      implicit double precision (a-h,o-z)
      character*(*) a
      character*1   char(12)
      data char     /'0','1','2','3','4','5','6','7','8','9',
     1               '+','-'/
c
      i1=1
      i2=la
c
c     ----- sign -----
c
      isign=1
      if(a(i1:i1).eq.char(12)) then
         isign=-1
         i1=i1+1
      elseif(a(i1:i1).eq.char(11)) then
         isign= 1
         i1=i1+1
      else
         isign= 1
      endif
      na=i2-i1+1
c
      n=0
      do ia=i1,i2
         ib=i2-ia
         do i=1,10
            if(a(ia:ia).eq.char(i)) then
               n=n+(i-1)*10**(ib)
            endif
         enddo
      enddo
      n=n*isign
c
      return
      end
      subroutine hnd_dparsr(a,la,x)
      implicit double precision (a-h,o-z)
      logical rep
      character*(*) a
      character*1   char(17)
      data char     /'0','1','2','3','4','5','6','7','8','9',
     1               '+','-','.','e','d','E','D'/
      data zero     /0.0d+00/
      data ten      /1.0d+01/
c     
      i1 =1            
      i2 =la
c
c     ----- sign -----
c
      isign=1
      if(a(i1:i1).eq.char(12)) then
         isign=-1
         i1=i1+1
      elseif(a(i1:i1).eq.char(11)) then
         isign= 1
         i1=i1+1
      else
         isign= 1
      endif
c
c     ----- exponent -----
c
      ie2=i2
      do ie=i1,i2
         if(a(ie:ie).eq.char(14).or.a(ie:ie).eq.char(15).or. 
     1      a(ie:ie).eq.char(16).or.a(ie:ie).eq.char(17)) go to 10
      enddo
      iexp=0
      go to 50
c
 10   ie2=i2
      ie1=ie+1
      i2 =ie-1
      iexp=1
      if(a(ie1:ie1).eq.char(12)                           ) then
         iexp=-1
      endif
      if(a(ie1:ie1).eq.char(12).or.a(ie1:ie1).eq.char(11)) then
         ie1=ie1+1
      endif
      itmp=0
      do i=ie1,ie2
         do j=1,10
            if(a(i:i).eq.char(j)) go to 30
         enddo
         go to 100
 30      itmp=itmp*10+j-1
      enddo
      iexp=iexp*itmp
c
c     ----- the number itself -----
c
 50   continue
      rep=.false.
      dum=zero
      do i=i1,i2
         if(a(i:i).ne.char(13)) then
            do j=1,10
               if(a(i:i).eq.char(j)) go to 70
            enddo
            go to 100
 70         dum=dum*ten+dble(j-1)
         else
            if(rep) go to 100
            iexp=iexp+i-i2
            rep=.true.
         endif
      enddo
      dum=dum*dble(isign)*ten**iexp
      x  =dum
      return
c
 100  continue       
      return
      end
      SUBROUTINE HND_HNDERR(LERR,ERRMSG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8  ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION    ERRMSG(LERR)
C
      WRITE(IW,9999)
      WRITE(IW,9998) ERRMSG
      WRITE(IW,9997)
C
      CALL ERRQUIT(' CALLS IT QUIT FROM HND_HNDERR ',0)
 9999 FORMAT(/,1X,72(1H-),/,31X,11HJOB STOPPED)
 9998 FORMAT(/,1X,9A8)
 9997 FORMAT(/,1X,72(1H-))
      END
      SUBROUTINE HND_GEOCLS(NFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8  ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION    ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-GEOCLS-'/
      CLOSE(UNIT=NFT,STATUS='KEEP',ERR=10)
      RETURN
   10 WRITE(IW,9999) NFT
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(/,' ----- ERROR CLOSING UNIT ',I3,' IN -GEOCLS- . STOP .')
      END
      SUBROUTINE HND_GEOOPN(NFT,GEOFIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 GEOFIL
      CHARACTER*8  ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION    ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-GEOOPN-'/
      OPEN(UNIT=NFT,FILE=GEOFIL,STATUS='OLD',
     1     ACCESS='SEQUENTIAL',FORM=  'FORMATTED',ERR=10)
      RETURN
   10 WRITE(IW,9999) NFT
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(/,' ----- ERROR OPENING FILE -',A80,'-',/,
     1         '       AS UNIT ',I3,' IN -GEOOPN- . STOP.')
      END
      SUBROUTINE HND_ZDAT(IAT,NAT,NVAR,ZVAL,ZLST,
     1                ZMT,NUMZMT,PRSZMT,FLGZMT,ZMTCHR,
     2                NUMVAR,PRSVAR,VARCHR,FRZVAR,FRZVAL,LST,
     3                IZMAT,IZ,IZFRZ,MXIZMT,NZMOD,DBUG,
     $     zvarname, zvarsign)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     THIS ROUTINE TAKES THE CHARACTER STRING VALUES FOR THE
C     -Z- MATRIX INPUT FOR THE INDICES, AND THE VALUES AND TRANSFORMS
C     THEM INTO INTEGER AND REAL VALUES. THE CHARACTER STRING
C     VALUES CAN BE EITHER NUMERIC (INTEGER AND REAL), OR A VARIABLE
C     THAT GETS REPLACED BY A VALUE FOR THE VARIABLE.
C
C     INPUT PARAMETERS:
C        IAT     - INTEGER    -  ATOM NUMBER BEING CONSIDERED
C        NAT     - INTEGER    -  TOTAL NUMBER OF ATOMS
C        NVAR    - INTEGER    -  TOTAL NUMBER OF VARIABLES
C        NUMZMT  - INTEGER    -  ARRAY OF NUMBERS OF WORDS PER LINE
C        PRSZMT  - CHARACTER*80 -  ARRAY OF PARSED -Z- DATA
C        FLGZMT  - INTEGER    -  ARRAY OF TYPES FOR -Z- DATA
C                                  0 ... BLANK
C                                  1 ... FLOATING
C                                  2 ... INTEGER
C                                  3 ... ALPHANUMERIC
C        ZMTCHR  - INTEGER    -  ARRAY OF LENGTHS FOR -Z- DATA
C        NUMVAR  - INTEGER    -  ARRAY OF NUM. OF PIECES OF INFO
C                                  FOR EACH VARIABLE FOR SUBSTITUTION
C        PRSVAR  - CHARACTER*80 -  ARRAY OF VARIABLES FOR SUBSTITUTION
C        VARCHR  - INTEGER    -  ARRAY OF LENGTHS FOR VARIABLES
C
C     OUTPUT PARAMETERS:
C        ZVAL    - REAL*8       -  ARRAY OF -Z- MATRIX VALUES
C        ZMT     - INTEGER    -  ARRAY OF -Z- MATRIX INDICES
c        zvarname- char*8     -  array of -z- matrix variable names
c     .                          or blank if numeric value given.
c     .  zvarsign- real*8     -  array of -z-matrix variables signs
c     .                          (+1 or -1 as input, 1 if numeric given)
C
C     S. CHIN: 11/08/90 - IBM KINGSTON, NY
C
      PARAMETER     (MXATOM=500)
      PARAMETER     (MAXGEO=MXATOM+1,MAXZMT=40,MAXVAR=256)
      LOGICAL       DBUG
      LOGICAL       LST
      LOGICAL       READY
      INTEGER     FLGZMT
      CHARACTER*80  PRSZMT
      INTEGER     ZMTCHR
      CHARACTER*80  PRSVAR
      INTEGER     VARCHR
      INTEGER     ZMT
      LOGICAL       FRZVAL
      LOGICAL       FRZVAR
      LOGICAL       CART
      LOGICAL       IZFRZ
      CHARACTER*8   ATNAME
      CHARACTER*1   PLUS
      CHARACTER*1   MINUS
      CHARACTER*8   ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_XYZGEO/XX(MAXGEO),YY(MAXGEO),ZZ(MAXGEO),
     1                  ATNAME(MAXGEO),ATNUM(MAXGEO),CART(MAXGEO)
      DIMENSION     ERRMSG(3)
      DIMENSION     NUMZMT(       MAXGEO)
      DIMENSION     PRSZMT(MAXZMT,MAXGEO)
      DIMENSION     FLGZMT(MAXZMT,MAXGEO)
      DIMENSION     ZMTCHR(MAXZMT,MAXGEO)
      DIMENSION     PRSVAR(MAXZMT,MAXVAR)
      DIMENSION     VARCHR(MAXZMT,MAXVAR)
      DIMENSION     FRZVAR(       MAXVAR)
      DIMENSION     NUMVAR(       MAXVAR)
      DIMENSION     FRZVAL(     3,MAXGEO)
      DIMENSION       ZVAL(     3,MAXGEO)
      DIMENSION        ZMT(     5,MAXGEO)
      DIMENSION       ZLST(     3,MAXGEO)
      DIMENSION      IZMAT(MXIZMT)
      DIMENSION      IZFRZ(MXIZMT)

      character*(*) zvarname(*)
      double precision zvarsign(*)

      DATA ERRMSG /'PROGRAM ','STOP IN ','- ZDAT -'/
      DATA ZERO   /0.0D+00/
      DATA TWO    /2.0D+00/
      DATA PIDEG  /180.0D+00/
      DATA PLUS   /'+'/
      DATA MINUS  /'-'/
c
      do i = nzmod+1,nzmod+3
         zvarname(i) = ' '       ! Default for numeric parameter
         zvarsign(i) = 1d0       ! Default for numeric parameter
      end do
C
C     ---- THIS ATOM IS BEING INPUTED WITH CARTESIAN COORDINATES ----
C
      IF(FLGZMT(3,IAT).EQ.1.AND.FLGZMT(4,IAT).EQ.1) THEN
             IF(NUMZMT(  IAT).EQ.4.AND.FLGZMT(1,IAT).EQ.3.AND.
     1          FLGZMT(2,IAT).EQ.1                            ) THEN
c           READ(PRSZMT(1,IAT),*) ATNAME(IAT)
c           READ(PRSZMT(2,IAT),*) XX(IAT)
c           READ(PRSZMT(3,IAT),*) YY(IAT)
c           READ(PRSZMT(4,IAT),*) ZZ(IAT)
            call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
            call hnd_dparsr(prszmt(2,iat),zmtchr(2,iat),xx(iat)      )
            call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),yy(iat)      )
            call hnd_dparsr(prszmt(4,iat),zmtchr(4,iat),zz(iat)      )
                 ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
         ELSEIF(NUMZMT(  IAT).EQ.5.AND.FLGZMT(1,IAT).EQ.3.AND.
     1          FLGZMT(2,IAT).EQ.2.AND.FLGZMT(5,IAT).EQ.1     ) THEN
c           READ(PRSZMT(1,IAT),*) ATNAME(IAT)
c           READ(PRSZMT(3,IAT),*) XX(IAT)
c           READ(PRSZMT(4,IAT),*) YY(IAT)
c           READ(PRSZMT(5,IAT),*) ZZ(IAT)
            call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
            call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),xx(iat)      )
            call hnd_dparsr(prszmt(4,iat),zmtchr(4,iat),yy(iat)      )
            call hnd_dparsr(prszmt(5,iat),zmtchr(5,iat),zz(iat)      )
                 ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
         ELSE
            WRITE(IW,9994) IAT
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         CART(IAT) = .TRUE.
         RETURN
      ENDIF
C
C     ----- CONSTRUCT THE -Z- MATRIX -----
C
      IF(IAT.EQ.1.AND.(NUMZMT(IAT).NE.1.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9989)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IAT.EQ.2.AND.(NUMZMT(IAT).NE.3.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9988)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IAT.EQ.3.AND.(NUMZMT(IAT).NE.5.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9987)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- ANCHOR ATOM              -----
C
      ZMT(1,IAT) = IAT
      IF(FLGZMT(1,IAT) .EQ. 3) THEN
c        READ (PRSZMT(1,IAT),*) ATNAME(IAT)
         call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
               ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
      ELSE
         WRITE(IW,9999)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 1) GO TO 100
C
C     ----- BOND LENGTH -----
C
      IF(FLGZMT(2,IAT) .EQ. 2) THEN
c        READ (PRSZMT(2,IAT),*) ZMT(2,IAT)
         call hnd_dparsi(prszmt(2,iat),zmtchr(2,iat),zmt(2,iat)   )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(2,IAT) (1:8) .EQ. PRSZMT(1,JAT) (1:8)) THEN
               ZMT(2,IAT)=JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(3,IAT) .EQ. 1) THEN
c        READ (PRSZMT(3,IAT),*) ZVAL(1,IAT)
         call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),zval(1,iat)  )
                                ZLST(1,IAT)=ZVAL(1,IAT)
      ELSEIF(FLGZMT(3,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(3, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(3, IAT) (1:ZMTCHR(3, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))   ) THEN
               zvarname(nzmod+1) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
c              READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                       ZLST(1,IAT)=ZVAL(1,IAT)
               FRZVAL(1,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(3, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(3, IAT) (3:ZMTCHR(3, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))   ) THEN
               zvarname(nzmod+1) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
               IF(PRSZMT(3, IAT) (2:2).EQ.PLUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                          ZLST(1,IAT)=ZVAL(1,IAT)
                  FRZVAL(1,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(3, IAT) (2:2).EQ.MINUS) THEN
                  zvarsign(nzmod+1) = -1d0
c                 READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                          ZLST(1,IAT)=ZVAL(1,IAT)
                  FRZVAL(1,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
                  ENDIF
                  ZVAL(1,IAT)=-ZVAL(1,IAT)
                  ZLST(1,IAT)=-ZLST(1,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9998) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9998) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 2) GO TO 100
C
C     ----- BOND ANGLE -----
C
      IF(FLGZMT(4,IAT) .EQ. 2) THEN
c        READ (PRSZMT(4,IAT),*) ZMT(3,IAT)
         call hnd_dparsi(prszmt(4,iat),zmtchr(4,iat),zmt(3,iat)  )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(4,IAT) (1:8).EQ.PRSZMT(1,JAT) (1:8)) THEN
               ZMT(3,IAT) = JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(5,IAT) .EQ. 1) THEN
c        READ (PRSZMT(5,IAT),*) ZVAL(2,IAT)
         call hnd_dparsr(prszmt(5,iat),zmtchr(5,iat),zval(2,iat)  )
                                ZLST(2,IAT)=ZVAL(2,IAT)
      ELSEIF(FLGZMT(5,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(5, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(5, IAT) (1:ZMTCHR(5, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))    ) THEN
               zvarname(nzmod+2) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
c              READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                       ZLST(2,IAT)=ZVAL(2,IAT)
               FRZVAL(2,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(5, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(5, IAT) (3:ZMTCHR(5, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))    ) THEN
               zvarname(nzmod+2) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
               IF(PRSZMT(5, IAT) (2:2).EQ.PLUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                          ZLST(2,IAT)=ZVAL(2,IAT)
                  FRZVAL(2,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(5, IAT) (2:2).EQ.MINUS) THEN
                  zvarsign(nzmod+2) = -1d0
c                 READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                          ZLST(2,IAT)=ZVAL(2,IAT)
                  FRZVAL(2,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
                  ENDIF
                  ZVAL(2,IAT)=-ZVAL(2,IAT)
                  ZLST(2,IAT)=-ZLST(2,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9997) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9997) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(ZVAL(2,IAT).LT.ZERO.OR.ZVAL(2,IAT).GT.PIDEG) THEN
         WRITE(IW,9986) ZVAL(2,IAT),IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 3) GO TO 100
C
C     ----- TORSION -----
C
      IF(FLGZMT(6,IAT) .EQ. 2) THEN
c        READ (PRSZMT(6,IAT),*) ZMT(4,IAT)
         call hnd_dparsi(prszmt(6,iat),zmtchr(6,iat),zmt(4,iat)  )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(6,IAT) (1:8).EQ.PRSZMT(1,JAT) (1:8)) THEN
               ZMT(4,IAT)=JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(7,IAT) .EQ. 1) THEN
c        READ (PRSZMT(7,IAT),*) ZVAL(3,IAT)
         call hnd_dparsr(prszmt(7,iat),zmtchr(7,iat),zval(3,iat)  )
                                ZLST(3,IAT)=ZVAL(3,IAT)
      ELSEIF(FLGZMT(7,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(7, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(7, IAT) (1:ZMTCHR(7, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))    ) THEN
c              READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
               zvarname(nzmod+3) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                       ZLST(3,IAT)=ZVAL(3,IAT)
               FRZVAL(3,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(7, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(7, IAT) (3:ZMTCHR(7, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))    ) THEN
               zvarname(nzmod+3) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
               IF(PRSZMT(7, IAT) (2:2).EQ.PLUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                          ZLST(3,IAT)=ZVAL(3,IAT)
                  FRZVAL(3,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(7, IAT) (2:2).EQ.MINUS) THEN
                  zvarsign(nzmod+3) = -1d0
c                 READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                          ZLST(3,IAT)=ZVAL(3,IAT)
                  FRZVAL(3,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
                  ENDIF
                  ZVAL(3,IAT)=-ZVAL(3,IAT)
                  ZLST(3,IAT)=-ZLST(3,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9996) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9996) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- LAST INTEGER PARAMETER -----
C
      IF(NUMZMT(IAT).EQ.8) THEN
         IF(FLGZMT(8,IAT).EQ.2) THEN
c           READ(PRSZMT(8,IAT),*) ZMT(5,IAT)
         call hnd_dparsi(prszmt(8,iat),zmtchr(8,iat),zmt(5,iat)  )
         ELSE
            WRITE(IW,9993) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
C     ----- POSSIBLE ERRORS -----
C
      IF(ZMT(5,IAT).EQ.0) THEN
         IF(ZVAL(3,IAT).LT.-TWO*PIDEG.OR.
     1      ZVAL(3,IAT).GT. TWO*PIDEG) THEN
            WRITE(IW,9985) ZVAL(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ZVAL(3,IAT).LT.-PIDEG) ZVAL(3,IAT)=ZVAL(3,IAT)+TWO*PIDEG
         IF(ZVAL(3,IAT).GT. PIDEG) ZVAL(3,IAT)=ZVAL(3,IAT)-TWO*PIDEG
C
         IF(ZLST(3,IAT).LT.-TWO*PIDEG.OR.
     1      ZLST(3,IAT).GT. TWO*PIDEG) THEN
            WRITE(IW,9985) ZLST(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZLST(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ZLST(3,IAT).LT.-PIDEG) ZLST(3,IAT)=ZLST(3,IAT)+TWO*PIDEG
         IF(ZLST(3,IAT).GT. PIDEG) ZLST(3,IAT)=ZLST(3,IAT)-TWO*PIDEG
      ELSE
         IF(ZVAL(3,IAT).LT. ZERO .OR.
     1      ZVAL(3,IAT).GT.PIDEG) THEN
            WRITE(IW,9985) ZVAL(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
C
         IF(ZLST(3,IAT).LT. ZERO .OR.
     1      ZLST(3,IAT).GT.PIDEG) THEN
            WRITE(IW,9985) ZLST(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZLST(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
  100 CONTINUE
      IF(DBUG) THEN
         WRITE(IW,9984) IAT,ATNAME(IAT)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(  ZVAL(I,IAT),I=1,3),
     1                                         (FRZVAL(I,IAT),I=1,3)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(  ZLST(I,IAT),I=1,3),
     1                                         (FRZVAL(I,IAT),I=1,3)
      ENDIF
C
C     ----- CREATE -IZMAT- FOR -HONDO- INPUT -----
C
      IF(IAT.GT.1) THEN
         NZMOD=NZMOD+1
         IZ=IZ+1
         IZMAT(IZ)=1
         IZ=IZ+1
         IZMAT(IZ)=ZMT(1,IAT)
         IZ=IZ+1
         IZMAT(IZ)=ZMT(2,IAT)
         IZFRZ(NZMOD)=FRZVAL(1,IAT)
C
         IF(IAT.GT.2) THEN
            NZMOD=NZMOD+1
            IZ=IZ+1
            IZMAT(IZ)=2
            IZ=IZ+1
            IZMAT(IZ)=ZMT(1,IAT)
            IZ=IZ+1
            IZMAT(IZ)=ZMT(2,IAT)
            IZ=IZ+1
            IZMAT(IZ)=ZMT(3,IAT)
            IZFRZ(NZMOD)=FRZVAL(2,IAT)
C
            IF(IAT.GT.3) THEN
               IF(ZMT(5,IAT).EQ.0) THEN
                  NZMOD=NZMOD+1
                  IZ=IZ+1
                  IZMAT(IZ)=3
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(1,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(2,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(3,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(4,IAT)
                  IZFRZ(NZMOD)=FRZVAL(3,IAT)
               ELSE
                  NZMOD=NZMOD+1
                  IZ=IZ+1
                  IZMAT(IZ)=7   ! RJH was 2 but discarded zmt(5,iat)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(1,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(2,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(4,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(5,IAT)
                  IZFRZ(NZMOD)=FRZVAL(3,IAT)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(IZ.GT.MXIZMT) THEN
         WRITE(IW,9991) IZ,MXIZMT,IAT
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      RETURN
 9999 FORMAT(' FIRST PARSED -Z- DATA IS NOT A CHARACTER STRING. STOP')
 9998 FORMAT(' THE 3-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9997 FORMAT(' THE 5-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9996 FORMAT(' THE 7-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9995 FORMAT(' IAT=',I5,' ZMAT=',5I3,3F12.5,3(2X,L4))
 9994 FORMAT(' SOMETHING IS WRONG WITH THE CARTESIAN COORDINATES INPUT',
     1       ' DATA FOR THIS ATOM, IAT = ',I3,' . STOP')
 9993 FORMAT(' THE 8-TH (INTEGER) PARAMETER FOR ATOM = ',I4,
     1       ' IS OF THE WRONG TYPE. STOP')
 9992 FORMAT(' ATOM REFERED TO FOR IAT = ',I3,' IS INVALID. STOP')
 9991 FORMAT(' TOO MANY -IZMAT- ENTRIES, IZ,MXIZMT = ',2I5,
     1       ' FOR -IAT- = ',I5,' STOP.')
 9989 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -1- IS INCORRECT. STOP')
 9988 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -2- IS INCORRECT. STOP')
 9987 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -3- IS INCORRECT. STOP')
 9986 FORMAT(' BOND ANGLE = ',F10.2,' FOR ATOM = ',I4,' IS OUT OF',
     1       ' ALLOWED RANGE ... 0.0 TO 180.0 . STOP')
 9985 FORMAT(' TORSION ANGLE OR BOND ANGLE = ',F10.2,' FOR ATOM = ',I4,
     1       ' IS OUT OF ALLOWED RANGE ... -360.0 TO 360.0 . STOP')
 9984 FORMAT(' IAT=',I5,' ATNAME=',1X,A8)
      END
      SUBROUTINE HND_ZXYZ(NAT,NZMT,ZVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ROUTINE COMPUTES :
C
C     THE CARTESIAN COORDINATES OF A MOLECULE FROM BONDING COORDINATES
C     (BOND LENGTHS, BOND AND DIHEDRAL ANGLES).
C
C     ATOMS NA, NB, NC, HAVE KNOWN COORDINATES AND ARE NOT COLLINEAR.
C     RCD   =ZVAL(1,I) : C-D BOND LENGTH (IN ANGSTROM)
C     THBCD =ZVAL(2,I) : B-C-D BOND ANGLE (IN DEGREES)
C     PHABCD=ZVAL(3,I) : DIHEDRAL ANGLE (IN DEGREES) OF C-D
C                        RELATIVE TO A-B MEASURED CLOCKWISE
C                        ( + ) ALONG THE DIRECTION B TO C
C
C     MEANING OF PARAMETER NZMT(5,I):
C          =0     ANGLES FROM DATA
C          =+/- 1 ZVAL(3,I) IS BOND ANGLE A-C-D
C          =2     BOND ANGLE B-C-D IS TETRAHEDRAL
C          =3     ATOMS B-C-D ARE COLLINEAR
C
      PARAMETER     (MXATOM=500)
      PARAMETER     (MAXGEO=MXATOM+1,MAXWRD=40,MAXVAR=256)
      LOGICAL       DBUG
      LOGICAL       CART
      CHARACTER*8   ATNAME
      CHARACTER*8   ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_XYZGEO/XX(MAXGEO),YY(MAXGEO),ZZ(MAXGEO),
     1                  ATNAME(MAXGEO),ATNUM(MAXGEO),CART(MAXGEO)
      DIMENSION   NZMT(5,MAXGEO)
      DIMENSION   ZVAL(3,MAXGEO)
      DIMENSION   ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','- ZXYZ -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TWO,THREE /2.0D+00,3.0D+00/
      DATA PT1,MAX   /0.1D+00,6/
C
      DBUG=.FALSE.
C
      PIFAC=3.1415926536D+00/180.0D+00
      COSTD=-ONE/THREE
      SINTD= TWO/THREE* SQRT(TWO)
      NA   =0
      NB   =0
      NC   =0
      ND   =0
      RCD  =ZERO
C
C     ----- ATOM -1- AND -2- -----
C
      IAT=1
      IF(CART(IAT)) THEN
         IF(DBUG) THEN
            WRITE(IW,9997) IAT
         ENDIF
      ELSE
         ND=NZMT(1,IAT)
         XX(ND)=ZERO
         YY(ND)=ZERO
         ZZ(ND)=ZERO
         CART(IAT)=.TRUE.
         IF(DBUG) THEN
            WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
         ENDIF
      ENDIF
C
      IF(NAT.GT.1) THEN
         IAT=2
         IF(CART(IAT)) THEN
            IF(DBUG) THEN
               WRITE(IW,9997) IAT
            ENDIF
         ELSE
            ND=NZMT(1,IAT)
            NC=NZMT(2,IAT)
            XX(ND)=ZERO
            YY(ND)=ZERO
            ZZ(ND)=ZVAL(1,IAT)
            RCD   =ZVAL(1,IAT)
            CART(IAT)=.TRUE.
            IF(DBUG) THEN
               WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
            ENDIF
         ENDIF
C
         IF(NAT.GT.2) THEN
C
C     ----- ATOM -3- -----
C
            IAT=3
            IF(CART(IAT)) THEN
               IF(DBUG) THEN
                  WRITE(IW,9997) IAT
               ENDIF
            ELSE
                   IF(NZMT(5,IAT).EQ.3) THEN
                  CCOS=-ONE
                  SSIN= ZERO
               ELSEIF(NZMT(5,IAT).EQ.2) THEN
                  CCOS=COSTD
                  SSIN=SINTD
               ELSEIF(NZMT(5,IAT).EQ.0) THEN
                  THBCD=ZVAL(2,IAT)*PIFAC
                  CCOS= COS(THBCD)
                  SSIN= SIN(THBCD)
               ELSE
                  WRITE(IW,9999) IAT
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               ND =NZMT(1,IAT)
               NC =NZMT(2,IAT)
               NB =NZMT(3,IAT)
               RCD=ZVAL(1,IAT)
               IF(DBUG) THEN
                  WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
               ENDIF
               XX(ND)=         RCD*SSIN
               YY(ND)=ZERO
               ZZ(ND)=ZZ(NC) - RCD*CCOS
               IF(ZZ(NC).LT.ZZ(NB)) THEN
                  ZZ(ND)=ZZ(NC) + RCD*CCOS
               ENDIF
               CART(IAT) = .TRUE.
            ENDIF
C
C     ----- ATOM -4- AND HIGHER -----
C
            IF(NAT.GT.3) THEN
               DO 20 IAT=4,NAT
                  IF(CART(IAT)) THEN
                     IF(DBUG) THEN
                        WRITE(IW,9997) IAT
                     ENDIF
                     GO TO 20
                  ENDIF
                  NA =NZMT(4,IAT)
                  NB =NZMT(3,IAT)
                  NC =NZMT(2,IAT)
                  ND =NZMT(1,IAT)
                  RCD=ZVAL(1,IAT)
                  IF(DBUG) THEN
                     WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
                  ENDIF
C
C     ----- LINEAR CONNECTION FOR THIS ATOM : D...C...B -----
C
                  IF(NZMT(5,IAT).EQ.3) THEN
                     RCB= SQRT((XX(NB)-XX(NC))**2+
     1                         (YY(NB)-YY(NC))**2+(ZZ(NB)-ZZ(NC))**2)
                     XX(ND)=XX(NC)+RCD*(XX(NC)-XX(NB))/RCB
                     YY(ND)=YY(NC)+RCD*(YY(NC)-YY(NB))/RCB
                     ZZ(ND)=ZZ(NC)+RCD*(ZZ(NC)-ZZ(NB))/RCB
                  ELSE
C
C     ----- OTHER CONNECTIONS FOR THIS ATOM : DEFINE LOCAL FRAME -----
C
                     XB=XX(NB)-XX(NC)
                     YB=YY(NB)-YY(NC)
                     ZB=ZZ(NB)-ZZ(NC)
                     RCB= SQRT(XB*XB+YB*YB+ZB*ZB)
                     XB=XB/RCB
                     YB=YB/RCB
                     ZB=ZB/RCB
                     T11=XB
                     T21=YB
                     T31=ZB
                     XA=XX(NA)-XX(NC)
                     YA=YY(NA)-YY(NC)
                     ZA=ZZ(NA)-ZZ(NC)
                     RCA= SQRT(XA*XA+YA*YA+ZA*ZA)
                     XA=XA/RCA
                     YA=YA/RCA
                     ZA=ZA/RCA
                     DOT=XA*XB+YA*YB+ZA*ZB
                     XA=XA-DOT*XB
                     YA=YA-DOT*YB
                     ZA=ZA-DOT*ZB
                     DUM= SQRT(XA*XA+YA*YA+ZA*ZA)
                     XA=XA/DUM
                     YA=YA/DUM
                     ZA=ZA/DUM
                     T12=XA
                     T22=YA
                     T32=ZA
                     T13=T21*T32-T31*T22
                     T23=T31*T12-T11*T32
                     T33=T11*T22-T21*T12
C
C     ----- GET CARTESIAN COORDINATES IN LOCAL  FRAME -----
C
                         IF(NZMT(5,IAT).EQ.0) THEN
                        PHI  = ZVAL(2,IAT)*PIFAC
                        THETA=-ZVAL(3,IAT)*PIFAC
                        CCOSP= COS(PHI)
                        SSINP= SIN(PHI)
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSEIF(NZMT(5,IAT).EQ.-1.OR.
     1                      NZMT(5,IAT).EQ. 1    ) THEN
                        RAB  =SQRT((XX(NB)-XX(NA))**2+(YY(NB)-YY(NA))**2
     1                            +(ZZ(NB)-ZZ(NA))**2)
                        PHI=ZVAL(2,IAT)*PIFAC
                        BET=ZVAL(3,IAT)*PIFAC
                        ALP=PHI
                        GAM  =ACOS((RCB**2+RCA**2-RAB**2)/(TWO*RCB*RCA))
                        THETA=ACOS(( COS(BET)- COS(ALP)* COS(GAM))/
     1                             (           SIN(ALP)* SIN(GAM))  )
                        IF(NZMT(5,IAT).EQ.-1) THEN
                           THETA=-THETA
                        ENDIF
                        CCOSP= COS(PHI)
                        SSINP= SIN(PHI)
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSEIF(NZMT(5,IAT).EQ.2) THEN
                        CCOSP= COSTD
                        SSINP= SINTD
                        THETA=ZVAL(3,IAT)*PIFAC
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSE
                        WRITE(IW,9999) IAT
                        CALL HND_HNDERR(3,ERRMSG)
                     ENDIF
                     XXD=RCD*CCOSP
                     YYD=RCD*SSINP*CCOST
                     ZZD=RCD*SSINP*SSINT
C
C     ----- GET CARTESIAN COORDINATES IN MASTER FRAME -----
C
                     XX(ND)=T11*XXD+T12*YYD+T13*ZZD + XX(NC)
                     YY(ND)=T21*XXD+T22*YYD+T23*ZZD + YY(NC)
                     ZZ(ND)=T31*XXD+T32*YYD+T33*ZZD + ZZ(NC)
                  ENDIF
                  CART(IAT)=.TRUE.
   20             CONTINUE
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
 9999 FORMAT(' ILLEGAL 5-TH -ZMT- PARAMETER FOR -IAT- = ',I3,' STOP.')
 9998 FORMAT(' IAT,NA,NB,NC,ND,RCD = ',5I4,F12.5)
 9997 FORMAT(' IAT = ',I4,' IS ALREADY SPECIFIED IN CARTESIAN SPACE.')
      END
      SUBROUTINE HND_PRSQ(V,M,N,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_LISTNG/LIST
      DIMENSION V(NDIM,1)
      DIMENSION IC(5),C(5)
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      IF(LIST.EQ.0) WRITE (IW,9048) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      WRITE(IW,9328) J
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
 9328 FORMAT(/,' VEC =',I5,/,1X,10(1H-))
 9348 FORMAT(5(I5,F11.5))
      END
      SUBROUTINE HND_PREV(V,E,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- PRINT OUT E AND V-MATRICES
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_LISTNG/LIST
      DIMENSION V(NDIM,1),E(1)
      DIMENSION IC(5),C(5)
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9068) (E(I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9168) (E(I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9268) (E(I),I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      IF(LIST.EQ.0) WRITE (IW,9048) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      WRITE(IW,9328) J,E(J)
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(15X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9068 FORMAT(15X,10F11.5)
 9128 FORMAT(15X,7(6X,I3,6X))
 9148 FORMAT(10X,I5,7F15.10)
 9168 FORMAT(15X,7F15.10)
 9228 FORMAT(15X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
 9268 FORMAT(15X,7E15.8)
 9328 FORMAT(/,' VEC =',I5,' EIG/OCC =',F11.5,/,1X,31(1H-))
 9348 FORMAT(5(I5,F11.5))
      END
      SUBROUTINE HND_PRTR(D,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- PRINT OUT A TRIANGULAR MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW
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
      II = MAX0( I, J)
      JJ = MIN0( I, J)
      IJ = (II*(II-1))/2 + JJ
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
      SUBROUTINE HND_TFTR(H,F,Q,T,IA,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- H(M,M) = Q(DAGGER)(N,M) * F(N,N) * Q(N,M) -----
C
      DIMENSION H(1),F(1),Q(NDIM,1),T(1),IA(1)
      DATA SMALL /1.0D-11/
      DATA ZERO  /0.0D+00/
      IJ = 0
      DO 180 J = 1,M
      IK = 0
      DO 140 I = 1,N
      DUM = ZERO
      QIJ = Q(I,J)
      MAX = I-1
      IF (MAX .EQ. 0) GO TO 120
      DO 100 K = 1,MAX
      IK = IK+1
      T(K) = T(K)+ F(IK)*QIJ
      DUM  = DUM + F(IK)*Q(K,J)
  100 CONTINUE
  120 IK = IK+1
      T(I) = DUM + F(IK)*QIJ
  140 CONTINUE
      DO 160 I = 1,J
      IJ = IJ+1
      HIJ = ddot(N,Q(1,I),1,T,1)
      IF( ABS(HIJ).LT.SMALL) HIJ=ZERO
      H(IJ)=HIJ
  160 CONTINUE
  180 CONTINUE
      RETURN
      END
      SUBROUTINE HND_DIAGIV(A,VEC,EIG,IA,NVEC,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- GENERAL CALLING  ROUTINE FOR -----
C               -GIVEIS- OR -GIVDIA-
C
#include "mafdecls.fh"
#include "global.fh"
      logical status
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DIAOPT/MODDIA
      DIMENSION A(1),VEC(NDIM,1),EIG(1),IA(1)
      dimension x(1),ix(1)
      equivalence (dbl_mb(1),x(1)),(int_mb(1),ix(1))
c
c     ----- memory -----
c
      status=ma_push_get(MT_DBL,4*ndim,'mem_dbl dia',i_idbl,idbl)
      status=ma_push_get(MT_INT,3*ndim,'mem_int dia',i_iint,iint)
c
c     ----- call -----
c
      CALL HND_GIVDIA(A,VEC,EIG,IA,N,NDIM,
     1     x(idbl), x(idbl+  ndim), x(idbl+2*ndim), x(idbl+3*ndim),
     2                              x(idbl+2*ndim), x(idbl+2*ndim),
     3    ix(iint),ix(iint+  ndim),ix(iint+2*ndim)                )
c
c     ----- memory -----
c
      status=ma_pop_stack(i_iint)
      if(.not.status) then
         call errquit('hnd_diagiv,ma_pop_stack of -iint- failed',911)
      endif
      status=ma_pop_stack(i_idbl)
      if(.not.status) then
         call errquit('hnd_diagiv,ma_pop_stack of -idbl- failed',911)
      endif
c
      RETURN
      END
      SUBROUTINE HND_GIVDIA(A,VEC,EIG,IA,N,NDIM,
     1           w,gamma,beta,betasq,p,q,iposv,ivpos,iord)
      IMPLICIT REAL*8 (A-H,O-Z)
C        -----  A GIVENS HOUSHOLDER MATRIX DIAGONALIZATION  -----
C        -----  ROUTINE SAME AS EIGEN BUT WORKS WITH A      -----
C        -----  LINEAR ARRAY.                               -----
      DIMENSION A(1),VEC(NDIM,1),EIG(1),IA(1)
      DIMENSION W(*),GAMMA(*),BETA(*),BETASQ(*)
      DIMENSION P(*),Q(*),IPOSV(*),IVPOS(*),IORD(*)
      DATA ZERO,PT5,ONE,TWO /0.0D+00,0.5D+00,1.0D+00,2.0D+00/
      DATA RHOSQ /1.0D-22/
C
      IF(N.EQ.0) GO TO 560
      N1=N-1
      N2=N-2
      GAMMA(1)=A(1)
      IF(N2) 280,270,120
  120 DO 260 NR=1,N2
      IK=IA(NR+1)+NR
      B=A(IK)
      S=ZERO
      DO 130 I=NR,N2
      IJ=IA(I+2)+NR
  130 S=S+A(IJ)**2
C
C        -----  PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION ----
C
      A(IK)=ZERO
      IF(S.LE.ZERO) GO TO 250
      S=S+B*B
      SGN=+ONE
      IF(B.GE.ZERO) GO TO 160
      SGN=-ONE
  160 SQRTS= SQRT(S)
      D=SGN/(SQRTS+SQRTS)
      TEMP= SQRT(PT5+B*D)
      W(NR)=TEMP
      A(IK)=TEMP
      D=D/TEMP
      B=-SGN*SQRTS
C
C        -----  -D- IS FACTOR OF PROPORTIONALITY. NOW       -----
C        -----  COMPUTE AND SAVE -W- VECTOR. EXTRA SINGLY   -----
C        -----  SUBSCRIPTED -W- VECTOR FOR SPEED.           -----
C
      DO 170 I=NR,N2
      IJ=IA(I+2)+NR
      TEMP=D*A(IJ)
      W(I+1)=TEMP
  170 A(IJ)=TEMP
C
C        -----  PREMULTIPLY VECTOR -W- BY MATRIX -A- TO     -----
C        -----  OBTAIN -P- VECTOR. SIMULTANEOUSLY ACCUMULATE ----
C        -----  DOT PRODUCT -WP- -- SCALR -K-.              -----
C
      WTAW=ZERO
      DO 220 I=NR,N1
      SUM=ZERO
      II=IA(I+1)
      DO 180 J=NR,I
      IJ=II+J+1
  180 SUM=SUM+A(IJ)*W(J)
      I1=I+1
      IF(N1.LT.I1) GO TO 210
      DO 200 J=I1,N1
      IJ=IA(J+1)+I+1
  200 SUM=SUM+A(IJ)*W(J)
  210 P(I)=SUM
  220 WTAW=WTAW+SUM*W(I)
      DO 230 I=NR,N1
  230 Q(I)=P(I)-WTAW*W(I)
C
C        -----  NOW FORM -PAP- MATRIX, REQUIRED PART        -----
C
      DO 240 J=NR,N1
      QJ=Q(J)
      WJ=W(J)
      JJ=J+1
      DO 240 I=J,N1
      IJ=IA(I+1)+JJ
  240 A(IJ)=A(IJ)-TWO*(W(I)*QJ+WJ*Q(I))
  250 BETA(NR)=B
      BETASQ(NR)=B*B
      IL=IK+1
  260 GAMMA(NR+1)=A(IL)
  270 IJ=IA(N)+N-1
      B=A(IJ)
      BETA(N-1)=B
      BETASQ(N-1)=B*B
      IJ=IJ+1
      GAMMA(N)=A(IJ)
  280 BETASQ(N)=ZERO
C
C        -----  ADJOIN AN IDENTYTY MATRIX TO BE POST-       -----
C        -----  MULTIPLIED BY ROTATIONS                     -----
C
      DO 300 I=1,N
      DO 299 J=1,N
  299 VEC(I,J)=ZERO
  300 VEC(I,I)=ONE
      M=N
      SUM=ZERO
      NPAS=1
      GO TO 400
  310 SUM=SUM+SHIFT
      COSA=ONE
      G=GAMMA(1)-SHIFT
      PP=G
      PPBS=PP*PP+BETASQ(1)
      PPBR= SQRT(PPBS)
      DO 370 J=1,M
      COSAP=COSA
      IF(PPBS.NE.ZERO) GO TO 320
      SINA=ZERO
      SINA2=ZERO
      COSA=ONE
      GO TO 350
  320 SINA=BETA(J)/PPBR
      SINA2=BETASQ(J)/PPBS
      COSA=PP/PPBR
C
C        -----  POSTMULTIPLY IDENTITY BY -P- TRANSPOSE .    -----
C
      NT=J+NPAS
      IF(NT.LT.N) GO TO 330
      NT=N
  330 CONTINUE
      DO 340 I=1,NT
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)
  340 VEC(I,J)=TEMP
  350 DIA=GAMMA(J+1)-SHIFT
      U=SINA2*(G+DIA)
      GAMMA(J)=G+U
      G=DIA-U
      PP=DIA*COSA-SINA*COSAP*BETA(J)
      IF(J.NE.M) GO TO 360
      BETA(J)=SINA*PP
      BETASQ(J)=SINA2*PP*PP
      GO TO 380
  360 PPBS=PP*PP+BETASQ(J+1)
      PPBR= SQRT(PPBS)
      BETA(J)=SINA*PPBR
  370 BETASQ(J)=SINA2*PPBS
  380 GAMMA(M+1)=G
C
C        -----  TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT ----
C
      NPAS=NPAS+1
      IF(BETASQ(M).GT.RHOSQ) GO TO 410
  390 EIG(M+1)=GAMMA(M+1)+SUM
  400 BETA(M)=ZERO
      BETASQ(M)=ZERO
      M=M-1
      IF(M.EQ.0) GO TO 430
      IF(BETASQ(M).LE.RHOSQ) GO TO 390
C
C        -----  TAKE ROOT OF CORMER 2 BY 2 NEAREST TO       -----
C        -----  LOWER DIAGONAL IN VALUE AS ESTIMATE OF      -----
C        -----  EIGENVALUE TO USE FOR SHIFT                 -----
C
  410 A2=GAMMA(M+1)
      R2=PT5*A2
      R1=PT5*GAMMA(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP= SQRT(DIF*DIF+BETASQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF=  ABS(A2-R1)-  ABS(A2-R2)
      IF(DIF.LT.ZERO) GO TO 420
      SHIFT=R2
      GO TO 310
  420 SHIFT=R1
      GO TO 310
  430 EIG(1)=GAMMA(1)+SUM
      DO 440 J=1,N
      IPOSV(J)=J
      IVPOS(J)=J
  440  IORD(J)=J
      M=N
      GO TO 470
  450 DO 460 J=1,M
      IF(EIG(J).LE.EIG(J+1)) GO TO 460
      TEMP=EIG(J)
      EIG(J)=EIG(J+1)
      EIG(J+1)=TEMP
      ITEMP=IORD(J)
      IORD(J)=IORD(J+1)
      IORD(J+1)=ITEMP
  460 CONTINUE
  470 M=M-1
      IF(M.NE.0) GO TO 450
      IF(N1.EQ.0) GO TO 500
      DO 490 L=1,N1
      NV=IORD(L)
      NP=IPOSV(NV)
      IF(NP.EQ.L) GO TO 490
      LV=IVPOS(L)
      IVPOS(NP)=LV
      IPOSV(LV)=NP
      DO 480 I=1,N
      TEMP=VEC(I,L)
      VEC(I,L)=VEC(I,NP)
  480 VEC(I,NP)=TEMP
  490 CONTINUE
C
C        -----  BACK TRANSFORM THE VECTORS OF THE TRIPLE    -----
C        -----  DIAGONAL MATRIX.                            -----
C
  500 DO 550 NRR=1,N
      K=N1
  510 K=K-1
      IF(K.LE.0) GO TO 550
      SUM=ZERO
      DO 520 I=K,N1
      IJ=IA(I+1)+K
  520 SUM=SUM+VEC(I+1,NRR)*A(IJ)
      SUM=SUM+SUM
      DO 530 I=K,N1
      IJ=IA(I+1)+K
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(IJ)
      GO TO 510
  550 CONTINUE
  560 CONTINUE
      RETURN
      END
