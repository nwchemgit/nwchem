      SUBROUTINE HND_GIAOX(RTDB,BASIS,GEOM)
* $Id: prp.F,v 1.38 1999-04-10 00:15:45 windus Exp $
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "schwarz.fh"
#include "cfock.fh"
      double precision jfac(1),kfac(1)
      integer     ga_create_atom_blocked
      external    ga_create_atom_blocked
      logical     movecs_read
      external    movecs_read
      logical     oskel
      integer     g_d_1     
      integer     g_h_1     
      integer     g_f_1(1)
      INTEGER     RTDB
      INTEGER     BASIS	
      INTEGER     GEOM
      LOGICAL     STATUS
      LOGICAL     PROP_GIAO_PUT_INT
      EXTERNAL    PROP_GIAO_PUT_INT
      LOGICAL     PROP_GIAO_GET_INT
      EXTERNAL    PROP_GIAO_GET_INT
      LOGICAL     PROP_GIAO_PUT_D
      EXTERNAL    PROP_GIAO_PUT_D
      LOGICAL     PROP_GIAO_GET_D
      EXTERNAL    PROP_GIAO_GET_D
      LOGICAL     PROP_GIAO_PUT_U
      EXTERNAL    PROP_GIAO_PUT_U
      LOGICAL     PROP_GIAO_GET_U
      EXTERNAL    PROP_GIAO_GET_U
      LOGICAL     PROP_GIAO_PUT_D1
      EXTERNAL    PROP_GIAO_PUT_D1
      LOGICAL     PROP_GIAO_GET_D1
      EXTERNAL    PROP_GIAO_GET_D1
      LOGICAL     PROP_GIAO_PUT_U1
      EXTERNAL    PROP_GIAO_PUT_U1
      LOGICAL     PROP_GIAO_GET_U1
      EXTERNAL    PROP_GIAO_GET_U1
      PARAMETER   (MXIODA=255)
      PARAMETER   (MXATOM=500)
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP
      CHARACTER*8 SCFTYP
      CHARACTER*8 MEMTYP
      CHARACTER*8 HNDTYP
      CHARACTER*8 DSKTYP
      CHARACTER*8 FILTYP
      CHARACTER*8 PCKTYP
      CHARACTER*8 STNDRD
      CHARACTER*8 DIRECT
      CHARACTER*8 SEMID 
      CHARACTER*8 SCF   
      CHARACTER*8 RHF   
      LOGICAL     SOME
      LOGICAL     OUT
      LOGICAL     DBUG    
      LOGICAL     DBUG0
      LOGICAL     DISK
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_SQFILE/IJK,IPK
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
C
      COMMON/HND_CVGGIA/CNVA,CNVB,DAMP
      COMMON/HND_PARGIA/VNEW,VOLD,MAXIT
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
C
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
C
      COMMON/HND_QMTOPT/TOLQMT,NQMT
      COMMON/HND_INTOPT/NOPK,NOK,NOSQUR
      COMMON/HND_HNDOPT/HNDTYP,DSKTYP,FILTYP,PCKTYP
      COMMON/HND_INTFIL/NINTMX
      COMMON/HND_PCKLAB/LABSIZ
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      COMMON/HND_SCFMEM/IFCK,IDUM,IDEN,JDUM(3),IBUF,KDUM(2),NEED
      COMMON/HND_DIRCAS/I2CASE
      COMMON/HND_DIRMEM/LDUM(2),JDEN,JFCK,MDUM(2),NDIM
      COMMON/HND_DIRGIA/IFCK0,IFCKX,IFCKY,IFCKZ,IDENS,IDIM
C
      DIMENSION    X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION    B20(6)
      DIMENSION    CM(3)
      DIMENSION    ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','- GIAOX-'/
      DATA MEMTYP  /'IN CORE '/
      DATA STNDRD  /'STANDARD'/
      DATA DIRECT  /'DIRECT  '/
      DATA SEMID   /'SEMI-DIR'/
      DATA SCF     /'SCF     '/
      DATA RHF     /'RHF     '/
      DATA ZERO    /0.0D+00/
      DATA ONE     /1.0D+00/
      DATA TWO     /2.0D+00/
      DATA CVEL    /137.0359895D+00/
      DATA PPM     /26.62566914D+00/
C
      DBUG0=.FALSE. 
      DBUG =.FALSE.
      OUT  =.FALSE.
      OUT  =OUT.OR.DBUG
      SOME =.TRUE.
      SOME =SOME.OR.OUT
C
      DISK  = DSKTYP.NE.MEMTYP
C
      if(ga_nodeid().eq.0) then
      IF(SOME) THEN
         WRITE(IW,9999)
      ENDIF
      CALL HND_FIFLSH(IW)
      endif
C
      WFNTYP=SCF        
      SCFTYP=RHF
      NSETO =0
      IF(WFNTYP.NE.SCF) THEN 
         IF(SCFTYP.NE.RHF) THEN
            IF(NSETO.NE.0) THEN
               WRITE(IW,9998)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ENDIF
      ENDIF
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM*NUM
c---
      l0 =num
      noc=num
      nintmx=0
      labsiz=1
c---
C
      IFT1=11
      IFT2=12
      IFT3=13
      IFT4=14
      CALL HND_REWFIL(IFT1)
      CALL HND_REWFIL(IFT2)
      CALL HND_REWFIL(IFT3)
      CALL HND_REWFIL(IFT4)
C
C     ----- INPUT -----
C
      CALL HND_GIAINP
      CALL HND_FIFLSH(IW)
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      NEED =L3+L1+NAT*3+NAT*9+NAT*9+NAT*9
      NEEDI=L2+L2+L2+L2+L3+L3+L3+L3+L3+L1+L1
      NEEDJ=L3+L3+L3+L3+L3+L3+L3+L2
      NEED =NEED+MAX(NEEDI,NEEDJ)+L3+L1+L3+3*NOC
      if(dbug) then
         write(iw,*) 'in hnd_giaox, num = ',num
         write(iw,*) '              nat = ',nat
         write(iw,*) '               l1 = ',l1 
         write(iw,*) '               l2 = ',l2 
         write(iw,*) '               l3 = ',l3 
         write(iw,*) '            needi = ',needi
         write(iw,*) '            needj = ',needj
         write(iw,*) '            need  = ',need 
         write(iw,*) '           nintmx = ',nintmx
         write(iw,*) '           labsiz = ',labsiz
      endif
      STATUS=MA_PUSH_GET(MT_DBL,NEED,'MEM INIT',I_INIT,INIT)
      IF(.NOT.STATUS) THEN
         CALL ERRQUIT('HND_GIAO, MALLOC OF INIT FAILED',911)
      ENDIF
C
C     ----- SET POINTERS FOR PARTITIONING OF CORE -----
C
      I00=INIT
      I01=I00+L3
      I02=I01+L1
      I03=I02+NAT*3
      I04=I03+NAT*9
      I05=I04+NAT*9
      I06=I05+NAT*9
C
C     ----- MORE MEMORY -----
C
      I10=I06
      I20=I10+L2
      I30=I20+L2
      I40=I30+L2
      I50=I40+L2
      I60=I50+L3
      I70=I60+L3
      I80=I70+L3
      I90=I80+L3
      I100=I90+L3
      I110=I100+L1
      I120=I110+L1
      I130=I110+NINTMX*(1+LABSIZ)
      LASTI=I130-1
C
      J10=I10
      J11=J10+L3
      J12=J11+L3
      J20=J12+L3
      J21=J20+L3
      J22=J21+L3
      J30=J22+L3
      J40=J30+L3
      J50=J40+L2
      J60=J50+MAX(L1,NINTMX*(1+LABSIZ))
      LASTJ=J60-1
C
      LAST =MAX(LASTI,LASTJ)
      NEXT =LAST+1
      NEED=L3+l1+L3+3*NOC+(LAST-INIT)
      IF(DBUG) THEN
         WRITE(IW,9994) NEED
      ENDIF
C
C     ----- GET SCF ORBITALS -----
C
      CALL HND_DAREAD(IDAF,IODA,X(I00),L3,15)
      CALL HND_DAREAD(IDAF,IODA,X(I01),L1,17)
      
      call hnd_giao_get_vecs(rtdb,basis,geom,
     1                       x(i50),x(i60),x(i70),
     2                       x(i00),x(i01),
     3                       scftyp,
     4                       nclosed ,nopen ,nvirt ,
     5                       nclosedu,nopenu,nvirtu)

      na  =nclosed
      nqmt=nclosed+nopen+nvirt
      l0  = nqmt
      m1  = nqmt
      m2  =(nqmt*(nqmt+1))/2
      m3  = nqmt* nqmt
      noc = na

c--
c--   do i=1,l1
c--      x(i+(2-1)*l1+i00-1)=-x(i+(2-1)*l1+i00-1)
c--      x(i+(5-1)*l1+i00-1)=-x(i+(5-1)*l1+i00-1)
c--   enddo
c--
      IF(OUT) THEN
         WRITE(IW,*) 'C(0)'
         CALL HND_PREVL(X(I00),X(I01),L0,L1,L1)
      ENDIF
C
c--
      out=.false.
c--
      IF(OUT) THEN
C
C     ----- THE -DPL- CPHF EQUATIONS AS A CHECK -----
C
      WRITE(IW,*) '-DPL- CPHF EQUATIONS AS A CHECK'
      CALL HND_FIFLSH(IW)
C
      DO I=1,3
         CM(I)=ZERO
      ENDDO
      CALL HND_DIPINT(X(I10),X(I20),X(I30),CM,l1)
      CALL HND_DAWRIT(IDAF,IODA,X(I10),L2,33,NAV)
      CALL HND_DAWRIT(IDAF,IODA,X(I20),L2,34,NAV)
      CALL HND_DAWRIT(IDAF,IODA,X(I30),L2,35,NAV)
      IFLD=1
      STATUS=          PROP_GIAO_PUT_INT(X(I10),L1,L2,IFLD)
      IFLD=2
      STATUS=STATUS.OR.PROP_GIAO_PUT_INT(X(I20),L1,L2,IFLD)
      IFLD=3
      STATUS=STATUS.OR.PROP_GIAO_PUT_INT(X(I30),L1,L2,IFLD)
      IF(.NOT.STATUS) THEN
         CALL ERRQUIT('HND_GIAO, PROP_GIAO_PUT_INT FAILED',911)
      ENDIF
C
      MFLD=1
      NFLD=1
      DO IFLD=MFLD,NFLD
         IF(OUT) THEN
            WRITE(IW,9996)
            WRITE(IW,*) 'FIELD DIRECTION',IFLD
         ENDIF
C
C     ----- GET PERTURBATION -----
C
         IF(IFLD.EQ.1) THEN
            IDAR=33
         ELSEIF(IFLD.EQ.2) THEN
            IDAR=34
         ELSEIF(IFLD.EQ.3) THEN
            IDAR=35
         ENDIF
         CALL HND_DAREAD(IDAF,IODA,X(I10),L2,IDAR)
         STATUS=PROP_GIAO_GET_INT(X(I10),L1,L2,IFLD)
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIAO, PROP_GIAO_GET_INT FAILED',911)
         ENDIF
C
C     ----- INITIAL U(1) -----  
C
         CALL HND_GIAU_DPL(X(I01),X(I00),X(I10),X(I60),X(I70),
     1                 NUM,NA,L0,L1)
C
         DMP=DAMP
         ITER=-1
  100    ITER=ITER+1
         IF(DBUG0) THEN
            WRITE(IW,*) 'ITERATION',ITER,MAXIT
         ENDIF
C
C     ----- D(1) -----
C
         CALL HND_GIAD_DPL(X(I00),X(I60),X(I70),X(I20),
     1                 NUM,NA,L0,L1)
C
C      ----- CHECK CONVERGENCE -----
C
         CALL HND_GIACVG_DPL(X(I20),X(I10),X(I60),X(I70),
     1                   NUM,NA,L0,L1,IEND,ITER,SOME)
         IF((IEND.NE.1).OR.(IEND.EQ.1.AND.DMP.NE.ZERO)) THEN      
C
C     ----- DAMPING -----
C
            IF(ITER.EQ.0) THEN
               CALL HND_DAWRIT(IDAF,IODA,X(I20),L2,31,NAV)
               STATUS=PROP_GIAO_PUT_D(X(I20),L1,L2)
               IF(.NOT.STATUS) THEN
               CALL ERRQUIT('HND_GIAO, PROP_GIAO_PUT_D FAILED',911)
               ENDIF
            ELSE
               IF(IEND.EQ.1) THEN
                  DMP=ZERO
               ENDIF
               IF(DMP.NE.ZERO) THEN
                  CALL HND_DAREAD(IDAF,IODA,X(I30),L2,31)
                  STATUS=PROP_GIAO_GET_D(X(I30),L1,L2)
                  IF(.NOT.STATUS) THEN
                  CALL ERRQUIT('HND_GIAO, PROP_GIAO_GET_D FAILED',911)
                  ENDIF
                  DO I=1,L2
                     X(I-1+I20)=DMP*X(I-1+I30)+(ONE-DMP)*X(I-1+I20)
                  ENDDO
                  CALL HND_DAWRIT(IDAF,IODA,X(I20),L2,31,NAV)
                  STATUS=PROP_GIAO_PUT_D(X(I20),L1,L2)
                  IF(.NOT.STATUS) THEN
                  CALL ERRQUIT('HND_GIAO, PROP_GIAO_PUT_D FAILED',911)
                  ENDIF
               ELSE
                  CALL HND_DAWRIT(IDAF,IODA,X(I20),L2,31,NAV)
                  STATUS=PROP_GIAO_PUT_D(X(I20),L1,L2)
                  IF(.NOT.STATUS) THEN
                  CALL ERRQUIT('HND_GIAO, PROP_GIAO_PUT_D FAILED',911)
                  ENDIF
               ENDIF
            ENDIF
            IF(DBUG0) THEN
               WRITE(IW,*) 'D(1)'
               CALL HND_PRTR(X(I20),L1)
            ENDIF
c
c     ----- make square matrix with -nw- format -----
c
      do j=1,l1
         do i=1,l1
            ii=max0(i,j)
            jj=min0(i,j)
            ij=(ii*(ii-1))/2 + jj
            x(i+l1*(j-1)+i50-1)=x(ij+i20-1)
         enddo
      enddo
      call hnd_hndnw_tran(x(i50),x(i60),l1)
      if(dbug0) then
         call hnd_prsq(x(i50),l1,l1,l1)
         call hnd_prsq(x(i60),l1,l1,l1)
      endif
C
C     ----- F(1) -----
C

c     ------ create ga block -----

      g_d_1    =ga_create_atom_blocked (geom, basis,'d1 in ga' )
      g_h_1    =ga_create_atom_blocked (geom, basis,'h1 in ga' )
      g_f_1(1) =ga_create_atom_blocked (geom, basis,'f1 in ga')

      call ga_zero(g_d_1)
      call ga_zero(g_h_1)
      call ga_zero(g_f_1)

      if(ga_nodeid().eq.0) then
         call ga_put(g_d_1,1,l1,1,l1,x(i60),l1)
         if(dbug0) then
            call ga_print(g_d_1)
         endif
      endif

c     ----- jfac and kfac for fock operator -----

      nfock   =  1
      jfac(1) =  1.0
      kfac(1) = -0.5

c     ----- get tol2e from scf_get_info.F -----

      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) then
          tol2e = 1d-9
      endif

c     ----- load symmetry information from scf_get_info.F -----

      If (.not.rtdb_get(rtdb, 'scf:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif

      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call scf_get_fock_param(rtdb,tol2e)
      call fock_2e(geom,basis,nfock,jfac,kfac,tol2e,oskel,
     &             g_d_1,g_f_1(1))
      call fock_2e_tidy(rtdb)
      call int_terminate()
      call schwarz_tidy()

      if(dbug0) then
         call ga_print(g_f_1(1))
      endif

c     ----- download g_f_1 to local memory:x(i50) -----

      call ga_get(g_f_1(1),1,num,1,num,x(i60),num)
      call hnd_nwhnd_tran(x(i60),x(i50),l1)
      if(dbug0) then
         call hnd_prsq(x(i50),l1,l1,l1)
      endif
      do j=1,l1
         do i=1,l1
            ii=max0(i,j)
            jj=min0(i,j)
            ij=(ii*(ii-1))/2 + jj
            x(ij+i30-1)=x(i+l1*(j-1)+i50-1)
         enddo
      enddo
      if(dbug0) then
         call hnd_prtr(x(i30),l1)
      endif

c     ----- release memory pointer -----

      status = ga_destroy ( g_d_1      )
      status = ga_destroy ( g_h_1      )
      status = ga_destroy ( g_f_1(1)   )

C
C     ----- ADD H(1) TO F(1) -----
C
            DO I=1,L2
               X(I-1+I30)=X(I-1+I30)-(-X(I-1+I10))
            ENDDO
            IF(DBUG0) THEN
               WRITE(IW,*) 'F(1)'
               CALL HND_PRTR(X(I30),L1)
            ENDIF
C
C     ----- NEW U(1) -----
C
            CALL HND_GIAU_DPL(X(I01),X(I00),X(I30),X(I60),X(I70),
     1                    NUM,NA,L0,L1)
C
            IF(DBUG0) THEN
               WRITE(IW,*) 'END OF ITERATION',ITER,MAXIT
            ENDIF
            IF(ITER.GT.MAXIT) THEN      
               WRITE(IW,9997)   
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            GO TO 100
         ELSE 
C
C     ----- CONVERGED ; SAVE U(1) AND D(1) -----
C
            WRITE(IW,*) '-DPL- CPHF CONVERGED'
         ENDIF
C
      ENDDO
         IF(OUT) THEN
C
C     ----- FOR DEBUGGING GET -F00- -----
C
            CALL HND_DAREAD(IDAF,IODA,X(I20),L2,16)

            call hnd_giao_get_dens(rtdb,basis,geom,
     1                             x(i20),x(i30),x(i50),x(i60),
     2                             scftyp)

            IF(DBUG) THEN
               WRITE(IW,*) 'D(0)'
               CALL HND_PRTR(X(I20),L1)
            ENDIF

c
c     ----- make square matrix with -nw- format -----
c
      do j=1,l1
         do i=1,l1
            ii=max0(i,j)
            jj=min0(i,j)
            ij=(ii*(ii-1))/2 + jj
            x(i+l1*(j-1)+i50-1)=x(ij+i20-1)
         enddo
      enddo
      call hnd_hndnw_tran(x(i50),x(i60),l1)
      if(dbug0) then
         call hnd_prsq(x(i50),l1,l1,l1)
         call hnd_prsq(x(i60),l1,l1,l1)
      endif


c     ------ create ga block -----

      g_d_1    =ga_create_atom_blocked (geom, basis,'d1 in ga' )
      g_h_1    =ga_create_atom_blocked (geom, basis,'h1 in ga' )
      g_f_1(1) =ga_create_atom_blocked (geom, basis,'f1 in ga')

      call ga_zero(g_d_1)
      call ga_zero(g_h_1)
      call ga_zero(g_f_1)

      if(ga_nodeid().eq.0) then
         call ga_put(g_d_1,1,l1,1,l1,x(i60),l1)
         if(dbug0) then
            call ga_print(g_d_1)
         endif
      endif

c     ----- jfac and kfac for fock operator -----

      nfock   =  1
      jfac(1) =  1.0
      kfac(1) = -0.5

c     ----- get tol2e from scf_get_info.F -----

      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) then
          tol2e = 1d-9
      endif

c     ----- load symmetry information from scf_get_info.F -----

      If (.not.rtdb_get(rtdb, 'scf:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif

      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call scf_get_fock_param(rtdb,tol2e)
      call fock_2e(geom,basis,nfock,jfac,kfac,tol2e,oskel,
     &             g_d_1,g_f_1(1))
      call fock_2e_tidy(rtdb)
      call int_terminate()
      call schwarz_tidy()

      if(dbug0) then
         call ga_print(g_f_1(1))
      endif


c     ----- download g_f_1 to local memory:x(i50) -----

      call ga_get(g_f_1(1),1,num,1,num,x(i60),num)
      call hnd_nwhnd_tran(x(i60),x(i50),l1)
      if(dbug0) then
         call hnd_prsq(x(i50),l1,l1,l1)
      endif
      do j=1,l1
         do i=1,l1
            ii=max0(i,j)
            jj=min0(i,j)
            ij=(ii*(ii-1))/2 + jj
            x(ij+i30-1)=x(i+l1*(j-1)+i50-1)
         enddo
      enddo
      if(dbug0) then
         call hnd_prtr(x(i30),l1)
      endif

c     ----- release memory pointer -----

      status = ga_destroy ( g_d_1      )
      status = ga_destroy ( g_h_1      )
      status = ga_destroy ( g_f_1(1)   )

            IF(DBUG) THEN
               WRITE(IW,*) 'F(0)'
               CALL HND_PRTR(X(I30),L1)
            ENDIF
C
         ENDIF
C
      ENDIF
      CALL HND_FIFLSH(IW)
C
C     ----- THE -GIAO- CPHF EQUATIONS NOW -----
C
      IF(OUT) THEN
         WRITE(IW,*) '-GIAO- CPHF EQUATIONS NOW'
      ENDIF
      CALL HND_FIFLSH(IW)
      DO I=1,3*NAT
         X(I+I02-1)=ZERO
      ENDDO
      DO I=1,9*NAT
         X(I+I03-1)=ZERO
         X(I+I04-1)=ZERO
         X(I+I05-1)=ZERO
      ENDDO
C
C     ----- CALCULATE -GIAO- PERTURBATION OPERATORS -----
C           -S10-     ON     -IFT1-
C           -H10-     ON     -IFT2-
C           -H01-     ON     -IFT3-
C
      IF(OUT) THEN
         WRITE(IW,*) '-H20- DIAMAGNETIC SUSCEPTIBILITY'
         WRITE(IW,*) 'L1,NOC= ',L1,NOC
         CALL HND_FIFLSH(IW)
      ENDIF
      CALL HND_GIAH20(X(J10),X(I00),B20,X(J40),L1,NOC,L1)   
            if(ga_nodeid().eq.0) then
      WRITE(IW,9993)
      WRITE(IW,9992) B20
            endif
      CALL HND_FIFLSH(IW)
C
      IF(OUT) THEN
         WRITE(IW,*) '-H01-'
         CALL HND_FIFLSH(IW)
      ENDIF
      CALL HND_GIAH01(X(J10),X(I00),X(J30),X(J40),L1,NOC,L0,L1)
      CALL HND_FIFLSH(IW)
C
      IF(OUT) THEN
         WRITE(IW,*) '-S10-'
         CALL HND_FIFLSH(IW)
      ENDIF
      CALL HND_GIAS10(X(J10),X(I00),X(J30),X(J40),L1,NOC,L0,L1)   
      CALL HND_FIFLSH(IW)
C
      IF(OUT) THEN
         WRITE(IW,*) '-H10-'
         CALL HND_FIFLSH(IW)
      ENDIF
C

      call hnd_giao_get_dens(rtdb,basis,geom,
     1                       x(i40),x(i50),x(i60),x(i70),
     2                       scftyp)
      do i=1,l2
         x(i+j40-1)=x(i+i40-1)
      enddo
      if(out) then
         write(iw,*) '-doo- in hnd_giaox'
         call hnd_prtr(x(j40),l1)
      endif

      IFCK0=J30
      IFCKX=J20
      IFCKY=J21
      IFCKZ=J22
      IDENS=J40
      IDIM =L1
      CALL HND_GIAH10(X(J10),X(I00),X(J30),X(J40),X(J50),X(J30),
     1            L1,NOC,L0,L1,NEXT)   
      CALL HND_FIFLSH(IW)
C
C     ----- READY FOR THE -CPHF- ITERATIONS -----
C
      MFLD=1
      NFLD=3
      DO IFLD=MFLD,NFLD
         IF(OUT) THEN
            WRITE(IW,9996)
            WRITE(IW,*) 'FIELD DIRECTION',IFLD
            CALL HND_FIFLSH(IW)
         ENDIF
C
C     ----- INITIAL -U(1)- FROM -H(1)- AND -S(1)- -----
C           -U(1)-  .....  RECORD -32- OF -IDAF-
C           -D(1)-  .....  RECORD -31- OF -IDAF-
C
         DO I=1,L3
            X(I+J10-1)=ZERO
         ENDDO
         CALL HND_GIAFIL(IFLD)
         CALL HND_GIAU1(X(I01),X(J30),X(J11),X(J10),X(J20),
     1              NUM,NA,L0,L1)
C
         DMP=DAMP
         ITER=-1
  200    ITER=ITER+1
         IF(DBUG) THEN
            WRITE(IW,*) 'ITERATION',ITER,MAXIT
            CALL HND_FIFLSH(IW)
         ENDIF
C
C     ----- D(1) -----
C
         CALL HND_GIAV1(X(I00),X(J30),X(J10),NUM,NA,L0,L1)
         CALL HND_GIAD1(X(I00),X(J10),X(J20),NUM,NA,L1)
C
C      ----- CHECK CONVERGENCE -----
C
         CALL HND_GIACVG(X(J20),X(J21),X(J30),X(J21),
     1               NUM,NA,L0,L1,IEND,ITER,SOME)
         IF((IEND.NE.1).OR.(IEND.EQ.1.AND.DMP.NE.ZERO)) THEN
C
C     ----- DAMPING -----
C
            IF(ITER.EQ.0) THEN
            if(ga_nodeid().eq.0) then
               CALL HND_DAWRIT(IDAF,IODA,X(J20),L3,31,NAV)
               STATUS=PROP_GIAO_PUT_D1(X(J20),L1,L3)
               IF(.NOT.STATUS) THEN
               CALL ERRQUIT('HND_GIAO, PROP_GIAO_PUT_D1 FAILED',911)
               ENDIF
            endif
            call ga_sync()
            ELSE
               IF(IEND.EQ.1) THEN
                  DMP=ZERO
               ENDIF
               IF(DMP.NE.ZERO) THEN
            if(ga_nodeid().eq.0) then
                  CALL HND_DAREAD(IDAF,IODA,X(J21),L3,31)
                  STATUS=PROP_GIAO_GET_D1(X(J21),L1,L3)
                  IF(.NOT.STATUS) THEN
                  CALL ERRQUIT('HND_GIAO, PROP_GIAO_GET_D1 FAILED',911)
                  ENDIF
            endif
            call ga_brdcst(788+1+msgdbl,x(j21),mdtob(l3),0)
            call ga_sync()
                  DO I=1,L3
                     X(I-1+J20)=DMP*X(I-1+J21)+(ONE-DMP)*X(I-1+J20)
                  ENDDO
            if(ga_nodeid().eq.0) then
                  CALL HND_DAWRIT(IDAF,IODA,X(J20),L3,31,NAV)
                  STATUS=PROP_GIAO_PUT_D1(X(J20),L1,L3)
                  IF(.NOT.STATUS) THEN
                  CALL ERRQUIT('HND_GIAO, PROP_GIAO_PUT_D1 FAILED',911)
                  ENDIF
            endif
            call ga_sync()
               ELSE
            if(ga_nodeid().eq.0) then
                  CALL HND_DAWRIT(IDAF,IODA,X(J20),L3,31,NAV)
                  STATUS=PROP_GIAO_PUT_D1(X(J20),L1,L3)
                  IF(.NOT.STATUS) THEN
                  CALL ERRQUIT('HND_GIAO, PROP_GIAO_PUT_D1 FAILED',911)
                  ENDIF
               ENDIF
            endif
            call ga_sync()
            ENDIF
C
C     ----- K(1) -----
C
            IF(DBUG) THEN
               WRITE(IW,*) 'D(1)'
               CALL HND_PRSQ(X(J20),L1,L1,L1)
            ENDIF

c     ------ create ga block -----

      g_d_1    =ga_create_atom_blocked (geom, basis,'d1 in ga' )
      g_h_1    =ga_create_atom_blocked (geom, basis,'h1 in ga' )
      g_f_1(1) =ga_create_atom_blocked (geom, basis,'f1 in ga')

      call ga_zero(g_d_1)
      call ga_zero(g_h_1)
      call ga_zero(g_f_1)

      if(ga_nodeid().eq.0) then
         call ga_put(g_d_1,1,l1,1,l1,x(j20),l1)
         if(dbug0) then
            call ga_print(g_f_1)
            call ga_print(g_d_1)
         endif
      endif

c     ----- jfac and kfac for fock operator -----

      nfock   =  1
      jfac(1) =  0.0
      kfac(1) = -1.0

c     ----- get tol2e from scf_get_info.F -----

      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) then
          tol2e = 1d-9
      endif

c     ----- load symmetry information from scf_get_info.F -----

      If (.not.rtdb_get(rtdb, 'scf:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif

      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call scf_get_fock_param(rtdb,tol2e)
      call giao_2e(geom,basis,nfock,jfac,kfac,tol2e,oskel,
     &             g_d_1,g_f_1(1))
      call fock_2e_tidy(rtdb)
      call int_terminate()
      call schwarz_tidy()

      if(dbug0) then
         call ga_print(g_f_1(1))
      endif

c     ----- download g_f_1 to local memory:x(i50) -----

      call ga_get(g_f_1(1),1,num,1,num,x(j11),num)

      if(dbug0) then
         write(iw,*) 'K10(1) in ao basis'
         call hnd_prsq(x(j11),l1,l1,l1)
      endif

c     ----- release memory pointer -----

      status = ga_destroy ( g_d_1      )
      status = ga_destroy ( g_h_1      )
      status = ga_destroy ( g_f_1(1)   )

C
            CALL HND_TFSQ(X(J10),X(J11),X(I00),X(J40),L0,L1,L1)
            IF(DBUG0) THEN
               WRITE(IW,*) 'K10(1) IN CANONICAL MO BASIS'
               CALL HND_PRSQ(X(J10),L1,L1,L1)
            ENDIF
C
C     ----- NEW U(1) -----
C
            CALL HND_GIAFIL(IFLD)
            CALL HND_GIAU1(X(I01),X(J30),X(J11),X(J10),X(J20),
     1                 NUM,NA,L0,L1)
C
            IF(DBUG) THEN
               WRITE(IW,*) 'END OF ITERATION',ITER,MAXIT
            ENDIF
            IF(ITER.GT.MAXIT) THEN
               WRITE(IW,9997)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            GO TO 200
         ELSE
C
C     ----- CONVERGED ; SAVE U(1) AND D(1) -----
C           -D(1)-    TO    -31-
C           -U(1)-    TO    -32-
C
         if(ga_nodeid().eq.0) then
            WRITE(IW,*) '-GIAO- CPHF CONVERGED'
            CALL HND_DAWRIT(IDAF,IODA,X(J20),L3,31,NAV)
            STATUS=PROP_GIAO_GET_D1(X(J20),L1,L3)
            IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIAO, PROP_GIAO_GET_D1 FAILED',911)
            ENDIF
         endif
         call ga_brdcst(789+3+msgdbl,x(j20),mdtob(l3),0)
         call ga_sync()
         if(ga_nodeid().eq.0) then
            CALL HND_DAWRIT(IDAF,IODA,X(J30),L3,32,NAV)
            STATUS=PROP_GIAO_GET_U1(X(J30),L1,L3)
         endif
         call ga_brdcst(789+3+msgdbl,x(j30),mdtob(l3),0)
         call ga_sync()
            IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIAO, PROP_GIAO_GET_U1 FAILED',911)
            ENDIF
            IF(OUT) THEN
               WRITE(IW,*) 'CONVERGED -U1-'
               CALL HND_PRSQ(X(J30),L0,L0,L1)
               WRITE(IW,*) 'CONVERGED -D1-'
               CALL HND_PRSQ(X(J20),L1,L1,L1)
            ENDIF
         ENDIF
C
C     ----- -H01- CONTRIBUTION TO MAGNETIC SHIELDING -----
C           -D(1)-    FROM  -31-
C           -D(0)-    FROM  -16-
C
         CALL HND_DAREAD(IDAF,IODA,X(J40),L2,16)
         call hnd_giao_get_dens(rtdb,basis,geom,
     1                          x(j40),x(i50),x(i60),x(i70),
     2                          scftyp)
         if(ga_nodeid().eq.0) then
         CALL HND_DAREAD(IDAF,IODA,X(J30),L3,31)
         STATUS=PROP_GIAO_GET_D1(X(J30),L1,L3)
         IF(.NOT.STATUS) THEN
         CALL ERRQUIT('HND_GIAO, PROP_GIAO_GET_D1 FAILED',911)
         ENDIF
         endif
         call ga_brdcst(789+3+msgdbl,x(j30),mdtob(l3),0)
         call ga_sync()
         CALL HND_GIAE11(X(I02),X(I03),X(J10),X(J30),X(J40),L1,IFLD)
C
      ENDDO
C
C     ----- -H11- CONTRIBUTIONS TO MAGNETIC SHIELDING -----
C
      CALL HND_DAREAD(IDAF,IODA,X(J40),L2,16)
      call hnd_giao_get_dens(rtdb,basis,geom,
     1                       x(j40),x(i50),x(i60),x(i70),
     2                       scftyp)
      CALL HND_GIAH11(X(I04),X(I05),X(J40))
C
      UNITS =PPM/TWO    
      DO I=1,9*NAT
         X(I+I03-1)= X(I+I03-1)*UNITS
         X(I+I04-1)= X(I+I04-1)*UNITS
         X(I+I05-1)= X(I+I05-1)*UNITS
      ENDDO
      if(ga_nodeid().eq.0) then
      CALL HND_GIAPRT(X(I03),NAT,1)
      IF(OUT) THEN
         CALL HND_GIAPRT(X(I04),NAT,2)
         CALL HND_GIAPRT(X(I05),NAT,3)
      ENDIF
      endif
C
      DO I=1,9*NAT
         X(I+I04-1)=(X(I+I04-1)+X(I+I05-1))
      ENDDO
      if(ga_nodeid().eq.0) then
      CALL HND_GIAPRT(X(I04),NAT,4)
      endif
C
      DO I=1,9*NAT
         X(I+I03-1)=(X(I+I03-1)+X(I+I04-1))
      ENDDO
      if(ga_nodeid().eq.0) then
      CALL HND_GIAPRT(X(I03),NAT,5)
      endif
C
C     ----- ALL DONE -----
C
      if(ga_nodeid().eq.0) then
      WRITE(IW,9996)
      CALL HND_FIFLSH(IW)
      endif
C
      RETURN
 9999 FORMAT(
     1 10X,57(1H-),/,
     2 10X,'-giao- magnetic susceptibility/chemical shielding tensors',/,
     3 10X,57(1H-))
 9998 FORMAT(' -giao- implemented for closed shell systems only. stop')
 9997 FORMAT(' too many iterations in the -giao- cphf . stop. ')
 9996 FORMAT(/)
 9995 FORMAT(I7,3F15.8)
 9994 FORMAT(' memory needed in -giaox- . need = ',i10)
 9993 FORMAT(/,' diamagnetic susceptibility tensor -b20- (a.u.) = ',/,
     1         ' ------------------------------------------------ ',/)
 9992 FORMAT(' b20(xx) = ',F20.10,/,
     1       ' b20(yy) = ',F20.10,/,
     2       ' b20(zz) = ',F20.10,/,
     3       ' b20(xy) = ',F20.10,/,
     4       ' b20(xz) = ',F20.10,/,
     5       ' b20(yz) = ',F20.10,/)
      END
      subroutine hnd_giao_get_vecs(rtdb,basis,geom,
     1                             v1nw,v2nw,e1nw,
     2                             v1hnd,e1hnd,
     3                             scftyp,
     4                             nclosed ,nopen ,nvirt ,
     5                             nclosedu,nopenu,nvirtu)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "schwarz.fh"
#include "util.fh"
c
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      double precision      v1nw,v2nw,e1nw
      double precision      v1hnd,e1hnd
      integer               nclosed ,nopen ,nvirt 
      integer               nclosedu,nopenu,nvirtu
c
      logical  int_normalize
      external int_normalize

      integer       rtdb              ! [input] database handle
      integer       geom, basis       ! handles
      logical       status
      character*255 title
      character*255 movecs_in, movecs_out
      character*8   scftyp 
      dimension     v1nw(*),v2nw(*),e1nw(*)
      dimension     v1hnd(*),e1hnd(*)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'prop:input vectors',  1, 
     $     movecs_in )) then
           call util_file_name('movecs', .false.,.false.,movecs_in )
      endif
      call util_file_name_resolve(movecs_in , .false.)
      if (.not. rtdb_cget(rtdb, 'prop:output vectors',  1,
     $     movecs_out)) then
           call util_file_name('movecs', .false.,.false.,movecs_out)
      endif
      call util_file_name_resolve(movecs_out, .false.)
c
      call hnd_giao_vecs_read(rtdb,basis,geom,movecs_in,
     1                        v1nw,v2nw,e1nw,
     2                        v1hnd,e1hnd,
     3                        scftyp,
     4                        nclosed ,nopen ,nvirt ,
     5                        nclosedu,nopenu,nvirtu)
c
      return
      end
      subroutine hnd_giao_vecs_read(rtdb,basis,geom,movecs,
     1                              v1nw,v2nw,e1nw,
     2                              v1hnd,e1hnd,
     3                              scftyp,
     4                              nclosed ,nopen ,nvirt ,
     5                              nclosedu,nopenu,nvirtu)
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"

      logical       status
      integer       rtdb, basis, geom
      character*(*) movecs
      character*4   scftype
      character*8   scftyp 
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      logical  movecs_read_header, movecs_read
      external movecs_read_header, movecs_read

      integer g_molec
  
      integer       max_sh_bf, max_at_bf 
      integer       nat, nbf, nopen , nclosed , nvirt 
      integer                 nopenu, nclosedu, nvirtu  
      integer       ndens, nbf2, nocc, i, iset
      character*255 title_molec, basis_name_molec
      character*20  scftype_molec
      integer       nbf_molec, nsets_molec, nmo_molec(2)
      integer                               nmo
c
      integer   l_eval, k_eval, l_occ, k_occ
c
      double precision v1nw(*),v2nw(*),e1nw(*)
      double precision v1hnd(*),e1hnd(*)
c
      logical dbug
      dbug=.false.
c
c     get information about basis set
c
      status = geom_ncent(geom,nat)
      if ( .not.status )
     & call errquit('vecs_read: could not get natoms',0)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) 
     & call errquit('vecs_read: could not get nbf',0)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_sh_bf',0)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_at_bf',0)
c
c     read vector file header
c
      if (.not. movecs_read_header(movecs, title_molec,
     &    basis_name_molec, scftype_molec, nbf_molec,
     &    nsets_molec, nmo_molec, 2))
     &    call errquit('vecs_read: failed to read header', 911)
      if(dbug) then
         write(6,9999)
     1     movecs
     2   ,title_molec
     3   ,basis_name_molec
     4   ,scftype_molec
     5   ,nbf_molec
     6   ,nsets_molec
     7   ,nmo_molec
      endif
c
      if (nsets_molec.eq.1) then
         scftype = 'RHF'
         scftyp  = 'RHF     '
      else
         scftype = 'UHF'
         scftyp  = 'UHF     '
      endif
      call inp_ucase(scftype)

C     # of eigenvalues and density matrices
      if ( scftype .eq. 'UHF' ) then
        nbf2  = 2 * nbf
        ndens = 3
        nmo   = max(nmo_molec(1),nmo_molec(2))
      else if ( scftype .eq. 'RHF' ) then
        nbf2  = nbf
        ndens = 1
        nmo   =     nmo_molec(1)
      endif

C     allocate and initialize global and local memory

C     occupation numbers
      status = MA_push_get(MT_DBL, nbf2,'occ. numbers', l_occ, 
     $                                                  k_occ)
      if (.not.status) 
     $     call errquit('property: could not allocate l_occ',911)
C     eigenvalues
      status = MA_push_get(MT_DBL, nbf2,'MO eval', l_eval, k_eval)
      if (.not.status) 
     $     call errquit('property: could not allocate l_eval',911)
      
      status = ga_create(MT_DBL, nbf_molec, nmo,' molec MOs',
     $                           nbf_molec, nmo,g_molec)
      if (.not.status)
     $     call errquit('property: GA for molecular MOs',911)
c
c     read ( alpha ) orbitals
c
      iset=1
      status = movecs_read(movecs,iset,
     $         dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
      if (.not.status)
     $  call errquit('property: could not read molecular MOs',911)
       
      nclosed = 0
      nopen   = 0
      nvirt   = 0
      do i=1,nmo_molec(1)
         if(dbl_mb(k_occ+i-1).eq.2.0d0) then
            nclosed=nclosed+1
         elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
            nopen=nopen+1
         else
            nvirt=nvirt+1
         endif
      enddo
      if(dbug) then
         write(6,9998) nclosed , nopen , nvirt
      endif
c
c     ----- stop if not closed-shell scf ! -----
c
      if(scftype_molec.ne.'scf'.or.
     $   nsets_molec.ne.1      .or.
     $   nopen.ne.0               ) then
         call errquit('giao: not a closed-shell wfn. stop',911)
      else
         if(dbug) then
            write(6,*) 'this is a closed-shell wavefunction'
         endif
      endif
c
c     read ( beta ) orbitals
c
      if( nsets_molec.gt.1 ) then
            iset=2
         status = movecs_read(movecs,iset,
     $            dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
         if (.not.status)
     $     call errquit('property: could not read molecular MOs',911)

         nclosedu= 0
         nopenu  = 0
         nvirtu  = 0
         do i=1,nmo_molec(2)
            if(dbl_mb(k_occ+i-1).eq.2.0d0) then
               nclosedu=nclosedu+1
            elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
               nopenu=nopenu+1
            else
               nvirtu=nvirtu+1
            endif
         enddo
         if(dbug) then
            write(6,9997) nclosedu, nopenu, nvirtu
         endif
      endif

      status = ga_destroy ( g_molec )

      call hnd_giao_vecs_make (geom, basis, 
     $     dbl_mb(k_eval), dbl_mb(k_occ), nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     v1nw, v2nw, e1nw, v1hnd, e1hnd )

c
c     free temporary arrays
c
      status = MA_pop_stack ( l_eval)
      status = MA_pop_stack ( l_occ )

      return
 9999 format(
     1   /,'movecs             = ---',a24,'---',
     2   /,'title_molec        = ---',a24,'---',
     3   /,'basis_name_molec   = ---',a24,'---',
     4   /,'scftype_molec      = ---',a24,'---',
     5   /,'nbf_molec          = ',i12,
     6   /,'nsets_molec        = ',i12,
     7   /,'nmo_molec          = ',i12,i12 
     8      )
 9998 format(' nclosed , nopen , nvirt  = ',3i10)
 9997 format(' nclosedu, nopenu, nvirtu = ',3i10)
      end
      subroutine hnd_giao_vecs_make (geom, basis,  
     $     evals, occ, nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     v1nw, v2nw, e1nw, v1hnd, e1hnd )

c     read mo vecs 

      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer geom, basis, ndens
      integer nbf, nbf2, nclosed, nopen, nclosedu, nopenu
      double precision evals, occ
      double precision da, db, d1nw, d2nw
      double precision v1nw ,v2nw , e1nw
      double precision v1hnd,e1hnd
      character*4 scftype
      character*(*) movecs

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  movecs_read, file_read_ga
      external movecs_read, file_read_ga

      dimension evals(nbf2), occ(nbf2)
      dimension v1nw(*), v2nw(*), e1nw(*)
      dimension v1hnd(*), e1hnd(*)

C     local variables
      integer nocc, l_ind, k_ind, l_ind2, k_ind2, g_vecs, g_vecs2
      integer i, j, ij
      integer iopen
      logical status
c
      logical dbug
      dbug=.false.

      if(dbug) then
         write(6,*) 'in -hnd_giao_vecs_make'
      endif
C  get MO vectors from file 

      g_vecs = ga_create_atom_blocked ( geom, basis, 'mo vectors' )
      if ( scftype.eq. 'UHF' ) then
        g_vecs2 = ga_create_atom_blocked (geom, basis, 'uhf mo vectors')
      endif

      status = movecs_read (movecs, 1, occ, evals, g_vecs)
      if (.not.status) then
        call errquit('giao_vecs_make: could not read mo vectors', 110)
      endif

      if (scftype .eq. 'UHF') then
        if (.not. movecs_read(movecs, 2, occ(nbf+1), evals(nbf+1),
     $       g_vecs2)) then
        call errquit('giao_vecs_make: could not read uhf vectors', 110)
        endif
      endif
      
      if(dbug) then
         write(6,*) 'passing by here..... do not stop for now'
      endif
c

      if (scftype .eq. 'RHF') then
 
        call ga_get(g_vecs,1,nbf,1,nbf,v1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran_vecs(v1nw,v1hnd,nbf)
        do i=1,nbf
           e1hnd(i)=evals(i)
        enddo

        if(dbug) then
           write(6,*) 'alpha vectors'
           call hnd_prev(v1hnd,evals,nbf,nbf,nbf)
        endif

      else  ! UHF

        call ga_get(g_vecs2,1,nbf,1,nbf,v2nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran_vecs(v2nw,v1nw,nbf)

        if(dbug) then
           write(6,*) 'beta  vectors'
           call hnd_prev(v1nw,evals(1+nbf),nbf,nbf,nbf)
        endif

      endif

c

      status = ga_destroy ( g_vecs )
      if ( scftype .eq. 'UHF' ) then
        status = ga_destroy ( g_vecs2 )
      endif

      return
      end
      subroutine hnd_giao_get_dens(rtdb,basis,geom,
     1                             d1hnd,d2hnd,   
     2                             d1nw,d2nw,   
     3                             scftyp)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "schwarz.fh"
#include "util.fh"
c
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      double precision      d1nw,d2nw
      double precision      d1hnd,d2hnd
c
      integer       rtdb              ! [input] database handle
      integer       geom, basis       ! handles
      logical       status
      character*255 title
      character*255 movecs_in, movecs_out
      character*8   scftyp 
      dimension     d1nw(*),d2nw(*)
      dimension     d1hnd(*),d2hnd(*)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'prop:input vectors',  1, 
     $     movecs_in )) then
           call util_file_name('movecs', .false.,.false.,movecs_in )
      endif
      call util_file_name_resolve(movecs_in , .false.)
      if (.not. rtdb_cget(rtdb, 'prop:output vectors',  1,
     $     movecs_out)) then
           call util_file_name('movecs', .false.,.false.,movecs_out)
      endif
      call util_file_name_resolve(movecs_out, .false.)
c
      call hnd_giao_dens_read(rtdb,basis,geom,movecs_in,
     1                        d1nw,d2nw,
     2                        d1hnd,d2hnd,
     3                        scftyp)
c
      return
      end
      subroutine hnd_giao_dens_read(rtdb,basis,geom,movecs,
     1                              d1nw,d2nw,
     2                              d1hnd,d2hnd,
     3                              scftyp)
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"

      logical       status
      integer       rtdb, basis, geom
      character*(*) movecs
      character*4   scftype
      character*8   scftyp 
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      logical  movecs_read_header, movecs_read
      external movecs_read_header, movecs_read

      integer g_molec
  
      integer       max_sh_bf, max_at_bf 
      integer       nat, nbf, nopen , nclosed , nvirt 
      integer                 nopenu, nclosedu, nvirtu  
      integer       ndens, nbf2, nocc, i, iset
      character*255 title_molec, basis_name_molec
      character*20  scftype_molec
      integer       nbf_molec, nsets_molec, nmo_molec(2)
      integer                               nmo
c
      integer   l_eval, k_eval, l_occ, k_occ
c
      double precision d1nw(*),d2nw(*)
      double precision d1hnd(*),d2hnd(*)
c
      logical dbug
      dbug=.false.
c
c     get information about basis set
c
      status = geom_ncent(geom,nat)
      if ( .not.status )
     & call errquit('vecs_read: could not get natoms',0)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) 
     & call errquit('vecs_read: could not get nbf',0)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_sh_bf',0)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status )
     & call errquit('vecs_read: could not get max_at_bf',0)
c
c     read vector file header
c
      if (.not. movecs_read_header(movecs, title_molec,
     &    basis_name_molec, scftype_molec, nbf_molec,
     &    nsets_molec, nmo_molec, 2))
     &    call errquit('vecs_read: failed to read header', 911)
      if(dbug) then
         write(6,9999)
     1     movecs
     2   ,title_molec
     3   ,basis_name_molec
     4   ,scftype_molec
     5   ,nbf_molec
     6   ,nsets_molec
     7   ,nmo_molec
      endif
c
      if (nsets_molec.eq.1) then
         scftype = 'RHF'
         scftyp  = 'RHF     '
      else
         scftype = 'UHF'
         scftyp  = 'UHF     '
      endif
      call inp_ucase(scftype)

C     # of eigenvalues and density matrices
      if ( scftype .eq. 'UHF' ) then
        nbf2  = 2 * nbf
        ndens = 3
        nmo   = max(nmo_molec(1),nmo_molec(2))
      else if ( scftype .eq. 'RHF' ) then
        nbf2  = nbf
        ndens = 1
        nmo   =     nmo_molec(1)
      endif

C     allocate and initialize global and local memory

C     occupation numbers
      status = MA_push_get(MT_DBL, nbf2,'occ. numbers', l_occ, 
     $                                                  k_occ)
      if (.not.status) 
     $     call errquit('property: could not allocate l_occ',911)
C     eigenvalues
      status = MA_push_get(MT_DBL, nbf2,'MO eval', l_eval, k_eval)
      if (.not.status) 
     $     call errquit('property: could not allocate l_eval',911)
      
      status = ga_create(MT_DBL, nbf_molec, nmo,' molec MOs',
     $                           nbf_molec, nmo,g_molec)
      if (.not.status)
     $     call errquit('property: GA for molecular MOs',911)
c
c     read ( alpha ) orbitals
c
      iset=1
      status = movecs_read(movecs,iset,
     $         dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
      if (.not.status)
     $  call errquit('property: could not read molecular MOs',911)
       
      nclosed = 0
      nopen   = 0
      nvirt   = 0
      do i=1,nmo_molec(1)
         if(dbl_mb(k_occ+i-1).eq.2.0d0) then
            nclosed=nclosed+1
         elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
            nopen=nopen+1
         else
            nvirt=nvirt+1
         endif
      enddo
      if(dbug) then
         write(6,9998) nclosed , nopen , nvirt
      endif
c
c     read ( beta ) orbitals
c
      if( nsets_molec.gt.1 ) then
            iset=2
         status = movecs_read(movecs,iset,
     $            dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
         if (.not.status)
     $     call errquit('property: could not read molecular MOs',911)

         nclosedu= 0
         nopenu  = 0
         nvirtu  = 0
         do i=1,nmo_molec(2)
            if(dbl_mb(k_occ+i-1).eq.2.0d0) then
               nclosedu=nclosedu+1
            elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
               nopenu=nopenu+1
            else
               nvirtu=nvirtu+1
            endif
         enddo
         if(dbug) then
            write(6,9997) nclosedu, nopenu, nvirtu
         endif
      endif

      status = ga_destroy ( g_molec )

      call hnd_giao_dens_make (geom, basis, 
     $     dbl_mb(k_eval), dbl_mb(k_occ), nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     d1nw, d2nw, d1hnd, d2hnd )

c
c     free temporary arrays
c
      status = MA_pop_stack ( l_eval)
      status = MA_pop_stack ( l_occ )

      return
 9999 format(
     1   /,'movecs             = ',a24,
     2   /,'title_molec        = ',a24,
     3   /,'basis_name_molec   = ',a24,
     4   /,'scftype_molec      = ',a24,
     5   /,'nbf_molec          = ',i12,
     6   /,'nsets_molec        = ',i12,
     7   /,'nmo_molec          = ',i12,i12 
     8      )
 9998 format(' nclosed , nopen , nvirt  = ',3i10)
 9997 format(' nclosedu, nopenu, nvirtu = ',3i10)
      end
      subroutine hnd_giao_dens_make (geom, basis,  
     $     evals, occ, nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     d1nw, d2nw, d1hnd, d2hnd )

c     read mo vecs 

      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer geom, basis, ndens
      integer nbf, nbf2, nclosed, nopen, nclosedu, nopenu
      double precision evals, occ
      double precision d1nw ,d2nw 
      double precision d1hnd,d2hnd
      character*4 scftype
      character*(*) movecs

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  movecs_read, file_read_ga
      external movecs_read, file_read_ga

      dimension evals(nbf2), occ(nbf2)
      dimension d1nw(*), d2nw(*)
      dimension d1hnd(*), d2hnd(*)

C     local variables
      integer g_vecs, g_vecs2
      integer g_dens, g_dens2
      integer nocc, l_ind, k_ind, l_ind2, k_ind2
      integer i, j, ij
      integer iopen
      logical status
c
      logical dbug
      dbug=.false.

      if(dbug) then
         write(6,*) 'in -hnd_giao_dens_make'
      endif

C  get MO vectors from file 

      g_vecs = ga_create_atom_blocked ( geom, basis, 'mo vectors' )
      g_dens = ga_create_atom_blocked ( geom, basis, 'mo density' )
      if ( scftype.eq. 'UHF' ) then
        g_vecs2 = ga_create_atom_blocked (geom, basis, 'uhf mo vectors')
        g_dens2 = ga_create_atom_blocked (geom, basis, 'uhf mo density')
      endif

      status = movecs_read(movecs,1,occ,evals,g_vecs)
      if (.not.status) then
        call errquit('giao_dens_make: could not read mo vectors', 110)
      endif

      if (scftype .eq. 'UHF') then
        status = movecs_read(movecs,2,occ(nbf+1),evals(nbf+1),g_vecs2)
        if (.not.status) then
        call errquit('giao_dens_make: could not read uhf vectors', 110)
        endif
      endif
      
      if(dbug) then
         write(6,*) 'passing by here..... do not stop for now'
      endif
c

      if (scftype .eq. 'RHF') then

        if(nclosed.eq.0) then
          call ga_zero(g_dens)
        else
          call ga_matmul_patch('n','t',2.0d+00,0.0d+00,
     1         g_vecs,1,nbf,1,nclosed,
     2         g_vecs,1,nclosed,1,nbf,
     3         g_dens,1,nbf,1,nbf)
        endif
        if(nopen.gt.0) then
          do iopen=1,nopen
            call ga_matmul_patch('n','t',occ(nclosed+iopen),1.0d+00,
     1           g_vecs,1,nbf,nclosed+iopen,nclosed+iopen,
     2           g_vecs,nclosed+iopen,nclosed+iopen,1,nbf,
     3           g_dens,1,nbf,1,nbf)
          enddo
        endif
        call ga_symmetrize(g_dens)
        call ga_get(g_dens,1,nbf,1,nbf,d1nw,nbf)
 
c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d1hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'alpha density'
           call hnd_prtr(d1hnd,nbf) 
        endif

      else  ! UHF

        call ga_zero(g_dens )
        call ga_zero(g_dens2)

        if(nopen.gt.0) then
          do iopen=1,nopen
            call ga_matmul_patch('n','t',occ(iopen),1.0d+00,
     1           g_vecs,1,nbf,iopen,iopen,
     2           g_vecs,iopen,iopen,1,nbf,
     3           g_dens,1,nbf,1,nbf)
          enddo
        endif
        if(nopenu.gt.0) then
          do iopen=1,nopenu
            call ga_matmul_patch('n','t',occ(nbf+iopen),1.0d+00,
     1           g_vecs2,1,nbf,iopen,iopen, 
     2           g_vecs2,iopen,iopen,1,nbf,
     3           g_dens2,1,nbf,1,nbf)
          enddo
        endif

        call ga_symmetrize(g_dens)
        call ga_get(g_dens,1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d2hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'alpha density'
           call hnd_prtr(d1hnd,nbf) 
        endif

c

        call ga_symmetrize(g_dens2)
        call ga_get(g_dens2,1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       triangularize

        ij=0
        do i=1,nbf
          do j=1,i
            ij=ij+1
            d2hnd(ij)=d2nw(i+nbf*(j-1))
          enddo
        enddo

        if(dbug) then
           write(6,*) 'beta  density'
           call hnd_prtr(d2hnd,nbf) 
        endif

      endif

c

      status = ga_destroy ( g_vecs )
      status = ga_destroy ( g_dens )
      if ( scftype .eq. 'UHF' ) then
        status = ga_destroy ( g_vecs2 )
        status = ga_destroy ( g_dens2 )
      endif

      return
      end
      logical function prop_giao_put_int(gint,int1,int2,ig)    
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer int1
      integer int2
      integer ig
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision gint(int2)     
      character*255 giao_int_fil
      data ift /31/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.int',
     1     .false.,.false.,giao_int_fil)
      open(unit=ift, file=giao_int_fil, form='unformatted',
     1     access='sequential', status='unknown', err=911)
c
      if(out) then
         write(iw,*) 'giao_put_int = ', giao_int_fil
         write(iw,*) 'int1,int2,ig = ',int1,int2,ig
      endif
      if(dbug) then
         call hnd_prtr(gint,int1)
      endif
c
      rewind ift
      if(ig.gt.1) then
         do i=1,ig-1
            read(ift)
            if(dbug) then
               write(iw,*) 'position i,ig = ',i,ig
            endif
         enddo
      endif
      write(ift) gint 
      rewind ift
c
      close(ift,status='keep')
c
      prop_giao_put_int=.true.
      return
c
 911  call errquit('prop_giao_put_int: open failed',0)
c
      end
      logical function prop_giao_get_int(gint,int1,int2,ig)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer int1
      integer int2
      integer ig
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      integer ios
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision gint(int2)
      character*255 giao_int_fil
      data ift /31/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.int',
     1     .false.,.false.,giao_int_fil)
      open(unit=ift, file=giao_int_fil, form='unformatted',
     1     access='sequential', status='old', err=911)
      if(out) then
         write(iw,*) 'giao_get_int = ', giao_int_fil
         write(iw,*) 'int1,int2,ig = ',int1,int2,ig
         do i=1,3
            if(dbug) then
               write(iw,*) 'position i,ig = ',i,ig
            endif
            read(ift,iostat=ios,end=912,err=912)
         enddo
      endif
c
      rewind ift
      if(ig.gt.1) then
         do i=1,ig-1
            read(ift)
            if(dbug) then
               write(iw,*) 'position i,ig = ',i,ig
            endif
         enddo
      endif
      read(ift,end=910,err=910,iostat=ios) gint
c
      close(ift,status='keep')
c
      if(dbug) then
         call hnd_prtr(gint,int1) 
      endif
c
      prop_giao_get_int=.true.
      return
c
c     ----- unable to complete the read ---
c
 910  close(ift,status='keep')
      if(out) then
         write(iw,*) 'iostat = ',ios
         write(iw,*) 'unable to complete the read in _giao_get_int 910'
      endif
      prop_giao_get_int=.false.
      return
c
c     ----- unable to open ; stop ... -----
c
 911  call errquit('prop_giao_get_int: open failed',911)
      return
c
c     ----- unable to complete the read ---
c
 912  close(ift,status='keep')
      if(out) then
         write(iw,*) 'position i,ig = ',i,ig
         write(iw,*) 'iostat = ',ios
         write(iw,*) 'unable to complete the read in _giao_get_int 912'
      endif
      prop_giao_get_int=.false.
      return
c
      end
      logical function prop_giao_put_d(d,l1,l2)    
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l2
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision d(l2)     
      character*255 giao_d_fil
      data ift /32/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.d',
     1     .false.,.false.,giao_d_fil)
      open(unit=ift, file=giao_d_fil, form='unformatted',
     1     access='sequential', status='unknown', err=911)
c
      if(out) then
         write(iw,*) 'giao_put_d = ',giao_d_fil
         write(iw,*) 'l1,l2      = ',l1,l2
      endif
      if(dbug) then
         call hnd_prtr(d,l1)
      endif
c
      rewind ift
      write(ift) d 
      rewind ift
c
      close(ift,status='keep')
c
      prop_giao_put_d=.true.
      return
c
 911  call errquit('prop_giao_put_d: open failed',0)
c
      end
      logical function prop_giao_get_d(d,l1,l2)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l2
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      integer ios
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision d(l2)
      character*255 giao_d_fil
      data ift /32/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.d',
     1     .false.,.false.,giao_d_fil)
      open(unit=ift, file=giao_d_fil, form='unformatted',
     1     access='sequential', status='old', err=911)
      if(out) then
         write(iw,*) 'giao_get_d = ',giao_d_fil
         write(iw,*) 'l1,l2      = ',l1,l2
      endif
c
      rewind ift
      read(ift,end=910,err=910,iostat=ios) d
c
      close(ift,status='keep')
c
      if(dbug) then
         call hnd_prtr(d,l1) 
      endif
c
      prop_giao_get_d=.true.
      return
c
c     ----- unable to complete the read ---
c
 910  close(ift,status='keep')
      if(out) then
         write(iw,*) 'iostat = ',ios
         write(iw,*) 'unable to complete the read in _giao_get_d 910'
      endif
      prop_giao_get_d=.false.
      return
c
c     ----- unable to open ; stop ... -----
c
 911  call errquit('prop_giao_get_d: open failed',911)
      return
c
      end
      logical function prop_giao_put_d1(d1,l1,l3)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l3
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision d1(l3)
      character*255 giao_d_fil
      data ift /32/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.d',
     1     .false.,.false.,giao_d_fil)
      open(unit=ift, file=giao_d_fil, form='unformatted',
     1     access='sequential', status='unknown', err=911)
c
      if(out) then
         write(iw,*) 'giao_put_d = ',giao_d_fil
         write(iw,*) 'l1,l3      = ',l1,l3
      endif
      if(dbug) then
         call hnd_prsq(d1,l1,l1,l1)
      endif
c
      rewind ift
      write(ift) d1
      rewind ift
c
      close(ift,status='keep')
c
      prop_giao_put_d1=.true.
      return
c
 911  call errquit('prop_giao_put_d1: open failed',0)
c
      end
      logical function prop_giao_get_d1(d1,l1,l3)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l3
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      integer ios
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision d1(l3)
      character*255 giao_d_fil
      data ift /32/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.d',
     1     .false.,.false.,giao_d_fil)
      open(unit=ift, file=giao_d_fil, form='unformatted',
     1     access='sequential', status='old', err=911)
      if(out) then
         write(iw,*) 'giao_get_d = ',giao_d_fil
         write(iw,*) 'l1,l3      = ',l1,l3
      endif
c
      rewind ift
      read(ift,end=910,err=910,iostat=ios) d1
c
      close(ift,status='keep')
c
      if(dbug) then
         call hnd_prsq(d1,l1,l1,l1)
      endif
c
      prop_giao_get_d1=.true.
      return
c
c     ----- unable to complete the read ---
c
 910  close(ift,status='keep')
      if(out) then
         write(iw,*) 'iostat = ',ios
         write(iw,*) 'unable to complete the read in _giao_get_d1 910'
      endif
      prop_giao_get_d1=.false.
      return
c
c     ----- unable to open ; stop ... -----
c
 911  call errquit('prop_giao_get_d1: open failed',911)
      return
c
      end
      logical function prop_giao_put_u(u,l1,l3)    
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l3
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision u(l3)     
      character*255 giao_u_fil
      data ift /33/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.u',
     1     .false.,.false.,giao_u_fil)
      open(unit=ift, file=giao_u_fil, form='unformatted',
     1     access='sequential', status='unknown', err=911)
c
      if(out) then
         write(iw,*) 'giao_put_u = ',giao_u_fil
         write(iw,*) 'l1,l3      = ',l1,l3
      endif
      if(dbug) then
         call hnd_prsq(u,l1,l1,l1)
      endif
c
      rewind ift
      write(ift) u 
      rewind ift
c
      close(ift,status='keep')
c
      prop_giao_put_u=.true.
      return
c
 911  call errquit('prop_giao_put_u: open failed',0)
c
      end
      logical function prop_giao_get_u(u,l1,l3)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l3
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      integer ios
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision u(l3)
      character*255 giao_u_fil
      data ift /33/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.u',
     1     .false.,.false.,giao_u_fil)
      open(unit=ift, file=giao_u_fil, form='unformatted',
     1     access='sequential', status='old', err=911)
      if(out) then
         write(iw,*) 'giao_get_u = ',giao_u_fil
         write(iw,*) 'l1,l3      = ',l1,l3
      endif
c
      rewind ift
      read(ift,end=910,err=910,iostat=ios) u
c
      close(ift,status='keep')
c
      if(dbug) then
         call hnd_prsq(u,l1,l1,l1) 
      endif
c
      prop_giao_get_u=.true.
      return
c
c     ----- unable to complete the read ---
c
 910  close(ift,status='keep')
      if(out) then
         write(iw,*) 'iostat = ',ios
         write(iw,*) 'unable to complete the read in _giao_get_u 910'
      endif
      prop_giao_get_u=.false.
      return
c
c     ----- unable to open ; stop ... -----
c
 911  call errquit('prop_giao_get_u: open failed',911)
      return
c
      end
      logical function prop_giao_put_u1(u,l1,l3)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l3
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision u(l3)
      character*255 giao_u_fil
      data ift /33/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.u',
     1     .false.,.false.,giao_u_fil)
      open(unit=ift, file=giao_u_fil, form='unformatted',
     1     access='sequential', status='unknown', err=911)
c
      if(out) then
         write(iw,*) 'giao_put_u1 = ',giao_u_fil
         write(iw,*) 'l1,l3       = ',l1,l3
      endif
      if(dbug) then
         call hnd_prsq(u,l1,l1,l1)
      endif
c
      rewind ift
      write(ift) u
      rewind ift
c
      close(ift,status='keep')
c
      prop_giao_put_u1=.true.
      return
c
 911  call errquit('prop_giao_put_u1: open failed',0)
c
      end
      logical function prop_giao_get_u1(u,l1,l3)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer l1
      integer l3
      integer i
      integer j
      integer ir
      integer iw
      integer ift
      integer ios
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision u(l3)
      character*255 giao_u_fil
      data ift /33/
c
      dbug=.false.
      dbug=dbug.or.util_print('giao_fil', print_debug)
      out =.false.
      out =out.or.dbug
c
      call util_file_name('giao.u',
     1     .false.,.false.,giao_u_fil)
      open(unit=ift, file=giao_u_fil, form='unformatted',
     1     access='sequential', status='old', err=911)
      if(out) then
         write(iw,*) 'giao_get_u1 = ',giao_u_fil
         write(iw,*) 'l1,l3       = ',l1,l3
      endif
c
      rewind ift
      read(ift,end=910,err=910,iostat=ios) u
c
      close(ift,status='keep')
c
      if(dbug) then
         call hnd_prsq(u,l1,l1,l1)
      endif
c
      prop_giao_get_u1=.true.
      return
c
c     ----- unable to complete the read ---
c
 910  close(ift,status='keep')
      if(out) then
         write(iw,*) 'iostat = ',ios
         write(iw,*) 'unable to complete the read in _giao_get_u1 910'
      endif
      prop_giao_get_u1=.false.
      return
c
c     ----- unable to open ; stop ... -----
c
 911  call errquit('prop_giao_get_u1: open failed',911)
      return
c
      end
      SUBROUTINE HND_GIAPRT(CS,NAT,IFLG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXBFN=3072)
      PARAMETER (TWO  =2.0D+00)
      PARAMETER (THREE=3.0D+00)
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      DIMENSION CS(3,3,*)
      DIMENSION A(6),AXS(3,3),EIG(3)
      DIMENSION IA(3)
C
      IA(1)=0
      IA(2)=1
      IA(3)=3
C
      IF(IFLG.EQ.1) THEN
C
C     ----- PARAMAGNETIC -----
C
         WRITE(IW,9995) 
      ELSEIF(IFLG.EQ.2) THEN
C
C     ----- DIAMAGNETIC -----
C
         WRITE(IW,9996) 
      ELSEIF(IFLG.EQ.3) THEN
C
C     ----- PARAMAGNETIC -----
C
         WRITE(IW,9997) 
      ELSEIF(IFLG.EQ.4) THEN
C
C     -----              -----
C
         WRITE(IW,9998) 
      ELSEIF(IFLG.EQ.5) THEN
C
C     ----- PARAMAGNETIC + DIAMAGNETIC -----
C
         WRITE(IW,9999) 
      ENDIF
C
C     ----- LOOP OVER ATOMS -----
C
      DO IAT=1,NAT
         WRITE(IW,9989) IAT,ANAM(IAT) 
C
         WRITE(IW,9994)
         WRITE(IW,9993) CS(1,1,IAT),CS(1,2,IAT),CS(1,3,IAT)
         WRITE(IW,9992) CS(2,1,IAT),CS(2,2,IAT),CS(2,3,IAT)
         WRITE(IW,9991) CS(3,1,IAT),CS(3,2,IAT),CS(3,3,IAT)
C
         IF(IFLG.NE.5) THEN
            AVERAG=(CS(1,1,IAT)+CS(2,2,IAT)+CS(3,3,IAT))/THREE
            WRITE(IW,9988) AVERAG
         ENDIF
C
         IF(IFLG.EQ.5) THEN
            A(1)= CS(1,1,IAT)
            A(2)=(CS(2,1,IAT)+CS(1,2,IAT))/TWO
            A(3)= CS(2,2,IAT)
            A(4)=(CS(3,1,IAT)+CS(1,3,IAT))/TWO
            A(5)=(CS(3,2,IAT)+CS(2,3,IAT))/TWO
            A(6)= CS(3,3,IAT)
            CALL HND_DIAAXS(A,AXS,EIG,IA,3,3,3)
            AVERAG=(EIG(1) + EIG(2) + EIG(3))/THREE
            ANISO = EIG(1) -(EIG(2) + EIG(3))/TWO
            WRITE(IW,9987) AVERAG,ANISO
            WRITE(IW,9986) (    J ,J=1,3)
            WRITE(IW,9985) (EIG(J),J=1,3)
            WRITE(IW,9984)
            DO I=1,3
               WRITE(IW,9983) I,(AXS(I,J),J=1,3)
            ENDDO
         ENDIF
C
      ENDDO
C
      RETURN
 9999 FORMAT(//,' SHIELDING TENSOR (PPM) = ',
     1          ' ( PARAMAGNETIC + DIAMAGNETIC ) ',
     2        /,' ------------------------ ')
 9998 FORMAT(//,' SHIELDING TENSOR (PPM) = ',
     1          ' (              P00 * H11 ) ',
     2        /,' ------------------------ ')
 9997 FORMAT(//,' SHIELDING TENSOR (PPM) = ',
     1          ' ( PARAMAGNETIC P00 * H11 ) ',
     2        /,' ------------------------ ')
 9996 FORMAT(//,' SHIELDING TENSOR (PPM) = ',
     1          ' (  DIAMAGNETIC P00 * H11 ) ',
     2        /,' ------------------------ ')
 9995 FORMAT(//,' SHIELDING TENSOR (PPM) = ',
     1          ' ( PARAMAGNETIC P10 * H01 ) ',
     2        /,' ------------------------ ')
 9994 FORMAT(/,7X,7X,'Bx',3X,7X,'By',3X,7X,'Bz')
 9993 FORMAT('      x',3F12.4)
 9992 FORMAT('      y',3F12.4)
 9991 FORMAT('      z',3F12.4)
 9990 FORMAT(/)
 9989 FORMAT(/,1X,I3,2X,A8,/,1X,'-------------')               
 9988 FORMAT(12X,' ISOTROPIC = ',F12.4)
 9987 FORMAT(12X,' ISOTROPIC = ',F12.4,/,
     1       12X,'ANISOTROPY = ',F12.4)
 9986 FORMAT(/,7X,3(7X,I1,4X))
 9985 FORMAT(' EIG = ',3F12.4)
 9984 FORMAT(' VEC = ') 
 9983 FORMAT(6X,I1,3F12.4)
      END
      SUBROUTINE HND_GIAINP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
C
C     ----- THIS ROUTINE READS THE -$GIA- NAMELIST INFORMATION -----
C           AND INITIALIZES VARIABLES CONTROLLING THE CPHF RUN.
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_CVGGIA/CNVA,CNVB,DAMPD
      COMMON/HND_PARGIA/VNEW,VOLD,ITMAX
      DATA CNV1  /1.0D-03/
      DATA CNV2  /1.0D-03/
      DATA DAMP  /0.2D+00/
      DATA MAXIT /50/
C
      NAMELIST /GIA/ MAXIT,CNV1,CNV2,DAMP
C
c--   REWIND IR
c--   READ(IR,GIA,END=10,ERR=10)
c--10 CONTINUE
C
      ITMAX =MAXIT
      CNVA  =CNV1
      CNVB  =CNV2
      DAMPD =DAMP
C
      if(ga_nodeid().eq.0) then
      WRITE(IW,9999) ITMAX,CNVA,CNVB,DAMPD
      endif
C
      RETURN
 9999 FORMAT(/,' ----- $GIA -----',/,
     1       /,' MAXIT  = ',I10  ,/,' CNVA   = ',F10.5,
     2       /,' CNVB   = ',F10.5,/,' DAMP   = ',F10.5)
      END
      SUBROUTINE HND_GIACVG_DPL(D,H,U,U0,NUM,NOC,NORB,NDIM,
     1                      IEND,ITER,SOME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- THIS ROUTINE CHECKS FOR CONVERGENCE OF THE CPHF -----
C           CALCULATION.  CONVERGENCE IS DETERMINED BY TWO
C           CRITERIA: THE CONVERGENCE OF THE TENSOR ELEMENT
C           FOR THE XX, YY, OR ZZ DIRECTION AND THE
C           CONVERGENCE OF THE INDIVIDUAL ELEMENTS OF THE U
C           MATRIX.
C
      LOGICAL  STATUS
      LOGICAL  PROP_GIAO_PUT_U
      EXTERNAL PROP_GIAO_PUT_U
      LOGICAL  PROP_GIAO_GET_U
      EXTERNAL PROP_GIAO_GET_U
      LOGICAL SOME
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (MXIODA=255)   
      COMMON/HND_CVGGIA/CNVA,CNVB,dummy_damp
      COMMON/HND_PARGIA/VNEW,VOLD,idummy_maxit_or_itmax
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      DIMENSION D(*),H(*),U(NDIM,*),U0(NDIM,*)
C
      IF(ITER.GT.0) THEN    
C
C     ----- ITER > 1 -----
C
         VOLD  = VNEW
         VNEW  =-HND_PRPDOT(D,H,NUM)
         DIFF  = VOLD - VNEW
         UDELM = ZERO
         UVAL  = ZERO
         IVAL  = 0
         JVAL  = 0
C
         NUM3=NUM*NUM
         CALL HND_DAREAD(IDAF,IODA,U0,NUM3,32)
         STATUS=PROP_GIAO_GET_U(U0,NUM,NUM3)
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIACVG_DPL PROP_GIAO_GET_U FAILED ',911)
         ENDIF
C
         ICON=0
         DO J=NOC+1,NORB
            DO I=1,NOC
               UDEL = U(I,J) - U0(I,J)
               IF(ABS(UDEL).GE.CNVB) ICON = ICON + 1
               IF(ABS(UDEL).GT.ABS(UDELM)) THEN
                  UDELM = UDEL
                  UVAL = U(I,J)
                  IVAL = I
                  JVAL = J
               ENDIF
            ENDDO
         ENDDO
         NUM3=NUM*NUM
         CALL HND_DAWRIT(IDAF,IODA,U,NUM3,32,NAV)
         STATUS=PROP_GIAO_PUT_U(U,NUM,NUM3)
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIACVG_DPL PROP_GIAO_GET_U FAILED ',911)
         ENDIF
C
         IF(ABS(DIFF).LT.CNVA.AND.ICON.EQ.0) THEN 
            IEND = 1
         ELSE
            IEND = 0
         ENDIF
         IF(SOME) THEN
            WRITE(IW,9999) ITER,VNEW,DIFF,UDELM,UVAL,IVAL,JVAL,
     1                     ICON,IEND
         ENDIF
         RETURN
C
      ELSE
C
C     ----- ITER = 1 -----
C
         VOLD  = ZERO
         VNEW  =-HND_PRPDOT(D,H,NUM)
         DIFF  = VOLD - VNEW
         UDELM = ZERO
         UVAL  = ZERO
         IVAL  = 0
         JVAL  = 0
C
         NUM3=NUM*NUM
         CALL HND_DAWRIT(IDAF,IODA,U,NUM3,32,NAV)
         STATUS=PROP_GIAO_PUT_U(U,NUM,NUM3)
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIACVG_DPL PROP_GIAO_GET_U FAILED ',911)
         ENDIF
C
         ICON = NOC * (NORB - NOC)
         IEND=0
         IF(SOME) THEN
            WRITE(IW,9999) ITER,VNEW,DIFF,UDELM,UVAL,IVAL,JVAL,
     1                     ICON,IEND
         ENDIF
      ENDIF
      RETURN
 9999 FORMAT(1X,I3,F20.9,3F16.9,2I4,I10,I3)
      END
      SUBROUTINE HND_GIAD_DPL(C,U,X,D,NUM,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- CONSTRUCT THE DERIVATIVE DENSITY MATRIX: -----
C                  D' = 2*((C'C+) + (CC'+))
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE =1.0D+00)
      PARAMETER (TWO =2.0D+00)
      DIMENSION C(NDIM,*),U(NDIM,*),X(NDIM,*),D(*)
C
      NUM2=(NUM*(NUM+1))/2
      DO I=1,NUM2
         D(I)=ZERO
      ENDDO
      DO J=1,NORB
         DO I=1,NUM
            DUM=ZERO
            DO K=1,NORB
               DUM=DUM+C(I,K)*U(K,J)
            ENDDO
            X(I,J)=DUM
         ENDDO
      ENDDO
      DO I=1,NOC
         IC=1
         DO J=1,NUM
            DO K=1,J
               D(IC) = D(IC)+C(J,I)*X(K,I)+X(J,I)*C(K,I)
               IC=IC+1
            ENDDO
         ENDDO
      ENDDO
      DO I=1,NUM2
         D(I)=D(I)*TWO
      ENDDO
      RETURN
      END
      SUBROUTINE HND_GIAU_DPL(E,C,F,U,FC,NUM,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- CONSTRUCT THE TRANSFORMATION MATRIX -U- USED -----
C           TO GET THE FIRST ORDER WAVEFUNCTION BY
C           PROJECTING OUT THE 'IJ' ELEMENT OF THE
C           DERIVATIVE FOCK MATRIX AND WEIGHTING IT
C           BY THE ENERGY DENOMINATOR.
C
      PARAMETER (ZERO=0.0D+00)
      LOGICAL   DBUG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION E(*),C(NDIM,*),F(*),U(NDIM,*),FC(NDIM,*)
C
      DBUG=.FALSE.
C
      DO J=1,NORB
         DO IS=1,NUM
            FC(IS,J)=ZERO
         ENDDO
      ENDDO
      DO J=NOC+1,NORB
         IST = 0
         DO IS=1,NUM
            CSJ = C(IS,J)
            YY = ZERO
            DO IT=1,IS
               IST = IST+1
               FC(IT,J) = FC(IT,J) + F(IST)*CSJ
               XX = F(IST)*C(IT,J)
               YY = YY + XX
            ENDDO
            FC(IS,J) = FC(IS,J) + YY - XX
         ENDDO
      ENDDO
C
C     ----- F(1) .... IF DEBUG -----
C
      IF(DBUG) THEN
         DO J=1,NORB
            DO I=1,NORB
               U(I,J)=ZERO
            ENDDO
         ENDDO
         DO I=1,NOC
            DO J=NOC+1,NORB
               YY = ZERO
               DO IS=1,NUM
                  YY = YY + C(IS,I)*FC(IS,J)
               ENDDO
               U(I,J) =  YY
               U(J,I) =  YY
            ENDDO
         ENDDO
         IF(DBUG) THEN
            WRITE(IW,*) 'F(1) IN -GIAU-'
            CALL HND_PRSQ(U,NORB,NORB,NDIM)
         ENDIF
      ENDIF
C
C     ----- U(1) -----
C
      DO J=1,NORB
         DO I=1,NORB
            U(I,J)=ZERO
         ENDDO
      ENDDO
      DO I=1,NOC
         DO J=NOC+1,NORB
            YY = ZERO
            DO IS=1,NUM
               YY = YY + C(IS,I)*FC(IS,J)
            ENDDO
            DIFF = E(I) - E(J)
            YY = YY / DIFF
            U(I,J) = -YY
            U(J,I) =  YY
         ENDDO
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,*) 'U(1) IN -GIAU-'
         CALL HND_PRSQ(U,NORB,NORB,NDIM)
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_GIAU1(E,U,S,F,H,NUM,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- CONSTRUCT THE TRANSFORMATION MATRIX -U- USED -----
C           TO GET THE FIRST ORDER WAVEFUNCTION BY
C           PROJECTING OUT THE 'IJ' ELEMENT OF THE
C           DERIVATIVE FOCK MATRIX AND WEIGHTING IT
C           BY THE ENERGY DENOMINATOR.
C
      PARAMETER (ZERO=0.0D+00)
      LOGICAL   DBUG
      LOGICAL   OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
      DIMENSION E(*)
      DIMENSION F(NDIM,*),H(NDIM,*)
      DIMENSION S(NDIM,*),U(NDIM,*)
      DATA TWO /2.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
C     ----- U(1) -----
C
      DO J=1,NORB
         DO K=1,NORB
            U(K,J)=ZERO
         ENDDO
      ENDDO
C
C     ----- READ OVERLAP DERIVATIVES -----
C
      CALL HND_GIARD(IFT1,S,NDIM*NDIM)
C
C     -----    OVERLAP      TERM FOR THE ( OCC,OCC) BLOCK -----
C
      DO J=1,NOC
         DO K=1,NOC
            U(K,J)= - S(K,J) / TWO 
         ENDDO
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,*) 'U(1) ( OCC,OCC) BLOCK IN -GIAU1-'
         CALL HND_PRSQ(U,NORB,NORB,NDIM)
         CALL HND_FIFLSH(IW)
      ENDIF
      IF(OUT) THEN
         WRITE(IW,*) '-EIG- in -GIAU1-'
         CALL HND_PRSQ(E,1,NORB,NORB)
      ENDIF
C
C     ----- H(1) -----
C
      CALL HND_GIARD(IFT2,H,NDIM*NDIM)
C
C     ----- U CONTAINS THE OVERLAP DERIVATIVE TERM -----
C
      DO J=1,NOC
         DO K=NOC+1,NORB
            DUM    = (   E(J)   - E(K)   )
            DUM1   = (   H(K,J) - F(K,J) ) / DUM
            DUM2   = ( - S(K,J) * E(J)   ) / DUM
            U(K,J) =     DUM1   + DUM2
            IF(DBUG) THEN
               WRITE(IW,9999) K,J,DUM
               WRITE(IW,9998) H(K,J),F(K,J),DUM1
               WRITE(IW,9997) S(K,J),E(J)  ,DUM2,U(K,J)
            ENDIF
         ENDDO
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,*) 'U(1) ( OCC,OCC) BLOCK IN -GIAU1-'
         WRITE(IW,*) 'U(1) (VIRT,OCC) BLOCK IN -GIAU1-'
         CALL HND_PRSQ(U,NORB,NORB,NDIM)
         CALL HND_FIFLSH(IW)
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -GIAU1- K,J = ',2I5,' E(J)-E(K) = ',F15.8)
 9998 FORMAT(' H(K,J),F(K,J),DUM1 = ',3F15.8)
 9997 FORMAT(' S(K,J),E(J)  ,DUM2 = ',3F15.8,/,
     1       '               DUM  = ',30X,F15.8)
      END
      SUBROUTINE HND_GIAV1(V0,U,V1,NUM,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      LOGICAL   OUT
      COMMON/HND_IOFILE/IR,IW
      DIMENSION V0(NDIM,*),U(NDIM,*),V1(NDIM,*)
C
      OUT=.FALSE.
C
      DO J=1,NOC
         DO I=1,NUM
            DUM=ZERO
            DO K=1,NORB
               DUM=DUM+V0(I,K)*U(K,J)
            ENDDO
            V1(I,J)=DUM
         ENDDO
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,*) '-V(1)- IN -GIAV1-'
         CALL HND_PRSQ(V1,NORB,NUM,NDIM)
      ENDIF
      RETURN
      END
      SUBROUTINE HND_GIAD1(V0,V1,D1,NUM,NOC,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- CONSTRUCT THE DERIVATIVE DENSITY MATRIX -----
C                  D' = 2*((V'V+) - (VV'+))
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (TWO =2.0D+00)
      LOGICAL   DBUG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION V0(NDIM,*),V1(NDIM,*),D1(NDIM,*)
C
      DBUG=.FALSE.
C
      DO J=1,NUM
         DO I=1,NUM
            D1(I,J)=ZERO
         ENDDO
      ENDDO
      DO IOC=1,NOC
         DO J=1,NUM
            DO I=1,NUM
               D1(I,J) = D1(I,J)
     1                 + (V0(I,IOC)*V1(J,IOC)-V1(I,IOC)*V0(J,IOC))
            ENDDO
         ENDDO
      ENDDO
      DO J=1,NUM
         DO I=1,NUM
            D1(I,J)=D1(I,J)*TWO
         ENDDO
      ENDDO
C
      IF(DBUG) THEN
         WRITE(IW,*) 'D(1) IN -GIAD1-'
         CALL HND_PRSQ(D1,NUM,NUM,NDIM)
      ENDIF
      RETURN
      END
      SUBROUTINE HND_GIAWT(IFT,H,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(NDIM)
      WRITE(IFT) H
      RETURN
      END
      SUBROUTINE HND_GIARD(IFT,H,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(NDIM)
      READ(IFT) H
      RETURN
      END
      SUBROUTINE HND_GIAFIL(IFLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- POSITION THE -GIA- FILES -----
C
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
C
      DBUG=.FALSE.
C
      CALL HND_REWFIL(IFT1)
      GO TO (30,20,10),IFLD
   10 READ(IFT1)
   20 READ(IFT1)
   30 CONTINUE
      IF(DBUG) THEN
         WRITE(IW,*) '-IFT1- POSITIONED TO -IFLD- = ',IFLD
         CALL HND_FIFLSH(IW)
      ENDIF
C
      CALL HND_REWFIL(IFT2)
      GO TO (130,120,110),IFLD
  110 READ(IFT2)
  120 READ(IFT2)
  130 CONTINUE
      IF(DBUG) THEN
         WRITE(IW,*) '-IFT2- POSITIONED TO -IFLD- = ',IFLD
         CALL HND_FIFLSH(IW)
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_GIAK00(GAUGE,F10,F00,W,T,NUM,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
      DIMENSION F10(NDIM,NDIM,*)
      DIMENSION F00(NDIM,NDIM,*)
      DIMENSION   W(NDIM,NDIM)
      DIMENSION   T(*)
      DIMENSION GAUGE(3,NOC)
      DATA ZERO    /0.0D+00/
      DATA ONE     /1.0D+00/
      DATA TWO     /2.0D+00/
      DATA CVEL    /137.0359895D+00/
C
      IA(I)=(I*(I-1))/2
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
C     ----- GAUGE DEPENDENT TERMS FIRST -----
C                  KX , KY , KZ 
C
      IF(OUT) THEN
         WRITE(IW,*) '-K00- GAUGE DEPENDENT TERMS'      
      ENDIF
      CALL HND_REWFIL(IFT2)
      DO IOP=1,3
         IF(OUT) THEN
            WRITE(IW,*) 'READ     RX, RY, RZ ... = ',IOP   
         ENDIF
         NORB2=(NORB*(NORB+1))/2
         CALL HND_GIARD(IFT2,T,NORB2)
         DO IORB=1,NORB
            DO JORB=1,IORB
               IJ=IA(IORB)+JORB
               F10(IORB,JORB,IOP)=T(IJ)
               F10(JORB,IORB,IOP)=T(IJ)
            ENDDO
         ENDDO
      ENDDO
C
      IEXCH=1
      IF(OUT) THEN
         WRITE(IW,*) 'READ K0             ... = ',IEXCH
      ENDIF
      CALL HND_REWFIL(IFT3)
      NORB2=(NORB*(NORB+1))/2
      CALL HND_GIARD(IFT3,T,NORB2)
      DO IORB=1,NORB
         DO JORB=1,IORB
            IJ=IA(IORB)+JORB
              W(IORB,JORB)=T(IJ)
              W(JORB,IORB)=T(IJ)
         ENDDO
      ENDDO
C
C     ----- WARNING ....DEACTIVATED !!!!!!
C
c     CALL HND_MXM(W(1,1),NDIM,F10(1,1,1),NDIM,F00(1,1,1),NDIM,
c    1        NORB,NORB,NORB)
c     CALL HND_MXM(W(1,1),NDIM,F10(1,1,2),NDIM,F00(1,1,2),NDIM,
c    1        NORB,NORB,NORB)
c     CALL HND_MXM(W(1,1),NDIM,F10(1,1,3),NDIM,F00(1,1,3),NDIM,
c    1        NORB,NORB,NORB)
C
C     ----- WE ARE READY ... -----
C
      CALL HND_REWFIL(IFT1)
      CALL HND_REWFIL(IFT4)
      DO IOC=1,NOC
         GX=GAUGE(1,IOC)
         GY=GAUGE(2,IOC)
         GZ=GAUGE(3,IOC)
         IF(OUT) THEN
            WRITE(IW,9999) IOC,GX,GY,GZ
         ENDIF
         CALL HND_GIARD(IFT1,F10(1,1,1),NDIM*NDIM)
         CALL HND_GIARD(IFT1,F10(1,1,2),NDIM*NDIM)
         CALL HND_GIARD(IFT1,F10(1,1,3),NDIM*NDIM)
         IF(DBUG) THEN
            WRITE(IW,*) 'B(X) PERTURBATION BEFORE -K00-'
            CALL HND_PRSQ(F10(1,1,1),NORB,NORB,NDIM)
            WRITE(IW,*) 'B(Y) PERTURBATION BEFORE -K00-'
            CALL HND_PRSQ(F10(1,1,2),NORB,NORB,NDIM)
            WRITE(IW,*) 'B(Z) PERTURBATION BEFORE -K00-'
            CALL HND_PRSQ(F10(1,1,3),NORB,NORB,NDIM)
         ENDIF
         DO IFLD=1,3
            IF(IFLD.EQ.1) THEN
               C1=GZ
               C2=GY
               IOP1=2
               IOP2=3
            ELSEIF(IFLD.EQ.2) THEN
               C1=GX
               C2=GZ
               IOP1=3
               IOP2=1
            ELSEIF(IFLD.EQ.3) THEN
               C1=GY
               C2=GX
               IOP1=1
               IOP2=3
            ENDIF 
            VELFAC=ONE/(TWO*CVEL)
            DO JORB=1,NORB
               DO IORB=1,NORB
                  F10(IORB,JORB,IFLD)=F10(IORB,JORB,IFLD)
     1      -VELFAC*( C1*(F00(IORB,JORB,IOP1)-F00(JORB,IORB,IOP1))
     2               -C2*(F00(IORB,JORB,IOP2)-F00(JORB,IORB,IOP2)) )
               ENDDO
            ENDDO
         ENDDO
         CALL HND_GIAWT(IFT4,F10(1,1,1),NDIM*NDIM)
         CALL HND_GIAWT(IFT4,F10(1,1,2),NDIM*NDIM)
         CALL HND_GIAWT(IFT4,F10(1,1,3),NDIM*NDIM)
         IF(DBUG) THEN
            WRITE(IW,*) 'B(X) PERTURBATION WITH -K00- GAUGE DEPENDENT'
            CALL HND_PRSQ(F10(1,1,1),NORB,NORB,NDIM)
            WRITE(IW,*) 'B(Y) PERTURBATION WITH -K00- GAUGE DEPENDENT'
            CALL HND_PRSQ(F10(1,1,2),NORB,NORB,NDIM)
            WRITE(IW,*) 'B(Z) PERTURBATION WITH -K00- GAUGE DEPENDENT'
            CALL HND_PRSQ(F10(1,1,3),NORB,NORB,NDIM)
         ENDIF
      ENDDO
      CALL HND_REWFIL(IFT1)
      CALL HND_REWFIL(IFT4)
C
C     ----- GAUGE INDEPENDENT TERMS NOW -----
C
      IF(OUT) THEN
         WRITE(IW,*) '-K00- GAUGE INDEPENDENT TERMS'      
      ENDIF
      CALL HND_REWFIL(IFT3)
      READ(IFT3)
      DO IOP=2,4
         IF(OUT) THEN
            WRITE(IW,*) 'READ     KX, KY, KZ ... = ',IOP  
         ENDIF
         NORB2=(NORB*(NORB+1))/2
         CALL HND_GIARD(IFT3,T,NORB2)
         DO IORB=1,NORB
            DO JORB=1,IORB
               IJ=IA(IORB)+JORB
               F00(IORB,JORB,IOP)=T(IJ)
               F00(JORB,IORB,IOP)=T(IJ)
            ENDDO
         ENDDO
      ENDDO
C
C     -----  X ,  Y ,  Z -----
C
      CALL HND_REWFIL(IFT2)
      DO IOP=1,3
         IF(OUT) THEN
            WRITE(IW,*) 'READ     RX, RY, RZ ... = ',IOP   
         ENDIF
         NORB2=(NORB*(NORB+1))/2
         CALL HND_GIARD(IFT2,T,NORB2)
         DO IORB=1,NORB
            DO JORB=1,IORB
               IJ=IA(IORB)+JORB
               F10(IORB,JORB,IOP)=T(IJ)
               F10(JORB,IORB,IOP)=T(IJ)
            ENDDO
         ENDDO
      ENDDO
C
C     ----- WE ARE READY ... -----
C
      CALL HND_REWFIL(IFT1)
      CALL HND_REWFIL(IFT4)
      DO IOC=1,NOC
         GX=GAUGE(1,IOC)
         GY=GAUGE(2,IOC)
         GZ=GAUGE(3,IOC)
         IF(OUT) THEN
            WRITE(IW,9999) IOC,GX,GY,GZ
         ENDIF
         DO IFLD=1,3
            IF(IFLD.EQ.1) THEN
               IOP1=3
               IOP2=2
            ELSEIF(IFLD.EQ.2) THEN
               IOP1=1
               IOP2=3
            ELSEIF(IFLD.EQ.3) THEN
               IOP1=2
               IOP2=1
            ENDIF
            CALL HND_GIARD(IFT4,W(1,1),NDIM*NDIM)
            DO JORB=1,NORB
               DO IORB=1,NORB
                  DUM=ZERO
                  DO LORB=1,NORB
                     DUM=DUM+F00(IORB,LORB,IOP1)*F10(LORB,JORB,IOP2)
     1                      -F00(JORB,LORB,IOP1)*F10(LORB,IORB,IOP2)
     2                      -F00(IORB,LORB,IOP2)*F10(LORB,JORB,IOP1)
     3                      +F00(JORB,LORB,IOP2)*F10(LORB,IORB,IOP1)
                  ENDDO
                  VELFAC=ONE/(TWO*CVEL)
                  W(IORB,JORB)=W(IORB,JORB)+VELFAC*DUM
               ENDDO
            ENDDO
            CALL HND_GIAWT(IFT1,W(1,1),NDIM*NDIM)
            IF(DBUG) THEN
               IF(IFLD.EQ.1) THEN
            WRITE(IW,*) 'B(X) PERTURBATION WITH -K00- GAUGE INDEPEND.'
               ELSEIF(IFLD.EQ.2) THEN
            WRITE(IW,*) 'B(Y) PERTURBATION WITH -K00- GAUGE INDEPEND.'
               ELSEIF(IFLD.EQ.3) THEN
            WRITE(IW,*) 'B(Z) PERTURBATION WITH -K00- GAUGE INDEPEND.'
               ENDIF
               CALL HND_PRSQ(W(1,1),NORB,NORB,NDIM)
            ENDIF
         ENDDO
      ENDDO
      CALL HND_REWFIL(IFT1)
      CALL HND_REWFIL(IFT4)
C
      RETURN
 9999 FORMAT(' IN -GIAK00- , GAUGE FOR -IOC- = ',I5,
     1       ' GX,GY,GZ = ',3F15.9)
      END
      SUBROUTINE HND_GIAS10(S10,C,S00,T,NUM,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
      DIMENSION S10(NDIM*NDIM,*)
      DIMENSION C(NDIM,*)
      DIMENSION S00(*)
      DIMENSION T(*)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
      DATA TWO  /2.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      IF(OUT) THEN
         WRITE(IW,*) 'START -GIAS10-'
         CALL HND_FIFLSH(IW)
      ENDIF
      DO IOP=1,3
         DO IJ=1,NDIM*NDIM
            S10(IJ,IOP)=ZERO
         ENDDO
      ENDDO
C
C     ----- OVERLAP DERIVATIVES -S10- -----
C
      CALL HND_S10INT(S10(1,1),S10(1,2),S10(1,3),S00,NDIM)
C
C     ----- TRANSFORM TO CANONICAL BASIS -----
C
      DO IOP=1,3
         IF(DBUG) THEN
            WRITE(IW,*) '-S10- OPERATORS IN -AO- BASIS'
            CALL HND_PRSQ(S10(1,IOP),NUM,NUM,NDIM)
         ENDIF
         CALL HND_TFSQ(S10(1,IOP),S10(1,IOP),C,T,NORB,NUM,NDIM)
         IF(OUT) THEN
            WRITE(IW,*) '-S10- OPERATORS IN -MO- BASIS'
            CALL HND_PRSQ(S10(1,IOP),NUM,NUM,NDIM)
         ENDIF
      ENDDO
C
C     ----- WRITE OUT -----
C
      CALL HND_REWFIL(IFT1)
      CALL HND_GIAWT(IFT1,S10(1,1),NDIM*NDIM)
      CALL HND_GIAWT(IFT1,S10(1,2),NDIM*NDIM)
      CALL HND_GIAWT(IFT1,S10(1,3),NDIM*NDIM)
      CALL HND_REWFIL(IFT1)
C
      IF(OUT) THEN
         WRITE(IW,*) 'END -GIAS10-'
         CALL HND_FIFLSH(IW)
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_S10INT(S10X,S10Y,S10Z,S00,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXPRIM=2048)
      LOGICAL NORM
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ,CX,CY,CZ
      DIMENSION S10X(NDIM,*),S10Y(NDIM,*),S10Z(NDIM,*)
      DIMENSION  S00(NDIM,*)
      DIMENSION SX(225),SY(225),SZ(225)
      DIMENSION  S(225)
      DIMENSION DIJ(225)
      DIMENSION  XS(6,5), YS(6,5), ZS(6,5)
      DIMENSION XXS(6,5),YYS(6,5),ZZS(6,5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0E+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     ----- ISHELL -----
C
      DO II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITMOD=LIT+1
C
C     ----- JSHELL -----
C
      DO JJ=1,NSHELL
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      QIJX=YI*ZJ-ZI*YJ
      QIJY=ZI*XJ-XI*ZJ
      QIJZ=XI*YJ-YI*XJ
      TIJX=XI-XJ
      TIJY=YI-YJ
      TIJZ=ZI-ZJ
C
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            S (IJ)=ZERO
            SX(IJ)=ZERO
            SY(IJ)=ZERO
            SZ(IJ)=ZERO
         ENDDO    
      ENDDO    
C
C     ----- I PRIMITIVE -----
C
      DO IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.LE.TOL) THEN       
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 DIJ(IJ)=DUM2
C
C     ----- X , Y , Z INTEGRALS -----
C
      T = SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO J=1,LJT
         NJ=J
         DO I=1,LITMOD
            NI=I
            CALL HND_SXYZ
             XS(I,J)=XINT*T
             YS(I,J)=YINT*T 
             ZS(I,J)=ZINT*T 
         ENDDO
      ENDDO
C
      DO J=1,LJT
         DO I=1,LIT
            XXS(I,J)=XS(I+1,J)
            YYS(I,J)=YS(I+1,J)
            ZZS(I,J)=ZS(I+1,J)
         ENDDO
      ENDDO
C
      IJ=0
      DO I=MINI,MAXI
         IX=IJX(I)
         IY=IJY(I)
         IZ=IJZ(I)
         DO J=MINJ,MAXJ
            JX=IJX(J)
            JY=IJY(J)
            JZ=IJZ(J)
            IJ=IJ+1
            DUM = XS(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
            DUMX=XXS(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
            DUMY= XS(IX,JX)*YYS(IY,JY)* ZS(IZ,JZ)
            DUMZ= XS(IX,JX)* YS(IY,JY)*ZZS(IZ,JZ)
         S (IJ)=S (IJ)+DIJ(IJ)*(     DUM                    )
         SX(IJ)=SX(IJ)+DIJ(IJ)*(QIJX*DUM+TIJY*DUMZ-TIJZ*DUMY)
         SY(IJ)=SY(IJ)+DIJ(IJ)*(QIJY*DUM+TIJZ*DUMX-TIJX*DUMZ)
         SZ(IJ)=SZ(IJ)+DIJ(IJ)*(QIJZ*DUM+TIJX*DUMY-TIJY*DUMX)
         ENDDO
      ENDDO
C
            ENDIF
         ENDDO
      ENDDO
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            ILOC=LOCI+I
            JLOC=LOCJ+J
            S00 (ILOC,JLOC)=S (IJ)
            S10X(ILOC,JLOC)=SX(IJ)
            S10Y(ILOC,JLOC)=SY(IJ)
            S10Z(ILOC,JLOC)=SZ(IJ)
         ENDDO
      ENDDO
C
         ENDDO
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,*) '-S00 - IN AO BASIS'
         CALL HND_PRSQ(S00 ,NUM,NUM,NDIM)
         WRITE(IW,*) '-S10X- IN AO BASIS'
         CALL HND_PRSQ(S10X,NUM,NUM,NDIM)
         WRITE(IW,*) '-S10Y- IN AO BASIS'
         CALL HND_PRSQ(S10Y,NUM,NUM,NDIM)
         WRITE(IW,*) '-S10Z- IN AO BASIS'
         CALL HND_PRSQ(S10Z,NUM,NUM,NDIM)
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_GIAH10(H10,C00,F00,D00,T,H00,
     1                  NUM,NOC,NORB,NDIM,NEXT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
      DIMENSION H10(NDIM*NDIM,*)
      DIMENSION C00(NDIM,*)
      DIMENSION F00(*)
      DIMENSION D00(*)
      DIMENSION T(*)
      DIMENSION H00(*)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
      DATA TWO  /2.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      IF(OUT) THEN
         WRITE(IW,*) 'START -GIAH10-'
         IF(DBUG) THEN
            WRITE(IW,*) 'D00...'
            CALL HND_PRTR(D00,NUM)
         ENDIF
         CALL HND_FIFLSH(IW)
      ENDIF
      DO IOP=1,6
         DO IJ=1,NDIM*NDIM
            H10(IJ,IOP)=ZERO
         ENDDO
      ENDDO
C
C     ----- L OPERATOR -----
C
      CALL HND_GIAL10(H10(1,1),H10(1,2),H10(1,3),H00,NDIM)
      CALL HND_FIFLSH(IW)
C
C     ----- T + V OPERATORS -----
C
      CALL HND_GITV10(H10(1,4),H10(1,5),H10(1,6),H00,NDIM)
C
C     ----- SUM UP INTO 1-E OPERATORS -----
C
      DO IOP=1,3
         DO IJ=1,NDIM*NDIM
            H10(IJ,IOP)=H10(IJ,IOP)+H10(IJ,IOP+3)
         ENDDO
      ENDDO
      IF(DBUG) THEN
         DO IOP=1,3
            WRITE(IW,*) '1-E       OPERATORS IN -AO- BASIS'
            CALL HND_PRSQ(H10(1,IOP),NUM,NUM,NDIM)
            CALL HND_FIFLSH(IW)
         ENDDO
      ENDIF
C
C     ----- WRITE OUT -----
C
      CALL HND_REWFIL(IFT2)
      CALL HND_GIAWT(IFT2,H10(1,1),NDIM*NDIM)
      CALL HND_GIAWT(IFT2,H10(1,2),NDIM*NDIM)
      CALL HND_GIAWT(IFT2,H10(1,3),NDIM*NDIM)
      CALL HND_REWFIL(IFT2)
C
C     ----- NOW 2-E OPERATORS -----
C
      IF(DBUG) THEN
         WRITE(IW,*) 'D00...'
         CALL HND_PRTR(D00,NUM)
      ENDIF
      CALL HND_GIAG10(H10(1,4),F00,D00,NDIM,NEXT)
C
C     ----- GLOBAL SUM FOR PARALLEL PROCESSING -----
C
      DO IOP=4,6
c--      CALL HND_LCPADD(H10(1,IOP),H00,NDIM*NDIM)
         CALL GA_DGOP(2468,H10(1,IOP),NDIM*NDIM,'+')
      ENDDO
      IF(DBUG) THEN
         DO IOP=4,6
            WRITE(IW,*) '      2-E OPERATORS IN -AO- BASIS'
            CALL HND_PRSQ(H10(1,IOP),NUM,NUM,NDIM)
            CALL HND_FIFLSH(IW)
         ENDDO
      ENDIF
C
C     ----- TRANSFORM TO CANONICAL BASIS -----
C
      DO IOP=1,3
         DO IJ=1,NDIM*NDIM
            H10(IJ,IOP)=(H10(IJ,IOP)+H10(IJ,IOP+3))
         ENDDO
         IF(OUT) THEN
            WRITE(IW,*) '1-E + 2-E OPERATORS IN -AO- BASIS'
            CALL HND_PRSQ(H10(1,IOP),NUM,NUM,NDIM)
            CALL HND_FIFLSH(IW)
         ENDIF
         CALL HND_TFSQ(H10(1,IOP),H10(1,IOP),C00,T,NORB,NUM,NDIM)
         IF(OUT) THEN
            WRITE(IW,*) '1-E + 2-E OPERATORS IN -MO- BASIS'
            CALL HND_PRSQ(H10(1,IOP),NUM,NUM,NDIM)
            CALL HND_FIFLSH(IW)
         ENDIF
      ENDDO
C
C     ----- WRITE OUT -----
C
      CALL HND_REWFIL(IFT2)
      CALL HND_GIAWT(IFT2,H10(1,1),NDIM*NDIM)
      CALL HND_GIAWT(IFT2,H10(1,2),NDIM*NDIM)
      CALL HND_GIAWT(IFT2,H10(1,3),NDIM*NDIM)
      CALL HND_REWFIL(IFT2)
      CALL HND_FIFLSH(IW)
C
      IF(OUT) THEN
         WRITE(IW,*) 'END -GIAH10-'
         CALL HND_FIFLSH(IW)
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_GIAL10(H10X,H10Y,H10Z,S00,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXPRIM=2048)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      DIMENSION H10X(NDIM,*),H10Y(NDIM,*),H10Z(NDIM,*)
      DIMENSION S00 (NDIM,*)
      DIMENSION XL(225),YL(225),ZL(225)
      DIMENSION SS(225)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0E+00/
C
      OUT=.FALSE.
      IF(OUT) THEN
         WRITE(IW,*) 'START -GIAL10-'
         CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- ISHELL -----
C
      DO II=1,NSHELL
         MINI=KMIN(II)
         MAXI=KMAX(II)
         LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
         DO JJ=1,NSHELL
            MINJ=KMIN(JJ)
            MAXJ=KMAX(JJ)
            LOCJ=KLOC(JJ)-MINJ
C
C     ----- L INTEGRALS -----
C
            CALL HND_L10INT(II,JJ,XL,YL,ZL,SS)
C
C     ----- L OPERATORS -----
C
            IJ=0
            DO I=MINI,MAXI
               DO J=MINJ,MAXJ 
                  ILOC=LOCI+I
                  JLOC=LOCJ+J
                  IJ=IJ+1
                  S00 (ILOC,JLOC)= SS(IJ)
                  H10X(ILOC,JLOC)=-XL(IJ)
                  H10Y(ILOC,JLOC)=-YL(IJ)
                  H10Z(ILOC,JLOC)=-ZL(IJ)
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
      IF(OUT) THEN        
         WRITE(IW,*) '- S00    - OPERATOR IN AO BASIS'
         CALL HND_PRSQ(S00 ,NUM,NUM,NDIM)
         WRITE(IW,*) '- L10 -X - OPERATOR IN AO BASIS'
         CALL HND_PRSQ(H10X,NUM,NUM,NDIM)
         CALL HND_FIFLSH(IW)
         WRITE(IW,*) '- L10 -Y - OPERATOR IN AO BASIS'
         CALL HND_PRSQ(H10Y,NUM,NUM,NDIM)
         CALL HND_FIFLSH(IW)
         WRITE(IW,*) '- L10 -Z - OPERATOR IN AO BASIS'
         CALL HND_PRSQ(H10Z,NUM,NUM,NDIM)
         CALL HND_FIFLSH(IW)
      ENDIF
      IF(OUT) THEN
         WRITE(IW,*) 'END   -GIAL10-'
         CALL HND_FIFLSH(IW)
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_L10INT(II,JJ,XL,YL,ZL,SS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXPRIM=2048)
      LOGICAL DBUG
      LOGICAL OUT
      LOGICAL NORM
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ,CX,CY,CZ
      DIMENSION XL(*),YL(*),ZL(*),SS(*)
      DIMENSION CIJ(225)
      DIMENSION  XS(5,6), YS(5,6), ZS(5,6)
      DIMENSION XSX(5,6),YSY(5,6),ZSZ(5,6)
      DIMENSION XSD(5,6),YSD(5,6),ZSD(5,6)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0E+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(OUT) THEN
         WRITE(IW,*) 'IN    -H10INT-',II,JJ
         CALL HND_FIFLSH(IW)
      ENDIF
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     ----- ISHELL -----
C
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      LJTMOD=LJT+1
C
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            SS(IJ)=ZERO
            XL(IJ)=ZERO
            YL(IJ)=ZERO
            ZL(IJ)=ZERO
         ENDDO
      ENDDO
C
C     ----- I PRIMITIVE -----
C
      DO IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.LE.TOL) THEN        
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 CIJ(IJ)=DUM2
C
C     ----- INTEGRAL COMPONENTS -----
C
      T = SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO J=1,LJTMOD
         NJ=J
         DO I=1,LIT
            NI=I
            CALL HND_SXYZ
             XS(I,J)=XINT*T
             YS(I,J)=YINT*T
             ZS(I,J)=ZINT*T
         ENDDO
      ENDDO
C
      DO J=1,LJT
         DO I=1,LIT
            XSX(I,J)=XS(I,J+1)
            YSY(I,J)=YS(I,J+1)
            ZSZ(I,J)=ZS(I,J+1)
         ENDDO
      ENDDO
C
      CALL HND_L10XYZ(XSD,YSD,ZSD,XS,YS,ZS,LIT,LJT,AJ)
C
C     ----- L INTEGRALS -----
C
      IJ=0
      DO I=MINI,MAXI
         IX=IJX(I)
         IY=IJY(I)
         IZ=IJZ(I)
         DO J=MINJ,MAXJ
            JX=IJX(J)
            JY=IJY(J)
            JZ=IJZ(J)
            IJ=IJ+1
            DUM =  XS(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
            DUMX=  XS(IX,JX)*YSY(IY,JY)*ZSD(IZ,JZ)
     1           - XS(IX,JX)*YSD(IY,JY)*ZSZ(IZ,JZ)
            DUMY= XSD(IX,JX)* YS(IY,JY)*ZSZ(IZ,JZ)
     1           -XSX(IX,JX)* YS(IY,JY)*ZSD(IZ,JZ)
            DUMZ= XSX(IX,JX)*YSD(IY,JY)* ZS(IZ,JZ)
     1           -XSD(IX,JX)*YSY(IY,JY)* ZS(IZ,JZ)
            SS(IJ)=SS(IJ)+CIJ(IJ)*DUM 
            XL(IJ)=XL(IJ)+CIJ(IJ)*DUMX
            YL(IJ)=YL(IJ)+CIJ(IJ)*DUMY
            ZL(IJ)=ZL(IJ)+CIJ(IJ)*DUMZ
         ENDDO
      ENDDO
C
            ENDIF
         ENDDO
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,*) 'OUT   -H10INT-'
         CALL HND_FIFLSH(IW)
      ENDIF
      RETURN
      END
      SUBROUTINE HND_L10XYZ(DX,DY,DZ,X,Y,Z,NI,NJ,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  X(5,*), Y(5,*), Z(5,*)
      DIMENSION DX(5,*),DY(5,*),DZ(5,*)
C
      DO I=1,NI
         DX(I,1)= (-(AJ+AJ)*X(I,2))
         DY(I,1)= (-(AJ+AJ)*Y(I,2))
         DZ(I,1)= (-(AJ+AJ)*Z(I,2))
      ENDDO
C
      IF(NJ.EQ.1) RETURN
C
      DO J=2,NJ
         DO I=1,NI
            DX(I,J)= (DBLE(J-1)*X(I,J-1)-(AJ+AJ)*X(I,J+1))
            DY(I,J)= (DBLE(J-1)*Y(I,J-1)-(AJ+AJ)*Y(I,J+1))
            DZ(I,J)= (DBLE(J-1)*Z(I,J-1)-(AJ+AJ)*Z(I,J+1))
         ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_GITV10(H10X,H10Y,H10Z,H00,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXPRIM=2048)
      CHARACTER*8 ERRMSG
      LOGICAL NORM
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ,CCX,CCY,CCZ
      COMMON/HND_RYS/XX,U(5),W(5),NROOTS
      DIMENSION H10X(NDIM,*),H10Y(NDIM,*),H10Z(NDIM,*)
      DIMENSION H00 (NDIM,*)
      DIMENSION CIJ(225)
      DIMENSION G(225)
      DIMENSION V(225)
      DIMENSION  GX(225),GY(225),GZ(225)
      DIMENSION  VX(225),VY(225),VZ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION  XS(6,7)  , YS(6,7)  , ZS(6,7)
      DIMENSION  XG(6,5)  , YG(6,5)  , ZG(6,5)
      DIMENSION  XV(6,5,5), YV(6,5,5), ZV(6,5,5)
      DIMENSION XXS(6,7)  ,YYS(6,7)  ,ZZS(6,7)
      DIMENSION XXG(6,5)  ,YYG(6,5)  ,ZZG(6,5)
      DIMENSION XXV(6,5,5),YYV(6,5,5),ZZV(6,5,5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-GITV10-'/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(OUT) THEN
         WRITE(IW,*) 'START -GITV10-'
         CALL HND_FIFLSH(IW)
      ENDIF
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     ----- ISHELL -----
C
      DO II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITMOD=LIT+1
C
C     ----- JSHELL -----
C
      DO JJ=1,NSHELL
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      LJTMOD=LJT+2
C
      QIJX=YI*ZJ-ZI*YJ
      QIJY=ZI*XJ-XI*ZJ
      QIJZ=XI*YJ-YI*XJ
      TIJX=XI-XJ
      TIJY=YI-YJ
      TIJZ=ZI-ZJ
C
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
      NROOTS=(1+LIT+LJT-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            G (IJ)=ZERO
            GX(IJ)=ZERO
            GY(IJ)=ZERO
            GZ(IJ)=ZERO
            V (IJ)=ZERO
            VX(IJ)=ZERO
            VY(IJ)=ZERO
            VZ(IJ)=ZERO
         ENDDO
      ENDDO
C
C     ----- I PRIMITIVE -----
C
      DO IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.LE.TOL) THEN       
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 CIJ(IJ)=DUM2
C
C     ----- KINETIC ENERGY -----
C
      T = SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      DO J=1,LJTMOD
         NJ=J
         DO I=1,LITMOD
            NI=I
            CALL HND_SXYZ
            XS(I,J)=XINT*T
            YS(I,J)=YINT*T
            ZS(I,J)=ZINT*T
         ENDDO
      ENDDO
C
      DO J=1,LJTMOD
         DO I=1,LIT
            XXS(I,J)=XS(I+1,J)
            YYS(I,J)=YS(I+1,J)
            ZZS(I,J)=ZS(I+1,J)
         ENDDO
      ENDDO
C
      CALL HND_T10XYZ( XG, YG, ZG, XS, YS, ZS,LIT,LJT,AJ)
      CALL HND_T10XYZ(XXG,YYG,ZZG,XXS,YYS,ZZS,LIT,LJT,AJ)
C
      IJ=0
      DO I=MINI,MAXI
         IX=IJX(I)
         IY=IJY(I)
         IZ=IJZ(I)
         DO J=MINJ,MAXJ
            JX=IJX(J)
            JY=IJY(J)
            JZ=IJZ(J)
            IJ=IJ+1
            DUM = XG(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
     1          + XS(IX,JX)* YG(IY,JY)* ZS(IZ,JZ)
     2          + XS(IX,JX)* YS(IY,JY)* ZG(IZ,JZ)
            DUMX=XXG(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
     1          +XXS(IX,JX)* YG(IY,JY)* ZS(IZ,JZ)
     2          +XXS(IX,JX)* YS(IY,JY)* ZG(IZ,JZ)
            DUMY= XG(IX,JX)*YYS(IY,JY)* ZS(IZ,JZ)
     1          + XS(IX,JX)*YYG(IY,JY)* ZS(IZ,JZ)
     2          + XS(IX,JX)*YYS(IY,JY)* ZG(IZ,JZ)
            DUMZ= XG(IX,JX)* YS(IY,JY)*ZZS(IZ,JZ)
     1          + XS(IX,JX)* YG(IY,JY)*ZZS(IZ,JZ)
     2          + XS(IX,JX)* YS(IY,JY)*ZZG(IZ,JZ)
      G (IJ)=G (IJ)+CIJ(IJ)*(     DUM                    )
      GX(IJ)=GX(IJ)+CIJ(IJ)*(QIJX*DUM+TIJY*DUMZ-TIJZ*DUMY)
      GY(IJ)=GY(IJ)+CIJ(IJ)*(QIJY*DUM+TIJZ*DUMX-TIJX*DUMZ)
      GZ(IJ)=GZ(IJ)+CIJ(IJ)*(QIJZ*DUM+TIJX*DUMY-TIJY*DUMX)
         ENDDO
      ENDDO
C
C     ----- NUCLEAR ATTRACTION -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO IC=1,NAT
      ZNUC=-ZAN(IC)
      CX=C(1,IC)
      CY=C(2,IC)
      CZ=C(3,IC)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO IROOT=1,NROOTS
         UU=U(IROOT)*AA
         WW=W(IROOT)*ZNUC
         TT=ONE/(AA+UU)
         T = SQRT(TT)
         X0=(AAX+UU*CX)*TT
         Y0=(AAY+UU*CY)*TT
         Z0=(AAZ+UU*CZ)*TT
         DO J=1,LJT
            NJ=J
            DO I=1,LITMOD
               NI=I
               CALL HND_SXYZ
               XV(I,J,IROOT)=XINT
               YV(I,J,IROOT)=YINT
               ZV(I,J,IROOT)=ZINT*WW
            ENDDO
         ENDDO
      ENDDO
C
      DO IROOT=1,NROOTS
         DO J=1,LJT
            DO I=1,LIT
               XXV(I,J,IROOT)=XV(I+1,J,IROOT)
               YYV(I,J,IROOT)=YV(I+1,J,IROOT)
               ZZV(I,J,IROOT)=ZV(I+1,J,IROOT)
            ENDDO
         ENDDO
      ENDDO
C
      IJ=0
      DO I=MINI,MAXI
         IX=IJX(I)
         IY=IJY(I)
         IZ=IJZ(I)
         DO J=MINJ,MAXJ
            JX=IJX(J)
            JY=IJY(J)
            JZ=IJZ(J)
            DUM =ZERO
            DUMX=ZERO
            DUMY=ZERO
            DUMZ=ZERO
            DO IROOT=1,NROOTS
      DUM =DUM + XV(IX,JX,IROOT)* YV(IY,JY,IROOT)* ZV(IZ,JZ,IROOT)
      DUMX=DUMX+XXV(IX,JX,IROOT)* YV(IY,JY,IROOT)* ZV(IZ,JZ,IROOT)
      DUMY=DUMY+ XV(IX,JX,IROOT)*YYV(IY,JY,IROOT)* ZV(IZ,JZ,IROOT)
      DUMZ=DUMZ+ XV(IX,JX,IROOT)* YV(IY,JY,IROOT)*ZZV(IZ,JZ,IROOT)
            ENDDO
            DUM =DUM *(AA1*PI212)
            DUMX=DUMX*(AA1*PI212)
            DUMY=DUMY*(AA1*PI212)
            DUMZ=DUMZ*(AA1*PI212)
            IJ=IJ+1
      V (IJ)=V (IJ)+CIJ(IJ)*(     DUM                    )
      VX(IJ)=VX(IJ)+CIJ(IJ)*(QIJX*DUM+TIJY*DUMZ-TIJZ*DUMY)
      VY(IJ)=VY(IJ)+CIJ(IJ)*(QIJY*DUM+TIJZ*DUMX-TIJX*DUMZ)
      VZ(IJ)=VZ(IJ)+CIJ(IJ)*(QIJZ*DUM+TIJX*DUMY-TIJY*DUMX)
         ENDDO
      ENDDO
C
C     ----- END LOOP OVER ATOMS -----
C
      ENDDO
C
C     ----- END LOOP OVER PRIMITIVES -----
C
            ENDIF
         ENDDO
      ENDDO
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            ILOC=LOCI+I
            JLOC=LOCJ+J
            H00 (ILOC,JLOC)=G (IJ)+V (IJ)
            H10X(ILOC,JLOC)=GX(IJ)+VX(IJ)
            H10Y(ILOC,JLOC)=GY(IJ)+VY(IJ)
            H10Z(ILOC,JLOC)=GZ(IJ)+VZ(IJ)
         ENDDO
      ENDDO
C
C     ----- END LOOP OVER SHELLS -----
C
         ENDDO
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,*) '- H00    - IN AO BASIS'
         CALL HND_PRSQ(H00 ,NUM,NUM,NDIM)
         WRITE(IW,*) '- TV10-X - IN AO BASIS'
         CALL HND_PRSQ(H10X,NUM,NUM,NDIM)
         WRITE(IW,*) '- TV10-Y - IN AO BASIS'
         CALL HND_PRSQ(H10Y,NUM,NUM,NDIM)
         WRITE(IW,*) '- TV10-Z - IN AO BASIS'
         CALL HND_PRSQ(H10Z,NUM,NUM,NDIM)
         CALL HND_FIFLSH(IW)
         WRITE(IW,*) 'END   -GITV10-'
         CALL HND_FIFLSH(IW)
      ENDIF
      RETURN
 9997 FORMAT(' IN -GITV10- THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I2,/,' LIT,LJT,NROOTS= ',3I3)
      END
      SUBROUTINE HND_T10XYZ(XT,YT,ZT,XS,YS,ZS,NI,NJ,AJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XT(6,*),YT(6,*),ZT(6,*)
      DIMENSION XS(6,*),YS(6,*),ZS(6,*)
C
      DO I=1,NI
         XT(I,1)=(XS(I,1  )            -XS(I,3  )*(AJ+AJ))*AJ
         YT(I,1)=(YS(I,1  )            -YS(I,3  )*(AJ+AJ))*AJ
         ZT(I,1)=(ZS(I,1  )            -ZS(I,3  )*(AJ+AJ))*AJ
      ENDDO     
C
      IF(NJ.EQ.1) RETURN
C
      DO I=1,NI
         XT(I,2)=(XS(I,2  )*DBLE(2+2-1)-XS(I,4  )*(AJ+AJ))*AJ
         YT(I,2)=(YS(I,2  )*DBLE(2+2-1)-YS(I,4  )*(AJ+AJ))*AJ
         ZT(I,2)=(ZS(I,2  )*DBLE(2+2-1)-ZS(I,4  )*(AJ+AJ))*AJ
      ENDDO
C
      IF(NJ.EQ.2) RETURN
C
      DO J=3,NJ
         DO I=1,NI
         XT(I,J)=(XS(I,J  )*DBLE(J+J-1)-XS(I,J+2)*(AJ+AJ))*AJ
     1           -XS(I,J-2)*DBLE(((J-1)*(J-2))/2)
         YT(I,J)=(YS(I,J  )*DBLE(J+J-1)-YS(I,J+2)*(AJ+AJ))*AJ
     1           -YS(I,J-2)*DBLE(((J-1)*(J-2))/2)
         ZT(I,J)=(ZS(I,J  )*DBLE(J+J-1)-ZS(I,J+2)*(AJ+AJ))*AJ
     1           -ZS(I,J-2)*DBLE(((J-1)*(J-2))/2)
         ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_GIAH20(H20,C,B20,T,NUM,NOC,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      DIMENSION H20(NDIM*NDIM,*)
      DIMENSION C(NDIM,*)
      DIMENSION T(*)
      DIMENSION CM(3)     
      DIMENSION RR(6)     
      DIMENSION B20(6)
      DATA CVEL  /137.0359895D+00/
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA TWO   /2.0D+00/
      DATA EIGHT /8.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      VELFAC=ONE/(EIGHT*CVEL*CVEL)
C
C     ----- R**2 OPERATORS : XX YY ZZ XY XZ YZ -----
C         FOR SIMPLICITY, RECALCULATE FOR EACH MO
C
      DO IOP=1,6
         RR(IOP)=ZERO 
      ENDDO
C
      IF(OUT) THEN
          WRITE(IW,*) 'IN -GIAH20-'
          CALL HND_PRSQ(C,NOC,NUM,NDIM)
      ENDIF
C
      DO IOC=1,NOC
C
C     ----- R**2 INTEGRALS : XX YY ZZ XY XZ YZ -----
C
         IF(OUT) THEN
            WRITE(IW,*) 'CALLING -QDPINT-'
         ENDIF
         CM(1)=ZERO
         CM(2)=ZERO
         CM(3)=ZERO
         do j=1,6     
            do i=1,ndim*ndim
               h20(i,j)=zero
            enddo
         enddo
         CALL HND_QDPINT(H20(1,1),H20(1,2),H20(1,3),
     1                   H20(1,4),H20(1,5),H20(1,6),CM,NDIM)
         IF(OUT.AND.IOC.EQ.1) THEN
            DO IOP=1,6
               WRITE(IW,*) 'R**2 OPERATOR, IOP,IOC = ',IOP,IOC
               CALL HND_PRTR(H20(1,IOP),NUM)
            ENDDO
         ENDIF
C
C     ----- SUSCEPTIBILITY CONTRIBUTION -----
C
         DO IOP=1,6
            DO I=1,NUM
               DUM=ZERO
               DO J=1,NUM
                  II=MAX(I,J)
                  JJ=MIN(I,J)
                  IJ=(II*(II-1))/2+JJ
                  DUM=DUM+H20(IJ,IOP)*C(J,IOC)
               ENDDO
               T(I)=DUM
            ENDDO
            DUM=ZERO
            DO I=1,NUM
               DUM=DUM+C(I,IOC)*T(I)
            ENDDO
            RR(IOP)=RR(IOP)+DUM 
         ENDDO
      ENDDO
C
      DO IOP=1,6
         RR(IOP)=RR(IOP)*TWO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (RR(IOP),IOP=1,6)
      ENDIF
C
C     ----- DIAMAGNETIC SUSCEPTIBILITY TENSOR -----
C
      B20(1)=RR(2)+RR(3)
      B20(2)=RR(1)+RR(3)
      B20(3)=RR(1)+RR(2)
      B20(4)=-TWO*RR(4)
      B20(5)=-TWO*RR(5)
      B20(6)=-TWO*RR(6)
      IF(OUT) THEN
         WRITE(IW,9996)
         WRITE(IW,9997) (B20(IOP),IOP=1,6)
      ENDIF
C
      RETURN
 9998 FORMAT(' IN -GIAH20- , -R2- TENSOR = ')
 9997 FORMAT('    XX = ',F20.10,/,
     1       '    YY = ',F20.10,/,
     2       '    ZZ = ',F20.10,/,
     3       '    XY = ',F20.10,/,
     4       '    XZ = ',F20.10,/,
     5       '    YZ = ',F20.10,/)
 9996 FORMAT(' IN -GIAH20- , -B20- TENSOR = ')
      END
      SUBROUTINE HND_GIAE11(E01,E11,H01,D10,D00,NDIM,IFLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (TWO =2.0D+00)
      LOGICAL   OUT
      LOGICAL   DBUG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      DIMENSION E01(3,*),E11(3,3,*)            
      DIMENSION H01(NDIM,*)
      DIMENSION D10(NDIM,*)
      DIMENSION D00(*)
C
      OUT=.FALSE.
C
      CALL HND_REWFIL(IFT3)
      DO KAT=1,NAT
         DO KXYZ=1,3
            CALL HND_GIARD(IFT3,H01,NDIM*NDIM)
            DUM=ZERO
            DO J=1,NUM
               DO I=1,NUM
                  DUM=DUM+D10(I,J)*H01(I,J)
               ENDDO
            ENDDO
            DUM=DUM*TWO
            E11(KXYZ,IFLD,KAT)=DUM
         ENDDO
         IF(OUT) THEN
            WRITE(IW,9999) IFLD 
            CALL HND_PRSQ(E11(1,1,KAT),3,3,3)
         ENDIF
      ENDDO
C
      RETURN
 9999 FORMAT(' IN -GIAE11- -E11- FOR -IFLD- = ',I5)
      END
      SUBROUTINE HND_GIAH01(H01,C00,UUU,TTT,NBF,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXPRIM=2048)
      CHARACTER*8 ERRMSG
      LOGICAL NORM
      LOGICAL SOME
      LOGICAL OUT
      LOGICAL DBUG
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_RYS/YY,U5(5),W5(5),NROOTS
      COMMON/HND_HFK/XX,U9(9),W9(9),MROOTS
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ
     1                              ,CX,CY,CZ
      DIMENSION H01(NDIM,NDIM,*)
      DIMENSION C00(*)
      DIMENSION TTT(*)
      DIMENSION UUU(*)
      DIMENSION CIJ(225)
      DIMENSION H01X(225),H01Y(225),H01Z(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION   XV(5,6,6),  YV(5,6,6),  ZV(5,6,6)
      DIMENSION  DXV(5,6,6), DYV(5,6,6), DZV(5,6,6)
      DIMENSION  XVX(5,6,6), YVY(5,6,6), ZVZ(5,6,6)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-GIAH01-'/
      DATA MAXRYS /6/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      NPWR=1+1
C
C     ----- KAT -----
C
      CALL HND_REWFIL(IFT3)
      DO KAT=1,NAT
C
C     ----- ISHELL -----
C
      DO II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO JJ=1,NSHELL
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      LJTMOD=LJT+1
C
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NPWR-2)/2+1
      MROOTS=NROOTS
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9999) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            IJ=IJ+1
            H01X(IJ)=ZERO
            H01Y(IJ)=ZERO
            H01Z(IJ)=ZERO
         ENDDO
      ENDDO
C
C     ----- I PRIMITIVE -----
C
      DO IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.LE.TOL) THEN       
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ 
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
      CIJ(IJ)=DUM2
  360 CONTINUE
C
C     ----- -H01- TERM -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
C
         ZNUC=ONE
         CX=C(1,KAT)
         CY=C(2,KAT)
         CZ=C(3,KAT)
         XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
         YY=XX
         IF(NROOTS.LE.3) CALL HND_RT123
         IF(NROOTS.EQ.4) CALL HND_ROOT4
         IF(NROOTS.EQ.5) CALL HND_ROOT5
         IF(NROOTS.EQ.6) THEN
            CALL HND_DROOT
         ELSE
            DO IROOT=1,NROOTS
               U9(IROOT)=U5(IROOT)
               W9(IROOT)=W5(IROOT)
            ENDDO
         ENDIF
         DO IROOT=1,NROOTS
            UU=U9(IROOT)*AA
            WW=W9(IROOT)*ZNUC
            WW=WW*(UU+UU)
            TT=ONE/(AA+UU)
            T = SQRT(TT)
            X0=(AAX+UU*CX)*TT
            Y0=(AAY+UU*CY)*TT
            Z0=(AAZ+UU*CZ)*TT
            DO J=1,LJTMOD
               NJ=J
               DO I=1,LIT
                  NI=I
                  CALL HND_DSXYZ
                   XV(I,J,IROOT)=XINT
                   YV(I,J,IROOT)=YINT
                   ZV(I,J,IROOT)=ZINT*WW
                  CALL HND_DVXYZ
                  DXV(I,J,IROOT)=XINT
                  DYV(I,J,IROOT)=YINT
                  DZV(I,J,IROOT)=ZINT*WW
               ENDDO
            ENDDO
         ENDDO
C
         DO IROOT=1,NROOTS
            CALL HND_H01XYZ(
     1            XVX(1,1,IROOT),YVY(1,1,IROOT),ZVZ(1,1,IROOT),
     2             XV(1,1,IROOT), YV(1,1,IROOT), ZV(1,1,IROOT),
     3           LIT,LJT,AJ)
         ENDDO
C
         IJ=0
         DO I=MINI,MAXI
            IX=IJX(I)
            IY=IJY(I)
            IZ=IJZ(I)
            DO J=MINJ,MAXJ
               JX=IJX(J)
               JY=IJY(J)
               JZ=IJZ(J)
               DUMX=ZERO
               DUMY=ZERO
               DUMZ=ZERO
               DO IROOT=1,NROOTS
                  DUMX=DUMX+
     1   XV(IX,JX,IROOT)*DYV(IY,JY,IROOT)*ZVZ(IZ,JZ,IROOT)
     2 - XV(IX,JX,IROOT)*YVY(IY,JY,IROOT)*DZV(IZ,JZ,IROOT)
                  DUMY=DUMY+
     1  XVX(IX,JX,IROOT)* YV(IY,JY,IROOT)*DZV(IZ,JZ,IROOT)
     2 -DXV(IX,JX,IROOT)* YV(IY,JY,IROOT)*ZVZ(IZ,JZ,IROOT)
                  DUMZ=DUMZ+
     1  DXV(IX,JX,IROOT)*YVY(IY,JY,IROOT)* ZV(IZ,JZ,IROOT)
     2 -XVX(IX,JX,IROOT)*DYV(IY,JY,IROOT)* ZV(IZ,JZ,IROOT)
               ENDDO     
               IJ=IJ+1
               DUM=PI212*AA1*CIJ(IJ)
               IF(DBUG) THEN
                  WRITE(IW,9998) II,JJ,I,J,DUMX,DUMY,DUMZ,DUM
               ENDIF
               DUMX =DUMX*DUM
               DUMY =DUMY*DUM
               DUMZ =DUMZ*DUM
               H01X(IJ)=H01X(IJ)+DUMX
               H01Y(IJ)=H01Y(IJ)+DUMY
               H01Z(IJ)=H01Z(IJ)+DUMZ
            ENDDO     
         ENDDO     
C
            ENDIF
         ENDDO      
      ENDDO      
C
      IJ=0
      DO I=MINI,MAXI
         DO J=MINJ,MAXJ
            ILOC=LOCI+I
            JLOC=LOCJ+J
            IJ=IJ+1
            H01(ILOC,JLOC,1)= H01X(IJ)
            H01(ILOC,JLOC,2)= H01Y(IJ)
            H01(ILOC,JLOC,3)= H01Z(IJ)
         ENDDO
      ENDDO
      IF(DBUG) THEN
         IJ=0
         DO I=MINI,MAXI
            DO J=MINJ,MAXJ
               IJ=IJ+1
               WRITE(IW,9998) II,JJ,I,J,H01X(IJ),H01Y(IJ),
     1                                           H01Z(IJ)
            ENDDO
         ENDDO
      ENDIF
C
         ENDDO      
      ENDDO      
C
C     ----- PRINT -----
C
      IF(OUT) THEN
         WRITE(IW,9997) KAT
         CALL HND_PRSQ(H01(1,1,1),NUM,NUM,NDIM)
         WRITE(IW,9996) KAT
         CALL HND_PRSQ(H01(1,1,2),NUM,NUM,NDIM)
         WRITE(IW,9995) KAT
         CALL HND_PRSQ(H01(1,1,3),NUM,NUM,NDIM)
      ENDIF
C
C     ----- WRITE TO -IFT3- -----
C
      CALL HND_GIAWT(IFT3,H01(1,1,1),NDIM*NDIM)
      CALL HND_GIAWT(IFT3,H01(1,1,2),NDIM*NDIM)
      CALL HND_GIAWT(IFT3,H01(1,1,3),NDIM*NDIM)
C
      IF(OUT) THEN
         CALL HND_TFSQ(UUU,H01(1,1,1),C00,TTT,NORB,NUM,NDIM)
         WRITE(IW,9994) KAT
         CALL HND_PRSQ(UUU,NORB,NORB,NDIM)
         CALL HND_TFSQ(UUU,H01(1,1,2),C00,TTT,NORB,NUM,NDIM)
         WRITE(IW,9993) KAT
         CALL HND_PRSQ(UUU,NORB,NORB,NDIM)
         CALL HND_TFSQ(UUU,H01(1,1,3),C00,TTT,NORB,NUM,NDIM)
         WRITE(IW,9992) KAT
         CALL HND_PRSQ(UUU,NORB,NORB,NDIM)
      ENDIF
C
      ENDDO
C
      CALL HND_REWFIL(IFT3)
C
      RETURN
 9999 FORMAT(' IN -GIAH01- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,' LIT,LJT,NROOTS= ',3I3)
 9998 FORMAT(' II,JJ,I,J,H01X,Y,Z ... = ',4I4,4F15.10)
 9997 FORMAT(' -H01 X- IN -AO- FOR ATOM = ',I5)
 9996 FORMAT(' -H01 Y- IN _AO- FOR ATOM = ',I5)
 9995 FORMAT(' -H01 Z- IN -AO- FOR ATOM = ',I5)
 9994 FORMAT(' -H01 X- IN -MO- FOR ATOM = ',I5)
 9993 FORMAT(' -H01 Y- IN _MO- FOR ATOM = ',I5)
 9992 FORMAT(' -H01 Z- IN -MO- FOR ATOM = ',I5)
      END
      SUBROUTINE HND_H01XYZ(DX,DY,DZ,X,Y,Z,NI,NJ,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  X(5,*), Y(5,*), Z(5,*)
      DIMENSION DX(5,*),DY(5,*),DZ(5,*)
C
      DO I=1,NI
         DX(I,1)= (-(AJ+AJ)*X(I,2))
         DY(I,1)= (-(AJ+AJ)*Y(I,2))
         DZ(I,1)= (-(AJ+AJ)*Z(I,2))
      ENDDO
C
      IF(NJ.EQ.1) RETURN
C
      DO J=2,NJ
         DO I=1,NI
            DX(I,J)= (DBLE(J-1)*X(I,J-1)-(AJ+AJ)*X(I,J+1))
            DY(I,J)= (DBLE(J-1)*Y(I,J-1)-(AJ+AJ)*Y(I,J+1))
            DZ(I,J)= (DBLE(J-1)*Z(I,J-1)-(AJ+AJ)*Z(I,J+1))
         ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_GIAH11(E11D,E11P,D0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXPRIM=2048)
      CHARACTER*8 ERRMSG
      LOGICAL NORM
      LOGICAL SOME
      LOGICAL OUT
      LOGICAL MORE
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_NSHEL/
     1EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     1KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),KNG(MXSHEL),
     2KLOC(MXSHEL),KMIN(MXSHEL),KMAX(MXSHEL),NSHELL
      COMMON/HND_BASPAR/NORMF,NORMP,ITOL
      COMMON/HND_RYS/YY,U5(5),W5(5),NROOTS
      COMMON/HND_HFK/XX,U9(9),W9(9),MROOTS
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ
     1                              ,CX,CY,CZ
      DIMENSION E11D(3,3,*),E11P(3,3,*)            
      DIMENSION D0(*)
      DIMENSION CIJ(225)
      DIMENSION DIJ(225)
      DIMENSION V00(4,225)
      DIMENSION H01(3,225)
      DIMENSION H11D(9,225)
      DIMENSION H11P(9,225)
      DIMENSION H11T(9,225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION   VX0 (6,6,6),  VY0 (6,6,6),  VZ0 (6,6,6)
      DIMENSION  DVX0 (6,6,6), DVY0 (6,6,6), DVZ0 (6,6,6)
      DIMENSION   VX0D(6,6,6),  VY0D(6,6,6),  VZ0D(6,6,6)
      DIMENSION   VX1 (6,6,6),  VY1 (6,6,6),  VZ1 (6,6,6)
      DIMENSION  DVX1 (6,6,6), DVY1 (6,6,6), DVZ1 (6,6,6)
      DIMENSION   VX1X(6,6,6),  VY1Y(6,6,6),  VZ1Z(6,6,6)
      DIMENSION  DVX1X(6,6,6), DVY1Y(6,6,6), DVZ1Z(6,6,6)
      DIMENSION   VX1D(6,6,6),  VY1D(6,6,6),  VZ1D(6,6,6)
      DIMENSION  XVX1 (6,6,6), YVY1 (6,6,6), ZVZ1 (6,6,6)
      DIMENSION XDVX1 (6,6,6),YDVY1 (6,6,6),ZDVZ1 (6,6,6)
      DIMENSION  XVX1D(6,6,6), YVY1D(6,6,6), ZVZ1D(6,6,6)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-GIAH11-'/
      DATA MAXRYS /6/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
      IF(OUT) THEN
         WRITE(IW,*) 'IN -GIAH11-'
         CALL HND_PRTR(D0,NUM)
      ENDIF
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      DO KAT=1,NAT
         DO J=1,3
            DO I=1,3
               E11D(I,J,KAT)=ZERO
               E11P(I,J,KAT)=ZERO
            ENDDO
         ENDDO
      ENDDO
C
C     ----- ISHELL -----
C
      DO II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
      LITMOD=LIT+1
C
C     ----- JSHELL -----
C
      DO JJ=1,NSHELL
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      LJTMOD=LJT+1
C
      QIJX=YI*ZJ-ZI*YJ
      QIJY=ZI*XJ-XI*ZJ
      QIJZ=XI*YJ-YI*XJ
      TIJX=XI   -   XJ
      TIJY=YI   -   YJ
      TIJZ=ZI   -   ZJ
C
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
      NPWR=1+1+1
      NROOTS=(LIT+LJT+NPWR-2)/2+1
      MROOTS=NROOTS
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9999) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(OUT) THEN
         IJ=0
         DO I=MINI,MAXI
            IX=IJX(I)
            IY=IJY(I)
            IZ=IJZ(I)
            DO J=MINJ,MAXJ
               JX=IJX(J)
               JY=IJY(J)
               JZ=IJZ(J)
               IJ=IJ+1
               DO IX =1,4
                  V00(IX ,IJ)=ZERO
               ENDDO
               DO IX =1,3
                  H01(IX ,IJ)=ZERO
               ENDDO
               DO IXY=1,9
                  H11D(IXY,IJ)=ZERO
                  H11P(IXY,IJ)=ZERO
                  H11T(IXY,IJ)=ZERO
               ENDDO
            ENDDO
         ENDDO
         WRITE(IW,*) 'IN -GIAH11- II,JJ,XI,XJ,QIJX,TIJX ...'
         WRITE(IW,*)  II , JJ
         WRITE(IW,*) XI,YI,ZI
         WRITE(IW,*) XJ,YJ,ZJ
         WRITE(IW,*) QIJX,QIJY,QIJZ
         WRITE(IW,*) TIJX,TIJY,TIJZ
      ENDIF
C
C     ----- I PRIMITIVE -----
C
      DO IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      DO JG=J1,J2
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.LE.TOL) THEN       
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      DO 360 J=MINJ,MAXJ 
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      GO TO 350
  260 DUM2=DUM1*CDJ
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      ILOC=MAX(LOCI+I,LOCJ+J)
      JLOC=MIN(LOCI+I,LOCJ+J)
      NN=(ILOC*(ILOC-1))/2+JLOC
      DEN0=D0(NN)*TWO
      IJ=IJ+1
      CIJ(IJ)=DUM2
      DIJ(IJ)=DEN0
  360 CONTINUE          
C
C     ----- -H11- DIA + PARAMAGNETIC TERMS -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO KAT=1,NAT
C
         MORE=OUT.AND.KAT.EQ.1
C
         ZNUC=ONE
         CX=C(1,KAT)
         CY=C(2,KAT)
         CZ=C(3,KAT)
         XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
         YY=XX
         IF(NROOTS.LE.3) CALL HND_RT123
         IF(NROOTS.EQ.4) CALL HND_ROOT4
         IF(NROOTS.EQ.5) CALL HND_ROOT5
         IF(NROOTS.EQ.6) THEN
            CALL HND_DROOT
         ELSE
            DO IROOT=1,NROOTS
               U9(IROOT)=U5(IROOT)
               W9(IROOT)=W5(IROOT)
            ENDDO
         ENDIF
         DO IROOT=1,NROOTS
            UU=U9(IROOT)*AA
            WW=W9(IROOT)*ZNUC
            VV=WW
            WW=WW*(UU+UU)
            TT=ONE/(AA+UU)
            T = SQRT(TT)
            X0=(AAX+UU*CX)*TT
            Y0=(AAY+UU*CY)*TT
            Z0=(AAZ+UU*CZ)*TT
            DO J=1,LJTMOD
               NJ=J
               DO I=1,LITMOD
                  NI=I
                  CALL HND_DSXYZ
C
C     ----- FOR 1/R -----
C
               IF(MORE) THEN
                   VX0(I,J,IROOT)=XINT
                   VY0(I,J,IROOT)=YINT
                   VZ0(I,J,IROOT)=ZINT*VV
               ENDIF
C
C     ----- FOR X/R**3 -----
C
                   VX1(I,J,IROOT)=XINT
                   VY1(I,J,IROOT)=YINT
                   VZ1(I,J,IROOT)=ZINT*WW
                  CALL HND_DVXYZ
                  DVX1(I,J,IROOT)=XINT
                  DVY1(I,J,IROOT)=YINT
                  DVZ1(I,J,IROOT)=ZINT*WW
               ENDDO
            ENDDO
         ENDDO
C
         DO IROOT=1,NROOTS
            DO J=1,LJT
               DO I=1,LIT
                    VX1X(I,J,IROOT)= VX1(I,J+1,IROOT)
                    VY1Y(I,J,IROOT)= VY1(I,J+1,IROOT)
                    VZ1Z(I,J,IROOT)= VZ1(I,J+1,IROOT)
                   DVX1X(I,J,IROOT)=DVX1(I,J+1,IROOT)
                   DVY1Y(I,J,IROOT)=DVY1(I,J+1,IROOT)
                   DVZ1Z(I,J,IROOT)=DVZ1(I,J+1,IROOT)
               ENDDO
            ENDDO
            DO J=1,LJTMOD
               DO I=1,LIT
                   XVX1 (I,J,IROOT)= VX1(I+1,J,IROOT)
                   YVY1 (I,J,IROOT)= VY1(I+1,J,IROOT)
                   ZVZ1 (I,J,IROOT)= VZ1(I+1,J,IROOT)
                  XDVX1 (I,J,IROOT)=DVX1(I+1,J,IROOT)
                  YDVY1 (I,J,IROOT)=DVY1(I+1,J,IROOT)
                  ZDVZ1 (I,J,IROOT)=DVZ1(I+1,J,IROOT)
               ENDDO
            ENDDO
         ENDDO
         DO IROOT=1,NROOTS
C
C     ----- FOR 1/R DERIVATIVES WRT. CENTERS -I- AND -J- -----
C
         IF(MORE) THEN
            CALL HND_V00XYZ(
     1           DVX0(1,1,IROOT),DVY0(1,1,IROOT),DVZ0(1,1,IROOT),
     2           VX0D(1,1,IROOT),VY0D(1,1,IROOT),VZ0D(1,1,IROOT),
     3            VX0(1,1,IROOT), VY0(1,1,IROOT), VZ0(1,1,IROOT),
     4           LIT,LJT,AI,AJ)
         ENDIF
C
C     ----- D/DX ... OPERATORS -----
C
            CALL HND_H11XYZ(
     1            VX1D(1,1,IROOT), VY1D(1,1,IROOT), VZ1D(1,1,IROOT),
     2            VX1 (1,1,IROOT), VY1 (1,1,IROOT), VZ1 (1,1,IROOT),
     3           LIT,LJT,AJ)
            CALL HND_H11XYZ(
     1           XVX1D(1,1,IROOT),YVY1D(1,1,IROOT),ZVZ1D(1,1,IROOT),
     2           XVX1 (1,1,IROOT),YVY1 (1,1,IROOT),ZVZ1 (1,1,IROOT),
     3           LIT,LJT,AJ)
         ENDDO
C
         IJ=0
         DO I=MINI,MAXI
            IX=IJX(I)
            IY=IJY(I)
            IZ=IJZ(I)
            DO J=MINJ,MAXJ
               JX=IJX(J)
               JY=IJY(J)
               JZ=IJZ(J)
               DUMV  =ZERO
               TRANSX=ZERO
               TRANSY=ZERO
               TRANSZ=ZERO
               DUMVX =ZERO
               DUMVY =ZERO
               DUMVZ =ZERO
               DUMXYZ=ZERO
               DUMXX =ZERO
               DUMXY =ZERO
               DUMXZ =ZERO
               DUMYX =ZERO
               DUMYY =ZERO
               DUMYZ =ZERO
               DUMZX =ZERO
               DUMZY =ZERO
               DUMZZ =ZERO
               DUM1  =ZERO
               DUM2  =ZERO
               DUM3  =ZERO
               DUMX1 =ZERO
               DUMX2 =ZERO
               DUMX3 =ZERO
               DUMY1 =ZERO
               DUMY2 =ZERO
               DUMY3 =ZERO
               DUMZ1 =ZERO
               DUMZ2 =ZERO
               DUMZ3 =ZERO
               DO IROOT=1,NROOTS
C
      IF(MORE) THEN
C
C     ----- 1/R -----
C
      DUMV  =DUMV  +
     1    VX0(IX,JX,IROOT)* VY0(IY,JY,IROOT)* VZ0(IZ,JZ,IROOT)
C
C     ----- TRANSLATION INVARIANCE FOR NUCLEAR DERIVATIVES OF 1/R -----
C
      TRANSX=TRANSX+
     1   DVX0 (IX,JX,IROOT)* VY0 (IY,JY,IROOT)* VZ0 (IZ,JZ,IROOT)
     2 + DVX1 (IX,JX,IROOT)* VY1 (IY,JY,IROOT)* VZ1 (IZ,JZ,IROOT)
     3 +  VX0D(IX,JX,IROOT)* VY0 (IY,JY,IROOT)* VZ0 (IZ,JZ,IROOT)
      TRANSY=TRANSY+
     1    VX0 (IX,JX,IROOT)*DVY0 (IY,JY,IROOT)* VZ0 (IZ,JZ,IROOT)
     2 +  VX1 (IX,JX,IROOT)*DVY1 (IY,JY,IROOT)* VZ1 (IZ,JZ,IROOT)
     3 +  VX0 (IX,JX,IROOT)* VY0D(IY,JY,IROOT)* VZ0 (IZ,JZ,IROOT)
      TRANSZ=TRANSZ+
     1    VX0 (IX,JX,IROOT)* VY0 (IY,JY,IROOT)*DVZ0 (IZ,JZ,IROOT)
     2 +  VX1 (IX,JX,IROOT)* VY1 (IY,JY,IROOT)*DVZ1 (IZ,JZ,IROOT)
     3 +  VX0 (IX,JX,IROOT)* VY0 (IY,JY,IROOT)* VZ0D(IZ,JZ,IROOT)
C
C     ----- X/R**3 -----
C
      DUMVX1=DUMVX1+
     1   DVX1(IX,JX,IROOT)* VY1(IY,JY,IROOT)* VZ1(IZ,JZ,IROOT)
      DUMVY1=DUMVY1+
     1    VX1(IX,JX,IROOT)*DVY1(IY,JY,IROOT)* VZ1(IZ,JZ,IROOT)
      DUMVZ1=DUMVZ1+
     1    VX1(IX,JX,IROOT)* VY1(IY,JY,IROOT)*DVZ1(IZ,JZ,IROOT)
C
      ENDIF
C
C     ----- FOR H(1,1) DIAMAGNETIC -----
C
      DUMXX =DUMXX +
     1  DVX1X(IX,JX,IROOT)* VY1 (IY,JY,IROOT)* VZ1 (IZ,JZ,IROOT)
      DUMXY =DUMXY +
     1   VX1X(IX,JX,IROOT)*DVY1 (IY,JY,IROOT)* VZ1 (IZ,JZ,IROOT)
      DUMXZ =DUMXZ +
     1   VX1X(IX,JX,IROOT)* VY1 (IY,JY,IROOT)*DVZ1 (IZ,JZ,IROOT)
      DUMYX =DUMYX +
     1  DVX1 (IX,JX,IROOT)* VY1Y(IY,JY,IROOT)* VZ1 (IZ,JZ,IROOT)
      DUMYY =DUMYY +
     1   VX1 (IX,JX,IROOT)*DVY1Y(IY,JY,IROOT)* VZ1 (IZ,JZ,IROOT)
      DUMYZ =DUMYZ +
     1   VX1 (IX,JX,IROOT)* VY1Y(IY,JY,IROOT)*DVZ1 (IZ,JZ,IROOT)
      DUMZX =DUMZX +
     1  DVX1 (IX,JX,IROOT)* VY1 (IY,JY,IROOT)* VZ1Z(IZ,JZ,IROOT)
      DUMZY =DUMZY +
     1   VX1 (IX,JX,IROOT)*DVY1 (IY,JY,IROOT)* VZ1Z(IZ,JZ,IROOT)
      DUMZZ =DUMZZ +
     1   VX1 (IX,JX,IROOT)* VY1 (IY,JY,IROOT)*DVZ1Z(IZ,JZ,IROOT)
C
C     ----- FOR H(1,1) PARAMAGNETIC -----
C
      DUM1=DUM1+
     1    VX1 (IX,JX,IROOT)* DVY1 (IY,JY,IROOT)*  VZ1D(IZ,JZ,IROOT)
     2 -  VX1 (IX,JX,IROOT)*  VY1D(IY,JY,IROOT)* DVZ1 (IZ,JZ,IROOT)
      DUM2=DUM2+
     1    VX1D(IX,JX,IROOT)*  VY1 (IY,JY,IROOT)* DVZ1 (IZ,JZ,IROOT)
     2 - DVX1 (IX,JX,IROOT)*  VY1 (IY,JY,IROOT)*  VZ1D(IZ,JZ,IROOT)
      DUM3=DUM3+
     1   DVX1 (IX,JX,IROOT)*  VY1D(IY,JY,IROOT)*  VZ1 (IZ,JZ,IROOT)
     2 -  VX1D(IX,JX,IROOT)* DVY1 (IY,JY,IROOT)*  VZ1 (IZ,JZ,IROOT)
C
      DUMX1=DUMX1+
     1   XVX1 (IX,JX,IROOT)* DVY1 (IY,JY,IROOT)*  VZ1D(IZ,JZ,IROOT)
     2 - XVX1 (IX,JX,IROOT)*  VY1D(IY,JY,IROOT)* DVZ1 (IZ,JZ,IROOT)
      DUMX2=DUMX2+
     1   XVX1D(IX,JX,IROOT)*  VY1 (IY,JY,IROOT)* DVZ1 (IZ,JZ,IROOT)
     2 -XDVX1 (IX,JX,IROOT)*  VY1 (IY,JY,IROOT)*  VZ1D(IZ,JZ,IROOT)
      DUMX3=DUMX3+
     1  XDVX1 (IX,JX,IROOT)*  VY1D(IY,JY,IROOT)*  VZ1 (IZ,JZ,IROOT)
     2 - XVX1D(IX,JX,IROOT)* DVY1 (IY,JY,IROOT)*  VZ1 (IZ,JZ,IROOT)
      DUMY1=DUMY1+
     1    VX1 (IX,JX,IROOT)*YDVY1 (IY,JY,IROOT)*  VZ1D(IZ,JZ,IROOT)
     2 -  VX1 (IX,JX,IROOT)* YVY1D(IY,JY,IROOT)* DVZ1 (IZ,JZ,IROOT)
      DUMY2=DUMY2+
     1    VX1D(IX,JX,IROOT)* YVY1 (IY,JY,IROOT)* DVZ1 (IZ,JZ,IROOT)
     2 - DVX1 (IX,JX,IROOT)* YVY1 (IY,JY,IROOT)*  VZ1D(IZ,JZ,IROOT)
      DUMY3=DUMY3+
     1   DVX1 (IX,JX,IROOT)* YVY1D(IY,JY,IROOT)*  VZ1 (IZ,JZ,IROOT)
     2 -  VX1D(IX,JX,IROOT)*YDVY1 (IY,JY,IROOT)*  VZ1 (IZ,JZ,IROOT)
      DUMZ1=DUMZ1+
     1    VX1 (IX,JX,IROOT)* DVY1 (IY,JY,IROOT)* ZVZ1D(IZ,JZ,IROOT)
     2 -  VX1 (IX,JX,IROOT)*  VY1D(IY,JY,IROOT)*ZDVZ1 (IZ,JZ,IROOT)
      DUMZ2=DUMZ2+
     1    VX1D(IX,JX,IROOT)*  VY1 (IY,JY,IROOT)*ZDVZ1 (IZ,JZ,IROOT)
     2 - DVX1 (IX,JX,IROOT)*  VY1 (IY,JY,IROOT)* ZVZ1D(IZ,JZ,IROOT)
      DUMZ3=DUMZ3+
     1   DVX1 (IX,JX,IROOT)*  VY1D(IY,JY,IROOT)* ZVZ1 (IZ,JZ,IROOT)
     2 -  VX1D(IX,JX,IROOT)* DVY1 (IY,JY,IROOT)* ZVZ1 (IZ,JZ,IROOT)
C
               ENDDO     
               IJ=IJ+1
               DUM=PI212*AA1*CIJ(IJ)
               DUMV  =DUMV  *DUM
               TRANSX=TRANSX*DUM
               TRANSY=TRANSY*DUM
               TRANSZ=TRANSZ*DUM
               DUMVX =DUMVX *DUM
               DUMVY =DUMVY *DUM
               DUMVZ =DUMVZ *DUM
               DUMXX =DUMXX *DUM
               DUMXY =DUMXY *DUM
               DUMXZ =DUMXZ *DUM
               DUMYX =DUMYX *DUM
               DUMYY =DUMYY *DUM
               DUMYZ =DUMYZ *DUM
               DUMZX =DUMZX *DUM
               DUMZY =DUMZY *DUM
               DUMZZ =DUMZZ *DUM
               DUM1  =DUM1  *DUM
               DUM2  =DUM2  *DUM
               DUM3  =DUM3  *DUM
               DUMX1 =DUMX1 *DUM
               DUMX2 =DUMX2 *DUM
               DUMX3 =DUMX3 *DUM
               DUMY1 =DUMY1 *DUM
               DUMY2 =DUMY2 *DUM
               DUMY3 =DUMY3 *DUM
               DUMZ1 =DUMZ1 *DUM
               DUMZ2 =DUMZ2 *DUM
               DUMZ3 =DUMZ3 *DUM
C
            IF(MORE) THEN
C
C     ----- 1/R AND TRANSLATION INVARIANCE -----
C
               V00(1,IJ)=V00(1,IJ)+DUMV
               V00(2,IJ)=V00(2,IJ)+DUMVX
               V00(3,IJ)=V00(3,IJ)+DUMVY
               V00(4,IJ)=V00(4,IJ)+DUMVZ
C
            ENDIF
C
C     ----- H(1,1) DIAMAGNETIC -----
C
               DUMXYZ= DUMXX +DUMYY +DUMZZ
               TMPXXD=(DUMXYZ-DUMXX)
               TMPXYD=(      -DUMXY)
               TMPXZD=(      -DUMXZ)
               TMPYXD=(      -DUMYX)
               TMPYYD=(DUMXYZ-DUMYY)
               TMPYZD=(      -DUMYZ)
               TMPZXD=(      -DUMZX)
               TMPZYD=(      -DUMZY)
               TMPZZD=(DUMXYZ-DUMZZ)
            IF(MORE) THEN
               H11D(1,IJ)=H11D(1,IJ)+TMPXXD
               H11D(2,IJ)=H11D(2,IJ)+TMPXYD
               H11D(3,IJ)=H11D(3,IJ)+TMPXZD
               H11D(4,IJ)=H11D(4,IJ)+TMPYXD
               H11D(5,IJ)=H11D(5,IJ)+TMPYYD
               H11D(6,IJ)=H11D(6,IJ)+TMPYZD
               H11D(7,IJ)=H11D(7,IJ)+TMPZXD
               H11D(8,IJ)=H11D(8,IJ)+TMPZYD
               H11D(9,IJ)=H11D(9,IJ)+TMPZZD
            ENDIF
               E11D(1,1,KAT)=E11D(1,1,KAT)+TMPXXD*DIJ(IJ)
               E11D(1,2,KAT)=E11D(1,2,KAT)+TMPXYD*DIJ(IJ)
               E11D(1,3,KAT)=E11D(1,3,KAT)+TMPXZD*DIJ(IJ)
               E11D(2,1,KAT)=E11D(2,1,KAT)+TMPYXD*DIJ(IJ)
               E11D(2,2,KAT)=E11D(2,2,KAT)+TMPYYD*DIJ(IJ)
               E11D(2,3,KAT)=E11D(2,3,KAT)+TMPYZD*DIJ(IJ)
               E11D(3,1,KAT)=E11D(3,1,KAT)+TMPZXD*DIJ(IJ)
               E11D(3,2,KAT)=E11D(3,2,KAT)+TMPZYD*DIJ(IJ)
               E11D(3,3,KAT)=E11D(3,3,KAT)+TMPZZD*DIJ(IJ)
C
C     ----- H(1,1) PARAMAGNETIC -----
C
            IF(MORE) THEN
               H01(1,IJ)=H01(1,IJ)+DUM1
               H01(2,IJ)=H01(2,IJ)+DUM2
               H01(3,IJ)=H01(3,IJ)+DUM3
            ENDIF
C  
C     ----- A BIT UNCLEAR ... YX ... INSTEAD OF ... XY ... -----
C
               TMPXXP=(QIJX*DUM1+TIJY*DUMZ1-TIJZ*DUMY1)
               TMPYXP=(QIJX*DUM2+TIJY*DUMZ2-TIJZ*DUMY2)
               TMPZXP=(QIJX*DUM3+TIJY*DUMZ3-TIJZ*DUMY3)
               TMPXYP=(QIJY*DUM1+TIJZ*DUMX1-TIJX*DUMZ1)
               TMPYYP=(QIJY*DUM2+TIJZ*DUMX2-TIJX*DUMZ2)
               TMPZYP=(QIJY*DUM3+TIJZ*DUMX3-TIJX*DUMZ3)
               TMPXZP=(QIJZ*DUM1+TIJX*DUMY1-TIJY*DUMX1)
               TMPYZP=(QIJZ*DUM2+TIJX*DUMY2-TIJY*DUMX2)
               TMPZZP=(QIJZ*DUM3+TIJX*DUMY3-TIJY*DUMX3)
C   
            IF(MORE) THEN
               H11P(1,IJ)=H11P(1,IJ)+TMPXXP
               H11P(2,IJ)=H11P(2,IJ)+TMPXYP
               H11P(3,IJ)=H11P(3,IJ)+TMPXZP
               H11P(4,IJ)=H11P(4,IJ)+TMPYXP
               H11P(5,IJ)=H11P(5,IJ)+TMPYYP
               H11P(6,IJ)=H11P(6,IJ)+TMPYZP
               H11P(7,IJ)=H11P(7,IJ)+TMPZXP
               H11P(8,IJ)=H11P(8,IJ)+TMPZYP
               H11P(9,IJ)=H11P(9,IJ)+TMPZZP
            ENDIF
               E11P(1,1,KAT)=E11P(1,1,KAT)+TMPXXP*DIJ(IJ)
               E11P(1,2,KAT)=E11P(1,2,KAT)+TMPXYP*DIJ(IJ)
               E11P(1,3,KAT)=E11P(1,3,KAT)+TMPXZP*DIJ(IJ)
               E11P(2,1,KAT)=E11P(2,1,KAT)+TMPYXP*DIJ(IJ)
               E11P(2,2,KAT)=E11P(2,2,KAT)+TMPYYP*DIJ(IJ)
               E11P(2,3,KAT)=E11P(2,3,KAT)+TMPYZP*DIJ(IJ)
               E11P(3,1,KAT)=E11P(3,1,KAT)+TMPZXP*DIJ(IJ)
               E11P(3,2,KAT)=E11P(3,2,KAT)+TMPZYP*DIJ(IJ)
               E11P(3,3,KAT)=E11P(3,3,KAT)+TMPZZP*DIJ(IJ)
C
C     ----- H(1,1) DIA + PARA MAGNETIC -----
C
               TMPXX =TMPXXD+TMPXXP
               TMPXY =TMPXYD+TMPXYP
               TMPXZ =TMPXZD+TMPXZP
               TMPYX =TMPYXD+TMPYXP
               TMPYY =TMPYYD+TMPYYP
               TMPYZ =TMPYZD+TMPYZP
               TMPZX =TMPZXD+TMPZXP
               TMPZY =TMPZYD+TMPZYP
               TMPZZ =TMPZZD+TMPZZP
            IF(MORE) THEN
               H11T(1,IJ)=H11T(1,IJ)+TMPXX 
               H11T(2,IJ)=H11T(2,IJ)+TMPXY 
               H11T(3,IJ)=H11T(3,IJ)+TMPXZ 
               H11T(4,IJ)=H11T(4,IJ)+TMPYX 
               H11T(5,IJ)=H11T(5,IJ)+TMPYY 
               H11T(6,IJ)=H11T(6,IJ)+TMPYZ 
               H11T(7,IJ)=H11T(7,IJ)+TMPZX 
               H11T(8,IJ)=H11T(8,IJ)+TMPZY 
               H11T(9,IJ)=H11T(9,IJ)+TMPZZ 
C
               IF(ABS(TRANSX).GE.1.0D-08.OR.ABS(TRANSY).GE.1.0D-08
     1                                  .OR.ABS(TRANSZ).GE.1.0D-08) THEN
                  WRITE(IW,9993)
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ENDIF
C
            ENDDO     
         ENDDO     
C
      ENDDO      
C
            ENDIF
         ENDDO      
      ENDDO      
C
      IF(OUT) THEN
         WRITE(IW,*) ' IN -GIAH11- II,JJ  = ',
     1                             II,JJ
         DO IX=1,4
            IJ=0
            DO I=MINI,MAXI
               DO J=MINJ,MAXJ
                  IJ=IJ+1
                  WRITE(IW,9997) IX,I,J,V00(IX,IJ),CIJ(IJ) 
               ENDDO
            ENDDO
         ENDDO
         DO IX=1,3
            IJ=0
            DO I=MINI,MAXI
               DO J=MINJ,MAXJ
                  IJ=IJ+1
                  WRITE(IW,9992) IX,I,J,H01(IX,IJ)
               ENDDO
            ENDDO
         ENDDO
         WRITE(IW,9996)
         DO IXY=1,9
            IJ=0
            DO I=MINI,MAXI
               DO J=MINJ,MAXJ
                  IJ=IJ+1
                  WRITE(IW,9995) IXY,I,J,H11D(IXY,IJ),H11P(IXY,IJ),
     1                                   H11T(IXY,IJ),DIJ(IJ)
               ENDDO
            ENDDO
         ENDDO
C 
         IF(OUT) THEN
            DO KAT=1,NAT
               WRITE(IW,*) ' IN -GIAH11- II,JJ,KAT = ',II,JJ,KAT
               CALL HND_PRSQ(E11D(1,1,KAT),3,3,3)
               CALL HND_PRSQ(E11P(1,1,KAT),3,3,3)
            ENDDO
         ENDIF
      ENDIF
C
         ENDDO      
      ENDDO      
C
C     ----- PRINT -----
C
      IF(SOME) THEN
         WRITE(IW,9998) 
         DO KAT=1,NAT
            CALL HND_PRSQ(E11D(1,1,KAT),3,3,3)
         ENDDO
         WRITE(IW,9994) 
         DO KAT=1,NAT
            CALL HND_PRSQ(E11P(1,1,KAT),3,3,3)
         ENDDO
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -GIH11P- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,' LIT,LJT,NROOTS= ',3I3)
 9998 FORMAT(/,' E11D = P00 * H11D = ',
     1       /,' ------------------- ')
 9997 FORMAT(' -V00(IX , )- = ',3I5,2F15.10)
 9996 FORMAT(/)
 9995 FORMAT(' -H11(IXY, )- = ',3I5,4F15.10)
 9994 FORMAT(/,' E11P = P00 * H11P = ',
     1       /,' ------------------- ')
 9993 FORMAT(' SOMETHING WRONG WITH TRANSLATIONAL',
     1       ' INVARIANCE IN -GIAH11- . STOP. ')
 9992 FORMAT(' -H01(IX , )- = ',3I5,F15.10)
      END
      SUBROUTINE HND_H11XYZ(DX,DY,DZ,X,Y,Z,NI,NJ,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  X(6,*), Y(6,*), Z(6,*)
      DIMENSION DX(6,*),DY(6,*),DZ(6,*)
C
      DO I=1,NI
         DX(I,1)= (-(AJ+AJ)*X(I,2))
         DY(I,1)= (-(AJ+AJ)*Y(I,2))
         DZ(I,1)= (-(AJ+AJ)*Z(I,2))
      ENDDO
C
      IF(NJ.EQ.1) RETURN
C
      DO J=2,NJ
         DO I=1,NI
            DX(I,J)= (DBLE(J-1)*X(I,J-1)-(AJ+AJ)*X(I,J+1))
            DY(I,J)= (DBLE(J-1)*Y(I,J-1)-(AJ+AJ)*Y(I,J+1))
            DZ(I,J)= (DBLE(J-1)*Z(I,J-1)-(AJ+AJ)*Z(I,J+1))
         ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_V00XYZ(DX,DY,DZ,XD,YD,ZD,X,Y,Z,NI,NJ,AI,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  X(6,*), Y(6,*), Z(6,*)
      DIMENSION DX(6,*),DY(6,*),DZ(6,*)
      DIMENSION XD(6,*),YD(6,*),ZD(6,*)
C
C     ----- DERIVATIVES WITH RESPECT TO XJ ... -----
C
      DO I=1,NI
         XD(I,1)=-(-(AJ+AJ)*X(I,2))
         YD(I,1)=-(-(AJ+AJ)*Y(I,2))
         ZD(I,1)=-(-(AJ+AJ)*Z(I,2))
      ENDDO
      IF(NJ.GT.1) THEN  
         DO J=2,NJ
            DO I=1,NI
               XD(I,J)=-(DBLE(J-1)*X(I,J-1)-(AJ+AJ)*X(I,J+1))
               YD(I,J)=-(DBLE(J-1)*Y(I,J-1)-(AJ+AJ)*Y(I,J+1))
               ZD(I,J)=-(DBLE(J-1)*Z(I,J-1)-(AJ+AJ)*Z(I,J+1))
            ENDDO
         ENDDO
      ENDIF
C
C     ----- DERIVATIVES WITH RESPECT TO XI ... -----
C
      DO J=1,NJ
         DX(1,J)=-(-(AI+AI)*X(2,J))
         DY(1,J)=-(-(AI+AI)*Y(2,J))
         DZ(1,J)=-(-(AI+AI)*Z(2,J))
      ENDDO
      IF(NI.GT.1) THEN  
         DO J=1,NJ
            DO I=2,NI
               DX(I,J)=-(DBLE(I-1)*X(I-1,J)-(AI+AI)*X(I+1,J))
               DY(I,J)=-(DBLE(I-1)*Y(I-1,J)-(AI+AI)*Y(I+1,J))
               DZ(I,J)=-(DBLE(I-1)*Z(I-1,J)-(AI+AI)*Z(I+1,J))
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_GIAKD(KFLG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE HND_GIAKX(D,F,XX,IX,NINTMX,IA,NOPK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE HND_GIAK(D,F,XX,IX,NINTMX,ID,NOPK)
c     SUBROUTINE HND_GIAK(D,F,XX,IX,NINTMX,IA,NOPK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- -GIAK- FORMS A SKELETON MATRIX -----
C                   F=( H* + H )/2
C
C              F(I,J)=(H**(I,J) + H**(J,I))/2
C
C     INDICES IN LABELS ARE IN STANDARD ORDER_
C      I.GE.J , K.GE.L , (IJ).GE.(KL)
C
C     ALL CONTRIBUTIONS ARE MADE INTO LOWER HALF OF
C     SKELETON MATRIX.
C     ONLY OFF-DIAGONAL ELEMENTS NEED BE DIVIDED BY TWO,
C     TO OBTAIN THE CORRECT F MATRIX.
C
      PARAMETER (MXATOM=500)
      CHARACTER*8 ERRMSG
      INTEGER SHIFTR
      LOGICAL INDEX2
      LOGICAL DBLLAB
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_SQFILE/IJK,IPK
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_HNDOPT/HNDTYP,DSKTYP,FILTYP,PCKTYP
      COMMON/HND_PCKLAB/LABSIZ
      DIMENSION D(1),F(1)
c     DIMENSION IA(1)
      DIMENSION XX(NINTMX),IX(NINTMX)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','- GIAK -'/
      DATA ZERO,PT5  /0.0D+00,0.5D+00/
*not used      DATA PCK2,PCK4 /'-2-LABEL','-4-LABEL'/
*not used      DATA DSK4,DSK6 /'DSK-P -O','DSK-PK-O'/
C
      SHIFTR(IARG,IBIT)= ISHFT(IARG,-IBIT)
      LAND(IARG1,IARG2)= IAND(IARG1,IARG2)
      MASK1(IARG)=2**IARG-1
      MASK2(IARG)=2**IARG-1
      IWORD1(IARG)=LAND(SHIFTR(IARG,24),MASK1(8))
      IWORD2(IARG)=LAND(SHIFTR(IARG,16),MASK1(8))
      IWORD3(IARG)=LAND(SHIFTR(IARG, 8),MASK1(8))
      IWORD4(IARG)=LAND(SHIFTR(IARG, 0),MASK1(8))
      JWORD1(IARG)=LAND(SHIFTR(IARG,16),MASK2(16))
      JWORD2(IARG)=LAND(SHIFTR(IARG, 0),MASK2(16))
C
      IA(I)=(I*(I-1))/2
C
      DBLLAB=LABSIZ.EQ.2
      NUM2=(NUM*(NUM+1))/2
C
      DO 50 M=1,NUM2
   50    F(M)=ZERO
C
      IF(NOPK.NE.1) GO TO 1000
C
C     ----- INTEGRALS ARE NOT IN SUPERMATRIX FORM (NOPK=.TRUE.) -----
C
  100 CALL HND_PREAD(IJK,XX,IX,NXX,NINTMX)
      IF(NXX.EQ.0) GO TO 350
      NINT=IABS(NXX)
      IF(NINT.GT.NINTMX) CALL HND_HNDERR(3,ERRMSG)
      IF(DBLLAB) THEN
         DO 200 M=1,NINT
            LABEL1=IX((M-1)*LABSIZ+1)
            LABEL2=IX((M-1)*LABSIZ+2)
            I=JWORD1(LABEL1)
            J=JWORD2(LABEL1)
            K=JWORD1(LABEL2)
            L=JWORD2(LABEL2)
            NIJ=IA(I)+J
            NKL=IA(K)+L
            NIK=IA(I)+K
            NIL=IA(I)+L
            NJK=IA(MAX(J,K))+MIN(J,K)
            NJL=IA(MAX(J,L))+MIN(J,L)
            VAL=XX(M)
            VAL4=(VAL+VAL)+(VAL+VAL)
            F(NIK)=F(NIK)+VAL*D(NJL)
            F(NIL)=F(NIL)+VAL*D(NJK)
            F(NJK)=F(NJK)+VAL*D(NIL)
            F(NJL)=F(NJL)+VAL*D(NIK)
  200       CONTINUE
      ELSE
         DO 300 M=1,NINT
            LABEL=IX(M)
            I=IWORD1(LABEL)
            J=IWORD2(LABEL)
            K=IWORD3(LABEL)
            L=IWORD4(LABEL)
            NIJ=IA(I)+J
            NKL=IA(K)+L
            NIK=IA(I)+K
            NIL=IA(I)+L
            NJK=IA(MAX(J,K))+MIN(J,K)
            NJL=IA(MAX(J,L))+MIN(J,L)
            VAL=XX(M)
            VAL4=(VAL+VAL)+(VAL+VAL)
            F(NIK)=F(NIK)+VAL*D(NJL)
            F(NIL)=F(NIL)+VAL*D(NJK)
            F(NJK)=F(NJK)+VAL*D(NIL)
            F(NJL)=F(NJL)+VAL*D(NIK)
  300       CONTINUE
      ENDIF
      IF(NXX.GT.0) GO TO 100
  350 CONTINUE
      IF(NUM.EQ.1) GO TO 500
      DO 400 M=2,NUM
         DO 400 N=1,M-1
            NIJ=IA(M)+N
  400       F(NIJ)=F(NIJ)*PT5
  500 CONTINUE
      CALL HND_REWFIL(IJK)
      RETURN
C
C     ----- INTEGRALS ARE IN SUPERMATRIX FORM (NOPK=.FALSE.) -----
C
 1000 CONTINUE
      WRITE(IW,9999)
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(' INTEGRALS MAY NOT BE IN SUPERMATRIX FORM',
     1       ' FOR THE -GIAO- -K+0- OPERATOR. STOP')
      END
      SUBROUTINE HND_GIAS(GAUGE,R00,S,T,NUM,NOC,NORB,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FILGIA/IFT1,IFT2,IFT3,IFT4
      DIMENSION R00(NDIM,NDIM,*)
      DIMENSION   S(NDIM,NDIM)
      DIMENSION   T(*)
      DIMENSION GAUGE(3,NOC)
      DATA ZERO    /0.0D+00/
      DATA ONE     /1.0D+00/
      DATA TWO     /2.0D+00/
      DATA CVEL    /137.0359895D+00/
C
      IA(I)=(I*(I-1))/2
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
C     ----- ORTHONORMALIZATION DERIVATIVE TERM -----
C
      IF(OUT) THEN
         WRITE(IW,*) 'ORTHONORMALIZATION DERIVATIVE TERM'      
      ENDIF
      CALL HND_REWFIL(IFT2)
      DO IOP=1,3
         IF(OUT) THEN
            WRITE(IW,*) 'READ     RX, RY, RZ ... = ',IOP   
         ENDIF
         NORB2=(NORB*(NORB+1))/2
         CALL HND_GIARD(IFT2,T,NORB2)
         DO IORB=1,NORB
            DO JORB=1,IORB
               IJ=IA(IORB)+JORB
               R00(IORB,JORB,IOP)=T(IJ)
               R00(JORB,IORB,IOP)=T(IJ)
            ENDDO
         ENDDO
         IF(DBUG) THEN
            CALL HND_PRSQ(R00(1,1,IOP),NORB,NORB,NDIM)
         ENDIF
      ENDDO
C
C     ----- WE ARE READY ... -----
C
      CALL HND_REWFIL(IFT2)
      DO IFLD=1,3
         IF(IFLD.EQ.1) THEN
            IOP1=2
            IOP2=3
            IF(OUT) THEN
               WRITE(IW,*) 'B(X) PERTURBATION'
            ENDIF
         ELSEIF(IFLD.EQ.2) THEN
            IOP1=3
            IOP2=1
            IF(OUT) THEN
               WRITE(IW,*) 'B(Y) PERTURBATION'
            ENDIF
         ELSEIF(IFLD.EQ.3) THEN
            IOP1=1
            IOP2=3
            IF(OUT) THEN
               WRITE(IW,*) 'B(Z) PERTURBATION'
            ENDIF
         ENDIF 
         DO KORB=1,NORB
            DO JORB=1,NORB
               S(JORB,KORB)=ZERO
            ENDDO
         ENDDO
         VELFAC=ONE/(TWO*CVEL)
         DO KOC=1,NORB
            IF(KOC.LE.NOC) THEN
               GXK=GAUGE(1,KOC)
               GYK=GAUGE(2,KOC)
               GZK=GAUGE(3,KOC)
               IF(DBUG) THEN
                  WRITE(IW,9999) KOC,GXK,GYK,GZK
               ENDIF
            ELSE
               GXK=ZERO
               GYK=ZERO
               GZK=ZERO
            ENDIF
            DO JOC=1,NORB
               IF(JOC.LE.NOC) THEN
                  GXJ=GAUGE(1,JOC)
                  GYJ=GAUGE(2,JOC)
                  GZJ=GAUGE(3,JOC)
                  IF(DBUG) THEN
                     WRITE(IW,9999) JOC,GXJ,GYJ,GZJ
                  ENDIF
               ELSE
                  GXK=ZERO
                  GYK=ZERO
                  GZK=ZERO
               ENDIF
               GX=GXJ-GXK
               GY=GYJ-GYK
               GZ=GZJ-GZK
               IF(IFLD.EQ.1) THEN
                  C1=GZ
                  C2=GY
               ELSEIF(IFLD.EQ.2) THEN
                  C1=GX
                  C2=GZ
               ELSEIF(IFLD.EQ.3) THEN
                  C1=GY
                  C2=GX
               ENDIF
               S(JOC,KOC)=( C1*R00(JOC,KOC,IOP1)
     1                     -C2*R00(JOC,KOC,IOP2) )*(VELFAC/TWO)
            ENDDO
         ENDDO
         CALL HND_GIAWT(IFT2,S,NDIM*NDIM)
         IF(DBUG) THEN
            CALL HND_PRSQ(S,NORB,NORB,NDIM)
         ENDIF
      ENDDO
      CALL HND_REWFIL(IFT2)
C
      RETURN
 9999 FORMAT(' IN -GIAS- , GAUGE FOR -OCC- = ',I5,
     1       ' GX,GY,GZ = ',3F15.9)
      END
      SUBROUTINE HND_GIACVG(D,H,U,U0,NUM,NOC,NORB,NDIM,
     1                  IEND,ITER,SOME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "util.fh"
C
C     ----- THIS ROUTINE CHECKS FOR CONVERGENCE OF THE CPHF -----
C           CALCULATION.  CONVERGENCE IS DETERMINED BY TWO
C           CRITERIA: THE CONVERGENCE OF THE TENSOR ELEMENT
C           FOR THE XX, YY, OR ZZ DIRECTION AND THE
C           CONVERGENCE OF THE INDIVIDUAL ELEMENTS OF THE U
C           MATRIX.
C
      LOGICAL  STATUS
      LOGICAL  PROP_GIAO_PUT_U1
      EXTERNAL PROP_GIAO_PUT_U1
      LOGICAL  PROP_GIAO_GET_U1
      EXTERNAL PROP_GIAO_GET_U1
      LOGICAL DBUG
      LOGICAL  OUT
      LOGICAL SOME
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (MXIODA=255)   
      COMMON/HND_CVGGIA/CNVA,CNVB,dummy_damp
      COMMON/HND_PARGIA/VNEW,VOLD,idummy_maxit_or_itmax
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      DIMENSION D(NDIM,*),H(NDIM,*),U(NDIM,*),U0(NDIM,*)
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      IF(DBUG) THEN
         WRITE(IW,*) '-D(1) IN -GIACVG-'
         CALL HND_PRSQ(D,NUM,NUM,NDIM)
      ENDIF
      IF(OUT) THEN
         WRITE(IW,*) '-U(1) IN -GIACVG-'
         CALL HND_PRSQ(U,NORB,NORB,NDIM)
      ENDIF
C
      IF(ITER.GT.0) THEN    
C
C     ----- ITER > 1 -----
C
         VOLD  = VNEW
         VNEW  = VOLD           
         DIFF  = VOLD - VNEW
         UDELM = ZERO
         UVAL  = ZERO
         IVAL  = 0
         JVAL  = 0
C
         NUM3=NUM*NUM
      if(ga_nodeid().eq.0) then
         CALL HND_DAREAD(IDAF,IODA,U0,NUM3,32)
         STATUS=PROP_GIAO_GET_U1(U0,NUM,NUM3)
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIACVG PROP_GIAO_GET_U1 FAILED ',911)
         ENDIF
      endif
      call ga_brdcst(789+1+msgdbl,u0,mdtob(num3),0)
      call ga_sync()
         IF(DBUG) THEN
            WRITE(IW,*) '-PREVIOUS U(1) IN -GIACVG-'
            CALL HND_PRSQ(U0,NORB,NORB,NDIM)
         ENDIF
C
         ICON=0
         DO J=1,NOC
            DO I=NOC+1,NORB
               UDEL = U(I,J) - U0(I,J)
               IF(ABS(UDEL).GE.CNVB) ICON = ICON + 1
               IF(ABS(UDEL).GT.ABS(UDELM)) THEN
                  UDELM = UDEL
                  UVAL = U(I,J)
                  IVAL = I
                  JVAL = J
               ENDIF
            ENDDO
         ENDDO
C
         NUM3=NUM*NUM
      if(ga_nodeid().eq.0) then
         CALL HND_DAWRIT(IDAF,IODA,U,NUM3,32,NAV)
         STATUS=PROP_GIAO_PUT_U1(U,NUM,NUM3)
      endif
      call ga_sync()
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIACVG PROP_GIAO_GET_U1 FAILED ',911)
         ENDIF
C
         IF(ABS(DIFF).LT.CNVA.AND.ICON.EQ.0) THEN 
            IEND = 1
         ELSE
            IEND = 0
         ENDIF
      if(ga_nodeid().eq.0) then
         IF(SOME) THEN
            WRITE(IW,9999) ITER,VNEW,DIFF,UDELM,UVAL,IVAL,JVAL,
     1                     ICON,IEND
         ENDIF
      endif
         RETURN
C
      ELSE
C
C     ----- ITER = 1 -----
C
         VOLD  = ZERO
         VNEW  = VOLD            
         DIFF  = VOLD - VNEW
         UDELM = ZERO
         UVAL  = ZERO
         IVAL  = 0
         JVAL  = 0
C
         NUM3=NUM*NUM
      if(ga_nodeid().eq.0) then
         CALL HND_DAWRIT(IDAF,IODA,U,NUM3,32,NAV)
         STATUS=PROP_GIAO_PUT_U1(U,NUM,NUM3)
      endif
      call ga_brdcst(789+2+msgdbl,u0,mdtob(num3),0)
      call ga_sync()
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('HND_GIACVG PROP_GIAO_GET_U1 FAILED ',911)
         ENDIF
C
         ICON = NOC * (NORB - NOC)
         IEND=0
      if(ga_nodeid().eq.0) then
         IF(SOME) THEN
            WRITE(IW,9999) ITER,VNEW,DIFF,UDELM,UVAL,IVAL,JVAL,
     1                     ICON,IEND
         ENDIF
      endif
      ENDIF
      RETURN
 9999 FORMAT(1X,I3,F20.9,3F16.9,2I4,I10,I3)
      END
      SUBROUTINE HND_GIAK10(D,F,XX,IX,NINTMX,NDIM,NOPK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- THIS ROUTINE FORMS (NUM*NUM) FOCK MATRIX. IT IS
C           ASSUMED THAT BOTH  D AND F MATRICES ARE UNSYMMETRIC.----
C
C     ----- AT THE PRESENT TIME  THE NOPK=.FALSE. OPTION IS NOT
C           INCLUDED.                                           -----
C     ----- FOR THE 2E INTEGRALS INDICES IN LABELS ARE IN
C           STANDARD ORDER , THAT IS,
C               I.GE.J , K.GE.L , (IJ).GE.(KL)
C
      PARAMETER   (MXATOM=500)
      PARAMETER   (ZERO=0.0D+00)
      CHARACTER*8 ERRMSG
      CHARACTER*8 HNDTYP
      CHARACTER*8 DSKTYP
      CHARACTER*8 FILTYP
      CHARACTER*8 PCKTYP
      CHARACTER*8 DISK
      LOGICAL     DBUG
      LOGICAL     OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_SQFILE/IJK,IPK
      COMMON/HND_HNDOPT/HNDTYP,DSKTYP,FILTYP,PCKTYP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION D(NDIM,*),F(NDIM,*),XX(*),IX(*)
      DIMENSION ERRMSG(3)
      DATA   DISK /'DISK IO '/
      DATA ERRMSG /'PROGRAM ','STOP IN ','-GIAK10-'/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      IF(OUT) THEN
         WRITE(IW,*) '-D(1)- IN -GIAK10-'
         CALL HND_PRSQ(D,NUM,NUM,NDIM)
      ENDIF
C
      DO M2=1,NUM
         DO M1=1,NUM
            F(M1,M2)=ZERO
         ENDDO
      ENDDO
C
      IF(NOPK.NE.1) THEN
C
C     ----- INTEGRALS ARE IN SUPERMATRIX FORM (NOPK=.FALSE.) -----
C     ----- NOT IMPLIMENTED IN FPL SYSTEMS AT THIS POINT, SO ----
C
         WRITE (IW,9999)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(DSKTYP.EQ.DISK) THEN
C
C     ----- INTEGRALS ARE NOT IN SUPERMATRIX FORM (NOPK=.TRUE.) -----
C
         CALL HND_REWFIL(IJK)
   10    CALL HND_PREAD(IJK,XX,IX,NXX,NINTMX)
         IF(NXX.NE.0) THEN     
            NINT=IABS(NXX)
            IF(NINT.GT.NINTMX) THEN 
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            CALL HND_K10GIA(F,D,NDIM,XX,IX,NINT)
            IF(NXX.GT.0) GO TO 10
         ENDIF
         CALL HND_REWFIL(IJK)
C
      ELSE
C
C     ----- IN-CORE METHOD -----
C
         CALL HND_K10GIX(F,D,NDIM)
      ENDIF
      RETURN
 9999 FORMAT(' FOR -GIA- CALCULATIONS, THE INTEGRALS MAY NOT BE',
     1       ' IN SUPERMATRIX FORMAT. STOP.',/,
     2       ' SET -NPKFIL- = 1 IN -$INTGRL- .')
 9998 FORMAT(' IN-CORE METHOD NOT YET IMPLEMENTED')
      END
      SUBROUTINE HND_K10GIA(F,D,NDIM,XX,IX,NINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PT5=0.5D+00)
      INTEGER SHIFTR
      COMMON/HND_PCKLAB/LABSIZ
      DIMENSION F(NDIM,*),D(NDIM,*),XX(*),IX(*)
C
      SHIFTR(IARG,IBIT)= ISHFT(IARG,-IBIT)
      LAND(IARG1,IARG2)= IAND(IARG1,IARG2)
      MASK1(IARG)=2**IARG-1
      MASK2(IARG)=2**IARG-1
      IWORD1(IARG)=LAND(SHIFTR(IARG,24),MASK1(8))
      IWORD2(IARG)=LAND(SHIFTR(IARG,16),MASK1(8))
      IWORD3(IARG)=LAND(SHIFTR(IARG, 8),MASK1(8))
      IWORD4(IARG)=LAND(SHIFTR(IARG, 0),MASK1(8))
      JWORD1(IARG)=LAND(SHIFTR(IARG,16),MASK1(16))
      JWORD2(IARG)=LAND(SHIFTR(IARG, 0),MASK1(16))
C
      IF(LABSIZ.EQ.2) THEN
         DO M=1,NINT
            LABEL1=IX((M-1)*LABSIZ+1)
            LABEL2=IX((M-1)*LABSIZ+2)
            I=JWORD1(LABEL1)
            J=JWORD2(LABEL1)
            K=JWORD1(LABEL2)
            L=JWORD2(LABEL2)
            VAL  = XX(M)
            VALH = VAL*PT5
            F(I,K) = F(I,K)-D(J,L)*VALH
            F(I,L) = F(I,L)-D(J,K)*VALH
            F(J,K) = F(J,K)-D(I,L)*VALH
            F(J,L) = F(J,L)-D(I,K)*VALH
            F(K,I) = F(K,I)-D(L,J)*VALH
            F(K,J) = F(K,J)-D(L,I)*VALH
            F(L,I) = F(L,I)-D(K,J)*VALH
            F(L,J) = F(L,J)-D(K,I)*VALH
         ENDDO       
      ELSE
         DO M=1,NINT
            LABEL=IX(M)
            I=IWORD1(LABEL)
            J=IWORD2(LABEL)
            K=IWORD3(LABEL)
            L=IWORD4(LABEL)
            VAL  = XX(M)
            VALH = VAL*PT5
            F(I,K) = F(I,K)-D(J,L)*VALH
            F(I,L) = F(I,L)-D(J,K)*VALH
            F(J,K) = F(J,K)-D(I,L)*VALH
            F(J,L) = F(J,L)-D(I,K)*VALH
            F(K,I) = F(K,I)-D(L,J)*VALH
            F(K,J) = F(K,J)-D(L,I)*VALH
            F(L,I) = F(L,I)-D(K,J)*VALH
            F(L,J) = F(L,J)-D(K,I)*VALH
         ENDDO       
      ENDIF
      RETURN
      END
      SUBROUTINE HND_K10GIX(F,D,NUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PT5=0.5D+00)
      INTEGER SHIFTR
      COMMON/HND_LCM   /Y(1)
      DIMENSION F(NUM,*),D(NUM,*)
      DIMENSION IY(1)
      EQUIVALENCE (IY(1),Y(1))
C
      SHIFTR(IARG,IBIT)= ISHFT(IARG,-IBIT)
      LAND(IARG1,IARG2)= IAND(IARG1,IARG2)
      MASK1(IARG)=2**IARG-1
      JWORD1(IARG)=LAND(SHIFTR(IARG,16),MASK1(16))
      JWORD2(IARG)=LAND(SHIFTR(IARG, 0),MASK1(16))
C
      NUM2=(NUM *(NUM +1))/2
      NUM4=(NUM2*(NUM2+1))/2
      NY=0
      IJ=0
      DO I=1,NUM
      DO J=1,I
         IJ=IJ+1
         MKL=IY(IJ)
         IF(MKL.NE.0) THEN      
         DO MY=1,MKL
            VAL = Y(NY+MY+NUM2+NUM4)
            KL  =IY(NY+MY+NUM2     )
            K=JWORD1(KL)
            L=JWORD2(KL)
            VALH = VAL*PT5
            F(I,K) = F(I,K)-D(J,L)*VALH
            F(I,L) = F(I,L)-D(J,K)*VALH
            F(J,K) = F(J,K)-D(I,L)*VALH
            F(J,L) = F(J,L)-D(I,K)*VALH
            F(K,I) = F(K,I)-D(L,J)*VALH
            F(K,J) = F(K,J)-D(L,I)*VALH
            F(L,I) = F(L,I)-D(K,J)*VALH
            F(L,J) = F(L,J)-D(K,I)*VALH
         ENDDO
         ENDIF
         NY = NY+IJ
      ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE HND_GIAG10(H10,F00,D00,NDIM,NEXT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- -GIAO- 2-E INTEGRALS -----
C
      PARAMETER (MXIODA=255)
      LOGICAL   DBUG
      LOGICAL   OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      DIMENSION H10(NDIM*NDIM,*)
      DIMENSION F00(NDIM*NDIM)
      DIMENSION D00(*)
      DATA ZERO /0.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(OUT) THEN
         WRITE(IW,*) 'IN -GIAG10-'
         CALL HND_FIFLSH(IW)
      ENDIF
C
      DO IOP=1,3
         DO IJ=1,NDIM*NDIM
            H10(IJ,IOP)=ZERO
         ENDDO
      ENDDO
      DO IJ=1,NDIM*NDIM
         F00(IJ)=ZERO
      ENDDO
C
C     ----- DENSITY MATRIX -----
C
      NDIM2=(NDIM*(NDIM+1))/2
      CALL HND_DAREAD(IDAF,IODA,D00,NDIM2,16)
C
C     ----- -GIAO- INTEGRALS -----
C
      CALL HND_GIA2EI(NEXT)
C
      IF(OUT) THEN
         WRITE(IW,*) 'D00'
         CALL HND_PRTR(D00,NDIM)
         WRITE(IW,*) 'G00'
         CALL HND_PRSQ(F00,NDIM,NDIM,NDIM)
         WRITE(IW,*) 'G10-X'
         CALL HND_PRSQ(H10(1,1),NDIM,NDIM,NDIM)
         WRITE(IW,*) 'G10-Y'
         CALL HND_PRSQ(H10(1,2),NDIM,NDIM,NDIM)
         WRITE(IW,*) 'G10-Z'
         CALL HND_PRSQ(H10(1,3),NDIM,NDIM,NDIM)
      ENDIF
C
      RETURN
      END
      SUBROUTINE HND_TFSQ(H,F,Q,T,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- H(M,M) = Q(DAGGER)(N,M) * F(N,N) * Q(N,M) -----
C
      DIMENSION H(NDIM,*),F(NDIM,*),Q(NDIM,*),T(*)
      DATA SMALL /1.0D-11/
      DATA ZERO  /0.0D+00/
C
      DO I=1,N
         DO J=1,M
            DUM=ZERO
            DO K=1,N
               DUM=DUM+F(I,K)*Q(K,J)
            ENDDO
            T(J)=DUM
         ENDDO
         DO J=1,M
            F(I,J)=T(J)
         ENDDO
      ENDDO
C
      DO J=1,M
         DO I=1,M
            DUM=ZERO
            DO K=1,N
               DUM=DUM+Q(K,I)*F(K,J)
            ENDDO
            T(I)=DUM
         ENDDO
         DO I=1,M
            H(I,J)=T(I)
         ENDDO
      ENDDO
C
      RETURN
      END
      subroutine hnd_stubs
      entry hnd_lcpsnc 
      entry hnd_hdir 
      entry hnd_hstar 
      entry hnd_xstar 
      entry hnd_lcpadd 
      write(6,*) 'pass in faked -lcpadd-'
      stop
      return
      entry hnd_fiflsh 
      entry hnd_texit 
      return
      entry hnd_pkread 
      write(6,*) 'pass in faked -pkread-'
      return
      entry hnd_pkwrit 
      write(6,*) 'pass in faked -pkwrit-'
      return
      entry hnd_pwrit 
      write(6,*) 'pass in faked -pwrit -'
      return
      entry hnd_advfil 
      write(6,*) 'pass in faked -advfil-'
      return
      entry hnd_initpn 
      write(6,*) 'pass in faked -initpn-'
      return
      entry hnd_initpk 
      write(6,*) 'pass in faked -initpk-'
      return
      entry hnd_exitpn 
      write(6,*) 'pass in faked -exitpn-'
      return
      entry hnd_exitpk 
      write(6,*) 'pass in faked -exitpk-'
      return
      entry hnd_ppack 
      write(6,*) 'pass in faked -ppack -'
      return
      entry hnd_pkpack 
      write(6,*) 'pass in faked -pkpack-'
      return
      entry hnd_jkpack 
      write(6,*) 'pass in faked -jkpack-'
      return
      entry hnd_dipkwr 
      write(6,*) 'pass in faked -dipkwr-'
      return
      entry hnd_dipwrt 
      write(6,*) 'pass in faked -dipwrt-'
      return
      end
