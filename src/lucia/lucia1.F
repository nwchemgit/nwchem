      SUBROUTINE ICOPVE2(IIN,IOFF,NDIM,IOUT)
*
* IOUT(I) = IIN(IOFF-1+I),I = 1, NDIM
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input 
      DIMENSION IIN(*)
*. Output
      DIMENSION IOUT(*)
*
      DO I = 1, NDIM
        IOUT(I) = IIN(IOFF-1+I)
      END DO
*
      RETURN
      END
      
      FUNCTION ISTRN(STRING,IGROUP)
*
* A string belonging to group IGROUP is given.
* find actual number
*
* Jeppe Olsen, September 1993
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. Specific input
      INTEGER STRING(*)
*
      INCLUDE 'strbas.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*./STRINP/
C     COMMON/STRINP/NSTTYP,MNRS1(MXPSTT),MXRS1(MXPSTT),
C    &              MNRS3(MXPSTT),MXRS3(MXPSTT),NELEC(MXPSTT),
C    &              IZORR(MXPSTT),IAZTP,IBZTP,IARTP(3,10),IBRTP(3,10),
C    &              NZSTTP,NRSTTP,ISTTP(MXPSTT)
C     COMMON/STRBAS/KSTINF,KOCSTR(MXPSTT),KNSTSO(MXPSTT),KISTSO(MXPSTT),
C    &              KSTSTM(MXPSTT,2),KZ(MXPSTT),
C    &              KSTREO(MXPSTT),KSTSM(MXPSTT),KSTCL(MXPSTT),
C    &              KEL1(MXPSTT),KEL3(MXPSTT),KACTP(MXPSTT),
C    &              KCOBSM,KNIFSJ,KIFSJ,KIFSJO,KSTSTX
C    &             ,KNDMAP(MXPSTT),KNUMAP(MXPSTT)
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPOBS,MXPR4T),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),
C    &              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),
C    &              ITOOBS(MXPOBS)
      
      NEL = NELEC(IGROUP) 
CJO-DEC  
C     ISTRN = ISTRNM(STRING,NACOB,NELEC,WORK(KZ(IGROUP)),
C    &               WORK(KSTREO(IGROUP)),1)
      ISTRN = ISTRNM(STRING,NACOB,NEL,int_mb(KZ(IGROUP)),
     &               int_mb(KSTREO(IGROUP)),1)
CJO-DEC

*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' String '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ' Actual address', ISTRN
      END IF
*
      RETURN
      END
       SUBROUTINE UNPCPC(IRC,IR,NR,IC)
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* A compound index IRC is given for a given element in
* a column packed lower half of a matrix with NR rows.
*
* Find the corrsponding row and column numbers, IR, IC
* 
* The relation to be fulfilled is
* IRC = (IC-1)*NR+IR -IC*(IC-1)/2 
*
*. Start by solving equation assuming IR = IC
* the corresponding second equation  is
* -XC**2/2 + XC(NR+1.5) -(NR +IRC)
*
* Well according to Mister Bechmann 
      A = -0.5D0
      B = FLOAT(NR) + 1.5D0
      C = -FLOAT(NR + IRC)
      XC = (-B)/(2.0D0*A) 
     &   + 1.0D0/(2.0D0*A)* SQRT( B ** 2 - 4.0D0*A*C)
*. Round down to get column number
      IC = XC
*. Row number
      IR = IRC - (IC-1)*NR + IC*(IC-1)/2
*. Check
      IF(IC.GT.IR  .OR. IR.GT.NR .OR. IR .LE. 0 ) THEN
        WRITE(6,*) ' Dear Sir '
        WRITE(6,*) ' I am a subroutine called UNPCPC '
        WRITE(6,*) ' trying to do my best, but my programmer '
        WRITE(6,*) ' ( Jeppe, if you know him ) '
        WRITE(6,*) ' has left me with some problems '
        WRITE(6,*) ' So I quit in a microsecond or two '
        WRITE(6,*)
        WRITE(6,*) ' IRC IR IC NR ', IRC,IR,IC,NR
        STOP 'UNPCPC '
      END IF
*
      RETURN
      END 
      Subroutine ExpStr(NstrI,Iocstr,Nel,NorbI,NorbO,IreOrb,
     &                  IltoaI,IltoaO,IarcwO,Ninob,                             
     &                  IreStr,Isign,Iscr,IPRT)
*                                                                               
* Mapping between strings in different Ras Spaces                               
*                                                                               
* Jeppe Olsen , October 1991                                                    
*                                                                               
      Implicit Double Precision (A-H,O-Z)                                       
* =======                                                                       
* Input :                                                                       
* =======                                                                       
*. Occupation of Input strings                                                  
      Integer Iocstr(Nel,NstrI)                                                 
*. Lexical to output order for Input strings                                    
      Integer IltoaI(*)                                                         
*. Arcweights for output strings                                                
      Integer IarcwO(NorbO,Nel)                                                 
*. Lexical to output order for Output strings                                   
      Integer IltoaO(*)                                                         
*. Reordering of Orbitals                                                       
      Integer IreOrb(*)                                                         
* ========                                                                      
* Output :                                                                      
* ========                                                                      
*                                                                               
*. Actual order of input strings to actual order of output strings              
      Integer IreStr(*)                                                         
*. Sign shift ( is maybe neccessary sometime, somewhere, somehow)               
      Integer Isign(*)                                                          
* =======                                                                       
* Scratch                                                                       
* =======                                                                       
      Integer Iscr(*)                                                           
*. Iscr should atleast be of length 2 * Nel                                     
*                                                                               
* Loop over lexical ordering of input string,                                   
* get corresponding actual string,                                              
* change  orbitals to output form and obtain lexcical output number             
* obtain corresponding actual number                                            
*                                                                               
      Ntest = 0
      NTEST = MAX(NTEST,IPRT)
      IF( NTEST .GE. 1 ) THEN
        Write(6,*) ' ================================'
        Write(6,*) ' Hi man , Expstr at your service '
        Write(6,*) ' ================================'
      END IF
*
      Kstr1 = 1                                                                 
      Kstr2 = Kstr1 + Nel                                                       
      Do 100 IstrLi = 1, NstrI                                                  
        IstrAi = IltoaI(IstrLi)                                                 
        Do 10 Iel = 1, Nel                                                      
          Iscr(Kstr1-1+Iel) = IreOrb(Iocstr(Iel,IstrAi)+NinOB)-Ninob
   10   Continue                                                                
        If( Ntest.Ge. 500 ) Then
          write(6,*) ' Input string,Lexical and actual number '
          Write(6,*)   IstrLi,IstrAi
          Write(6,*) ' Corresponding string,original orbital numbers'
          Call Iwrtma(Iocstr(1,IstrAi),1,Nel,1,Nel)
          Write(6,*) '                     ,modified orbital numbers'
          Call Iwrtma(Iscr(Kstr1),1,Nel,1,Nel)
        End if
*. Order to increasing ordering                                                 
C         ORDSTR(IINST,IOUTST,NELMNT,ISIGN )                                    
        If(Nel.Ne.0)
     &  Call Ordstr(Iscr(Kstr1),Iscr(Kstr2),Nel,Isgn,NTEST)
        If( Ntest .Ge. 500 ) Then
          Write(6,*) ' Ordered string ' 
          Call Iwrtma(Iscr(Kstr2),1,Nel,1,Nel)
        End if
*. Lexical number of output string                                              
C                 ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
        IstrAO =  ISTRNM(Iscr(Kstr2),NORBO,Nel,iArcwO,IltoaO,1)          
        Irestr(IstrAI) = IstrAO                                                 
        If( Ntest. Ge. 500 ) then
          write(6,*) ' Reordered string,actual number '               
          Write(6,*) IstrAO                                          
        End if
        Isign(IstrAi) = isgn                                                    
  100 Continue                                                                  
*                                                                               
      If(Ntest .ge. 5 ) then                                                    
        Write(6,*) ' ----------------------'                                    
        Write(6,*) ' Expstr at your service' 
        Write(6,*) ' ----------------------'                                    
        Write(6,*) ' Input actual order to Output actual order '                
        If(Ntest.Ge. 1000) Then
          Nprt = NstrI
        Else
          Maxprt = 50  
          Nprt = Min(NstrI,Maxprt)
        End if
*
        If(Nprt.ne.Nstri) Write(6,*) ' Only a part printed out '                
        Call Iwrtma(IreStr,1,Nprt,1,Nprt)                                       
        Write(6,*) ' Sign shift ( given in actual input order )'                
        Call Iwrtma(Isign,1,Nprt,1,Nprt)                                        
      End if                                                                    
*                                                                               
      Return                                                                    
      End                                                                       
      SUBROUTINE COMCIM(ISBDET,NSBDET,NVAR,ICMP,HOLD)
*
* Construct complete CI matrix spanned by determinants defined by ISBDET in lower half form
*
* Compare with H if ICMP .ne. 0 
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (MXNSB = 100 )
      PARAMETER (MXNVAR = 10 000 )
      DIMENSION H(MXNSB*(MXNSB+1)/2)
      DIMENSION VEC1(MXNVAR),VEC2(MXNVAR)
      DIMENSION ISBDET(*)
      DIMENSION HOLD(*)
*
      DO 100 I = 1, NSBDET
        CALL SETVEC(VEC1,0.0D0,NVAR)
        VEC1(ISBDET(I)) = 1.0D0
        CALL MV7(VEC1,VEC2,0,0,0,0)
        DO 90 J = 1, I
         H(I*(I-1)/2+J) = VEC2(ISBDET(J))
*
         write(6,*) ' Matrix under construction '
        CALL PRSYM(H,NSBDET)
   90   CONTINUE
  100 CONTINUE
*
*. For test with alpha-alpha loop
      IF(ICMP.NE.0) THEN
        WRITE(6,*) ' Comparison of H and HOLD '
        CALL CMP2VC(H,HOLD,NSBDET*(NSBDET+1)/2,1.0D-11)
      END IF
*
      NTEST = 1
      IF(NTEST.NE.0) THEN
*
        WRITE(6,*) ' ========================================'
        WRITE(6,*) ' Complete CI matrix delivered as  COMCIM'
        WRITE(6,*) ' ========================================'
*
        CALL PRSYM(H,NSBDET)
      END IF
*
      RETURN
      END 
      SUBROUTINE DIHDJ(IASTR,IBSTR,NIDET,
     &                 JASTR,JBSTR,NJDET,
     &                 NAEL,NBEL,
     &IWORK,LWORK,NORB,ONEBOD,HAMIL,ISYM,NINOB,ECORE,ICOMBI,PSIGN,
     & IASTRM,IBSTRM,JASTRM,JBSTRM,
     & IGENSG,IASGN,IBSGN,JASGN,JBSGN,LIA,LIB,NDIF0,NDIF1,NDIF2,
     & IPRT)
*
* A set of left hand side determinants defined by string numbers 
* IASTR and IBSTR and a set of right hand side determinants 
* defined by JASTR and JBSTR are given.
*
* Obtain Hamiltonian matrix  < IA IB ! H ! JA JB >
*
* If Icombi .NE. 0 Spin combinations are assumed  for alpha and 
* beta strings with different orbital configurations
*   1/SQRT(2) * ( !I1A I2B! + PSIGN * !I2A I1B! )
*
* If ISYM .EQ. 0 FULL Hamiltonian is constructed
* If ISYM .NE. 0 LOWER half of hamiltonian is constructed
*
* JEPPE OLSEN JANUARY 1989
*
*. Modifed to work with string numbers instead of strings
*. March 93
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IASTR(*),IBSTR(*)
      DIMENSION JASTR(*),JBSTR(*)
      DIMENSION IASTRM(NAEL,*),IBSTRM(NBEL,*)
      DIMENSION JASTRM(NAEL,*),JBSTRM(NBEL,*)
      DIMENSION IASGN(*),IBSGN(*),JASGN(*),JBSGN(*)
*
      DIMENSION IWORK(*), HAMIL(*), ONEBOD(NORB,NORB)
      DIMENSION LIA(NAEL),LIB(NBEL)
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRT)
      IF( NTEST .GE. 5 ) THEN
         WRITE(6,*) ' =========================='
         WRITE(6,*) ' DIHDJ reporting to service'
         WRITE(6,*) ' =========================='
      END IF
*
*. To get rid of annoying compiler warnings
      JAEQJB = 0
      IEL1 = 0
      JEL1 = 0
      IPERM = 0
      JPERM = 0
      SIGN = 0.0D0

      IF(NTEST.GE. 200) THEN
        WRITE(6,*) ' Input alpha strings,IASTR,JASTR '
        CALL IWRTMA(IASTR,1,NIDET,1,NIDET)
        CALL IWRTMA(JASTR,1,NJDET,1,NJDET)
        WRITE(6,*) ' Input beta  strings,IBSTR,JBSTR '
        CALL IWRTMA(IBSTR,1,NIDET,1,NIDET)
        CALL IWRTMA(JBSTR,1,NJDET,1,NJDET)
      END IF
*
*. Scratch space : 4 vectors of length NORB
      KLFREE = 1
      KLIAE  = KLFREE
      KLFREE = KLIAE + NORB
      KLIBE  = KLFREE
      KLFREE = KLIBE + NORB
*
      KLJAE = KLFREE
      KLFREE = KLJAE + NORB
      KLJBE = KLFREE
      KLFREE = KLJBE + NORB
*
      IF( ISYM .EQ. 0 ) THEN
        LHAMIL = NIDET*NJDET
      ELSE
        LHAMIL = NIDET*(NIDET+1) / 2
      END IF
      CALL SETVEC(HAMIL,0.0D0,LHAMIL)
*
      NTERMS= 0
      NDIF0 = 0
      NDIF1 = 0
      NDIF2 = 0
*
*. Loop over J determinants
*
      DO 1000 JDET = 1,NJDET
        IF( NTEST .GE. 100 ) WRITE(6,*) '  ****** 1000 JDET ', JDET
* Expand JDET
        JASTAC =JASTR(JDET)
        JBSTAC =JBSTR(JDET)
*
        IF(IGENSG .GT. 0 ) THEN
         JXSGN = JASGN(JASTAC)*JBSGN(JBSTAC)
        ELSE
         JXSGN = 1
        END IF
*
        CALL ISETVC(IWORK(KLJAE),0,NORB)
        CALL ISETVC(IWORK(KLJBE),0,NORB)
        DO 40 IAEL = 1, NAEL
          IWORK(KLJAE-1+JASTRM(IAEL,JASTAC) ) = 1
   40   CONTINUE
*
        DO 50 IBEL = 1, NBEL
          IWORK(KLJBE-1+JBSTRM(IBEL,JBSTAC) ) = 1
   50   CONTINUE
*
        IF( ICOMBI .NE. 0 ) THEN
          IF(JASTAC .EQ. JBSTAC) THEN 
             JAEQJB = 1
          ELSE 
             JAEQJB = 0
          END IF
        END IF
*
*
        IF( NTEST .GE. 100 ) THEN
          WRITE(6,*) ' LOOP 1000 JDET =  ',JDET
          WRITE(6,*) ' JASTAC AND JBSTAC ', JASTAC,JBSTAC
          WRITE(6,*) ' Expanded ALPHA and BETA string '
          CALL IWRTMA(IWORK(KLJAE),1,NORB,1,NORB)
          CALL IWRTMA(IWORK(KLJBE),1,NORB,1,NORB)
        END IF
*
        IF( ISYM .EQ. 0 ) THEN
          MINI = 1
        ELSE
          MINI = JDET
        END IF
*
* Loop over I determinants 
*
        DO 900 IDET = MINI, NIDET
          IF(NTEST.GE. 100 ) 
     &    WRITE(6,*) '   LOOP 900 IDET .... ',IDET
          IASTAC = IASTR(IDET)
          IBSTAC = IBSTR(IDET)
*
          IF(IGENSG .GT. 0 ) THEN
           IXSGN = IASGN(IASTAC)*IBSGN(IBSTAC)
          ELSE
           IXSGN = 1
          END IF
*
          IF(IASTAC.EQ.IBSTAC) THEN 
            IAEQIB = 1
          ELSE 
            IAEQIB = 0
          END IF

*
          IF(ICOMBI.EQ.1 .AND. IAEQIB+JAEQJB.EQ.0 ) THEN 
              NLOOP = 2
          ELSE
              NLOOP = 1
          END IF
C
          DO 899 ILOOP = 1, NLOOP
           NTERMS = NTERMS + 1
* For second part of spin combinations strings should be swopped         
           IF(ILOOP.EQ.1) THEN
             CALL ICOPVE(IASTRM(1,IASTAC),LIA,NAEL)
             CALL ICOPVE(IBSTRM(1,IBSTAC),LIB,NBEL)
           ELSE IF (ILOOP.EQ.2) THEN
             CALL ICOPVE(IASTRM(1,IASTAC),LIB,NAEL)
             CALL ICOPVE(IBSTRM(1,IBSTAC),LIA,NBEL)
           END IF
*
* ==============================
*. Number of orbital differences 
* ==============================
*
           NACM = 0
           DO 61 IAEL = 1, NAEL
             NACM = NACM + IWORK(KLJAE-1+LIA(IAEL))
   61      CONTINUE
C
           NBCM = 0
           DO 62 IBEL = 1, NBEL
             NBCM = NBCM + IWORK(KLJBE-1+LIB(IBEL))
   62      CONTINUE
C
           NADIF = NAEL-NACM
           NBDIF = NBEL-NBCM
           IF( NTEST .GE. 100 ) THEN
             WRITE(6,*) '  LOOP 900 IDET ',IDET
             WRITE(6,*) ' Comparison : NADIF ,NBDIF ', NADIF,NBDIF
           END IF
*
           IF(NADIF+NBDIF .GT. 2 ) GOTO 898
*. Factor for combinations
           IF( ICOMBI .EQ. 0 ) THEN
             CONST = 1.0D0
           ELSE
             IF((JAEQJB +IAEQIB) .EQ.2 ) THEN
               CONST = 1.0D0
             ELSE IF( (JAEQJB+IAEQIB) .EQ. 1 ) THEN
               CONST = 1.0D0/SQRT(2.0D0)*(1.0D0+PSIGN)
              ELSE IF( (JAEQJB+IAEQIB) .EQ. 0 ) THEN
               IF( ILOOP .EQ. 1)  THEN
                 CONST = 1.0D0
               ELSE
                 CONST = PSIGN
               END IF
             END IF
           END IF
*. External sign factor
           IF(IXSGN*JXSGN .EQ. -1 ) CONST = - CONST
           IF(NTEST.GE.100) WRITE(6,*) ' CONST ', CONST
*
* ==================================================
*.. Find differing orbitals and sign for permutation
* ==================================================
*
* Expand idet
           CALL ISETVC(IWORK(KLIAE),0,NORB)
           CALL ISETVC(IWORK(KLIBE),0,NORB)
*
           DO 42 IAEL = 1, NAEL
             IWORK(KLIAE-1+LIA(IAEL)) = 1
   42      CONTINUE
*
           DO 52 IBEL = 1, NBEL
             IWORK(KLIBE-1+LIB(IBEL) ) = 1
   52      CONTINUE
*
*. One pair of differing alpha electrons
*
           IF(NADIF .EQ. 1 ) THEN
             DO 120 IAEL = 1,NAEL
               IF(IWORK(KLJAE-1+LIA(IAEL)).EQ.0) THEN
                 IA = LIA(IAEL)          
                 IEL1 = IAEL
                 GOTO 121
               END IF
  120        CONTINUE
  121        CONTINUE
*
             DO 130 JAEL = 1,NAEL
               IF(IWORK(KLIAE-1+JASTRM(JAEL,JASTAC)).EQ.0) THEN
                 JA = JASTRM(JAEL,JASTAC)
                 JEL1 = JAEL
                 GOTO 131
                END IF
  130        CONTINUE
  131        CONTINUE
             SIGNA = (-1)**(JEL1+IEL1)
             IF(NTEST.GE.100)
     &       WRITE(6,*) ' IA JA SIGNA... ',IA,JA,SIGNA
           END IF
*
*. One pair of differing beta electrons
*
           IF(NBDIF .EQ. 1 ) THEN
             DO 220 IBEL = 1,NBEL
               IF(IWORK(KLJBE-1+LIB(IBEL) ).EQ.0) THEN
                 IB = LIB(IBEL)
                 IEL1 = IBEL
                 GOTO 221
                END IF
  220        CONTINUE
  221        CONTINUE
C
             DO 230 JBEL = 1,NBEL
               IF(IWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)).EQ.0) THEN
                 JB = JBSTRM(JBEL,JBSTAC)
                 JEL1 = JBEL
                 GOTO 231
                END IF
  230        CONTINUE
  231        CONTINUE
             SIGNB = (-1)**(JEL1+IEL1)
             IF(NTEST.GE.100) 
     &       WRITE(6,*) ' IB JB SIGNB... ',IB,JB,SIGNB
           END IF
*
*. Two pairs of differing alpha electrons
*
           IF(NADIF .EQ. 2 ) THEN
             IDIFF = 0
             DO 320 IAEL = 1,NAEL
               IF(IWORK(KLJAE-1+LIA(IAEL)          ).EQ.0) THEN
                 IF( IDIFF .EQ. 0 ) THEN
                   IDIFF = 1
                   I1 = LIA(IAEL)          
                   IPERM = IAEL
                 ELSE
                   I2 = LIA(IAEL)          
                   IPERM = IAEL + IPERM
                   GOTO 321
                 END IF
               END IF
  320        CONTINUE
  321        CONTINUE
*
             JDIFF = 0
             DO 330 JAEL = 1,NAEL
               IF(IWORK(KLIAE-1+JASTRM(JAEL,JASTAC)).EQ.0) THEN
                 IF( JDIFF .EQ. 0 ) THEN
                   JDIFF = 1
                   J1 = JASTRM(JAEL,JASTAC)
                   JPERM = JAEL
                 ELSE
                   J2 = JASTRM(JAEL,JASTAC)
                   JPERM = JAEL + JPERM
                   GOTO 331
                 END IF
               END IF
  330        CONTINUE
  331        CONTINUE
             SIGN = (-1)**(IPERM+JPERM)
           END IF
*
*. Two pairs of differing beta electrons
*
           IF(NBDIF .EQ. 2 ) THEN
             IDIFF = 0
             DO 420 IBEL = 1,NBEL
               IF(IWORK(KLJBE-1+LIB(IBEL)          ).EQ.0) THEN
                 IF( IDIFF .EQ. 0 ) THEN
                   IDIFF = 1
                   I1 = LIB(IBEL)          
                   IPERM = IBEL
                 ELSE
                   I2 = LIB(IBEL)          
                   IPERM = IBEL + IPERM
                   GOTO 421
                  END IF
               END IF
  420        CONTINUE
  421        CONTINUE
*
             JDIFF = 0
             DO 430 JBEL = 1,NBEL
               IF(IWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)).EQ.0) THEN
                 IF( JDIFF .EQ. 0 ) THEN
                   JDIFF = 1
                   J1 = JBSTRM(JBEL,JBSTAC)
                   JPERM = JBEL
                 ELSE
                   J2 = JBSTRM(JBEL,JBSTAC)
                   JPERM = JBEL + JPERM
                   GOTO 431
                 END IF
               END IF
  430        CONTINUE
  431        CONTINUE
             SIGN = (-1)**(IPERM+JPERM)
           END IF
*
* =======================
* Value of matrix element
* =======================
*
        IF( NADIF .EQ. 2 .OR. NBDIF .EQ. 2 ) THEN
* 2 differences in alpha or beta strings
          NDIF2 = NDIF2 + 1
* SIGN * (I1 J1 ! I2 J2 ) - ( I1 J2 ! I2 J1 )
          XVAL = SIGN*( GTIJKL(I1,J1,I2,J2)-GTIJKL(I1,J2,I2,J1) )
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 1 ) THEN
*. 1 difference in alpha strings and one difference in beta string
          NDIF2 = NDIF2 + 1
* SIGN * (IA JA ! IB JB )
          XVAL = SIGNA*SIGNB* GTIJKL(IA,JA,IB,JB)
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 0 .OR.
* 1 differences in alpha or beta strings
     &           NADIF .EQ. 0 .AND. NBDIF .EQ. 1 )THEN
          NDIF1 = NDIF1 + 1
* SIGN *
*(  H(I1 J1 ) +
*  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 ! JORB JORB )
* -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB ! JORB J1 ) )
          IF( NADIF .EQ. 1 ) THEN
            I1 = IA 
            J1 = JA
            SIGN = SIGNA
          ELSE
            I1 = IB
            J1 = JB
            SIGN = SIGNB
          END IF
*
          XVAL = GETH1I(I1,J1)
          DO 520 JAEL = 1, NAEL
            JORB = JASTRM(JAEL,JASTAC)
            XVAL = XVAL + GTIJKL(I1,J1,JORB,JORB)
C?     write(6,*) ' I1 J1 JORB JORB ', I1, J1, JORB, JORB
C?     write(6,*) ' Integral : ', GTIJKL(I1,J1,JORB,JORB) 
C?     write(6,*) 'XVAL ', xval
  520     CONTINUE
          DO 521 JBEL = 1, NBEL
            JORB = JBSTRM(JBEL,JBSTAC)
            XVAL = XVAL + GTIJKL(I1,J1,JORB,JORB)
C?     write(6,*) ' I1 J1 JORB JORB ', I1, J1, JORB, JORB
C?     write(6,*) ' Integral : ', GTIJKL(I1,J1,JORB,JORB) 
C?     write(6,*) 'XVAL ', xval
  521     CONTINUE
          IF( NADIF .EQ. 1 ) THEN
            DO 522 JAEL = 1, NAEL
              JORB = JASTRM(JAEL,JASTAC)
              XVAL = XVAL - GTIJKL(I1,JORB,JORB,J1)
C?     write(6,*) ' I1  JORB JORB J1', I1,  JORB, JORB, J1
C?     write(6,*) ' Integral : ', GTIJKL(I1,JORB,JORB,J1) 
C?     write(6,*) 'XVAL ', xval
  522       CONTINUE
          ELSE
            DO 523 JBEL = 1, NBEL
              JORB = JBSTRM(JBEL,JBSTAC)
              XVAL = XVAL - GTIJKL(I1,JORB,JORB,J1)
C?     write(6,*) ' I1  JORB JORB J1', I1,  JORB, JORB, J1
C?     write(6,*) ' Integral : ', GTIJKL(I1,JORB,JORB,J1) 
C?     write(6,*) 'XVAL ', xval
  523       CONTINUE
          END IF
          XVAL = XVAL * SIGN
        ELSE IF( NADIF .EQ. 0 .AND. NBDIF .EQ. 0 ) THEN
*. Diagonal elements
          NDIF0 = NDIF0 + 1
C SUM(I,J OF JDET) H(I,J) + (I I ! J J ) - (I J ! J I )
C
          XVAL = ECORE
          DO 650 IAB = 1, 2
            IF(IAB .EQ. 1 ) THEN
              NIABEL = NAEL
            ELSE
              NIABEL = NBEL
            END IF
            DO 640 JAB = 1, 2
              IF(JAB .EQ. 1 ) THEN
                NJABEL = NAEL
              ELSE
                NJABEL = NBEL
              END IF
              DO 630 IEL = 1, NIABEL
                IF( IAB .EQ. 1 ) THEN
                  IORB = LIA(IEL)          
                ELSE
                  IORB = LIB(IEL)          
                END IF
                IF(IAB .EQ. JAB ) XVAL = XVAL + GETH1I(IORB,IORB)
C?     write(6,*) ' XVAL ', XVAL 
                DO 620 JEL = 1, NJABEL
                  IF( JAB .EQ. 1 ) THEN
                    JORB = LIA(JEL) 
                  ELSE
                    JORB = LIB(JEL)           
                  END IF
                  XVAL = XVAL + 0.5D0*GTIJKL(IORB,IORB,JORB,JORB)
*. test
C?        FAC = GTIJKL(IORB,IORB,JORB,JORB)
C?        write(6,*) ' IORB JORB (IORB IORB ! JORB JORB )'
C?        WRITE(6,*) IORB,JORB,FAC
C?     write(6,*) ' XVAL ', XVAL 
    
                  IF( IAB . EQ. JAB )
     &            XVAL = XVAL - 0.5D0*GTIJKL(IORB,JORB,JORB,IORB)
*. test
          FAC = GTIJKL(IORB,JORB,JORB,IORB)
C?        write(6,*) ' IORB JORB (IORB JORB ! JORB IORB )'
C?        WRITE(6,*) IORB,JORB,FAC
C?     write(6,*) ' XVAL ', XVAL 
  620           CONTINUE
  630         CONTINUE
  640       CONTINUE
  650     CONTINUE
        END IF
 
        IF( NTEST .GE. 100 ) 
     &  WRITE(6,*) ' CONST XVAL  ', CONST ,XVAL
        IF( ISYM .EQ. 0 ) THEN
          HAMIL((JDET-1)*NIDET+IDET) =
     &    HAMIL((JDET-1)*NIDET+IDET) + CONST * XVAL
        ELSE
          HAMIL((IDET-1)*IDET/2 + JDET ) =
     &    HAMIL((IDET-1)*IDET/2 + JDET ) + CONST * XVAL
        END IF
  898 CONTINUE
  899 CONTINUE
  900 CONTINUE
 1000 CONTINUE
 
      IF( IPRT .GT. 0 ) THEN
      WRITE(6,*)
     &'  Number of elements differing by 0 excitation.. ',NDIF0
 
      WRITE(6,*)
     &'  Number of elements differing by 1 excitation.. ',NDIF1
 
      WRITE(6,*)
     &'  Number of elements differing by 2 excitation.. ',NDIF2
 
      WRITE(6,*)
     &'  Number of vanishing elments                    ',
     &   NTERMS - NDIF0 - NDIF1 - NDIF2
      END IF
      IF( IPRT .GE. 10) THEN
        WRITE(6,*) '  Hamiltonian matrix '
        IF( ISYM .EQ. 0 ) THEN
          CALL WRTMAT(HAMIL,NIDET,NJDET,NIDET,NJDET)
        ELSE
          CALL PRSYM(HAMIL,NIDET)
        END IF
      END IF
C
C!    STOP ' ENFORCED STOP AT END OF DIHDJ '
      RETURN
      END
      SUBROUTINE DIHDJ2(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,
     &IWORK,LWORK,NORB,ONEBOD,HAMIL,ISYM,NINOB,ECORE,ICOMBI,PSIGN,
     &IPRT,NTERMS,NDIF0,NDIF1,NDIF2)
C
C A SET OF DETERMINANTS IA DEFINED BY ALPHA AND BETASTRINGS
C IASTR,IBSTR AND ANOTHER SET OF DETERMINATS DEFINED BY STRINGS
C JASTR AND JBSTR ARE GIVEN . OBTAIN CORRESPONDING HAMILTONIAN MATRIX
C
C IF ICOMBI .NE. 0 COMBINATIONS ARE USED FOR ALPHA AND BETA STRING
C THAT DIFFERS :
C   1/SQRT(2) * ( !I1A I2B! + PSIGN * !I2A I1B! )
C
C IF ISYM .EQ. 0 FULL HAMILTONIAN IS CONSTRUCTED
C IF ISYM .NE. 0 LOWER HALF OF HAMILTONIAN IS CONSTRUCTED
C
C JEPPE OLSEN JANUARY 1989
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION JASTR(NAEL,*),JBSTR(NBEL,*)
C
      DIMENSION IWORK(*), HAMIL(*), ONEBOD(NORB,NORB)
C
      NTEST = 00
*. To eliminate compiler warnings
      KLIAB = 0
      IAEQIB = 0
      JAEQJB = 0
      IEL1 = 0
      JEL1 = 0
      SIGNA = 0.0D0
      SIGNB = 0.0D0
      IPERM = 0   
      JPERM = 0
      SIGN = 0.0D0


C?    write(6,*) ' NIDET NJDET NORB ', NIDET,NJDET,NORB
C     write(6,*) ' One body matrix '
C     CALL WRTMAT(ONEBOD,NORB,NORB,NORB,NORB)
C SCATCH SPACE
C
C .. 1 : EXPANSION OF ALPHA AND BETA STRINGS OF TYPE I
C
C?    WRITE(6,*) ' MISTER DIHDJ SPEAKING  ICOMBI = ', ICOMBI
      KLFREE = 1
      KLIAE  = KLFREE
      KLFREE = KLIAE + NORB
      KLIBE  = KLFREE
      KLFREE = KLIBE + NORB
C
      KLJAE = KLFREE
      KLFREE = KLJAE + NORB
      KLJBE = KLFREE
      KLFREE = KLJBE + NORB
      IF( ICOMBI .NE. 0 ) THEN
        KLIAB  = KLFREE
        KLFREE = KLFREE + NIDET
      END IF
C?    WRITE (6,*) ' PSIGN ', PSIGN
C
      IF( ICOMBI .NE. 0 ) THEN
C SET UP ARRAY COMBARING ALPHA AND BETA STRINGS IN IDET LIST
        DO 13 IDET = 1, NIDET
          IAEQIB = 1
          DO 14 IEL = 1, NAEL
            IF(IASTR(IEL,IDET) .NE. IBSTR(IEL,IDET))IAEQIB = 0
   14     CONTINUE
          IWORK(KLIAB-1+IDET) = IAEQIB
   13   CONTINUE
C?      WRITE(6,*) ' IAEQIB ARRAY : '
C?      CALL IWRTMA(IWORK(KLIAB),1,NIDET,1,NIDET)
      END IF
C
      IF( ISYM .EQ. 0 ) THEN
        LHAMIL = NIDET*NJDET
      ELSE
        LHAMIL = NIDET*(NIDET+1) / 2
      END IF
      CALL SETVEC(HAMIL,0.0D0,LHAMIL)
C
COLD  NTERMS= 0
COLD  NDIF0 = 0
COLD  NDIF1 = 0
COLD  NDIF2 = 0
C.. LOOP OVER J DETERMINANTS
C
      DO 1000 JDET = 1,NJDET
        IF( NTEST .GE. 20 ) WRITE(6,*) '  ****** 1000 JDET ', JDET
C
C EXPAND JDET
        CALL ISETVC(IWORK(KLJAE),0,NORB)
        CALL ISETVC(IWORK(KLJBE),0,NORB)
C
        IF( ICOMBI .NE. 0 ) THEN
          JAEQJB = 1
          DO 32 IEL = 1, NAEL
            IF(JASTR(IEL,JDET) .NE. JBSTR(IEL,JDET))JAEQJB = 0
   32     CONTINUE
C?        WRITE(6,*) ' JAEQJB ', JAEQJB
        END IF
C
        DO 40 IAEL = 1, NAEL
          IWORK(KLJAE-1+JASTR(IAEL,JDET) ) = 1
   40   CONTINUE
C
        DO 50 IBEL = 1, NBEL
          IWORK(KLJBE-1+JBSTR(IBEL,JDET) ) = 1
   50   CONTINUE
C
        IF( NTEST .GE. 10 ) THEN
          WRITE(6,*) ' LOOP 1000 JDET =  ',JDET
          WRITE(6,*) ' JASTR AND JBSTR '
          CALL IWRTMA(JASTR(1,JDET),1,NAEL,1,NAEL)
          CALL IWRTMA(JBSTR(1,JDET),1,NBEL,1,NBEL)
          WRITE(6,*) ' EXPANDED ALPHA AND BETA STRING '
          CALL IWRTMA(IWORK(KLJAE),1,NORB,1,NORB)
          CALL IWRTMA(IWORK(KLJBE),1,NORB,1,NORB)
        END IF
C
        IF( ISYM .EQ. 0 ) THEN
          MINI = 1
        ELSE
          MINI = JDET
        END IF
        DO 900 IDET = MINI, NIDET
C?      WRITE(6,*) '   LOOP 900 IDET .... ',IDET
        IF( ICOMBI .EQ. 0 ) THEN
            NLOOP = 1
        ELSE
          IAEQIB = IWORK(KLIAB-1+IDET)
          IF(IAEQIB+JAEQJB .EQ. 0 ) THEN
            NLOOP = 2
          ELSE
            NLOOP = 1
          END IF
        END IF
C
        DO 899 ILOOP = 1, NLOOP
         NTERMS = NTERMS + 1
C?       WRITE(6,*) '   899 : ILOOP ' , ILOOP
C
C.. COMPARE DETERMINANTS
C
C SWAP IA AND IB FOR SECOND PART OF COMBINATIONS
        IF( ILOOP .EQ. 2 )
     &  CALL ISWPVE(IASTR(1,IDET),IBSTR(1,IDET),NAEL)
C
        NACM = 0
        DO 61 IAEL = 1, NAEL
          NACM = NACM + IWORK(KLJAE-1+IASTR(IAEL,IDET))
   61   CONTINUE
C
        NBCM = 0
        DO 62 IBEL = 1, NBEL
          NBCM = NBCM + IWORK(KLJBE-1+IBSTR(IBEL,IDET))
   62   CONTINUE
C
        NADIF = NAEL-NACM
        NBDIF = NBEL-NBCM
        IF( NTEST .GE. 10 ) THEN
          WRITE(6,*) '  LOOP 900 IDET ',IDET
          WRITE(6,*) ' COMPARISON , NADIF , NBDIF ', NADIF,NBDIF
        END IF
C
        IF(NADIF+NBDIF .GT. 2 ) GOTO 899
C
C
C  FACTOR FOR COMBINATIONS
       IF( ICOMBI .EQ. 0 ) THEN
         CONST = 1.0D0
       ELSE
         IF((JAEQJB +IAEQIB) .EQ.2 ) THEN
           CONST = 1.0D0
         ELSE IF( (JAEQJB+IAEQIB) .EQ. 1 ) THEN
           CONST = 1.0D0/SQRT(2.0D0)*(1.0D0+PSIGN)
          ELSE IF( (JAEQJB+IAEQIB) .EQ. 0 ) THEN
           IF( ILOOP .EQ. 1)  THEN
             CONST = 1.0D0
           ELSE
             CONST = PSIGN
           END IF
          END IF
        END IF
C
C.. FIND DIFFERING ORBITALS AND SIGN FOR PERMUTATION
C
C EXPAND IDET
        CALL ISETVC(IWORK(KLIAE),0,NORB)
        CALL ISETVC(IWORK(KLIBE),0,NORB)
C
          DO 42 IAEL = 1, NAEL
            IWORK(KLIAE-1+IASTR(IAEL,IDET) ) = 1
   42     CONTINUE
C
          DO 52 IBEL = 1, NBEL
            IWORK(KLIBE-1+IBSTR(IBEL,IDET) ) = 1
   52     CONTINUE
C
        IF(NADIF .EQ. 1 ) THEN
          DO 120 IAEL = 1,NAEL
            IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
              IA = IASTR(IAEL,IDET)
              IEL1 = IAEL
              GOTO 121
             END IF
  120     CONTINUE
  121     CONTINUE
C
          DO 130 JAEL = 1,NAEL
            IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
              JA = JASTR(JAEL,JDET)
              JEL1 = JAEL
              GOTO 131
             END IF
  130     CONTINUE
  131     CONTINUE
          SIGNA = (-1)**(JEL1+IEL1)
C?        WRITE(6,*) ' IA JA SIGNA... ',IA,JA,SIGNA
 
 
C
        END IF
 
        IF(NBDIF .EQ. 1 ) THEN
          DO 220 IBEL = 1,NBEL
            IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
              IB = IBSTR(IBEL,IDET)
              IEL1 = IBEL
              GOTO 221
             END IF
  220     CONTINUE
  221     CONTINUE
C
          DO 230 JBEL = 1,NBEL
            IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
              JB = JBSTR(JBEL,JDET)
              JEL1 = JBEL
              GOTO 231
             END IF
  230     CONTINUE
  231     CONTINUE
          SIGNB = (-1)**(JEL1+IEL1)
C?        WRITE(6,*) ' IB JB SIGNB... ',IB,JB,SIGNB
C
        END IF
        IF(NADIF .EQ. 2 ) THEN
          IDIFF = 0
          DO 320 IAEL = 1,NAEL
            IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
              IF( IDIFF .EQ. 0 ) THEN
                IDIFF = 1
                I1 = IASTR(IAEL,IDET)
                IPERM = IAEL
              ELSE
                I2 = IASTR(IAEL,IDET)
                IPERM = IAEL + IPERM
                GOTO 321
              END IF
            END IF
  320     CONTINUE
  321     CONTINUE
C
          JDIFF = 0
          DO 330 JAEL = 1,NAEL
            IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
              IF( JDIFF .EQ. 0 ) THEN
                JDIFF = 1
                J1 = JASTR(JAEL,JDET)
                JPERM = JAEL
              ELSE
                J2 = JASTR(JAEL,JDET)
                JPERM = JAEL + JPERM
                GOTO 331
              END IF
            END IF
  330     CONTINUE
  331     CONTINUE
          SIGN = (-1)**(IPERM+JPERM)
C
        END IF
C
        IF(NBDIF .EQ. 2 ) THEN
          IDIFF = 0
          DO 420 IBEL = 1,NBEL
            IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
              IF( IDIFF .EQ. 0 ) THEN
                IDIFF = 1
                I1 = IBSTR(IBEL,IDET)
                IPERM = IBEL
              ELSE
                I2 = IBSTR(IBEL,IDET)
                IPERM = IBEL + IPERM
                GOTO 421
               END IF
            END IF
  420     CONTINUE
  421     CONTINUE
C
          JDIFF = 0
          DO 430 JBEL = 1,NBEL
            IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
              IF( JDIFF .EQ. 0 ) THEN
                JDIFF = 1
                J1 = JBSTR(JBEL,JDET)
                JPERM = JBEL
              ELSE
                J2 = JBSTR(JBEL,JDET)
                JPERM = JBEL + JPERM
                GOTO 431
              END IF
            END IF
  430     CONTINUE
  431     CONTINUE
          SIGN = (-1)**(IPERM+JPERM)
C
        END IF
C
C OBTAIN VALUE OF HAMILTONIAN ELEMENT
C
        IF( NADIF .EQ. 2 .OR. NBDIF .EQ. 2 ) THEN
          NDIF2 = NDIF2 + 1
C SIGN * (I1 J1 ! I2 J2 ) - ( I1 J2 ! I2 J1 )
          I1 = I1 + NINOB
          I2 = I2 + NINOB
          J1 = J1 + NINOB
          J2 = J2 + NINOB
          XVAL = SIGN*( GTIJKL(I1,J1,I2,J2)-GTIJKL(I1,J2,I2,J1) )
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 1 ) THEN
          NDIF2 = NDIF2 + 1
C SIGN * (IA JA ! IB JB )
          IA = IA + NINOB
          IB = IB + NINOB
          JA = JA + NINOB
          JB = JB + NINOB
C?        WRITE(6,*) ' IA IB JA JB ', IA,IB,JA,JB
          XVAL = SIGNA*SIGNB* GTIJKL(IA,JA,IB,JB)
C?        WRITE(6,*) ' SIGNA SIGNB XVAL ',SIGNA,SIGNB,XVAL
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 0 .OR.
     &           NADIF .EQ. 0 .AND. NBDIF .EQ. 1 )THEN
          NDIF1 = NDIF1 + 1
C SIGN *
C(  H(I1 J1 ) +
C  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 ! JORB JORB )
C -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB ! JORB J1 ) )
          IF( NADIF .EQ. 1 ) THEN
            I1 = IA + NINOB
            J1 = JA + NINOB
            SIGN = SIGNA
          ELSE
            I1 = IB + NINOB
            J1 = JB + NINOB
            SIGN = SIGNB
          END IF
C?        WRITE(6,*) ' ONE DIFF I1 J1 SIGN : ',I1,J1,SIGN
C
          XVAL = GETH1I(I1-NINOB,J1-NINOB)
          DO 520 JAEL = 1, NAEL
            JORB = JASTR(JAEL,JDET)+NINOB
            XVAL = XVAL + GTIJKL(I1,J1,JORB,JORB)
  520     CONTINUE
          DO 521 JBEL = 1, NBEL
            JORB = JBSTR(JBEL,JDET)+NINOB
            XVAL = XVAL + GTIJKL(I1,J1,JORB,JORB)
  521     CONTINUE
          IF( NADIF .EQ. 1 ) THEN
            DO 522 JAEL = 1, NAEL
              JORB = JASTR(JAEL,JDET)+NINOB
              XVAL = XVAL - GTIJKL(I1,JORB,JORB,J1)
  522       CONTINUE
          ELSE
            DO 523 JBEL = 1, NBEL
              JORB = JBSTR(JBEL,JDET)+NINOB
              XVAL = XVAL - GTIJKL(I1,JORB,JORB,J1)
  523       CONTINUE
          END IF
          XVAL = XVAL * SIGN
        ELSE IF( NADIF .EQ. 0 .AND. NBDIF .EQ. 0 ) THEN
          NDIF0 = NDIF0 + 1
C SUM(I,J OF JDET) H(I,J) + (I I ! J J ) - (I J ! J I )
C
          XVAL = ECORE
          DO 650 IAB = 1, 2
            IF(IAB .EQ. 1 ) THEN
              NIABEL = NAEL
            ELSE
              NIABEL = NBEL
            END IF
            DO 640 JAB = 1, 2
              IF(JAB .EQ. 1 ) THEN
                NJABEL = NAEL
              ELSE
                NJABEL = NBEL
              END IF
              DO 630 IEL = 1, NIABEL
                IF( IAB .EQ. 1 ) THEN
                  IORB = IASTR(IEL,IDET)
                ELSE
                  IORB = IBSTR(IEL,IDET)
                END IF
                IF(IAB .EQ. JAB ) XVAL = XVAL + GETH1I(IORB,IORB)
                write(6,*)  ' XVAL after one body term ', XVAL
                IORB = IORB + NINOB
                DO 620 JEL = 1, NJABEL
                  IF( JAB .EQ. 1 ) THEN
                    JORB = IASTR(JEL,IDET)+NINOB
                  ELSE
                    JORB = IBSTR(JEL,IDET)+NINOB
                  END IF
                  XVAL = XVAL + 0.5D0*GTIJKL(IORB,IORB,JORB,JORB)
        write(6,*) ' XVAL ', XVAL
                  IF( IAB . EQ. JAB )
     &            XVAL = XVAL - 0.5D0*GTIJKL(IORB,JORB,JORB,IORB)
        write(6,*) ' XVAL ', XVAL
  620           CONTINUE
  630         CONTINUE
  640       CONTINUE
  650     CONTINUE
        END IF
C
        WRITE(6,*) ' CONST XVAL  ', CONST ,XVAL
        IF( ISYM .EQ. 0 ) THEN
          HAMIL((JDET-1)*NIDET+IDET) =
     &    HAMIL((JDET-1)*NIDET+IDET) + CONST * XVAL
        ELSE
          HAMIL((IDET-1)*IDET/2 + JDET ) =
     &    HAMIL((IDET-1)*IDET/2 + JDET ) + CONST * XVAL
        END IF
C RESTORE ORDER !!!
        IF( ILOOP .EQ. 2 )
     &  CALL ISWPVE(IASTR(1,IDET),IBSTR(1,IDET),NAEL)
  899 CONTINUE
  900 CONTINUE
 1000 CONTINUE
 
C     IF( IPRT .GT. 0 ) THEN
C     WRITE(6,*)
C    &'  Number of elements differing by 0 excitation.. ',NDIF0
c
C     WRITE(6,*)
C    &'  Number of elements differing by 1 excitation.. ',NDIF1
c
C     WRITE(6,*)
C    &'  Number of elements differing by 2 excitation.. ',NDIF2
c
C     WRITE(6,*)
C    &'  Number of vanishing elments                    ',
C    &   NTERMS - NDIF0 - NDIF1 - NDIF2
C     END IF
      IF( IPRT .GE. 2 ) THEN
        WRITE(6,*) '  HAMILTONIAN MATRIX '
        IF( ISYM .EQ. 0 ) THEN
          CALL WRTMAT(HAMIL,NIDET,NJDET,NIDET,NJDET)
        ELSE
          CALL PRSYM(HAMIL,NIDET)
        END IF
      END IF
C
C!    STOP ' ENFORCED STOP AT END OF DIHDJ '
      RETURN
      END
      SUBROUTINE H0DIAG(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,NROOT,
     &           EIGVAL,EIGVEC,SCR,NTESTG,ECORE)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Obtain the lowest NROOT eigenvectors
*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* NROOT : Number of roots to be obtained
*
* ======
* Output
* ======
* EIGVAL(*) : Converged eigen values
* EIGVEC(IROOT,*) : Complete eigenvector IROOT
*
* Note : The NROOT eigenpairs to be obtained are assumed
*        to be ' slightly ' perturbed eigensolutions
*        of PHP.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL HPQTVM
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*)
* Output
      DIMENSION EIGVAL(*),EIGVEC(NP1DM+NP2DM+NQDM,*)
* Scratch
      DIMENSION SCR(*)
*.SCR Should atleast be dimensioned ??????
      LOGICAL CONVER
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
* There are two routes :2 Iterative diagonalization of complete matrix
*                       1 complete diagonalizations of partitioned
*                         matrices.
*.
* ========
* Route 1 :
* ========
* The Q-space can be partitioned into the P -space
* to give the effective eigenvalue equation
*
* (PHP - PHQ  (QHQ-E)**-1 QHP ) VP = E VP
*
* This leads to a simple iterative scheme

*
* Note : Only NROOT = 1 tested
*
      NTESTL = 1  
      NTEST = MAX(NTESTG,NTESTL)
      CALL LUCIAQENTER('H0DIA')
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0DIAG speaking '
        WRITE(6,*) ' =============== '
      END IF
*
C?    write(6,*) ' QHQ as delivered '
C?    call wrtmat(QHQ,1,NQDM,1,NQDM)
C?    write(6,*) ' PHP as delivered '
C?    call PRSYM(PHP,NP1DM)
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
*. A bit of memory
*
      I12 = 2
      IF(I12.EQ.1. OR. I12. EQ.3 ) THEN
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*
 
      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. A PQ matrix
      KLPQ = KLFREE
      KLFREE = KLFREE + NPDM * NQDM
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
*
*. Initial eigenvalues
      CALL COPVEC(PHP,SCR(KLPP1),NPDM*(NPDM+1)/2)
      CALL EIGENL(SCR(KLPP1),SCR(KLPP2),NPDM,0,1)
*. Extract eigenvalues
      CALL XTRCDI(SCR(KLPP1),EIGVAL,NROOT,1)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Initial set of eigenvalues '
        CALL WRTMAT(EIGVAL,1,NROOT,1,NROOT)
      END IF
*. Largest allowed number of iterations
      MAXIT = 5
      DO 1000 IROOT = 1, NROOT
        CONVER = .FALSE.
        EINI = EIGVAL(IROOT)
        DO 900 ITER = 1, MAXIT
          IF(NTEST.GE.2) WRITE(6,*) ' Info from iteration ', ITER
*. Current eigenvalue and eigenvector
          E = EIGVAL(IROOT)
          CALL COPVEC(SCR(KLPP2+(IROOT-1)*NPDM),SCR(KLV1),NPDM)
* ==============================
* HPP - PHQ (QHQ - E) **-1 * QHP
* ==============================
*. QHP in KLPQ
C         TRPMT3(XIN,NROW,NCOL,XOUT)
          CALL TRPMT3(PHQ,NP1DM,NQDM,SCR(KLPQ))
*.Multiply with (QHQ-E)**-1
          DO 30 IP1 = 1, NP1DM
C           DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
            IOFF = KLPQ + (IP1-1)*NQDM
            CALL DIAVC3(SCR(IOFF),SCR(IOFF),QHQ,-E,NQDM,XDUMMY)
   30     CONTINUE
*. Multiply with PHQ
          CALL MATML4(SCR(KLPP1),PHQ,SCR(KLPQ),NP1DM,NP1DM,
     &                NP1DM,NQDM,NQDM,NP1DM,0)
*.
C?      write(6,*) ' PHQ (QHQ-E)-1 QHP matrix '
C?      CALL WRTMAT(SCR(KLPP1),NP1DM,NP1DM,NP1DM,NP1DM)
 
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          CALL TRIPAK(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL PRSYM(SCR(KLPP1),NPDM)
          END IF
*.Diagonalize
          CALL EIGENL(SCR(KLPP1),SCR(KLPP2),NPDM,0,1)
*. Extract eigenvalues
           EIGVAL(IROOT) = SCR(KLPP1-1+IROOT*(IROOT+1)/2)
           IF(NTEST.GE.2)
     &     WRITE(6,*) ' Eigenvalue ', EIGVAL(IROOT)
           IF(NTEST.GE.10) THEN
             WRITE(6,*) ' P-space eigenvector '
             CALL WRTMAT(SCR(KLPP2+(IROOT-1)*NPDM),
     &            1,NPDM,1,NPDM)
           END IF
*. Check for convergence
           EVALDF = ABS(E-EIGVAL(IROOT))
           EVECOV = SQRT(INPROD(SCR(KLV1),SCR(KLPP2+(IROOT-1)*NPDM),
     &                          NPDM) )
           IF(EVALDF.LT.1.0D-7.AND.EVECOV.GT.0.999D0) CONVER = .TRUE.
           IF(CONVER) GOTO 901
  900    CONTINUE
  901    CONTINUE
*. P-part of eigenvector
         CALL COPVEC(SCR(KLPP2+(IROOT-1)*NPDM),EIGVEC(1,IROOT),
     &               NPDM)
*. Obtain Q part of eigenvector
*.    -(QHQ-E)**-1 HQP XP
 
         CALL MATML4(SCR(KLV1),PHQ,EIGVEC(1,IROOT),
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
*
C             DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3(EIGVEC(NPDM+1,IROOT),SCR(KLV1),QHQ,
     &               -EIGVAL(IROOT),NQDM,XDUMMY)
         CALL SCALVE(EIGVEC(NPDM+1,IROOT),-1.0D0,NQDM)
*. Normalize
         XNORM = INPROD(EIGVEC(1,IROOT),EIGVEC(1,IROOT),NPQDM)
         SCALE = 1.0D0/SQRT(XNORM)
         CALL SCALVE(EIGVEC(1,IROOT),SCALE,NPQDM)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Initial and final eigenvalue ',
     &     EINI,EIGVAL(IROOT)
           WRITE(6,*)
     &   ' Part of eigenvector in Q space',
     &     SQRT(ABS(XNORM-1.0D0)/XNORM)
C?         WRITE(6,*) ' Eigenvector in PQ space '
C?         CALL WRTMAT(EIGVEC(1,IROOT),1,NPQDM,1,NPQDM)
         END IF
*
 1000 CONTINUE
       END IF    
       IF( I12. EQ.2 .OR. I12 .EQ.3 ) THEN
*
*. Iterative scheme
*
*
*. Initial eigenvalues and eigenvectors
*
       KLPP1 = 1
       KLFREE = KLPP1 + NP1DM*(NP1DM+1)/2
       KLPP2 = KLFREE
       KLFREE = KLPP2 + NP1DM*NP1DM
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LU4 = 37
        LU5 = 38
        LUDIA = 39
        CALL COPVEC(PHP,SCR(KLPP1),NP1DM*(NP1DM+1)/2)
        CALL EIGENL(SCR(KLPP1),SCR(KLPP2),NP1DM,0,1)
*. Extract eigenvalues and eigenvectors on LU1
        CALL XTRCDI(SCR(KLPP1),EIGVAL,NROOT,1)
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Initial set of eigenvalues '
          CALL WRTMAT(EIGVAL,1,NROOT,1,NROOT)
        END IF
*. Eigenvectors on LU1
        CALL SETVEC(EIGVEC(1,1),0.0D0,NPQDM)
        CALL REWINE(LU1,-1)
        DO 510 IROOT = 1, NROOT
          CALL COPVEC(SCR(KLPP2+(IROOT-1)*NP1DM),EIGVEC(1,1),NP1DM)
          CALL TODSC(EIGVEC(1,1),NPQDM,-1,LU1)
  510   CONTINUE
*
*. Iterations
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LUDIA)
*. Davidson CI diagonalization
        MINST = 1
        NBLK = 1
        INICI = -1
        MAXCIT = 15
        IPRTCI = MAX(NTEST-2,0)
        MXVCCI = 3 * NROOT
        SHIFT = 0.0D0  
        IDIAG = 1
        ICISTR = 1
        IDIAG = 1
        THRES= 1.0D-8
        CALL CIEIG5(HPQTVM,INICI,EIGVAL,SCR(KLV1),SCR(KLV2),
     &            MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,NPQDM ,
     &            NBLK,NROOT,MXVCCI,MAXCIT,LU1,IPRTCI,
     &            DUMMY,0,DUMMY,IDUMMY,
     &            0,0,0,DUMMY,ECORE,ICISTR,NPQDM,IDIAG,DUMMY,THRES)
C       CALL CIEIG5(HPQTVM,INICI,ENOT,SCR(KLV1),SCR(KLV2),
C    &            MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,NPQDM ,
C    &            NBLK,NROOT,MXVCCI,MAXCIT,LU1,IPRTCI,
C    &            DUMMY,0,DUMMY,IDUMMY,
C    &            0,0,0,DUMMY,0.0D0,1,NPQDM,1,DUMMY)
        CALL REWINE(LU1,-1)
        DO 1286 JROOT = 1, NROOT
         CALL FRMDSC(EIGVEC(1,JROOT),NPQDM,-1,LU1,IMZERO,IAMPACK)
 1286   CONTINUE
      END IF
*
      CALL LUCIAQEXIT('H0DIA')
      CALL lUCIAQSTAT
*
      RETURN
      END
      SUBROUTINE H0MAT(INTSPC,NROOT,ONEBOD,H0,SBEVC,SBEVL,ISBDET,
     &                 ISBIA,ISBIB,ISBCNF,
     &                 LUHDIA,LBLK,
     &                 MXP1,MXP2,MXQ,
     &                 MP1CSF,MP2CSF,MQCSF,NOCOB,
     &                 NPRCIV,NOCSF,IREFSM,IPRT,IPROCC,
     &                 VEC1,VEC2,H0SCR,IDC,PSSIGN,ECORE)
* Obtain preconditioner space corresponding to internalt space INTSPC
* Obtain Hamiltonian matrices correponding to this subspacw
* Obtain the first Nroot eigensolutions.
*
* Construct Preconditioner blocks of Hamilton matrix
*
* ======
*.Output
* ======
*
* CSF : NP1CSF,NP2CSF,NQCSF : Number of CSF's in the 3 primary subspaces
*
* NPRCIV : Number of parameters in preconditioner space
*
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csfbas.inc'
      COMMON/SPINFO_OLD/MULTSP,MS2P,
     &              MINOP,MAXOP,NTYP,NDPCNT(MXPCTP),NCPCNT(MXPCTP),
     &              NCNATS(MXPCTP,MXPCSM),NDTASM(MXPCSM),NCSASM(MXPCSM),
     &              NCNASM(MXPCSM)
*
      DIMENSION ONEBOD(*)
*
      NTEST = 0
      NTEST = MAX(IPRT,NTEST)
*
      MXSBDT = MXP1+MXP2+MXQ
      IF(NTEST.NE.0.AND.MXSBDT.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================================== '
        WRITE(6,*) '    Information about CI preconditioner '
        WRITE(6,*) ' ========================================== '
        WRITE(6,*)
       END IF
       WRITE(6,*) ' H0MAT : ecore ', ECORE
* Info on actual internal subspace
      IATP = IASTFI(INTSPC)
      IBTP = IBSTFI(INTSPC)
      MNR1 = MNR1IC(INTSPC)
      MXR1 = MXR1IC(INTSPC)
      MNR3 = MNR3IC(INTSPC)
      MXR3 = MXR3IC(INTSPC)
      NAEL = NAELCI(INTSPC)
      NBEL = NBELCI(INTSPC)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Allowed combination of alpha and beta strings
      CALL MEMMAN(KIOCOC,NOCTPA*NOCTPB,'ADDL  ',2,'IOCOC ')
      STOP ' update call to IAIBCM_GAS '
      CALL IAIBCM_GAS(MNR1,MXR3,NOCTPA,NOCTPB,WORK(KEL1(IATP)),
     &            WORK(KEL3(IATP)),WORK(KEL1(IBTP)),WORK(KEL3(IBTP)),
     &            WORK(KIOCOC),NTEST)
*
      IF(IDC.EQ.1) THEN
        ICOMBI = 0
        PSIGN = 0.0D0
      ELSE
        PSIGN = PSSIGN
        ICOMBI = 1
      END IF
*
      IF( NOCSF .NE. 0) THEN
*.Combinations expansion, PQ preconditioner
*
        IHAMSM = 1
        IWAY = 1
* strings are unsigned
        ISTSGN = 0
        CALL H0SD(LUHDIA,LBLK,VEC1,IWAY,NSBDET,NAEL,NBEL,
     &            ISMOST(1,IREFSM),WORK(KIOCOC),
     &            IHAMSM,H0,ONEBOD,NOCOB,0,
     &            ECORE,ICOMBI,PSIGN,NPRCIV,SBEVC,
     &            SBEVL,1,NCIVAR,ISBDET,ISBIA,ISBIB,NROOT,
     &            MXP1,MXP2,MXQ,
     &            MP1CSF,MP2CSF,MQCSF,
     &            WORK(KOCSTR(IATP)),WORK(KOCSTR(IBTP)),
     &            ISTSGN,IDUMMY,IDUMMY,
     &            INTSPC,IPRT,IPROCC)
      END IF
*
      IF(IPRT.NE.0 .AND. NPRCIV.LT.NSBDET) THEN
        NSBDET = NPRCIV
        WRITE(6,*)
     &  '  Number of elements in primary space REDUCED to . ',NSBDET
      END IF
*
      IF(NTEST.NE.0.AND.MXSBDT.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================================== '
        WRITE(6,*) '    End of information about CI preconditioner '
        WRITE(6,*) ' ============================================== '
        WRITE(6,*)
       END IF
*
      RETURN
      END 
      SUBROUTINE H0SD(LUDIA,LBLK,CIDIA,IWAY,MXPRDT,NAEL,NBEL,
     &                ISMOST,IOCOC,
     &                IHSYM,HAMIL,ONEBOD,NORB,
     &                NINOB,ECORE,ICOMBI,PSIGN,NPRDET,EIGVEC,
     &                EIGVAL,IDODGN,NDET,IDET,IA,IB,NROOT,
     &                MXP1DM,MXP2DM,MXQDM,NP1DM,NP2DM,NQDM,
     &                IASTR,IBSTR,ISTSGN,IASGN,IBSGN,
     &                INTSPC,IPRT,IPROCC)
*
* Construct an P1 P2 Q preconditioner matrix in the SD basis
* and/or find and print the lowest IPROCC elements
*
*.. Subspace :
* ============
*     
*     IWAY = 1 : choose the lowest values of
*                the CI diagonal. The number of combinations is
*                choosen so no degenerate levels are splitted .
*                the number of combinations used ,NPRDET, can
*                therefore be lower than mxsbdt.
*     IWAY = 2 : Choose the first NPRDET combinations,
*                stupid, but convenient for testing
*
* IDET contains adresses of elements chosen
*
* IDODGN .GT. 0 : DIAGONALIZE CONSTRUCTED HAMILTON MATRIX.
*                 EIGVL CONTAINS EIGENVALUES ON RETURN
*                 EIGVC CONTAINS EIGENVECTORS(COLUMNS) ON RETURN
*
* Obtained from EXPHAM, January 1993 
*
* Put into LUCIA form June 1993
* Combinations enabled September 1993
*
*. Core energy : not added to matrix, added to resulting eigenvalues
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. Input
      DIMENSION ISMOST(*)
      DIMENSION ONEBOD(*),CIDIA(*)
      INTEGER IASTR(NAEL,*),IBSTR(NBEL,*)
*. Output
      DIMENSION HAMIL(*)
      DIMENSION IDET(MXPRDT), EIGVAL(NROOT), EIGVEC(NROOT * MXPRDT )
      DIMENSION IA(MXPRDT),IB(MXPRDT)
*
*
      CALL LUCIAQENTER('EXPHAM')
      NTEST = 0
      NTEST = MAX(NTEST,IPRT)
      MXDM1 = MXP1DM + MXP2DM + MXQDM
      MXDM = MAX(MXDM1,IPROCC)
      IF (NTEST .GE. 5 ) THEN 
        WRITE(6,*) ' MXDM ', MXDM
        write(6,*) ' MXP1DM MXP2DM MXQDM ', MXP1DM, MXP2DM,MXQDM 
        write(6,*) ' LUDIA and LBLK ', LUDIA,LBLK
        write(6,*) ' INTSPC and IPRT ', INTSPC,IPRT
        write(6,*) ' NROOT ', NROOT
        write(6,*) ' IPROCC ', IPROCC
      END IF
*
* ====================
* 0 : Select subspace
* ====================
*
      IF( IWAY .EQ. 1) THEN
*.      Find number of combinations less or equal to MXDM
*       that does not separate degenerate pairs .
*. Used Scratch space : 6 * MXDM ( a bit extravagant )
        CALL MEMMAN(KL1,MAX(NAEL,3*MXDM+1),'ADDL  ',1,'KL1   ')
        CALL MEMMAN(KL2,MAX(NBEL,2*MXDM+1),'ADDL  ',2,'KL2   ')
        CALL MEMMAN(KL3,2*MXDM+1,'ADDL  ',2,'KL3   ')
        CALL FNDMND(LUDIA,LBLK,CIDIA,MXDM,NPRDET,WORK(KL1),
     &              WORK(KL2),IDET,WORK(KL3),NTEST )
      ELSE IF ( IWAY .EQ. 2 ) THEN
        CALL ISTVC2(IDET,0,1,NPRDET)
      END IF
*. Check for degenerencies on the boundaries between P1, P2 and Q space 
*. P1 - P2
      IF( MXP1DM .GT. 0 ) THEN
        IF(MXP1DM.GE.NPRDET ) THEN
          NP1DM = NPRDET
        ELSE
          IIDET = MXP1DM
 101      CONTINUE
          IF( ABS(WORK(KL3-1+IIDET+1)-WORK(KL3-1+IIDET))
     &       .LE. 0.000001) THEN
             IIDET = IIDET - 1
             GOTO 101
          END IF
          NP1DM = IIDET
        END IF
      ELSE
        NP1DM = 0
      END IF
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 ) 
     & WRITE(6,*) ' Actual dimension of P1 Space ', NP1DM
*. P2 - Q space
      IF(MXP2DM.GT.0) THEN
        IF(MXP1DM+MXP2DM.GE.NPRDET) THEN
          NP2DM = NPRDET - NP1DM
        ELSE 
          IIDET = MXP1DM + MXP2DM
 102      CONTINUE
          IF( ABS(WORK(KL3-1+IIDET+1)-WORK(KL3-1+IIDET))
     &       .LE. 0.0000001) THEN
             IIDET = IIDET - 1
             GOTO 102
          END IF
          NP2DM = IIDET - NP1DM
        END IF
      ELSE
        NP2DM = 0
      END IF
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 ) 
     & WRITE(6,*) ' Actual dimension of P2 Space ', NP2DM
*. Q space 
      IF(MXQDM.NE.0) THEN
        NQDM = NPRDET - NP1DM - NP2DM
      ELSE
        NQDM = 0
      END IF
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 ) 
     & WRITE(6,*) ' Actual dimension of Q Space ', NQDM
*
      NPDM = NP1DM + NP2DM
      NPRDET = NP1DM + NP2DM + NQDM
*. Copy over to MX numbers ( not nice but .. )
      MXP1DM = NP1DM
      MXP2DM = NP2DM
      MXQDM = NQDM
*
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 ) 
     & WRITE(6,*) 
     & ' Total number of combinations in subspace : ', NPRDET
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' IDET IN EXPHAM '
        CALL IWRTMA(IDET,1,NPRDET,1,NPRDET)
       
      END IF
*
* ============================================================= 
* Alpha and beta strings corresponding to selected combinations
* ============================================================= 
*
*
*. Convert determinant numbers to string numbers 
*
*. P dets
      CALL STRFDT(INTSPC,ISMOST,IOCOC,NPDM,IDET,IA,IB,ICOMBI)
*. Q dets
      IF(NQDM.NE.0) 
     &CALL STRFDT(INTSPC,ISMOST,IOCOC,NQDM,
     &            IDET(1+NPDM),IA(1+NPDM),IB(1+NPDM),ICOMBI)
*. and remaining to be printed
      IF(NPDM+NQDM.LT.IPROCC) THEN
        LREST = IPROCC - NPDM - NQDM
        CALL STRFDT(INTSPC,ISMOST,IOCOC,LREST,
     &       IDET(1+NPDM+NQDM),IA(1+NPDM+NQDM),
     &       IB(1+NPDM+NQDM),ICOMBI)
       END IF
*
      IF(IPROCC.NE.0) THEN
*. Print occupation of lowest IPROCC lowest SD 's
        WRITE(6,*)
        WRITE(6,'(A)')
     &  ' ========================================================='
        WRITE(6,'(A,I4,A)') 
     &  ' Occupation and energy of lowest ', IPROCC, ' combinations'
        WRITE(6,'(A)')
     &  ' ========================================================='
        WRITE(6,*)
        DO I = 1, IPROCC
*
          WRITE(6,'(A,I8,A,F18.10)')
     &    '  Energy  of combination ',IDET(I),' is ', 
     &       WORK(KL3-1+I)+ECORE
          WRITE(6,'(A)')
     &    '  Corresponding alpha - and beta string '
            WRITE(6,'(4X,10I4)')
     &      (IASTR(IEL,IA(I)),IEL = 1, NAEL )
            WRITE(6,'(4X,10I4)')
     &      (IBSTR(IEL,IB(I)),IEL = 1, NBEL )
C         ELSE
C           WRITE(6,'(4X,10(1X,A6))')
C    &      (IOBLAB(IASTR(IEL,IA)),IEL = 1, NAEL )
C           WRITE(6,'(4X,10(1X,A6))')
C    &      (IOBLAB(IBSTR(IEL,IB)),IEL = 1, NBEL )
C         END IF
        END DO
      END IF
      MXDM = MXDM1
      IF(MXDM.EQ.0) GOTO 9999
*
* ==================================
* Obtain Hamiltonian matrix elements 
* ==================================
*
*. Pointers
      KLPHP = 1
      KLPHQ = KLPHP + NPDM*(NPDM+1)/2
      KLQHQ = KLPHQ + NP1DM * NQDM
*. PHP Hamiltonian
      CALL LUCIAQENTER('DIHDJ')
*
* Scratch space for DIHDJ
      LSCR = 4 * NORB
      CALL MEMMAN(KSCR,LSCR,'ADDL  ',1,' H0SCR')
      CALL MEMMAN(KLIAST,NAEL,'ADDL  ',1,'LIA   ')
      CALL MEMMAN(KLIBST,NBEL,'ADDL  ',1,'LIB   ')
*
      ISTSGN = 0
      ECOREP = 0.0D0
      CALL DIHDJ(IA,IB,NPDM,IA,IB,
     &           NPDM,NAEL,NBEL,WORK(KSCR),LSCR,
     &           NORB,ONEBOD,HAMIL(KLPHP),
     &           1,NINOB,ECOREP,ICOMBI,PSIGN,IASTR,IBSTR,IASTR,
     &           IBSTR,ISTSGN,IASGN,IBSGN,IASGN,IBSGN,WORK(KLIAST),
     &           WORK(KLIBST),NDIF0,NDIF1,NDIF2,NTEST)
*. PHQ Hamiltonian
      IF(NQDM.NE.0)
     &CALL DIHDJ(IA,IB,NP1DM,IA(1+NPDM),IB(1+NPDM),NQDM,      
     &           NAEL,NBEL,WORK(KSCR),LSCR,
     &           NORB,ONEBOD,HAMIL(KLPHQ),
     &           0,NINOB,ECOREP,ICOMBI,PSIGN,IASTR,IBSTR,IASTR,
     &           IBSTR,ISTSGN,IASGN,IBSGN,IASGN,IBSGN,
     &           WORK(KLIAST),
     &           WORK(KLIBST),NDIF0,NDIF1,NDIF2,NTEST)
      CALL LUCIAQEXIT('DIHDJ')
*. QHQ Hamiltonian
      IF(LUDIA.LE.0) THEN
        DO 0607 IIDET = NPDM + 1, NPDM+NQDM
          HAMIL(KLQHQ + IIDET - NPDM-1) =
     &    CIDIA(IDET(IIDET))
 0607   CONTINUE
      ELSE IF (LUDIA.GT. 0 ) THEN
C            GATVCD(LU   ,LBLK,NGAT,IGAT,XGAT,SEGMNT,IPRT)
        IZERO = 0
        CALL GATVCD(LUDIA,LBLK,NQDM,IDET(NPDM+1),
     &              HAMIL(KLQHQ),CIDIA,IZERO)
      END IF
*
      IF(NTEST .GE. 20 ) THEN
        IF(NQDM .NE. 0 ) THEN
          WRITE(6,*) ' PHP, PHQ and QHQ parts of H0 '
        END IF
        WRITE(6,*) ' PHP '
        WRITE(6,*) ' === '
        CALL PRSYM(HAMIL(KLPHP),NPDM)     
        IF(NQDM.NE.0) THEN
          WRITE(6,*) ' PHQ '
          WRITE(6,*) ' === '
          CALL WRTMAT(HAMIL(KLPHQ),NP1DM,NQDM,NP1DM,NQDM)     
          WRITE(6,*) ' QHQ '
          WRITE(6,*) ' === '
          CALL WRTMAT(HAMIL(KLQHQ),1,NQDM,1,NQDM)     
        END IF
      END IF
*
* =================================
* Diagonalize to obtain lowest roots
* =================================
*
      IF(IDODGN .GT. 0 ) THEN
        CALL MEMMAN(KFREE,IDUMMY,'FREE  ',IDUMMY,'CDUMMY')
        CALL H0DIAG(HAMIL(KLPHP),HAMIL(KLPHQ),HAMIL(KLQHQ),
     &       NP1DM,NP2DM,NQDM,NROOT,EIGVAL,EIGVEC,WORK(KFREE),
     &       NTEST,ECORE )
        IF( NTEST .GE. 2 ) THEN
          WRITE(6,*) ' Eigenvalues of subspace Hamiltonian '
C         CALL WRTMAT(EIGVAL,1,NROOT,1,NROOT)
          WRITE(6,'(5F18.10)') (EIGVAL(IROOT),IROOT=1,NROOT)
        END IF
        IF( NTEST .GE. 50 ) THEN
          WRITE(6,*) ' Eigenvectors of subspace Hamiltonian '
          CALL WRTMAT(EIGVEC,NPRDET,NROOT,NPRDET,NROOT)
*
          WRITE(6,*) ' Subspace determinants '
          CALL IWRTMA(IDET,1,NPRDET,1,NPRDET)
        END IF
      END IF
*
* ===============================
* Analyze the Root approximations                   
* ===============================
*
C     IF( IPRT .GT. 0 ) THEN
C     DO 1869 IROOT = 1, NROOT
C       IOFF = (IROOT-1)*NPRDET + 1
C       CALL SETVEC(CIDIA,0.0D0,NDET)
C       CALL SCAVEC(CIDIA,EIGVEC(IOFF),IDET,NPRDET)
C       WRITE(6,*)
C       WRITE(6,'(A,I3)') '  Information about root ... ',IROOT
C       WRITE(6,'(A)')    '  ******************************'
C       WRITE(6,*)
C       WRITE(6,'(A,E15.8)') '   Energy .... ',EIGVAL(IROOT)
C       CUTOFF = 0.1D0
CTOBE   CALL ANACI(CIDIA,ICSYM,CUTOFF,100)
C1869 CONTINUE
C     END IF
*
      write(6,*) ' Final memory check in H0mat '
      call memchk
      write(6,*) ' Memory check passed '
*
 9999 CONTINUE
*
      CALL LUCIAQEXIT('EXPHAM')
      RETURN
      END
       SUBROUTINE HPARTV(PHP,P1HQ,QHQ,NP1,NP2,NQ,VECIN,VECUT,
     &                   E,SCR)
*
* Multiply partitioned precondition matrix
*
* PHP - P1HQ  (QHQ-E)**-1 QHP1 
*
* with vector VECIN to give vector VECUT 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION PHP(*),P1HQ(NP1,NQ),QHQ(NQ)
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECUT(*)
*. Scratch
      DIMENSION SCR(*)
*. Scr should at least be of length NQ 
*
      KLQVEC = 1
      KLFREE = KLQVEC + NQ
*.    QHP1 * VECIN
C          MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(P1HQ,VECIN,SCR(KLQVEC),NP1,NQ,1)
*. (QHQ-E)-1 * QHP1 * VECIN
      SHIFT = - E
      CALL DIAVC3(SCR(KLQVEC),SCR(KLQVEC),QHQ,SHIFT,NQ,NDSV)
*.- PHQ (QHQ-E)-1 * QHP1 * VECIN
      CALL MATVCC(P1HQ,SCR(KLQVEC),VECUT,NP1,NQ,0)
      CALL SCALVE(VECUT,-1.0D0,NP1)
      CALL SETVEC(VECUT(1+NP1),0.0D0,NP2)
* + PHP * VECIN
      NP = NP1 + NP2
      IJ = 0
      DO 100 I = 1, NP
        DO 90 J = 1, I
          IJ = IJ + 1
          VECUT(I) = VECUT(I)+ PHP(IJ)*VECIN(J)
          VECUT(J) = VECUT(J)+ PHP(IJ)*VECIN(I)
  90    CONTINUE
        VECUT(I) = VECUT(I)- PHP(IJ)*VECIN(I)
  100 CONTINUE
*
      NTEST = 10
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' HPARTV Speaking, Input and output vectors '
        CALL WRTMAT(VECIN,1,NP,1,NP)
        CALL WRTMAT(VECIN,1,NP,1,NP)
      END IF
*
      RETURN
      END
 
      SUBROUTINE HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT)
*
* Multiply P1P2Q preconditioner with a vector
* Jeppe Olsen , July 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General Input
      DIMENSION PHP(*),PHQ(NP1,NQ),QHQ(NQ)
*. PHP is in lower triangular form
*. Specific input
      DIMENSION VECIN(*)
*.Output
      DIMENSION VECUT(*)
*
      NP = NP1+NP2
      NPQ = NP + NQ
*
      CALL SETVEC(VECUT,0.0D0,NPQ)
*. PHQ * VECIN
      CALL MATVCC(PHQ,VECIN(1+NP),VECUT,NP1,NQ,0)
*. PHP * VECIN
      IJ = 0
      DO 60 I = 1, NP
      DO 50 J = 1, I
        IJ = IJ + 1
        VECUT(I) = VECUT(I)+PHP(IJ)*VECIN(J)
        VECUT(J) = VECUT(J)+PHP(IJ)*VECIN(I)
   50 CONTINUE
      VECUT(I) = VECUT(I)-PHP(IJ)*VECIN(I)
   60 CONTINUE
* QHP * VECIN
      CALL MATVCC(PHQ,VECIN,VECUT(1+NP),NP1,NQ,1)
*. QHQ * VECIN
      DO 100 I = 1, NQ
        VECUT(I+NP) = VECUT(I+NP)+QHQ(I)*VECIN(NP+I)
  100 CONTINUE
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from HPQTV '
        CALL WRTMAT(VECIN,1,NPQ,1,NPQ)
        CALL WRTMAT(VECUT,1,NPQ,1,NPQ)
      END IF
      RETURN
      END
      SUBROUTINE HPQTVM(VECIN,VECUT,IDUM,JDUM)
* 
* Outer routine for Preconditioner times vector 
*
*. LUCIA version
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,
C    &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR
C     COMMON/CRUN/MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV,
C    &            ICISTR,NOCSF,IDIAG
      
      

      
      COMMON/SHFT/SHIFT
*
      NPDM = MXP1+MXP2    
      NQDM = MXQ   
      KLPHP = KH0
      KLPHQ = KH0 + NPDM*(NPDM+1)/2
      KLQHQ = KLPHQ + MXP1*MXQ   
*
C          HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT)
      CALL HPQTV(MXP1,MXP2,MXQ,      
     &           WORK(KLPHP),WORK(KLPHQ),WORK(KLQHQ),
     &           VECIN,VECUT)         
*
      IF(SHIFT.NE.0.0D0) THEN
       CALL VECSUM(VECUT,VECUT,VECIN,1.0D0,SHIFT,NPDM+NQDM)
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from HPQTVM'
        CALL WRTMAT(VECIN,1,NPDM+NQDM,1,NPDM+NQDM)
        CALL WRTMAT(VECUT,1,NPDM+NQDM,1,NPDM+NQDM)
      END IF
*
      RETURN
      END 
      SUBROUTINE HPRITV(LUDIA,LUC,LUS,NVEC,NP1,NP2,NQ,H0,IPNTR,
     &                  Enot,VEC1,VEC2,SCR,LUIN)
*
* Obtain (Hpre-Enot)**(-1) VECIN  where vecin is vector stored in LUIN
* final vector is returned in VEC1
*
* Where Hpre is the P1P2Q initial matrix , modified to
* be exact in the subsplace spanned by the vectors .
*
* Jeppe Olsen, July 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. subspace preconditioner
      DIMENSION H0(*)
*.two scratch vectors
      DIMENSION VEC1(*),VEC2(*)
*.Local scratch space
      DIMENSION SCR(*)
*
      NTEST = 0
      KLFREE = 1
*
* ==================================================
* 1 : Obtain diverse matrices projected onto subspace
* ==================================================
*
*. Subspace hessian
      KLSBH = KLFREE
      KLFREE = KLFREE + NVEC ** 2
*. Subspace approximate Hessian
      KLSBH0 = KLFREE
      KLFREE = KLFREE + NVEC ** 2
* (H0 - E)**(-1) between 2 C vectors
      KLCCHI = KLFREE
      KLFREE = KLFREE + NVEC ** 2
      KLSH0I = KLFREE
* (H0 - E)**(-1) between  C vectors  and sigma vectors
      KLCSHI = KLFREE
      KLFREE = KLFREE + NVEC ** 2
* (H0 - E)**(-1) between   and sigma vectors
      KLSSHI = KLFREE
      KLFREE = KLFREE + NVEC ** 2
*
*Exact Hamiltonian in subspace
*
      CALL REWINO(LUC)
      DO 100 I = 1, NVEC
        CALL FRMDSC(VEC1,NVAR,-1,LUC,IMZERO,IAMPACK)
        CALL REWINO(LUS)
        DO 50 J = 1, I
          JI = (I-1)*NVEC+J
          IJ = (J-1)*NVEC+I
          CALL FRMDSC(VEC2,NVAR,-1,LUS,IMZERO,IAMPACK)
          SCR(KLSBH-1+IJ) = INPROD(VEC1,VEC2,NVAR)
          SCR(KLSBH-1+JI) = SCR(KLSBH-1+IJ)
   50   CONTINUE
  100 CONTINUE
*
*.Approximate Hamiltonian in subspace
*
      CALL REWINO(LUC)
      DO 200 I = 1, NVEC
        CALL FRMDSC(VEC1,NVAR,-1,LUC,IMZERO,IAMPACK)
        CALL REWINO(LUDIA)
        CALL FRMDSC(VEC2,NVAR,-1,LUDIA,IMZERO,IAMPACK)
C            H0TV(VECIN,VECUT,DIAG,NVAR,NPQDM,IPNTR,H0,
C    &                WORK,NP1,NP2,NQ)
        CALL H0TV(VEC1,VEC2,VEC2,NVAR,NPQDM,IPNTR,H0,
     &            SCR (KLFREE),NP1,NP2,NQ)
        CALL REWINO(LUC)
        DO 150 J = 1, I
          CALL FRMDSC(VEC2,NVAR,-1,LUC,IMZERO,IAMPACK)
          JI = (I-1)*NVEC+J
          IJ = (J-1)*NVEC+I
          SCR(KLSBH0-1+IJ) = INPROD(VEC1,VEC2,NVAR)
          SCR(KLSBH0-1+JI) = SCR(KLSBH0-1+IJ)
  150   CONTINUE
  200 CONTINUE
*
*.Inverted Approximate Hamiltonian in subspace  and between sigma vecs
*
      CALL REWINO(LUC)
      DO 300 I = 1, NVEC
        CALL FRMDSC(VEC1,NVAR,-1,LUC,IMZERO,IAMPACK)
        CALL REWINO(LUDIA)
        CALL FRMDSC(VEC2,NVAR,-1,LUDIA,IMZERO,IAMPACK)
C            H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
C    &              H0,SHIFT,WORK,XH0PSX,
C    &              NP1,NP2,NQ)
        CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPQDM,IPNTR,
     &       H0,-ENOT,SCR(KLFREE),XH0PSX,
     &       NP1,NP2,NQ)
*. In subspace
        CALL REWINO(LUC)
        DO 250 J = 1, I
          CALL FRMDSC(VEC1,NVAR,-1,LUC,IMZERO,IAMPACK)
          JI = (I-1)*NVEC+J
          IJ = (J-1)*NVEC+I
          SCR(KLCCHI-1+IJ) = INPROD(VEC1,VEC2,NVAR)
          SCR(KLCCHI-1+JI) = SCR(KLCCHI-1+IJ)
  250   CONTINUE
*. between subspace and sigma vectors
        CALL REWINO(LUS)
        DO 280 J = 1, NVEC
          CALL FRMDSC(VEC1,NVAR,-1,LUS,IMZERO,IAMPACK)
          IJ = (J-1)*NVEC+I
          SCR(KLCSHI-1+IJ) = INPROD(VEC1,VEC2,NVAR)
  280   CONTINUE
  300 CONTINUE
*
*.Inverted Approximate Hamiltonian in subspace  of sigma vectors
*
      CALL REWINO(LUS)
      DO 400 I = 1, NVEC
        CALL FRMDSC(VEC1,NVAR,-1,LUS,IMZERO,IAMPACK)
        CALL REWINO(LUDIA)
        CALL FRMDSC(VEC2,NVAR,-1,LUDIA,IMZERO,IAMPACK)
C            H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
C    &              H0,SHIFT,WORK,XH0PSX,
C    &              NP1,NP2,NQ)
        CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPQDM,IPNTR,
     &       H0,-ENOT,SCR(KLFREE),XH0PSX,
     &       NP1,NP2,NQ)
        CALL REWINO(LUS)
        DO 350 J = 1, I
          CALL FRMDSC(VEC1,NVAR,-1,LUS,IMZERO,IAMPACK)
          JI = (I-1)*NVEC+J
          IJ = (J-1)*NVEC+I
          SCR(KLSSHI-1+IJ) = INPROD(VEC1,VEC2,NVAR)
          SCR(KLSSHI-1+JI) = SCR(KLSSHI-1+IJ)
  350   CONTINUE
  400 CONTINUE
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' Subspace matrices '
        WRITE(6,*) ' ================= '
        WRITE(6,*)
        WRITE(6,*) ' Exact Hamiltonian in subspace '
        CALL WRTMAT(SCR(KLSBH),NVEC,NVEC,NVEC,NVEC)
        WRITE(6,*)
        WRITE(6,*) ' Approximate Hamiltonian in subspace '
        CALL WRTMAT(SCR(KLSBH0),NVEC,NVEC,NVEC,NVEC)
        WRITE(6,*)
        WRITE(6,*) ' H0-Enot inverted between C vectors  '
        CALL WRTMAT(SCR(KLCCHI),NVEC,NVEC,NVEC,NVEC)
        WRITE(6,*)
        WRITE(6,*) ' H0-Enot between C and S vectors  '
        CALL WRTMAT(SCR(KLCSHI),NVEC,NVEC,NVEC,NVEC)
        WRITE(6,*)
        WRITE(6,*) ' H0-Enot inverted between S vectors  '
        CALL WRTMAT(SCR(KLCCHI),NVEC,NVEC,NVEC,NVEC)
      END IF
*
* =================================================
* 2 : Obtain low rank matrix modifying (Hnot-E)**-1
* =================================================
*
*                  ( 0    1   )
* Obtain matrix P =(          )
*                  ( 1   H-H0 )
 
      KLP = KLFREE
      KLFREE = KLFREE + 4 * NVEC ** 2
      CALL SETVEC(SCR(KLP),0.0D0,4*NVEC**2)
      DO 500 I = NVEC+1,2*NVEC
        J = I - NVEC
        SCR(KLP-1+(I-1)*2*NVEC+J) = 1.0D0
        SCR(KLP-1+(J-1)*2*NVEC+I) = 1.0D0
        DO 450 J = NVEC+1,2*NVEC
          IJN = (I-1)*NVEC+J
          SCR(KLP-1+(I-1)*2*NVEC+J)=
     &    SCR(KLSBH -1 +  IJN) - SCR(KLSBH0-1+IJN)
  450   CONTINUE
  500 CONTINUE
* Invert P
      KLPI = KLFREE
      KLFREE = KLFREE + 4 * NVEC ** 2
      CALL INVMAT(SCR(KLPI),SCR(KLFREE),2*NVEC,2*NVEC,ISING)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Inverted P matrix  '
        CALL WRTMAT(SCR(KLPI),2*NVEC,2*NVEC,2*NVEC,2*NVEC)
      END IF
*
*            ( DelSig(T)(H0-Enot)-1 Delsig  Delsig(T)(H0-Enot)-1 X  )
* Obtain Q = (                                                      )
*            ( X(T)(H0-Enot)-1 Delsig            X(T)(H0-Enot)-1 X  )
*
* where Delsig = sigma - sigma0
      KLQ = KLFREE
      KLFREE = KLFREE + 4 * NVEC ** 2
*
* ============================
* DelSig(T)(H0-Enot)-1 Delsig can be rewritten to
* ============================
*
*          Sigma(T) (H0-Enot)-1 Sigma
*-Enot *   Sigma(T) (H0-Enot)-1 C
*-Enot *   C(T)     (H0-Enot)-1 Sigma
*+Enot**2* C(T)     (H0-Enot)-1 C
* - 2Hsub -Enot * 1 + H0sub
      CALL TRPMT3(SCR(KLCSHI),NVEC,NVEC,SCR(KLFREE) )
      CALL VECSUM(SCR(KLFREE),SCR(KLFREE),SCR(KLCSHI),-ENOT,-ENOT,
     &            NVEC ** 2 )
      CALL VECSUM(SCR(KLFREE),SCR(KLFREE),SCR(KLSSHI),1.0D0,ENOT**2,
     &            NVEC ** 2 )
      CALL VECSUM(SCR(KLFREE),SCR(KLFREE),SCR(KLSBH),1.0D0,-2.0D0,
     &            NVEC ** 2 )
      CALL VECSUM(SCR(KLFREE),SCR(KLFREE),SCR(KLSBH),1.0D0,1.0D0,
     &            NVEC ** 2 )
      DO 700 I = 1, NVEC
        SCR(KLFREE-1+(I-1)*NVEC+I) = SCR(KLFREE-1+(I-1)*NVEC+I)-ENOT
  700 CONTINUE
* Transfer to q
       DO 705 I = 1, NVEC
       DO 705 J = 1, NVEC
         SCR(KLQ-1+(I-1)*2*NVEC+J) = SCR(KLFREE-1+(I-1)*NVEC+J)
  705  CONTINUE
 
*
* ============================
* 2 : DelSig(T)(H0-Enot)-1 X   can be rewritten to
* ============================
*
*          Sigma(T) (H0-Enot)-1 C
*-Enot *       C(T) (H0-Enot)-1 C
*         -       1
      CALL TRPMT3(SCR(KLCSHI),NVEC,NVEC,SCR(KLFREE))
      CALL VECSUM(SCR(KLFREE),SCR(KLFREE),SCR(KLCCHI),1.0D0,-ENOT,
     &            NVEC ** 2 )
      DO 800 I = 1, NVEC
        SCR(KLFREE-1+(I-1)*NVEC+I) = SCR(KLFREE-1+(I-1)*NVEC+I)-1.0D0
  800 CONTINUE
* Transfer to q
       DO 805 I = 1, NVEC
       DO 805 J = 1, NVEC
         SCR(KLQ-1+(I+NVEC-1)*2*NVEC+J) = SCR(KLFREE-1+(I-1)*NVEC+J)
         SCR(KLQ-1+(I-1)*2*NVEC+J+NVEC) = SCR(KLFREE-1+(J-1)*NVEC+I)
  805  CONTINUE
*
* ======================
* 3 : X(T)(H0-Enot)-1 X
* ======================
      DO 900 I = 1, NVEC
      DO 900 J = 1, NVEC
        SCR(KLQ-1+(I-1)*2*NVEC+J+NVEC) = SCR(KLCCHI-1+(I-1)*NVEC+J)
  900 CONTINUE
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' The Q Matrix '
        WRITE(6,*) ' ============='
        CALL WRTMAT(SCR (KLQ),2*NVEC,2*NVEC,2*NVEC,2*NVEC)
      END IF
* Obtain (P-1 + Q) ** - 1
      CALL VECSUM(SCR (KLPI),SCR (KLPI),SCR(KLQ),1.0D0,1.0D0,4*NVEC**2)
      CALL INVMAT(SCR (KLPI),SCR (KLQ),2*NVEC,2*NVEC,ISING)
*
* ===========================================
*                 ( Delsig(T)(H0-Enot)-1 V )
* 3 : Obtain V1 = (                        )
*                 ( X(T)     (H0-Enot)-1 V )
* ============================================
*
      CALL REWINO(LUIN)
      CALL FRMDSC(VEC1,NVAR,-1,LUIN,IMZERO,IAMPACK)
      CALL REWINO(LUC)
      KLV1 = KLFREE
      KLFREE = KLFREE + 2 * NVEC
      DO 1000 I = 1, NVEC
        CALL FRMDSC(VEC2,NVAR,-1,LUC,IMZERO,IAMPACK)
        SCR(KLV1-1+I) = - INPROD(VEC1,VEC2,NVAR)
 1000 CONTINUE
      CALL REWINO(LUDIA)
      CALL FRMDSC(VEC2,NVAR,-1,LUDIA,IMZERO,IAMPACK)
C            H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
C    &              H0,SHIFT,WORK,XH0PSX,
C    &              NP1,NP2,NQ)
      CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDM,IPNTR,
     &            H0,-ENOT,SCR (KLFREE),XH0PSX,NP1,NP2,NQ)
      CALL REWINO(LUC)
      CALL REWINO(LUS)
      DO 1100 I = 1, NVEC
        CALL FRMDSC(VEC1,NVAR,-1,LUC,IMZERO,IAMPACK)
        X = INPROD(VEC1,VEC2,NVAR)
        SCR(KLV1-1+I) = SCR(KLV1-1+I)- ENOT * X
        SCR(KLV1-1+ NVAR+ I) =  X
        CALL FRMDSC(VEC1,NVAR,-1,LUS,IMZERO,IAMPACK)
        X = INPROD(VEC1,VEC2,NVAR)
        SCR(KLV1-1+I) = SCR(KLV1-1+I) + X
 1100 CONTINUE
*
* Obtain vector (H-enot)**-1
* (C- sum(i)Vec1(i)delsig(i)+sum(i)Vec1(Nvec+i)*X(I))
* which happens to be what this subroutine is all about )
*
      CALL REWINO(LUIN)
      CALL FRMDSC(VEC1,NVAR,-1,LUIN,IMZERO,IAMPACK)
      CALL REWINO(LUS)
      CALL REWINO(LUC)
      DO 1300 I = 1, NVEC
        CALL FRMDSC(VEC2,NVAR,-1,LUS,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-SCR(KLV1-1+I),NVAR)
        CALL FRMDSC(VEC2,NVAR,-1,LUC,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-SCR(KLV1-1+NVEC+I),NVAR)
 1300 CONTINUE
      CALL REWINO(LUDIA)
      CALL FRMDSC(VEC2,NVAR,-1,LUDIA,IMZERO,IAMPACK)
C            H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
C    &              H0,SHIFT,WORK,XH0PSX,
C    &              NP1,NP2,NQ)
      CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDM,IPNTR,
     &            H0,-ENOT,SCR (KLFREE),XH0PSX,NP1,NP2,NQ)
      CALL REWINO(LUC)
      DO 1400 I = 1, NVEC
        CALL FRMDSC(VEC1,NVAR,-1,LUC,IMZERO,IAMPACK)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,SCR(KLV1-1+I),NVAR)
 1400 CONTINUE
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Sihifted inverse preconditioner times vector '
        CALL WRTMAT(VEC2,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
      SUBROUTINE MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
*
* ITRNS = 0 : VOUT(I) = A(I,J)*VIN(J)
* ITRNS = 1 : VOUT(I) = A(J,I)*VIN(J)
*
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(NROW,NCOL)
       DIMENSION VIN(*),VOUT(*)
*
      IF(ITRNS.EQ.0) THEN
        DO 10 I = 1, NROW
   10   VOUT(I) = 0.0D0
        DO 100 J = 1, NCOL
         VINJ = VIN(J)
         DO 60 I = 1, NROW
           VOUT(I) = VOUT(I) + A(I,J)*VINJ
   60    CONTINUE
  100   CONTINUE
      ELSE IF( ITRNS.EQ.1) THEN
        DO 200 I = 1, NCOL
          X = 0.0D0
          DO 160 J = 1, NROW
            X = X + A(J,I)*VIN(J)
  160     CONTINUE
          VOUT(I) = X
  200   CONTINUE
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(ITRNS.EQ.0) THEN
          WRITE(6,*) ' Vectorout = matrix * vectorin (MATVCC) '
          WRITE(6,*) ' Input and output vectors '
          CALL WRTMAT(VIN,1,NCOL,1,NCOL)
          CALL WRTMAT(VOUT,1,NROW,1,NROW)
          WRITE(6,*) ' Matrix '
          CALL WRTMAT(A,NROW,NCOL,NROW,NCOL)
        ELSE 
          WRITE(6,*) ' Vectorout = matrix(T) * vectorin (MATVCC) '
          WRITE(6,*) ' Input and output vectors '
          CALL WRTMAT(VIN,1,NROW,1,NROW)
          CALL WRTMAT(VOUT,1,NCOL,1,NCOL)
          WRITE(6,*) ' Matrix (untransposed)'
          CALL WRTMAT(A,NROW,NCOL,NROW,NCOL)
        END IF
      END IF

*
      RETURN
      END
      SUBROUTINE PRMBLK(IDC,ISGV,IASM,IBSM,IATP,IBTP,PS,PL,
     &                  JATP,JBTP,JASM,JBSM,ISGN,ITRP,NPERM)
*
* A block of CI coefficients defined by by IATP,IASM,IBTP,IBSM is given
*
* Obtain the number of other blocks that can be obtained by spin 
* and relection symmetry.
*
* Jeppe Olsen, July 1993
*
* =====
* Output
* =====
* JATP(I),JASM(I),JBTP(I),JBSM(I) Indeces for Block I
* NPERM : Number of blocks  that can be obtained
* ITRP(I) = 1 => block should     be transposed
*         = 0 => block should not be transposed
* ISGN   : Sign to multiply previous block with to getnew sign
*
*
* There are four types of permutations
*
*    operation   *      JASM  *      JBSM  * JATP * JBTP * Iperm * Sign *
*   *********************************************************************
*   * Identity   *      IASM  *      IBSM  * IATP * IBTP *   0   * 1    *
*   * Ml         * ISGV(IASM) * ISGV(IBSM) * IATP * IBTP *   0   * PL   *
*   * Ms         *      IBSM  *      IASM  * IBTP * IATP *   1   * PS   *
*   * Ms+Ml      * ISGV(IBSM) * ISGV(IASM) * IBTP * IATP *   1   * PS PL*
*   *********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*.Input
      DIMENSION ISGV(*)
*.Output
      DIMENSION JATP(4),JBTP(4),JASM(4),JBSM(4),ISGN(4),ITRP(4)
*
*. To eliminate some compiler warnings 
      KASM = 0
      KBSM = 0
      KATP = 0
      KBTP = 0
      KSIGN = 0
      KTRP = 0
      LSIGN = 0
      LTRP = 0
*
      NPERM = 0
      DO 100 IPERM = 1, 4
        ISET = 0
        IF(IPERM.EQ.1) THEN
*
* Identity operation
*
          KASM = IASM
          KBSM = IBSM
          KATP = IATP
          KBTP = IBTP
          KSIGN = 1
          KTRP = 0
          ISET = 1
        ELSE IF(IPERM.EQ.2.AND.(IDC.EQ.3.OR.IDC.EQ.4)) THEN
*
* Ml reflection
*
          KASM = ISGV(IASM)
          KBSM = ISGV(IBSM)
          KATP = IATP
          KBTP = IBTP
          IF(PL.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PL .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 0
          ISET = 1
        ELSE IF(IPERM.EQ.3.AND.(IDC.EQ.2.OR.IDC.EQ.4)) THEN
*
* Ms reflection
*
          KASM = IBSM
          KBSM = IASM
          KATP = IBTP
          KBTP = IATP
          IF(PS.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PS .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 1
          ISET = 1
        ELSE IF(IPERM.EQ.4 .AND. IDC.EQ.4) THEN
*
* Ms Ml  reflection
*
          KASM = ISGV(IBSM)
          KBSM = ISGV(IASM)
          KATP = IBTP
          KBTP = IATP
          IF(PS*PL.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PS .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 1
          ISET = 1
        END IF
*
        IF(ISET.EQ.1) THEN
*. A new permutation was found, check and see if it was obtained previously
          INEW = 1
          DO 50 LPERM = 1, NPERM 
            IF(JATP(LPERM).EQ.KATP  .AND. JASM(LPERM).EQ.KASM .AND. 
     &         JBTP(LPERM).EQ.KBTP  .AND. JBSM(LPERM).EQ.KBSM) INEW = 0
   50     CONTINUE
          IF(INEW.EQ.1) THEN
*. The permutation was new, add it to the list
            NPERM = NPERM + 1
            JASM(NPERM) = KASM
            JBSM(NPERM) = KBSM
            JATP(NPERM) = KATP
            JBTP(NPERM) = KBTP
            IF(NPERM.EQ.1. OR. (NPERM.GE.1.AND.KSIGN.EQ.LSIGN))THEN
              ISGN(NPERM) = 1
            ELSE 
              ISGN(NPERM) = -1
            END IF
            LSIGN = KSIGN
            IF(NPERM.EQ.1. OR. (NPERM.GE.1.AND.KTRP.EQ.LTRP))THEN
              ITRP(NPERM) = 0
            ELSE 
              ITRP(NPERM) = 1
            END IF
            LTRP = KTRP
          END IF
        END IF
  100 CONTINUE
*
*. Should the block be trnasposed or scaled to return to initial form
      ITRP(NPERM+1) = LTRP
      ISGN(NPERM+1) = LSIGN
      IFNSGN = LSIGN
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)') ' Blocks obtained from IASM IBSM IATP IBTP ',
     &  IASM,IBSM,IATP,IBTP
        WRITE(6,*)
        WRITE(6,'(A)') ' JASM JBSM JATP JBTP Isgn Itrp  '
        WRITE(6,*)
        DO 10 IPERM = 1, NPERM
          WRITE(6,'(2x,6I4)') JASM(IPERM),JBSM(IPERM),JATP(IPERM),
     &                        JBTP(IPERM),ISGN(IPERM),ITRP(IPERM)
   10   CONTINUE
      END IF
*
      RETURN
      END 
      SUBROUTINE SCLDIA(A,FACTOR,NDIM,IPACK)
*
* scale diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*                assumed packed columnwise !!!!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*)
*
      IF( IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          II = (I-1)*NDIM + I
          A(II) = A(II) * FACTOR
  100   CONTINUE
      ELSE
        II = 1
        DO 200 I = 1, NDIM
          A(II) = A(II) * FACTOR
          II = II + NDIM - I + 1
  200   CONTINUE
      END IF
*
      RETURN
      END 
      SUBROUTINE SCTOGT(ISCAT,IGAT,NSCAT,NGAT,IPRT)

*
* A scattering array ISCAT is given, obtain the corresponding
* gather array IGAT. Dimensions of arrays can be different,
* a zero in IGAT indicates that there is no element defined 
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*. Input
      DIMENSION ISCAT(NSCAT)
*. output 
      DIMENSION IGAT(NGAT)
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRT)
      CALL ISETVC(IGAT,0,NGAT)
      DO 100 I = 1, NSCAT
       IGAT(ISCAT(I)) = I
  100 CONTINUE
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)  ' SCTOGT : ISCAT(input) and IGAT(output)'
        CALL IWRTMA(ISCAT,1,NSCAT,1,NSCAT)
        WRITE(6,*)
        CALL IWRTMA(IGAT ,1,NGAT ,1,NGAT )
      END IF
*
      RETURN
      END
      SUBROUTINE SDCMRF(CSD,CCM,ICOFF,IWAY,IATP,IBTP,IASM,IBSM,NA,NB,
     &                  IDC,PS,PL,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
*
* Change a block of coefficients bwtween combination format and 
* Slater determinant format
*
*     IWAY = 1 : SD => Combinations
*     IWAY = 2 : Combinations => SD
*
* Input 
* =====
* CSD : Block in determinant form
* CCM : Block in combination  form
* IWAY : as above
* IATP,IBTP : type of alpha- and beta- string
* NA,NB : Number of alpha- and beta- strings
* IDC  : Combination type
* PS   : Spin combination sign
* PL   : Ml   combination sign
* ISGVST : Ml reflection of strings
*
*
* If ISCALE .EQ. 0, no overall scaling is performed,
*                   the overall scale factor is returned
*                   as SCLFAC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
CNW   DIMENSION CSD(*),CCM(*),ISGVST(*)
      DIMENSION CSD(*),ISGVST(*)
      integer ccm
*
      CALL LUCIAQENTER('SDCMR')
      NTEST = 000
      IF(NTEST.GE.100)  THEN
       WRITE(6,*) ' SDCMRF: NA, NB =', NA,NB
       WRITE(6,*) ' PS, PL = ', PS, PL
       WRITE(6,*) ' IDC ', IDC
       WRITE(6,*) ' IWAY ', IWAY      
      END IF
*
      SQRT2  = SQRT(2.0D0)
      SQRT2I = 1.0D0/SQRT2
*
*. Is combination array packed ?
*
      SCLFAC = 1.0D0
      IPACK = 0
      FACTOR = 1.0D0
*
      IF(IDC.EQ.2.OR.IDC.EQ.4) THEN
         SIGN = PS
         FACTOR = SQRT2 
         IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) IPACK = 1
      ELSE IF( IDC.EQ.4.AND.IASM.EQ.ISGVST(IBSM)) THEN
        IF(IATP.EQ.IBTP) IPACK = 1
        SIGN = PS*PL
        FACTOR = 2.0D0
      END IF
*
      LDET = NA * NB
      IF( IPACK.EQ.0) THEN
        LCOMB = LDET
      ELSE
        LCOMB = NA*(NA+1)/2
      END IF
      IF(IDC.EQ.4.AND.IPACK.EQ.0) FACTOR = SQRT2
      IF(IWAY.EQ.2) FACTOR = 1.0D0/FACTOR
*
*. SD => combination transformation
*
      IF(IWAY .EQ. 1 ) THEN
        IF(IPACK.EQ.1) THEN                    
*. Pack to triangular form
          CALL TRIPK3G(CSD,CCM,ICOFF,1,NA,NA,SIGN)
C              TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
        ELSE
          call ga_put(ccm,icoff,icoff+na*nb-1,1,1,csd,1)
c         CALL COPVEC(CSD,CCM,NA*NB)
        END IF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          IF(ISCALE.EQ.1) THEN
            SCLFAC = 1.0D0
            CALL SCALVE(CCM,FACTOR,LCOMB)
          ELSE
            SCLFAC = FACTOR
          END IF
          IF(IPACK.EQ.1 ) THEN
            CALL SCLDIA(CCM,SQRT2I,NA,1) 
          END IF
        END IF
      END IF
*
*. Combination => SD transformation
*
      IF(IWAY.EQ.2) THEN
        IF(IPACK.EQ.1) THEN
*. Unpack from triangular form
           CALL TRIPK3G(CSD,CCM,ICOFF,IWAY,NA,NA,SIGN)
        ELSE
           call ga_get(ccm,icoff,icoff+na*nb-1,1,1,csd,1)
c          CALL COPVEC(CCM,CSD,NA*NB)
        ENDIF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          IF(ISCALE.EQ.1) THEN
            SCLFAC = 1.0D0
            CALL SCALVE(CSD,FACTOR,LDET)
          ELSE
            SCLFAC = FACTOR
          END IF
          IF(IPACK.EQ.1) THEN
             CALL SCLDIA(CSD,SQRT2,NA,0)
          END IF
        END IF
      END IF
*
      IF(NTEST.NE.0) THEN
C     IF(NTEST.NE.0.AND.IWAY.EQ.1) THEN
        WRITE(6,*) ' Information from SDCMRF '
   
        WRITE(6,'(A,6I4)') ' IWAY IATP IBTP IASM IBSM IDC ',
     &                   IWAY,IATP,IBTP,IASM,IBSM,IDC
        WRITE(6,'(A,I4,3X,2E13.6)') ' IPACK FACTOR SIGN',
     &  IPACK,FACTOR,SIGN
        IF(NTEST.GE. 100000 ) THEN
          WRITE(6,*) ' Slater determinant block '
          CALL WRTMAT(CSD,NA,NB,NA,NB)
          WRITE(6,*)
          WRITE(6,*) ' Combination block '
          IF(IPACK.EQ.1) THEN
CNW         CALL PRSM2(CCM,NA)
            call ga_print(ccm)
          ELSE
CNW         CALL WRTMAT(CCM,NA,NB,NA,NB)
            call ga_print(ccm)
          END IF
        END IF
      END IF
      CALL LUCIAQEXIT('SDCMR')
*
      RETURN
      END
      SUBROUTINE STFDT3(NDET,IDET,IASTR2,IBSTR2,NSMST,NOCTPA,NOCTPB,
     &                  NSSOA,NSSOB,ISSOA,ISSOB,IOCOC,
     &                  ISMOST,ICOMBI)
*
* Obtain actual numbers for alpha- and beta strings corresponding to
* Given determinant numbers
*
* Jeppe Olsen, January 1989 
*
*. Modified March 1993 so only string numbers are returned
*
* Changed into LUCIA form , June 1993
* Combinations enabled, september 1993
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*. Specific Input
      DIMENSION IDET(NDET)
*. General input
      DIMENSION NSSOA(NOCTPA,NSMST),NSSOB(NOCTPB,NSMST)
      DIMENSION ISSOA(NOCTPA,NSMST),ISSOB(NOCTPB,NSMST)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ISMOST(*)
*. Output
      DIMENSION IASTR2(NDET),IBSTR2(NDET)
*
COLD  IF( ICOMBI.NE.0) THEN
COLD    WRITE(6,*) ' STFDT3 does not work for COMBINATIONS'  
COLD    WRITE(6,*) ' Enforced STOP'
COLD    STOP 'STFDT3'
COLD  END IF
*
      IDETF = 1
      DO 1000 IASM = 1, NSMST
        IBSM = ISMOST(IASM)                      
        IF(IBSM.LE.0) GOTO 1000
        IF(ICOMBI.NE.0. AND. IASM.LT.IBSM) GOTO 1000
        DO 999 IATP = 1,NOCTPA
          NIA = NSSOA(IATP,IASM)
          IAOFF = ISSOA(IATP,IASM)
          IF(ICOMBI.NE.0.AND.IASM.EQ.IBSM) THEN
            ISYM1 = 1
            MXBTP = IATP
          ELSE
            ISYM1 = 0
            MXBTP = NOCTPB
          END IF
          DO 901 IBTP = 1, MXBTP 
            IF(IATP.EQ.IBTP.AND.ISYM1.EQ.1) THEN
              ISYM = 1
            ELSE
              ISYM = 0
            END IF
*
            IF(IOCOC(IATP,IBTP).LE.0) GOTO 901
            NIB = NSSOB(IBTP,IBSM)
            IF(NIA*NIB.EQ.0) GOTO 901
            IBOFF = ISSOB(IBTP,IBSM)
            IF(ISYM.EQ.0) THEN
              IDETL = IDETF + NIA*NIB-1
            ELSE
              IDETL = IDETF + NIA*(NIA+1)/2 - 1
            END IF
            DO 500 I = 1, NDET
              IF(IDET(I).GE.IDETF.AND.IDET(I).LE.IDETL) THEN
                IREL = IDET(I)-IDETF+1
                IF(ISYM.EQ.0) THEN
*. Find IAREL, IBREL so (IBREL-1)*NIA + IAREL = IREL
                  IBREL = (IREL-1)/NIA+1
                  IAREL = IREL - NIA*(IBREL-1)
                ELSE IF (ISYM.EQ. 1 ) THEN
*. Find IAREL, IBREL so IREL =
*                   (IBREL-1)*NIA + IAREL  - IBREL*(IBREL-1)/2
                   CALL UNPCPC(IREL,IAREL,NIA,IBREL)
                END IF
                IASTR2(I) = IAREL + IAOFF -1
                IBSTR2(I) = IBREL + IBOFF -1
              END IF
  500       CONTINUE
            IDETF = IDETL+1
  901     CONTINUE
  999   CONTINUE
 1000 CONTINUE
 1001 CONTINUE
C
      NTEST = 00
      IF( NTEST .GE. 2 ) THEN
         IF( ICOMBI .NE. 0 ) THEN
           WRITE(6,*) ' Combinations selected '
         ELSE
           WRITE(6,*) ' Determinants selected '
         END IF
C
         CALL IWRTMA(IDET,1,NDET,1,NDET)
         WRITE(6,*) ' Selected alpha and beta strings from STFDT3 '
         CALL IWRTMA(IASTR2,1,NDET,1,NDET)
         WRITE(6,*)
         CALL IWRTMA(IBSTR2,1,NDET,1,NDET)
      END IF
C
      RETURN
      END
      SUBROUTINE STRFDT(INTSPC,JSMOST,IOCOC,NSBDET,IDET,IASTR2,IBSTR2,
     &            ICOMBI)
*
* A CI expansion from internal subspace INTSPC and symmetry defined
* by JSMOST is considered. In this CI expansion, a set of NSBDET combinations
* are given in terms of combination numbers IDET. Obtain the number of 
* the corresponding alpha and betastrings
*
*
*. Revised March 22 1993 so only stringnumbers are stored
*  September 1993 : combinations enabled
*
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      DIMENSION IDET(*)
*. General Input
      DIMENSION JSMOST(*),IOCOC(*)
*. Common blocks
      INCLUDE 'strbas.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
*. 
*. Output
      DIMENSION IASTR2(*),IBSTR2(*)
*
      IATP = IASTFI(INTSPC)
      IBTP = IBSTFI(INTSPC)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
       CALL STFDT3(NSBDET,IDET,IASTR2,IBSTR2,NSMST,NOCTPA,NOCTPB,
     &             int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &             int_mb(KISTSO(IATP)),int_mb(KISTSO(IBTP)),
     &             IOCOC,JSMOST,ICOMBI)
*
      RETURN
      END
      SUBROUTINE XTRCDI(AMAT,DIAG,NDIM,ISYM)
C
C EXTRACT DIAGONAL OF A MATRIX
C
C IF ISYM .LE. 0 MATRIX IS ASSUMED STORED IN COMPLETE FORM
C IF ISYM .GT. 0 MATRIX IS ASSUMED PACKED ROWWISE IN
C                LOWER TRIANGULAR FORM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AMAT(1),DIAG(1)
C
      DO 100 I = 1, NDIM
        IF ( ISYM .EQ. 0 ) THEN
          II = (I-1)*NDIM + I
        ELSE
          II = I*(I+1)/2
        END IF
        DIAG(I) = AMAT(II)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE SCDTC2(RASVEC,ISMOST,ICBLTP,NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,IOCOC,IDC,IWAY,IMMLST,IPRNT)
*
* Scale elements of a RAS vector to transfer between
* combinations and packed determinants
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
* Combination storage mode is defined BY IDC
*
* General symmetry version , Feb 1991
*. GAS form of NSASAm NSBSO, Aug. 95
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RASVEC(*),NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ISMOST(*),ICBLTP(*),IMMLST(*)
*
COLD  LOGICAL DIAGBL
*
      NTEST = 10
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from SCDTC2 '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
        CALL WRTRS2(RASVEC,ISMOST,ICBLTP,IOCOC,
     &              NOCTPA,NOCTPB,NSASO,NSBSO,NSMST)
      END IF
*
      SQ2 = SQRT(2.0D0)
      SQ2I = 1.0D0/SQ2
*
      IBASE = 1
      DO 200 IASM = 1, NSMST
*
        IBSM = ISMOST(IASM)
        IF(IBSM.EQ.0.OR.ICBLTP(IASM).EQ.0) GOTO 200
        DO  100 IATP = 1, NOCTPA
          IF(ICBLTP(IASM).EQ.2) THEN
            IBTPMX = IATP
          ELSE
            IBTPMX = NOCTPB
          END IF
          NIA   = NSASO(IASM,IATP)
          DO 50 IBTP = 1,IBTPMX
            IF(IOCOC(IATP,IBTP).EQ.0) GOTO   50
*. Number of elements in this block
            NIB = NSBSO(IBSM,IBTP)
            IF(ICBLTP(IASM).EQ.2.AND.IATP.EQ.IBTP) THEN
                NELMNT =  NIA*(NIA+1)/2
            ELSE
                NELMNT =  NIA*NIB
            END IF
*Ms combinations
          IF(IDC.EQ.2) THEN
            IF(IWAY.EQ.1) THEN
              FACTOR = SQ2
            ELSE
              FACTOR = SQ2I
            END IF
            CALL SCALVE(RASVEC(IBASE),FACTOR,NELMNT)
            IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
              FACTOR = 1.0D0/FACTOR
              CALL SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
            END IF
*Ml combinations
          ELSE IF(IDC.EQ.3.AND.IMMLST(IASM).NE.IASM) THEN
            IF(IWAY.EQ.1) THEN
              FACTOR = SQ2
            ELSE
              FACTOR = SQ2I
            END IF
            CALL SCALVE(RASVEC(IBASE),FACTOR,NELMNT)
*Ml Ms combinations
          ELSE IF(IDC.EQ.4) THEN
            IF(IWAY.EQ.1) THEN
              IF(IASM.EQ.IBSM) THEN
                FACTOR = SQ2
              ELSE
                FACTOR = 2.0D0
              END IF
            ELSE IF(IWAY.EQ.2) THEN
              IF(IASM.EQ.IBSM) THEN
                FACTOR = SQ2I
              ELSE
                FACTOR = 0.5D0
              END IF
            END IF
            CALL SCALVE(RASVEC(IBASE),FACTOR,NELMNT)
            IF(IATP.EQ.IBTP) THEN
              IF(IWAY.EQ.1) THEN
                FACTOR = SQ2I
              ELSE IF(IWAY.EQ.2) THEN
                FACTOR = SQ2
              END IF
              CALL SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
            END IF
          END IF
*
          IBASE = IBASE + NELMNT
  50      CONTINUE
 100    CONTINUE
 200  CONTINUE
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Scaled vector '
        CALL WRTRS2(RASVEC,ISMOST,ICBLTP,IOCOC,
     &              NOCTPA,NOCTPB,NSASO,NSBSO,NSMST)
      END IF
*
      RETURN
      END
      FUNCTION GETH1I(IORB,JORB)
*
* Obtain one -electron integral H(IORB,JOB)
*
* Interface from EXPHAM to LUCIA
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*
      ISM = ISMFTO(IORB)
      ITP = ITPFSO(IREOTS(IORB))
      IREL = IORB - IOBPTS(ITP,ISM) + 1
*
      JSM = ISMFTO(JORB)
      JTP = ITPFSO(IREOTS(JORB))
      JREL = JORB - IOBPTS(JTP,JSM) + 1
*
      GETH1I = GETH1E(IREL,ITP,ISM,JREL,JTP,JSM)
*
      NTEST = 00
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' GETH1I : IORB JORB ', IORB, JORB
        WRITE(6,*) ' ISM ITP IREL ', ISM,ITP,IREL
        WRITE(6,*) ' JSM JTP JREL ', JSM,JTP,JREL 
        WRITE(6,*) ' GETH1I = ', GETH1I
      END IF
*
      RETURN
      END 
      FUNCTION GETH1I_2(IORB,JORB)
*
* Obtain one -electron integral H(IORB,JOB)
*
* Interface from EXPHAM to LUCIA
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*
      ISM = ISMFTO(IORB)
      ITP = ITPFSO(IREOTS(IORB))
      IREL = IORB - IOBPTS(ITP,ISM) + 1
*
      JSM = ISMFTO(JORB)
      JTP = ITPFSO(IREOTS(JORB))
      JREL = JORB - IOBPTS(JTP,JSM) + 1
*
      IJSM = MULTD2H(ISM,JSM)
*
      IF (IJSM.EQ.1) THEN
        GETH1I_2 = GETH1E(IREL,ITP,ISM,JREL,JTP,JSM)
      ELSE
        GETH1I_2 = 0D0
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' GETH1I : IORB JORB ', IORB, JORB
        WRITE(6,*) ' ISM ITP IREL ', ISM,ITP,IREL
        WRITE(6,*) ' JSM JTP JREL ', JSM,JTP,JREL 
        WRITE(6,*) ' GETH1I = ', GETH1I
      END IF
*
      RETURN
      END 
c $Id$
