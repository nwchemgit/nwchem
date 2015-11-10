*
* Codes for Hessian and orbital optimization
* 
      SUBROUTINE LIN_RESP(V1,IV1SM,V2,IV2SM)
*
* Calculate static linear response due to one-body perturbations
* V1, V2
*
* Sloppy version written to test HF Hessian/HF response densities 
* - only programmed/tested for total symmetric operators.
*
* Jeppe Olsen, May 1999
*
* proceeds as follows :
* 1 : Calculate Response due to pert V1 :
*     a : Fock matrix for V1       
*     b : Response from this Fock matrix as a density matrix
* 2 : Multiply with integrals V2
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input : Integrals, not triangular packed
      DIMENSION V1(*),V2(*)
*
      REAL * 8 INPROD
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LIN_RE')
*. V1 integrals triangular packed
      CALL MEMMAN(KLV1P,NINT1,'ADDL  ',2,'V1PACK')
C  TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
      CALL TRIPAK_BLKM(V1,dbl_mb(KLV1P),1,NTOOBS,NSMOB)
*. Fock matrix corresponding to V1
      LFOCK = 0
      DO ISM = 1, NSMOB
        LFOCK = LFOCK + NTOOBS(ISM)*NTOOBS(ISM)
      END DO
      CALL MEMMAN(KLFV1,LFOCK,'ADDL  ',2,'V1FOCK')
      CALL SWAPVE(dbl_mb(KLV1P),KINT1,NINT1)
      CALL FOCK_MAT(dbl_mb(KLFV1),1)
      CALL SWAPVE(dbl_mb(KLV1P),KINT1,NINT1)
*. Response density corresponding to this Fock matrix
      CALL MEMMAN(KLR1RS,LFOCK,'ADDL  ',2,'R1RS  ')
      CALL RESPDEN_FROM_F(dbl_mb(KLFV1),dbl_mb(KLR1RS))
*. And polarizability as expectation value
      ALPHA = INPROD(V2,dbl_mb(KLR1RS),LFOCK)
      WRITE(6,*) ' Polarizability ', ALPHA
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'LIN_RE')
*
      RETURN
      END 
      SUBROUTINE RESPDEN_FROM_F(FOCK,RESPDEN)
*
* Obtain Response contribution to density for method based upon Hartree-Fock 
* orbitals using a general Fock matrix. 
*
* Output response contribution to density is delivered in RESPDEN 
* as symmetrypacked block diagonal matrix. Both upper and lower halfs 
* are included.
*     Jeppe Olsen, Spring of 99
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input 
      DIMENSION FOCK(*)
*. Output
      DIMENSION RESPDEN(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'RESPD_')
* 
*. Non-redundant orbital excitations
*
*. Nonredundant type-type excitations
      CALL MEMMAN(KLTTACT,NGAS*NGAS,'ADDL  ',1,'TTACT ')
      CALL NONRED_TT_EXC(int_mb(KLTTACT),1,0)
*. Nonredundant orbital excitations
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
      CALL MEMMAN(KLOOEXCC,2*NTOOB*NTOOB,'ADDL  ',1,'OOEXCC')
      CALL NONRED_OO_EXC(NOOEXC,int_mb(KLOOEXC),int_mb(KLOOEXCC),
     &                   1,int_mb(KLTTACT),2)
*. Find inactive/active rotations (e.g. frozen core approximations)
*. and sort them to the end of OOEXCC array; return number of these
*. rotations on NCANCON
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CANCON')
      CALL MEMMAN(KLICC,2*NOOEXC,'ADDL  ',1,'ICANCN')
      CALL GET_NCANCON(NOOEXC,int_mb(KLOOEXCC),int_mb(KLOOEXC),
     &                 NCANCON,int_mb(KLICC))
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'CANCON')

*. Gradient from Fock matrix 
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1    ')
      CALL F_TO_E1(FOCK,dbl_mb(KLE1),1,int_mb(KLOOEXCC),NOOEXC,0)
*. If additional canonical conditions of frozen orbital have to be
*. applied: calculate here there contribution to the response density
*. (divide by orbital denominator) and modify the gradient accordingly
      CALL MEMMAN(KLKAP,NOOEXC,'ADDL  ',2,'KAPPA ')
      IF (NCANCON.GT.0)
     &     CALL FROCON(WORK(KLKAP),WORK(KLE1),WORK(KFIZ),WORK(KLOOEXCC),
     &     NOOEXC,NCANCON)
*. Construct orbital Hessian
C          E2_FUSK(ORBHES,NOOEXC,IOOEXCC)  
      NOOEXC2 = NOOEXC-NCANCON
      IF (NOOEXC2.GT.0) THEN
        CALL MEMMAN(KLE2,NOOEXC2*(NOOEXC2+1)/2,'ADDL  ',2,'E2    ')
        CALL E2_FUSK(WORK(KLE2),NOOEXC2,WORK(KLOOEXCC),0) 
*. Find inverted Hessian 
        CALL MEMMAN(KLE2_EXP,NOOEXC2*NOOEXC2,'ADDL  ',2,'E2_EXP') 
        CALL MEMMAN(KLSCR   ,NOOEXC2*NOOEXC2,'ADDL  ',2,'E2_SCR') 
        CALL TRIPAK(WORK(KLE2_EXP),WORK(KLE2),2,NOOEXC2,NOOEXC2)
        CALL INVMAT(WORK(KLE2_EXP),WORK(KLSCR),NOOEXC2,NOOEXC2,ISING)
*. Kappa = -E[2]**(-1) V[1]
        CALL MATVCB(WORK(KLE2_EXP),WORK(KLE1),WORK(KLKAP),
     &            NOOEXC2,NOOEXC2,0)
        ONEM = -1.0D0
        CALL SCALVE(WORK(KLKAP),ONEM,NOOEXC2)
      END IF
*. Response density is now [Kappa(T),Den1] = [Den1,Kappa]
*
*. For Hartree-Fock we simply have that for each nonvanisking excitation 
* below diagonal we obtain an element +2*kappa
      CALL MEMMAN(KLKAP1,NTOOB**2,'ADDL  ',2,'KAP1  ')
*. Expand kappa to complete lower half form
C     REFRM_KAPPA(XKAP_CMP,XKAP_FULL,IWAY,IVSM,
C    &                       IOOEXC,IOOEXCC,NOOEXC)
C?    WRITE(6,*) ' WORK(KLKAP) before REFRM '
C?    CALL WRTMAT(WORK(KLKAP),1,NOOEXC,1,NOOEXC)
      CALL REFRM_KAPPA(WORK(KLKAP),WORK(KLKAP1),1,1,
     &                 WORK(KLOOEXC),WORK(KLOOEXCC),NOOEXC)
*. Expand to complete matrix
C         TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
      CALL TRIPAK_BLKM(RESPDEN,WORK(KLKAP1),2,NTOOBS,NSMOB)
*. And multiply with 2
      TWO = 2.0D0
      LENGTH = 0
      DO ISM = 1, NSMOB
        LENGTH = LENGTH + NTOOBS(ISM) ** 2
      END DO
      CALL SCALVE(RESPDEN,TWO,LENGTH)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*)
        WRITE(6,*) ' Orbital relaxation term to density '
        WRITE(6,*) ' ================================== '
        WRITE(6,*)
        CALL APRBLM2(RESPDEN,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'RESPD_')
      RETURN
      END
      SUBROUTINE GET_NCANCON(NOOEXC,IOOEXCC,IOOEXC,NCANCON,IBUFF)
*
*     In frozen core calculations, multipliers appear, that do not
*     correspond to Brillioun but to canonical conditions; as they
*     need separate care-taking, we sort them to the end of the
*     IOOEXCC-array
*
*     by Andreas, Jan. 2005
*

      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'

      INTEGER, PARAMETER ::
     &     NTEST = 00
      
      INTEGER, INTENT(IN) ::
     &     NOOEXC
      INTEGER, INTENT(INOUT) ::
     &     IOOEXCC(2,NOOEXC), IBUFF(2,NOOEXC), IOOEXC(NTOOB,NTOOB)
      INTEGER, INTENT(OUT) ::
     &     NCANCON

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'IOOEXCC array on entry:'
        DO IEX = 1, NOOEXC
          WRITE(6,*) '  ',IEX,' -> ',IOOEXCC(1,IEX),':',IOOEXCC(2,IEX)
        END DO
      END IF

      NCANCON = 0
      DO IEX = 1, NOOEXC
        IORB = IOOEXCC(1,IEX)
        JORB = IOOEXCC(2,IEX)
        ITP = IHPVGAS(ITPFTO(IORB))
        JTP = IHPVGAS(ITPFTO(JORB))
        IF (ITP.EQ.3.OR.JTP.EQ.3) THEN
          WRITE(6,*) 'Still somewhat indecisive about valence spaces'
          STOP 'valence spaces?'
        END IF
        IF (ITP.EQ.JTP) THEN
          NCANCON = NCANCON+1
          IBUFF(1,NCANCON) = IOOEXCC(1,IEX)
          IBUFF(2,NCANCON) = IOOEXCC(2,IEX)
        END IF

      END DO

      IF (NTEST.GE.100) WRITE(6,*) NCANCON,
     &     ' canonical conditions identified'

      IF (NCANCON.GT.0) THEN
        ICNT = 0
        DO IEX = 1, NOOEXC
          IORB = IOOEXCC(1,IEX)
          JORB = IOOEXCC(2,IEX)
          ITP = IHPVGAS(ITPFTO(IORB))
          JTP = IHPVGAS(ITPFTO(JORB))
          IF (ITP.NE.JTP) THEN
            ICNT = ICNT + 1
            IBUFF(1,NCANCON+ICNT) = IOOEXCC(1,IEX)
            IBUFF(2,NCANCON+ICNT) = IOOEXCC(2,IEX)
          END IF
          
        END DO
        DO IEX = 1, ICNT
          IORB = IBUFF(1,NCANCON+IEX)
          JORB = IBUFF(2,NCANCON+IEX)
          IOOEXCC(1,IEX) = IORB
          IOOEXCC(2,IEX) = JORB
          IOOEXC(IORB,JORB) = IEX
          IOOEXC(JORB,IORB) = -IEX
        END DO
        DO IEX = 1,NCANCON
          IORB = IBUFF(1,IEX)
          JORB = IBUFF(2,IEX)
          IOOEXCC(1,ICNT+IEX) = IORB
          IOOEXCC(2,ICNT+IEX) = JORB
          IOOEXC(IORB,JORB) = ICNT+IEX
          IOOEXC(JORB,IORB) = -ICNT-IEX
        END DO

        IF (NTEST.GE.100) THEN
          WRITE(6,*) 'IOOEXCC array on exit:'
          DO IEX = 1, NOOEXC
            WRITE(6,*) '  ',IEX,' -> ',IOOEXCC(1,IEX),':',IOOEXCC(2,IEX)
          END DO
        END IF

      END IF

      RETURN
      END
      SUBROUTINE FROCON(KAP,E1,FOCK,IOOEXCC,NOOEXC,NCANCON)
*
*     treat frozen core/deleted orbital contributions to orbital response
*     a canonical fock-matrix is expected
*
*     the orbital response contributions are then
*
*       kap(i,J) = -E1(i,J)/(eps(i)-eps(J))
*
*     and the gradient is modified as
*
*       E1(i,a) = E(i,a) + kap(j,J)A(i,a,j,J)
*
*     Andreas, Feb. 2005
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'

      INTEGER, PARAMETER ::
     &     NTEST = 0

      INTEGER, INTENT(IN) ::
     &     NOOEXC, NCANCON, IOOEXCC(2,NOOEXC)
      REAL(8), INTENT(IN) ::
     &     FOCK(*)
      REAL(8), INTENT(INOUT) ::
     &     E1(NOOEXC), KAP(NOOEXC)

      INTEGER ::
     &     ISYOFF(8)

      IF (NTEST.GT.0) THEN
        WRITE(6,*) ' entered FROCON'
        WRITE(6,*) ' =============='
      END IF

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'Fock matrix provided'
        CALL APRBLM2(FOCK,NTOOBS,NTOOBS,NSMOB,1)
      END IF

      IDX = 0
      DO ISYM = 1, NSMOB
        ISYOFF(ISYM) = IDX
        IDX = IDX + NTOOBS(ISYM)*(NTOOBS(ISYM)+1)/2
      END DO

      ICANCON = NOOEXC-NCANCON+1

      DO IEX = ICANCON, NOOEXC
        IORB = IOOEXCC(1,IEX)
        JORB = IOOEXCC(2,IEX)
        ISM  = ISMFTO(IORB)
        JSM  = ISMFTO(JORB)
        II   = IREOTS(IORB)-IBSO(ISM)+1
        JJ   = IREOTS(JORB)-IBSO(JSM)+1
        IF (ISM.NE.JSM) THEN
          WRITE(6,*) 'Unexpected event in FROCON: ISM!=JSM'
          WRITE(6,*) ' ISM = ',ISM,'   JSM = ',JSM
          STOP 'frocon: unexpected event'
        END IF
        DENOM = FOCK(ISYOFF(ISM) + II*(II+1)/2) -
     &          FOCK(ISYOFF(JSM) + JJ*(JJ+1)/2)
        KAP(IEX) = -E1(IEX)/(4d0*DENOM)

        DO IEXBJ = 1, ICANCON-1
          IB = IOOEXCC(1,IEXBJ)
          IJ = IOOEXCC(2,IEXBJ)

          E1(IEXBJ) = E1(IEXBJ)
     &              + 4d0*KAP(IEX)*(4D0*GTIJKL(IB,IJ,IORB,JORB)
     &                             -GTIJKL(IB,JORB,IORB,IJ)
     &                             -GTIJKL(IB,IORB,JORB,IJ))

        END DO
      END DO

      IF(NTEST.GE.100) THEN  
        WRITE(6,*) ' Modified orbital gradient '  
        WRITE(6,*) ' ========================= '
        CALL WRT_EXCVEC(E1,IOOEXCC,NOOEXC-NCANCON)
        WRITE(6,*) ' Frozen core part of kappa '  
        WRITE(6,*) ' ========================= '
        CALL WRT_EXCVEC(KAP(ICANCON),IOOEXCC(1,ICANCON),NCANCON)

      END IF

      RETURN
      END

      SUBROUTINE RESPDEN_FROM_DEN(DEN1,DEN2,RESPDEN)
*
* Obtain Response contribution to density for method based upon Hartree-Fock 
* orbitals. Densities of correlated wavefunction is DEN1,DEN2
*
* Output response contribution to density is delivered in RESPDEN 
* as symmetrypacked block diagonal matrix. Both upper and lower halfs 
* are included.
*     Jeppe Olsen, Spring of 99
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
*. Input 
      DIMENSION DEN1(*),DEN2(*)
*. Output
      DIMENSION RESPDEN(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'RESPD_')
* 
*. Non-redundant orbital excitations
*
*. Nonredundant type-type excitations
      CALL MEMMAN(KLTTACT,NGAS*NGAS,'ADDL  ',1,'TTACT ')
      CALL NONRED_TT_EXC(WORK(KLTTACT),1,0)
*. Nonredundant orbital excitations
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
      CALL MEMMAN(KLOOEXCC,2*NTOOB*NTOOB,'ADDL  ',1,'OOEXCC')
      CALL NONRED_OO_EXC(NOOEXC,WORK(KLOOEXC),WORK(KLOOEXCC),
     &                   1,WORK(KLTTACT),2)
*. Orbital gradient with densities  DEN1, DEN2
*. Save current one- and two-electron densities
      LRHO1 = NTOOB**2
      CALL MEMMAN(KLRHO1_SAVE,LRHO1,'ADDL  ',2,'RHO1_SA')
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
      CALL MEMMAN(KLRHO2_SAVE,LRHO2,'ADDL  ',2,'RHO2_SA')
      CALL COPVEC(WORK(KRHO1),WORK(KLRHO1_SAVE),LRHO1)
      CALL COPVEC(WORK(KRHO2),WORK(KLRHO2_SAVE),LRHO2)
      CALL COPVEC(DEN1,WORK(KRHO1),LRHO1)
      CALL COPVEC(DEN2,WORK(KRHO2),LRHO2)
      CALL MEMMAN(KLFMAT,NTOOB**2,'ADDL  ',2,'FMAT  ') 
*Gradient using DEN1 and DEN2
      CALL FOCK_MAT(WORK(KLFMAT),2) 
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1    ')
C          F_TO_E1(F,E1,IEXCSM,IOOEXCC,NOOEXC,IMODE)
      CALL F_TO_E1(WORK(KLFMAT),WORK(KLE1),1,WORK(KLOOEXCC),NOOEXC,0)
* Restore order - atleast densities
      CALL COPVEC(WORK(KLRHO1_SAVE),WORK(KRHO1),LRHO1)
      CALL COPVEC(WORK(KLRHO2_SAVE),WORK(KRHO2),LRHO2)
*. Construct orbital Hessian
C          E2_FUSK(ORBHES,NOOEXC,IOOEXCC)  
      CALL MEMMAN(KLE2,NOOEXC*(NOOEXC+1)/2,'ADDL  ',2,'E2    ')
      CALL E2_FUSK(WORK(KLE2),NOOEXC,WORK(KLOOEXCC),0) 
*. Find inverted Hessian 
      CALL MEMMAN(KLE2_EXP,NOOEXC*NOOEXC,'ADDL  ',2,'E2_EXP') 
      CALL MEMMAN(KLSCR   ,NOOEXC*NOOEXC,'ADDL  ',2,'E2_SCR') 
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(WORK(KLE2_EXP),WORK(KLE2),2,NOOEXC,NOOEXC)
C          INVMAT(A,B,MATDIM,NDIM,ISING)
      CALL INVMAT(WORK(KLE2_EXP),WORK(KLSCR),NOOEXC,NOOEXC,ISING)
*. Kappa = -E[2]**(-1) V[1]
      CALL MEMMAN(KLKAP,NOOEXC,'ADDL  ',2,'KAPPA ')
C          MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(WORK(KLE2_EXP),WORK(KLE1),WORK(KLKAP),
     &            NOOEXC,NOOEXC,0)
*. Response density is now [Kappa(T),Den1] = [Den1,Kappa]
*
*. For Hartree-Fock we simply have that for each nonvanisking excitation 
* below diagonal we obtain an element +2*kappa
      CALL MEMMAN(KLKAP1,LENGTH,'ADDL  ',2,'KAP1  ')
*. Expand kappa to complete lower half form
C     REFRM_KAPPA(XKAP_CMP,XKAP_FULL,IWAY,IVSM,
C    &                       IOOEXC,IOOEXCC,NOOEXC)
      WRITE(6,*) ' WORK(KLKAP) before REFRM '
      CALL WRTMAT(WORK(KLKAP),1,NOOEXC,1,NOOEXC)
      CALL REFRM_KAPPA(WORK(KLKAP),WORK(KLKAP1),1,1,
     &                 WORK(KLOOEXC),WORK(KLOOEXCC),NOOEXC)
*. Expand to complete matrix
C         TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
      CALL TRIPAK_BLKM(RESPDEN,WORK(KLKAP1),2,NTOOBS,NSMOB)
*. And multiply with 2
      TWO = 2.0D0
      LENGTH = 0
      DO ISM = 1, NSMOB
        LENGTH = LENGTH + NTOOBS(ISM) ** 2
      END DO
      CALL SCALVE(RESPDEN,TWO,LENGTH)
*
* For general, not tested 
*. Extract symmetry blocks from complete one-electron density 
C     CALL MEMMAN(KLRHO1S,NTOOB**2,'ADDL  ','RHO1  ')
C     I RHO1SM = 1
C     CALL REORHO1(WORK(KRHO1),WORK(KLRHO1S),IRHO1SM)
*. Expand kappa to complete matrix
C     LENGTH = 0
C     DO ISM = 1, NSMOB
C       LENGTH = LENGTH + NTOOBS(ISM) ** 2
C     END DO
C     CALL MEMMAN(KLBLM1,LENGTH,'ADDL  ',2,'KLM1  ')
C     CALL MEMMAN(KLBLM2,LENGTH,'ADDL  ',2,'KLM2  ')
C     CALL MULT_BLOC_MAT(WORK(KLBLM1),WORK(
C      MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*)
        WRITE(6,*) ' Orbital relaxation term to density '
        WRITE(6,*) ' ================================== '
        WRITE(6,*)
        CALL APRBLM2(RESPDEN,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'RESPD_')
      RETURN
      END
      SUBROUTINE LIN_RESP_TO_V(V,IVSM,IVOP,NVOP,XKAPPA,IV12) 
*
* Obtain linear response due to external perturbation V 
* IF IV12 program assumes that V is a two-electron perturbation
* and that the two-electron term is in place
*
* Jeppe Olsen, Spring of 1999
*
*. General input
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
*. Specific input
      INTEGER IVOP(2,*)
      DIMENSION V(*)
*. Output
      DIMENSION XKAPPA(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LIN_RE')
*
*. Set up right hand side = -<0![E,V]!0> 
*
      LENGTH = NTOOB**2
      CALL MEMMAN(KLFMAT,LENGTH,'ADDL  ',2,'LFMAT ')    
      CALL MEMMAN(KLINT_SAVE,LENGTH,'ADDL  ',2,'INT_SA')    
      CALL MEMMAN(KLE1      ,LENGTH,'ADDL  ',2,'E1    ')    
*
      LENGTH = 0
      DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          IF(ISM.EQ.JSM) THEN
            LENGTH = LENGTH + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
          ELSE IF(ISM.GT.JSM) THEN
            LENGTH = LENGTH + NTOOBS(ISM)*NTOOBS(JSM)
          END IF
        END DO
      END DO
      CALL COPVEC(WORK(KINT1),WORK(KLINT_SAVE),NINT1)
      CALL COPVEC(V,WORK(KINT1),LENGTH) 
      IF(IVSM.EQ.1) THEN
        CALL FOCK_MAT(WORK(KLFMAT),IV12) 
      ELSE
        CALL FOCK_MAT2(WORK(KLFMAT),IV12,IVSM)
      END IF
      CALL F_TO_E1(WORK(KLFMAT),WORK(KLE1),IVSM,IVOP,NVOP,0)
*. Clean up
      CALL COPVEC(WORK(KLINT_SAVE),WORK(KINT1),NINT1)
*. Obtain orbital Hessian
      LENGTH = NVOP*(NVOP+1)/2        
      CALL MEMMAN(KLE2,LENGTH,'ADDL  ',2,'E2    ')
C          E2_FUSK(ORBHES,NOOEXC,IOOEXCC)  
      CALL E2_FUSK(WORK(KLE2),NVOP,IVOP,0) 
*. Find inverted Hessian 
      CALL MEMMAN(KLE2_EXP,NVOP*NVOP    ,'ADDL  ',2,'E2_EXP') 
      CALL MEMMAN(KLSCR   ,NVOP*NVOP    ,'ADDL  ',2,'E2_SCR') 
*
C       TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(WORK(KLE2_EXP),WORK(KLE2),2,NVOP,NVOP)
C      INVMAT(A,B,MATDIM,NDIM,ISING)
*      HER ER NOGET ROD ...
C     CALL TRIPAK(WORK(KLE2_EXP),WORK(KLSCR),NVOP,NVOP)
*. -E[2]*V[1]
C       MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(WORK(KLE2_EXP),WORK(KLE1),XKAPPA,NVOP,NVOP,0)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' Linear Response vector '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
C     WRT_EXCVEC(VEC,IOOEXCC,NOOEXC)
        CALL WRT_EXCVEC(XKAPPA,IVOP,NVOP)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ', IDUM,'LIN_RE')
*
      RETURN
      END
      SUBROUTINE F_TO_E1(F,E1,IEXCSM,IOOEXCC,NOOEXC,IMODE)
*
* Obtain gradient from F-matrix E1(I,J) = 2(F(I,J)-F(J,I))
*
* Jeppe Olsen, Spring of 1999
*
* version with direct addressing and IMODE: Andreas, summer/autumn 2004
*
*  IMODE:  0 --- transfer all gradient contributions
*  IMODE:  1 --- transfer only p/h rotation gradients of active shells
*  IMODE: -1 --- transfer only elements complementary to IMODE==1
*
*. General input
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      DIMENSION F(*), IOOEXCC(2,*)
*. Output
      DIMENSION E1(*)
*. Scratch
      DIMENSION ISOFF(NSMOB)
      LOGICAL L_AHAP
*

      NTEST = 00

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'Fock matrix:'
        CALL APRBLM2(F,NTOOBS,NTOOBS,NSMOB,0)
      END IF

*. Symmetry offset
      IIDX = 0
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IEXCSM)
        ISOFF(ISM) = IIDX
        IIDX = IIDX + NTOOBS(ISM)*NTOOBS(JSM)
      END DO
*     
      DO IEXC = 1, NOOEXC
        IORB = IOOEXCC(1,IEXC) 
        JORB = IOOEXCC(2,IEXC)
        ISM = ISMFTO(IORB)
        JSM = MULTD2H(ISM,IEXCSM)
        II = IREOTS(IORB) - IBSO(ISM) + 1
        JJ = IREOTS(JORB) - IBSO(JSM) + 1
        IJDX = ISOFF(ISM) + (JJ-1)*NTOOBS(ISM) + II
        JIDX = ISOFF(JSM) + (II-1)*NTOOBS(JSM) + JJ
        IF (IMODE.NE.0) THEN
          ITP = ITPFTO(IORB)
          JTP = ITPFTO(JORB)
          ! is it a particle/hole rotation between active shells?
          L_AHAP = I_IADX(ITP).EQ.2.AND.I_IADX(JTP).EQ.2.AND.
     &             IHPVGAS(ITP).NE.IHPVGAS(JTP)

          IF (.NOT.L_AHAP.AND.IMODE.EQ.1) CYCLE
          IF (L_AHAP.AND.IMODE.EQ.-1) CYCLE
        END IF

        E1(IEXC) = 2.0D0*(F(IJDX)-F(JIDX))
      END DO

      IF(NTEST.GE.10) THEN  
        WRITE(6,*) ' Orbital gradient '  
        WRITE(6,*) ' ================ '
        CALL WRT_EXCVEC(E1,IOOEXCC,NOOEXC)
      END IF
*
      RETURN
      END 
      SUBROUTINE F_TO_E1_OLD(F,E1,IEXCSM,IOOEXCC,NOOEXC)
*
* Obtain gradient from F-matrix E1(I,J) = 2(F(I,J)-F(J,I))
*
* Jeppe Olsen, Spring of 1999
*
*. General input
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cintfo.inc'
*. Specific input
      DIMENSION F(*), IOOEXCC(2,*)
*. Output
      DIMENSION E1(*)
*
C     WRITE(6,*) ' WORK(KINT1) ARRAY 1 : '
C     CALL WRTMAT(WORK(KINT1),1,NINT1, 1, NINT1)
*
*. Fetch routines for one-electron integrals picks integrals 
*. in WORK(KINT1).  Copy one-electron integrals to this array
*. and tell routines that 1-electron integrals do not have any
* permutational symmetry.
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'F_TO_ ')
*. Number of integrals with symmetry IEXCSM
      LENGTH = 0
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IEXCSM)
        LENGTH = LENGTH + NTOOBS(ISM)*NTOOBS(JSM)
      END DO
*.
      CALL MEMMAN(KLINT_SAVE,NINT1,'ADDL  ',2,'INT_SA') 
      CALL COPVEC(WORK(KINT1),WORK(KLINT_SAVE),NINT1)
      CALL COPVEC(F,WORK(KINT1),LENGTH)
      IH1FORM = 2
*
      DO IEXC = 1, NOOEXC
        IORB = IOOEXCC(1,IEXC) 
        JORB = IOOEXCC(2,IEXC)
        E1(IEXC) = 2.0D0*(GETH1_B(IORB,JORB)-GETH1_B(JORB,IORB))
      END DO
*. Cleanup 
      IH1FORM = 1
      CALL COPVEC(WORK(KLINT_SAVE),WORK(KINT1),NINT1)
*
C     WRITE(6,*) ' WORK(KINT1) ARRAY 2 : '
C     CALL WRTMAT(WORK(KINT1),1,NINT1, 1, NINT1)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN  
        WRITE(6,*) ' Orbital gradient '  
        WRITE(6,*) ' ================ '
        CALL WRT_EXCVEC(E1,IOOEXCC,NOOEXC)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'F_TO_ ')
*
      RETURN
      END 
      SUBROUTINE WRT_EXCVEC(VEC,IOOEXCC,NOOEXC)
*
* Write vector of 1-electron excitations
*
* Jeppe Olsen, Spring of 99
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOOEXCC(2,*)
      DIMENSION VEC(*)
*
      WRITE(6,*) ' Iorb Jorb             value'
      WRITE(6,*) '----------------------------------------------'
      DO IEXC = 1, NOOEXC
       WRITE(6,'(2I5,3X,E22.8)') 
     & IOOEXCC(1,IEXC),IOOEXCC(2,IEXC),VEC(IEXC)
      END DO
*
      RETURN
      END
      
      SUBROUTINE REFRM_KAPPA(XKAP_CMP,XKAP_FULL,IWAY,IVSM,
     &                       IOOEXC,IOOEXCC,NOOEXC)
*
* Reform Kappa between compact form, being a vector of nonred parameters
* Kappa_{ij},  and the full form XKAP_FULL being a complete lower 
* half form.
*
* It is assumed that XKAP_CMP is on TS order, whereas      
* XKAP_FULL is in ST form
*
* IWAY = 1 : CMP =>  FULL
*      = 2 : FULL = CMP
*
* Jeppe Olsen, Winter of 99
*. General input
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
      INTEGER IOOEXC(NTOOB,NTOOB), IOOEXCC(2,*)
*. Input and output
      DIMENSION XKAP_CMP(*),XKAP_FULL(*)
*
      IF(IVSM.NE.1) THEN
        WRITE(6,*) ' STOP : REFRM works only for IVSM = 1'
                     STOP': REFRM works only for IVSM = 1'
      END IF
* 
C?    WRITE(6,*) ' First element of XKAP_CMP', XKAP_CMP(1)
*
      ZERO = 0.0D0
      IJ_FULL = 1
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IVSM)
        IF(ISM.GE.JSM) THEN
          NI = NTOOBS(ISM)
          IOFF = IBSO(ISM)
          NJ = NTOOBS(JSM)
          JOFF = IBSO(JSM)
          IF(ISM.EQ.JSM) THEN
*.Rectangular block, rowwise packed
            NIJ = NI*(NI+1)/2  
          ELSE 
            NIJ = NI*NJ
          END IF
          IF(IWAY.EQ.1) CALL SETVEC(XKAP_FULL(IJ_FULL),ZERO,NIJ)
*
          DO I = 1, NI
          IF(ISM.EQ.JSM) THEN
            JMAX = I
          ELSE
            JMAX = NJ
          END IF
          DO J = 1, JMAX 
*. Absolut numbers in type order
            IABS = IREOST(IOFF-1+I) 
            JABS = IREOST(JOFF-1+J)
            IF(IOOEXC(IABS,JABS).NE.0) THEN
              IF(IOOEXC(IABS,JABS).GT.0) THEN
                IJ_CMP = IOOEXC(IABS,JABS)
                SIGN = 1.0D0
              ELSE
                IJ_CMP = -IOOEXC(IABS,JABS)
                SIGN = -1.0D0
              END IF
*
              IF(IWAY.EQ.2) THEN
                XKAP_CMP(IJ_CMP)   = SIGN*XKAP_FULL(IJ_FULL) 
              ELSE
                XKAP_FULL(IJ_FULL) = SIGN*XKAP_CMP(IJ_CMP)
              END IF
            END IF
*           ^ End of active excitation
            IJ_FULL = IJ_FULL + 1
          END DO
          END DO
*         ^ End of loop over IJ
        END IF
*       ^ End of loop assuring that ISM.GE.JSM
      END DO
*     ^ End of loop over symmetries of orbitals
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Compressed => Full '
        ELSE
          WRITE(6,*) ' Full => Compressed  '
        END IF
*
        WRITE(6,*) ' Excitations in compact form '
        WRITE(6,*) ' =========================== '
        WRITE(6,*) 
        CALL WRT_EXCVEC(XKAP_CMP,IOOEXCC,NOOEXC)
C            WRT_EXCVEC(VEC,IOOEXCC,NOOEXC)
        WRITE(6,*)
        WRITE(6,*) ' Excitations in matrix form '
        WRITE(6,*) ' =========================== '
        WRITE(6,*) 
        CALL APRBLM2(XKAP_FULL,NOCOBS,NOCOBS,NSMOB,1)
       END IF
*
      RETURN
      END
      SUBROUTINE ORB_RELAX(V,IVSM,XKAPPA)
*
* Obtain orbital relaxation term to expectation of 
* one-electron operator V.
*
* The first order orbital relaxation Kappa is input
*
* Jeppe Olsen, Winter of 99
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
*. Specific input : Kappa and V in lower half symmetrypacked forms
C     I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      DIMENSION XKAPPA(*),V(*)
* 
       IDUM = 1
       CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ORB_RE')
*
       CALL MEMMAN(KLVF,NTOOB**2,'ADDL  ',IDUM,'VF    ')
*. Construct F-matrix with V integrals
       RETURN
       END
*
      SUBROUTINE E2_FUSK(ORBHES,NOOEXC,IOOEXCC,IONLYF)  
*
* Fusk calculation of HF orbital Hessian 
*
* Jeppe Olsen, April 99
* Last revision, Aug. 29 2012, Jeppe Olsen, GTIJKL => GTIJKL_GN
*
* Using Formulae E2(AI,BJ ) = 4(delta_{ij}F_ab - delta_{ab}F_ij + 
*                               4g_aibj-g_abij-gajib 
*
* If IONLYF = 1, then only the Fock terms are included 
      INCLUDE 'wrkspc.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INTEGER IOOEXCC(2,*)
*. Output
      DIMENSION ORBHES(*)
*
      INCLUDE 'cintfo.inc'
*
      NTEST = 000
*
*. In the following we will fetch FI elements through GETH1
*  so copy inactive Fock matrix from  zero order space
*. Commented out, normally FI is in KINT1
*.  ... but not always, we are better on the safe side .... 
C     CALL SWAPVE(WORK(KINT1),WORK(KFIZ),NINT1)
*. Well I hope FI is in space
*
      IADR = 0
      DO IAI_EXC = 1, NOOEXC
      DO IBJ_EXC = 1, IAI_EXC 
        IA = IOOEXCC(1,IAI_EXC)
        II = IOOEXCC(2,IAI_EXC)
        IB = IOOEXCC(1,IBJ_EXC)
        IJ = IOOEXCC(2,IBJ_EXC)
*
        XIAJB = 0.0D0
*. One-electron terms
        IF(IA.EQ.IB) XIAJB = XIAJB - GETH1_B(II,IJ)
C       IF(IAI_EXC.EQ.1.AND.IBJ_EXC.EQ.1) 
C    &  WRITE(6,*) ' XIAJB 1 ', XIAJB
        IF(II.EQ.IJ) XIAJB = XIAJB + GETH1_B(IA,IB)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IA II IB IJ = ', IA,II,IB,IJ
          WRITE(6,*) ' Fock contribution the E2 element ', XIAJB
        END IF
C?      IF(IAI_EXC.EQ.1.AND.IBJ_EXC.EQ.1) 
C?   &  WRITE(6,*) ' XIAJB 2 ', XIAJB
*. Two-electron terms
*  4g_aibj-g_abij-gajib 
        IF(IONLYF.NE.1) THEN
          XIAJB = XIAJB + 4*GTIJKL_GN(IA,II,IB,IJ) 
     &                  -   GTIJKL_GN(IA,IB,II,IJ)
     &                  -   GTIJKL_GN(IA,IJ,II,IB)
C?        IF(IAI_EXC.EQ.1.AND.IBJ_EXC.EQ.1) 
C?   &    WRITE(6,*) ' XIAJB 3 ', XIAJB
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IA II IB IJ = ', IA,II,IB,IJ
          WRITE(6,*) ' Fock contribution the E2 element ', XIAJB
        END IF
        END IF
        IADR = IADR + 1
        ORBHES(IADR) = 4*XIAJB
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' Final Hessian element ', ORBHES(IADR)
*
      END DO
      END DO
*. And clean up 
C?    CALL SWAPVE(WORK(KINT1),WORK(KFIZ),NINT1)
*
      IF(NTEST.GE.100) THEN  
        WRITE(6,*) ' Orbital Hessian '
        IF(IONLYF.EQ.1) WRITE(6,*) ' (Only F-terms included ) '
        CALL PRSYM(ORBHES,NOOEXC)
      END IF
*
      RETURN
      END
      FUNCTION GETH1_B(I,J)
*
* Get one-electron matrix elements H(I,J), where I and J are in 
* type order
*
* Last modification, Sept. 24, 2012, Jeppe Olsen, Allow use of ina/sec orbitals 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc' 
*
      ISM = ISMFTO(I)
      ITP = ITPFTO(I)
      IREL = I - IOBPTS_GN(ITP,ISM)+1
*
      JSM = ISMFTO(J)
      JTP = ITPFTO(J)
      JREL = J - IOBPTS_GN(JTP,JSM)+1
* 
      GETH1_B = GETH1E(IREL,ITP,ISM,JREL,JTP,JSM) 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'I, IREL, ITP, ISM = ', I,IREL,ITP,ISM
        WRITE(6,*) 'J, JREL, JTP, JSM = ', J,JREL,JTP,JSM
        WRITE(6,*) ' 1 el integral for I,J = ', I,J, ' is ', GETH1_B
      END IF
*
      RETURN
      END 
      FUNCTION GETH1_B2(I,J,H)
*
* Get one-electron matrix elements H(I,J), where I and J are in 
* type order
*
* Version where integrals are passed throgh input arguments 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc' 
*
      DIMENSION H(*)
*
      ISM = ISMFTO(I)
      ITP = ITPFTO(I)
      IREL = I - IOBPTS(ITP,ISM)+1
*
      JSM = ISMFTO(J)
      JTP = ITPFTO(J)
      JREL = J - IOBPTS(JTP,JSM)+1
* 
      GETH1_B2 = GETH1EX(IREL,ITP,ISM,JREL,JTP,JSM,H)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' 1 el integral for I,J = ', I,J, ' is ', GETH1_B2
      END IF
      RETURN
      END 
      SUBROUTINE TEST_E12
*
* Routine for testing gradient and Hessian routines
*
* Jeppe Olsen, Winter of 99
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      WRITE(6,*)
      WRITE(6,*) ' ======================================='
      WRITE(6,*) ' The world of TEST_E12 has been entered '
      WRITE(6,*) ' ======================================='
      WRITE(6,*)
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'TEST_E')
*
*. First an easy one : The number of nonredundant parameters
*
*. Nonredundant type-type excitations
      CALL MEMMAN(KLTTACT,NGAS*NGAS,'ADDL  ',1,'TTACT ')
      CALL NONRED_TT_EXC(WORK(KLTTACT),1,0)
*. Nonredundant orbital excitations
      CALL MEMMAN(KLOOACT,NTOOB*NTOOB,'ADDL  ',1,'OOACT ')
      CALL MEMMAN(KLOOACTC,2*NTOOB*NTOOB,'ADDL  ',1,'OOACTC')
      CALL NONRED_OO_EXC(NOOEXC,WORK(KLOOACT),WORK(KLOOACTC),
     &                   1,WORK(KLTTACT),2)
*
* Calculate orbital Gradient 
*
      WRITE(6,*) 'calling EN_FROM_DEN'
      CALL EN_FROM_DENS(ENERGY,2,0)
      WRITE(6,'(A,F24.12)') ' Energy from EN_FROM_DEN',ENERGY 
C     WRITE(6,*) ' Construct fock matrix '
      CALL MEMMAN(KLFOO,NTOOB**2,'ADDL  ',2,'FOO   ')
      CALL FOCK_MAT(WORK(KLFOO),2)
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1    ')
      CALL E1_FROM_F(WORK(KLE1),WORK(KLFOO),1,WORK(KLOOACT),
     &               NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
*
* Orbital Hessian
*

      LE2 = NOOEXC*(NOOEXC+1)/2
      CALL MEMMAN(KLE2,LE2   ,'ADDL  ',2,'E2    ')
      CALL ORBHES(WORK(KLE2),WORK(KLOOACT),NOOEXC,1,WORK(KLTTACT))
C          ORBHES(OOHES,IOOEXC,NOOEXC,IOOSM,ITTACT)
      WRITE(6,*) ' End of TEST_E12' 
      STOP       ' End of TEST_E12'
      END
* 
      SUBROUTINE E1_FROM_F(E1,F,IOPSM,IOOEXC,IOOEXCC,
     &           NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)                 
*
* Obtain gradient E1 from Fock matrix as
*
*    E1(ij) = 2F(ij)-2F(ji)
*
* Jeppe Olsen, Jan. 1999   
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER IOOEXC(NTOOB,NTOOB), IOOEXCC(2,NOOEXC)
      INTEGER NTOOBS(NSMOB),IBSO(NSMOB), IREOST(*)
      INCLUDE 'multd2h.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
CNW   DIMENSION F(*)
      integer F, ilo, ihi, ifp
*. Output
      DIMENSION E1(*)
*
      NTEST = 00
*
      call luciaqenter('e1frf')
      ZERO = 0.0D0
C?    WRITE(6,*) ' NOOEXC = ', NOOEXC
      CALL SETVEC(E1,ZERO,NOOEXC)
*
      call ga_distribution(F,ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(F, ilo, ihi, 1, 1, ifp, idum)
      IOFF = 1
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IOPSM)
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)
        IB = IBSO(ISM)
        JB = IBSO(JSM)
C?      WRITE(6,*) ' ISM,JSM,NI,NJ,IB,JB ', ISM,JSM,NI,NJ,IB,JB
*
        DO IORB = 1, NI
          DO JORB = 1, NJ
            IORBABS = IB + IORB - 1
            JORBABS = JB + JORB - 1
*. And in Type-order as expected in IOOEXC
            IORBABS_T = IREOST(IORBABS)
            JORBABS_T = IREOST(JORBABS)
            IJEXC = IOOEXC(IORBABS_T,JORBABS_T)
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' IORB,JORB,IJEXC', IORB,JORB,IJEXC
            IINDX = IOFF-1+(JORB-1)*NI+IORB
            if (IINDX.ge.ilo.and.IINDX.le.ihi) then
            IF(IJEXC.GT.0) THEN
CNW           E1(IJEXC) = E1(IJEXC) + 2.0D0*F(IOFF-1+(JORB-1)*NI+IORB)
              E1(IJEXC) = E1(IJEXC) + 2.0D0*dbl_mb(ifp+IINDX-ilo)       
C?            WRITE(6,*) ' Updated E1 = ',  E1(IJEXC)
            ELSE IF(IJEXC.LT.0) THEN
CNW           E1(-IJEXC) = E1(-IJEXC) - 2.0D0*F(IOFF-1+(JORB-1)*NI+IORB)
              E1(-IJEXC) = E1(-IJEXC) - 2.0D0*dbl_mb(ifp+IINDX-ilo)
C?            WRITE(6,*) ' Updated E1 = ',  E1(-IJEXC)
            END IF
            endif
          END DO
        END DO
        IOFF = IOFF + NI*NJ
      END DO
      call ga_release(F,ilo,ihi,1,1)
      endif
      call ga_sync()
      call ga_dgop(1,E1,NOOEXC,'+')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' E1 vector '
        WRITE(6,*) ' ========='
C       CALL WRTMAT(E1,1,NOOEXC,1,NOOEXC)
C            WRT_EXCVEC(VEC,IOOEXCC,NOOEXC)
        CALL WRT_EXCVEC(E1,IOOEXCC,NOOEXC)
      END IF
      call luciaqexit('e1frf')
*
      RETURN
      END
      SUBROUTINE ORB_RESP(V,IVSM,XKAPPA)
*
* A set of one-integrals V of symmetry IVSM is given
* Obtain corresponding orbital response
*
* Jeppe Olsen, Jan. 1999
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cintfo.inc'
*. Specific input
      DIMENSION V(*)
*. Output
      DIMENSION XKAPPA(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'ORB_RS')
*
*.1 : Number of non-redundant operators of this symmetry and
*.    indeces of non-redundant operators
*. Space for excitation index matrix and Type-type excitation vector
      CALL MEMMAN(KLTTACT,NGAS*NGAS,'ADDL  ',1,'TTACT ')
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
      CALL MEMMAN(KLOOEXCC,2*NTOOB*NTOOB,'ADDL  ',1,'OOEXCC')
*. Matrix giving nonredundant type-type excitations
      CALL NONRED_TT_EXC(WORK(KLTTACT),1,0)
*. Orbital excitation indeces and number of orbital excitations
C          NONRED_OO_EXC(NOOEXC,IOOEXC,IPHSM,ITTACT)
      CALL NONRED_OO_EXC(NOOEXC,WORK(KLOOEXC),WORK(KLOOEXCC),
     &                   IVSM,WORK(KLTTACT),2)
*. Space for E1 vector and E2 matrix
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1VEC ')
      CALL MEMMAN(KLE2,(NOOEXC+1)*NOOEXC/2,'ADDL  ',2,'E2MAT ')

*
*.2 :  Generate gradient like vector <0![E,V]!0>
*. Fock matrix with 
      LV =0
      DO ISM = 1, NSMOB
       JSM = MULTD2H(ISM,IVSM)
       LV = LV + NTOOBS(ISM)*NTOOBS(JSM)
      END DO
      CALL MEMMAN(KLINT1P,NINT1,'ADDL  ',2,'INT1P ')
      CALL COPVEC(WORK(KINT1),WORK(KLINT1P),NINT1)
      CALL COPVEC(V,WORK(KINT1),LV)
      CALL FOCK_MAT2(WORK(KINT1),1,IVSM)
      CALL COPVEC(WORK(KLINT1P),WORK(KINT1),NINT1)
*. F to E1
C      E1_FROM_F(E1,F,IOPSM,IOOEXC,NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO)                 
C      CALL E1_FROM_F(WORK(KLE1),
*.3 :  Construct orbital Hessian
*.4 :  Cholesky decompose Hessian
*.5 :  Solve lin-eq E2*XKappa = - V1
      RETURN
      END
      SUBROUTINE NONRED_TT_EXC(ITTACT,ICISPC_ACT,INTRA)
*
* Set up list of non-redundant excitations between orbital types ITTACT
*
* Inactive is type 0
* active are type 1 - ngas
* secondary are type ngas + 1
*
* Initial version based on division of active orbital spaces into 
* inactive, active, secondary ( I_IAD array)
* 
* Active-active excitations are analyzed further, assuming that 
* the CI-space in question is space ICISPC_ACT. Setting ICISPC_ACT
* to zero eliminates this analysis
*
* If INTRA = 1, then Intra-space excitations are flagged active
*
* Jeppe Olsen, Jan. 99
*              Aug 31., 2003 : ICISPC_ACT added
*              June 10: Explicit inclusion of inactive and secondary
*              August 12: Option INTRA added
*. Last Revision; Oct. 14, 2012; Jeppe Olsen; Local construction of IAD
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      INTEGER ITTACT(0:NGAS+1,0:NGAS+1)
*. Local scratch 
      INTEGER IOC(MXPNGAS+2), I_IADL(MXPNGAS)
*
      NTEST = 000
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'RED_TT')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from NONRED_TT '
      END IF
* 
      CALL GASSPC2(I_IADL,IGSOCCX(1,1,ICISPC_ACT))
C          GASSPC2(I_IADP,IOC)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' I_IADL array '
        CALL IWRTMA(I_IADL,1,NGAS,1,NGAS)
      END IF
*
      DO IGAS = 0, NGAS+1
*. Intra Gas assumed redundant 
        IF(INTRA.EQ.1) THEN
          ITTACT(IGAS,IGAS) = 1
        ELSE
          ITTACT(IGAS,IGAS) = 0
        END IF
*
        DO JGAS = 0, IGAS-1
C?        WRITE(6,*) ' IGAS, JGAS, I_IADL(IGAS),I_IADL(JGAS) ',
C?   &                 IGAS, JGAS, I_IADL(IGAS),I_IADL(JGAS)
*
          I_AM_ACTIVE = 0
*
          IF(IGAS.EQ.NGAS+1.AND.0.LT.JGAS.AND.JGAS.LE.NGAS) THEN
            IF(I_IADL(JGAS).NE.3) I_AM_ACTIVE = 1
          END IF
*
          IF(IGAS.EQ.NGAS+1.AND.JGAS.EQ.0) THEN
           I_AM_ACTIVE = 1
          END IF
*
          IF(JGAS.EQ.0.AND.0.LT.IGAS.AND.IGAS.LT.NGAS+1) THEN
            IF(I_IADL(IGAS).NE.1) I_AM_ACTIVE = 1
          END IF
*
          IF(0.LT.IGAS.AND.IGAS.LT.NGAS+1.AND.
     &       0.LT.JGAS.AND.JGAS.LT.NGAS+1     ) THEN
            IF(I_IADL(IGAS) .NE. I_IADL(JGAS) ) THEN
              I_AM_ACTIVE = 1
            END IF
          END IF
          IF(I_AM_ACTIVE.EQ.1) THEN
*. Non-redundant
            ITTACT(IGAS,JGAS) = 1
            ITTACT(JGAS,IGAS) = -1
          ELSE
*. Assumed redundant
            ITTACT(IGAS,JGAS) = 0
            ITTACT(JGAS,IGAS) = 0
          END IF
        END DO
      END DO
*. Analyze active-active rotations further - if required
      IF(ICISPC_ACT.NE.0) THEN
*. Number of occupation classes 
       IATP = 1
       IBTP = 2
       NAEL = NELEC(IATP)
       NBEL = NELEC(IBTP)
       NEL = NAEL + NBEL
*
C      OCCLSE(IWAY,NOCCLS,IOCCLS,NEL,ICISPC,
C    &        I_DO_BASSPC,IBASSPC,NOBPT)
       CALL OCCLSE(1,NOCCLS,IDUM,NEL,ICISPC_ACT,0,IDUM,NOBPT)
C
       CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'OCCLS ')
       CALL OCCLSE(2,NOCCLS,int_mb(KLOCCLS),NEL,ICISPC_ACT,0,IDUM,NOBPT)
*
       DO IGAS_C = 1, NGAS
        DO IGAS_A = 1, IGAS_C-1
         IF(I_IADL(IGAS_C).EQ.2.AND.I_IADL(IGAS_A).EQ.2) THEN
C                 IS_SX_OP_REDUNDANT(IGAS_C,IGAS_A,NOCCLS,IOCCLS,NGAS)
          IRED =  IS_SX_OP_REDUNDANT(IGAS_C,IGAS_A,NOCCLS,
     &            int_mb(KLOCCLS),NGAS)
          IF(IRED.EQ.0) THEN
             ITTACT(IGAS_C,IGAS_A) = 1
             ITTACT(IGAS_A,IGAS_C) = -1
          END IF
         END IF
        END DO
       END DO
*
      END IF
*     ^ End if active-active rotations should be checked for a given CI space
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' List of nonredundant excitation between types '
        WRITE(6,*) ' INTRA = ', INTRA
        WRITE(6,*) ' (Inactive, gas1, gas2, ..., gasn, secondary) '
        CALL IWRTMA(ITTACT,NGAS+2,NGAS+2,NGAS+2,NGAS+2)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'RED_TT')
      RETURN
      END
      SUBROUTINE NONRED_TT_EXC_CC(ITTACT,ITSPC_ACT,ISPIN)
*
*     version of nonred_tt_exc for the coupled-cluster case
*       uses IHPVGAS_AB in addition to I_IAD
*
*     Andreas, based on Jeppe's original for CI
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
*. Output
      INTEGER ITTACT(0:NGAS+1,0:NGAS+1)
*. Local scratch 
      INTEGER IOC(MXPNGAS)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'RED_TT')
*
C?    WRITE(6,*) ' IADX array '
C?    CALL IWRTMA(I_IADX,1,NGAS,1,NGAS)
*
      DO IGAS = 1, NGAS
        DO JGAS = 1, IGAS
C?        WRITE(6,*) ' IGAS, JGAS, IADX(IGAS),IADX(JGAS) ',
C?   &                 IGAS,JGAS,I_IADX(IGAS),I_IADX(JGAS)
          IF(I_IADX(IGAS) .NE. I_IADX(JGAS) .OR.
     &       IHPVGAS_AB(IGAS,ISPIN).NE.IHPVGAS_AB(JGAS,ISPIN) ) THEN
*. Non-redundant
            ITTACT(IGAS,JGAS) = 1
            ITTACT(JGAS,IGAS) = -1
          ELSE
*. Assumed redundant
            ITTACT(IGAS,JGAS) = 0
            ITTACT(JGAS,IGAS) = 0
          END IF
        END DO
      END DO

*. Analyze active-active rotations further - if required
      IF(ITSPC_ACT.NE.0) THEN
        IF (LCMBSPC(ITSPC_ACT).GT.1) THEN
          WRITE(6,*)
     &         'Sorry, did not consider combination spaces so far!'
          STOP 'NONRED_TT_EXC_CC'
        END IF
     
        DO IGAS_C = 1, NGAS
          DO IGAS_A = 1, IGAS_C-1

            IF(I_IADX(IGAS_C).EQ.2.AND.I_IADX(IGAS_A).EQ.2) THEN
              IRED = 1
              ! rot. between two hole spaces:
              IF(IHPVGAS_AB(IGAS_C,ISPIN).EQ.1.AND.
     &           IHPVGAS_AB(IGAS_A,ISPIN).EQ.1) THEN                
                ! max. number of holes after shell completion different?
                MAXH_A = IGSOCCX(IGAS_A,2,ITSPC_ACT)-
     &                   IGSOCCX(IGAS_A,1,ITSPC_ACT)
                MAXH_C = IGSOCCX(IGAS_C,2,ITSPC_ACT)-
     &                   IGSOCCX(IGAS_C,1,ITSPC_ACT)
                IF (MAXH_A.NE.MAXH_C) IRED=0
              ! rot. between two particle spaces:
              ELSE IF(IHPVGAS_AB(IGAS_C,ISPIN).EQ.2.AND.
     &           IHPVGAS_AB(IGAS_A,ISPIN).EQ.2) THEN
                IF (IGAS_A.LE.1.OR.IGAS_C.LE.1) THEN
                  WRITE(6,*) 'Unexpected: IGAS_A, IGAS_C = ',
     &                        IGAS_A,IGAS_C
                  STOP 'NONRED_TT_EXC_CC'
                END IF
                ! max. number of particles different?
                MAXP_A = IGSOCCX(IGAS_A-1,2,ITSPC_ACT)-
     &                   IGSOCCX(IGAS_A-1,1,ITSPC_ACT) 
                MAXP_C = IGSOCCX(IGAS_C-1,2,ITSPC_ACT)-
     &                   IGSOCCX(IGAS_C-1,1,ITSPC_ACT)
                IF (MAXP_A.NE.MAXP_C) IRED=0
              ! rot. between two valence spaces:
              ELSE IF(IHPVGAS_AB(IGAS_C,ISPIN).EQ.3.AND.
     &           IHPVGAS_AB(IGAS_A,ISPIN).EQ.3) THEN
                WRITE(6,*) 'Sorry, did not consider more than 1 '//
     &                     'valence space so far ...'
                STOP 'NONRED_TT_EXC_CC'                
              END IF
              IF(IRED.EQ.0) THEN
                ITTACT(IGAS_C,IGAS_A) = 1
                ITTACT(IGAS_A,IGAS_C) = -1
              END IF
            END IF

          END DO
        END DO
*
      END IF
*     ^ End if active-active rotations should be checked for a 
*       given T-operator space
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' List of nonredundant excitation between types '
        CALL IWRTMA(ITTACT,NGAS,NGAS,NGAS,NGAS)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'RED_TT')
      RETURN
      END
      SUBROUTINE NONRED_OO_EXC2(NOOEXC,IOOEXC,IOOEXCC,IPHSM,ITTACT,
     &           I_RESTRICT_SUPSYM,MO_SUPSYM,N_INTER_EXC,N_INTRA_EXC,
     &           IFLAG)
*
* Obtain number of non-redundant orbital excitations
* (a+i aj and a+j ai counts as one excitation )
*
* Construct matrix IOOEXC giving index of non-redundant single excitation
* of symmetry IPHSM 
*
* If IFLAG .EQ. 1, only the number of nonvanishing excitations are 
* obtained
*
* Version with supersymmetry
*
* IPHMAT(IORB,JORB) = 0 => excitation between Iorb and Jorb is redundant
* IPHMAT(IORB,JORB) = IJ => a+IORB A JORB is redundant and has number 
*                           ABS(IJ)
*
* Obtain compact list of orbital excitation operators IOOEXCC
*
* Orbital numbers refers to type ordered numbers 
*
* IJ > 0 => Excitation a+Iorb aJorb is excitation operator
*    < 0 => Excitation a+Iorb aJorb is deexcication operator
*
* Jeppe Olsen, Jan. 99
*
* Supersymmetry added May 2012
* Frozen orbitals added, June 2012
* Intra-space excitations added, June 2012
*
*. Info on frozen orbitals obtained from orbinp
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INTEGER ITTACT(0:NGAS+1,0:NGAS+1) 
      INTEGER MO_SUPSYM(NTOOB)
*. Output
      INTEGER IOOEXC(NTOOB,NTOOB), IOOEXCC(2,NTOOB*NTOOB)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from NONRED_OO_EXC2'
        WRITE(6,*) ' ========================'
        WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Supersymmetry of symmetry-ordered orbitals' 
        CALL IWRTMA3(MO_SUPSYM,1,NTOOB,1,NTOOB)
        IF(NFRZ_ORB_ACT.GT.0) THEN
        WRITE(6,*) 
          WRITE(6,*) ' Frozen orbitals: '
          CALL IWRTMA(IFRZ_ORB,1,NFRZ_ORB_ACT,1,NFRZ_ORB_ACT)
        END IF
      END IF 
*
      IF(IFLAG.NE.1) THEN
        IZERO = 0
        CALL ISETVC(IOOEXC,IZERO,NTOOB ** 2 )
      END IF
*
* The interspace excitations
*
      IJEXC = 0
      NSUPSYM_ELI = 0
*. Loop over I orbitals, in type-ordered form
      DO IORB = 1, NTOOB
*. Symmetry and type of IORB
        ISM = ISMFTO(IORB)
        ITP = ITPFTO(IORB)
        DO JORB = 1, IORB - 1
          IF(NTEST.GE.10000) WRITE(6,'(A,2(2X,I3))') ' IORB, JORB = ',
     &    IORB, JORB
          JSM = ISMFTO(JORB)
          JTP = ITPFTO(JORB)
*. Right symmetry ?
          IJSM = MULTD2H(ISM,JSM)
          IMOKAY1 = 0
          IF( IJSM.EQ.IPHSM .AND. ITTACT(ITP,JTP).EQ.1.
     &        AND.ITP.NE.JTP) IMOKAY1 = 1
*. Is one of the orbitals frozen
          DO IFRZ = 1, NFRZ_ORB_ACT
            IF(IFRZ_ORB(IFRZ).EQ.IORB.OR.IFRZ_ORB(IFRZ).EQ.JORB)
     &      IMOKAY1 = 0
          END DO
          IMOKAY2 = 1
          IF(I_RESTRICT_SUPSYM.EQ.1) THEN
            IF(MO_SUPSYM(IREOTS(IORB)).NE.MO_SUPSYM(IREOTS(JORB)))
     &      IMOKAY2 = 0
          END IF
          IF(IMOKAY1.EQ.1.AND.IMOKAY2.EQ.0) THEN
            NSUPSYM_ELI = NSUPSYM_ELI  + 1
            IF(NTEST.GE.200) THEN
            WRITE(6,'(A,2(1X,I3))') 
     &      ' I and J of excitation eliminated by supersym ', IORB,JORB
            WRITE(6,'(A,2(1X,I2))') 
     &      ' I and J in symmetry order ', IREOTS(IORB),IREOTS(JORB)
            END IF
          END IF
          IF(IMOKAY1.EQ.1.AND.IMOKAY2.EQ.1) THEN
*. Correct symmetry and nonredundant
            IJEXC = IJEXC + 1 
            IF(IFLAG.NE.1) THEN
              IOOEXC(IORB,JORB) = IJEXC
              IOOEXC(JORB,IORB) = -IJEXC
*
              IOOEXCC(1,IJEXC) = IORB
              IOOEXCC(2,IJEXC) = JORB
            END IF
          END IF
        END DO
      END DO
      N_INTER_EXC = IJEXC
*
* The intraspace excitations
*
*. Loop over I orbitals, in type-ordered form
      DO IORB = 1, NTOOB
*. Symmetry and type of IORB
        ISM = ISMFTO(IORB)
        ITP = ITPFTO(IORB)
        DO JORB = 1, IORB - 1
          IF(NTEST.GE.10000) WRITE(6,'(A,2(2X,I3))') ' IORB, JORB = ',
     &    IORB, JORB
          JSM = ISMFTO(JORB)
          JTP = ITPFTO(JORB)
*. Right symmetry ?
          IJSM = MULTD2H(ISM,JSM)
          IMOKAY1 = 0
          IF( IJSM.EQ.IPHSM .AND. ITTACT(ITP,JTP).EQ.1.
     &        AND.ITP.EQ.JTP) IMOKAY1 = 1
*. Is one of the orbitals frozen
          DO IFRZ = 1, NFRZ_ORB_ACT
            IF(IFRZ_ORB(IFRZ).EQ.IORB.OR.IFRZ_ORB(IFRZ).EQ.JORB)
     &      IMOKAY1 = 0
          END DO
          IMOKAY2 = 1
          IF(I_RESTRICT_SUPSYM.EQ.1) THEN
            IF(MO_SUPSYM(IREOTS(IORB)).NE.MO_SUPSYM(IREOTS(JORB)))
     &      IMOKAY2 = 0
          END IF
          IF(IMOKAY1.EQ.1.AND.IMOKAY2.EQ.0) THEN
            NSUPSYM_ELI = NSUPSYM_ELI  + 1
            IF(NTEST.GE.200) THEN
            WRITE(6,'(A,2(1X,I3))') 
     &      ' I and J of excitation eliminated by supersym ', IORB,JORB
            WRITE(6,'(A,2(1X,I2))') 
     &      ' I and J in symmetry order ', IREOTS(IORB),IREOTS(JORB)
            END IF
          END IF
          IF(IMOKAY1.EQ.1.AND.IMOKAY2.EQ.1) THEN
*. Correct symmetry and nonredundant
            IJEXC = IJEXC + 1 
            IF(IFLAG.NE.1) THEN
              IOOEXC(IORB,JORB) = IJEXC
              IOOEXC(JORB,IORB) = -IJEXC
*
              IOOEXCC(1,IJEXC) = IORB
              IOOEXCC(2,IJEXC) = JORB
            END IF
          END IF
        END DO
      END DO
*
      NOOEXC = IJEXC
      N_INTRA_EXC = NOOEXC - N_INTER_EXC 
*
      IF(NTEST.GE.10) THEN
        WRITE(6,'(A,I4)') 
     &  '  Number of non-redundant orbital excitations', NOOEXC
        WRITE(6,'(A,2I4)') 
     &  '  Number of inter- and intra-space excitations ',
     &    N_INTRA_EXC, N_INTER_EXC
      END IF
      IF(NTEST.GE.100) THEN
        IF(I_RESTRICT_SUPSYM.EQ.1)
     &  WRITE(6,'(A,I4)') 
     &  ' Number of excitations eliminated by supersym', NSUPSYM_ELI
      END IF
*
      IF(IFLAG.NE.1.AND. NTEST.GE.100) THEN
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' List of non-redundant orbital excitations '
          CALL IWRTMA(IOOEXC,NTOOB,NTOOB,NTOOB,NTOOB)
        END IF
*
        WRITE(6,*) ' Orbital excitations in compact form '
        WRITE(6,*) ' ===================================='
        WRITE(6,*)
        DO IEXC = 1, NOOEXC
          WRITE(6,'(2I5)') IOOEXCC(1,IEXC),IOOEXCC(2,IEXC)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE ORBHES(OOHES,IOOEXC,NOOEXC,IOOSM,ITTACT)
*
* Obtain orbital part of orbital-orbital part of Hessian
* for non-redundant operators of symmetry IOOSM
*
* It is assumed that the Fock matrix, inactive Fock matrix, 
* active matrices have been calculated and stored in 
* KF, KFI/KINT, KFA  has been calculated and is 
* stored in KF whereas the raw one-electron integrals are 
* in KINTO
*
* Jeppe Olsen, January 99
*              Updated with inactive orbitals, June 2010
*
*
* Argument list
*
* OOHES : Space for orbital-orbital hessian, triangular packed
* IOOEXC : Index matrix for indeces of orbital excitations
* NOOEXC : Total number of orbital excitations 
* IOOSM  : Required symmetry of orbital excitations
* ITTACT : Matrix giving nonredundant type-type excitations
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOOEXC(NTOOB,NTOOB)
      INTEGER ITTACT(0:NGAS+1,0:NGAS+1)
*. Output
      DIMENSION OOHES(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ------------------------'
        WRITE(6,*) ' Information from ORBHES '
        WRITE(6,*) ' ------------------------'
        WRITE(6,*)
      END IF
*
      IDUM = 0
      CALL LUCIAQENTER('ORBHE')
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'ORBHES')
*. Scratch : five four index blocks and the generalized Fock matrix 
*
      MXTSOB = 0
      DO ISM = 1, NSMOB
       DO IGAS = 0, NGAS+1
         MXTSOB = MAX(MXTSOB,NOBPTS_GN(IGAS,ISM))
       END DO
      END DO
      LEN4 = MXTSOB * MXTSOB * MXTSOB * MXTSOB
      CALL MEMMAN(KLSCR1,LEN4,'ADDL  ',2,'4IND_1')
      CALL MEMMAN(KLSCR2,LEN4,'ADDL  ',2,'4IND_2')
      CALL MEMMAN(KLSCR3,LEN4,'ADDL  ',2,'4IND_3')
      CALL MEMMAN(KLSCR4,LEN4,'ADDL  ',2,'4IND_4')
      CALL MEMMAN(KLSCR5,LEN4,'ADDL  ',2,'4IND_5')
*
      DO IPTP = 0, NGAS+1
      DO IQTP = 0, NGAS+1
        IPQIND = IQTP*(NGAS+2)+IPTP
C?      WRITE(6,'(A,3I4)') ' IPTP, IQTP, ITTACT( = ',
C?   &                       IPTP, IQTP, ITTACT(IPTP,IQTP)
        IF(ITTACT(IPTP,IQTP).EQ.1) THEN
          DO IPSM = 1, NSMOB
            IQSM = MULTD2H(IPSM,IOOSM)
            DO IRTP = 0, NGAS+1
            DO ISTP = 0, NGAS+1
              IRSIND = ISTP*(NGAS+2)+IRTP
              DO IRSM = 1, NSMOB
C?            WRITE(6,'(A,3I4)') ' IRTP, ISTP, ITTACT( = ',
C?   &                             IRTP, ISTP, ITTACT(IRTP,ISTP)
              IF(ITTACT(IRTP,ISTP).EQ.1.AND.
     &           (IPQIND.GT.IRSIND.OR.
     &            (IPQIND.EQ.IRSIND.AND.IPSM.GE.IRSM)) )THEN
                ISSM = MULTD2H(IRSM,IOOSM)
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Block of Hessian to be calculated : '
                  WRITE(6,'(A,8I4)') 
     &            ' Type and sym for P,Q,R,S (IPTP IPSM ..) ',  
     &              IPTP, IPSM, IQTP, IQSM, IRTP, IRSM, ISTP, ISSM   
                END IF
*. Obtain Hessian block E2(P,Q,R,S)
                NP = NOBPTS_GN(IPTP,IPSM)
                NQ = NOBPTS_GN(IQTP,IQSM)
                NR = NOBPTS_GN(IRTP,IRSM)
                NS = NOBPTS_GN(ISTP,ISSM)
                IF(NP*NQ*NR*NS.NE.0) THEN
                  CALL GET_E2BLK(dbl_mb(KLSCR5),NP,IPTP,IPSM,
     &            NQ,IQTP,IQSM,NR,IRTP,IRSM,NS,ISTP,ISSM,
     &            KF,KFI,KFA,KINT1O,
     &            dbl_mb(KLSCR1),dbl_mb(KLSCR2),dbl_mb(KLSCR3),
     &            dbl_mb(KLSCR4))
                END IF
*. Scatter out
C               IPOFF = IBSO(IPSM)
C               DO IPPTP = 1, IPTP-1
C                 IPOFF = IPOFF + NOBPTS(IPPTP,IPSM)
C               END DO
*
C               IQOFF = IBSO(IQSM)-1
C               DO IQQTP = 1, IQTP
C                 IQOFF = IQOFF + NOBPTS(IQQTP,IQSM)
C               END DO
*
C               IROFF = IBSO(IRSM)
C               DO IRRTP = 1, IRTP-1
C                 IROFF = IROFF + NOBPTS(IRRTP,IRSM)
C               END DO
*
C               ISOFF = IBSO(ISSM)-1
C               DO ISSTP = 1, ISTP
C                 ISOFF = ISOFF + NOBPTS(ISSTP,ISSM)
C               END DO
                IPOFF = IOBPTS_GN(IPTP,IPSM)
                IQOFF = IOBPTS_GN(IQTP,IQSM)
                IROFF = IOBPTS_GN(IRTP,IRSM)
                ISOFF = IOBPTS_GN(ISTP,ISSM)
C?              WRITE(6,'(A,4I3)') ' IPOFF, IQOFF, IROFF, ISOFF',
C?   &                               IPOFF, IQOFF, IROFF, ISOFF
*
                DO IP = 1, NP
                DO IQ = 1, NQ
                DO IR = 1, NR
                DO IS = 1, NS
                  IPABS = IP + IPOFF-1
                  IQABS = IQ + IQOFF-1
                  IRABS = IR + IROFF-1
                  ISABS = IS + ISOFF-1
                  IPQ = IOOEXC(IPABS,IQABS)
                  IRS = IOOEXC(IRABS,ISABS)
                  IF(IPQ.GT.0.AND.IRS.GT.0) THEN
C?                WRITE(6,'(A,4I3)') ' IPABS, IQABS, IRABS, ISABS ',
C?   &            IPABS, IQABS, IRABS, ISABS
C?                WRITE(6,*) ' IPQ, IRS = ', IPQ, IRS
                  IPQRS = MAX(IPQ,IRS)*(MAX(IPQ,IRS)-1)/2+
     &                    MIN(IPQ,IRS)
                  I2PQRS = ((IS-1)*NR+IR-1)*NP*NQ+(IQ-1)*NP + IP
C?                WRITE(6,*) ' IPQRS, I2PQRS', IPQRS,I2PQRS
                  OOHES(IPQRS) = dbl_mb(KLSCR5-1+I2PQRS)
                  END IF
                END DO
                END DO
                END DO
                END DO
*               ^ End of loop over p,q,r,s
              END IF
*              ^ End if block should be calculated
              END DO
*             ^ End of loop over irsm
            END DO
            END DO
*           ^ End of loop over irtp, istp
          END DO
*         ^ End of loop over IPSM
        END IF
*       ^ End if pq op is nonredundant
      END DO
      END DO
*     ^ End of loop over ipsm, iqsm
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================'
        WRITE(6,*) ' Orbital-Orbital Hessian '
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        CALL PRSYM(OOHES,NOOEXC)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ORBHES')
      CALL LUCIAQEXIT('ORBHE')
*
      RETURN
      END
      SUBROUTINE GET_E2BLK_NORT(E2BLK,
     &           NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &           NPP,IPPTP,IPPSM,NQP,IQPTP,IQPSM,
     &                         F, FI, FA, H, T,
     &                         SCR1,SCR2,SCR3,SCR4)
*
* Obtain block of orbital-orbital Hessian as
*
* E2(p,q,p',q') = 0.5*(1+P(pq,p'q'))(1-P(p,q))(1-P(p',q'))
*                [2D(pp')h(qq')-(F(pp')+F(p'p))*T(qq') + 2Y(pqp'q')]
*
* Where T is the transformation matrix from the bio to the normal
* basis, i.e. S
*
* Y(pqp'q') = Sum_(rs:ina+act) 
*    d(qsrq') (ps!rp')
* +  d(qq'rs) (pp'!rs)
*+1/2d(qrq's) (pr!p's)
*+1/2d(rqsq') (rp!sp')
*
*
* Modifications arises from the reduction of Y to inactive/active
* orbitals
* 
* Jeppe Olsen, June 2012, Corrected ans finished, June 2013
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
*. To access blocks of one-electron integrals
      INCLUDE 'intform.inc'
*. The one-electron matrices are supposed to be in a form
*. where all indeces are in the original basis
      DIMENSION F(*), FI(*), FA(*), H(*)
*. Output 
      DIMENSION E2BLK(NP*NQ*NPP*NQP)
      INCLUDE 'cintfo.inc'
*. Scratch : each SCR matrix should be able to hold largest 4 index matrix
*            of orbitals with given type-symmetry 
*
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Hessian block E2(pq),(rs) will be calculated for '
       WRITE(6,'(A,2I4)') ' Type and symmetry of P ' ,IPTP, IPSM
       WRITE(6,'(A,2I4)') ' Type and symmetry of Q ' ,IQTP, IQSM
       WRITE(6,'(A,2I4)') ' Type and symmetry of P'' ' ,IPPTP, IPPSM
       WRITE(6,'(A,2I4)') ' Type and symmetry of Q'' ' ,IQPTP, IQPSM
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E2BLK ')
*
      ZERO = 0.0D0
      CALL SETVEC(E2BLK,ZERO,NP*NQ*NPP*NQP)
*
      LF = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
*
*. Loop over permutations P(pq),P(p'q')
      DO IPQ = 1, 2
       IF(IPQ .EQ. 1 ) THEN
         NPX = NP
         IPXSM = IPSM
         IPXTP = IPTP
         NQX = NQ
         IQXSM = IQSM
         IQXTP = IQTP 
         XPQ = 1.0D0
       ELSE
         NPX = NQ
         IPXSM = IQSM
         IPXTP = IQTP
         NQX = NP
         IQXSM = IPSM
         IQXTP = IPTP
         XPQ = -1.0D0
       END IF
       DO IPPQP = 1, 2
        IF(IPPQP.EQ.1) THEN
          NPPX = NPP
          IPPXSM = IPPSM
          IPPXTP = IPPTP
          NQPX= NQP
          IQPXSM = IQPSM
          IQPXTP = IQPTP
          XPPQP = 1.0D0
        ELSE
          NRR = NQP
          IRRSM = IQPSM
          IRRTP = IQPTP
          NSS = NPP  
          IPQXSM = IPPSM
          IPQXTP = IPPTP
          XPPQP = -1.0D0
        END IF
*
* F terms: -(F(pp')+F(p'p))*delta(qq')
*
        IF(IQXSM.EQ.IQPXSM.AND.IQXTP.EQ.IQPXTP) THEN
*. Obtain F(p,p'), F(p',p) blocks 
          IH1FORM = 2
          CALL COPVEC(F,WORK(KINT1),LF)
*. F(pp')
          CALL GETH1(SCR2,IPXSM,IPXTP,IPPXSM,IPPXTP)
*. F(rp)
          CALL GETH1(SCR3,IPPXSM,IPPXTP,IPXSM,IPXTP)
*. Cleanup
*
          DO IP = 1, NP
          DO IPP = 1, NPP
          DO IQ = 1, NQ
          DO IQP = 1, NQP
            IF(IPQ.EQ.1) THEN
              IPX = IP
              IQX = IQ
            ELSE
              IPX = IQ
              IQX = IP
            END IF
            IF(IPPQP.EQ.1) THEN
              IPPX = IPP
              IQPX = IQP
            ELSE
              IPPX = IQP
              IQPX = IPP
            END IF
            IF(IQX.EQ.IQPX) THEN
              IPQPQ = ((IQP-1)*NPP+IPP-1)*NQ*NP + (IQ-1)*NP+IP  
              E2BLK(IPQPQ) = E2BLK(IPQRS)
     &      - XPQ*XRS*(SCR3((IPX-1)*NPPX+IQPX)+SCR2((IPPX-1)*NPX+IPX))
            END IF
          END DO
          END DO
          END DO
          END DO
*         ^ End of loop  over p,q,p',q'
        END IF
*       ^ End if F- terms are active
*
* h-terms: 2D(pp')h(qq')
*
*. Nonvanishing if indeces PX and PPX corresponds to occupied orbitals
*
        IF(IPXSM.EQ.IPPXSM.AND.
     &     ((1.LE.IPXTP.AND.IPXTP.LE.NGAS.AND.
     &       1.LE.IPPXTP.AND.IPPXTP.LE.NGAS).OR.
     &        IPXTP.EQ.0.AND.IPPXTP.EQ.0       ) ) THEN 
*.h(qx,q'x)
          CALL COPVEC(H,WORK(KINT1),NINT1)
          CALL GETH1(SCR2,IQQSM,IQQTP,ISSSM,ISSTP)
*. D(px,p'x)
          IF(IPXTP.EQ.0.AND.IPPXTP.EQ.0) THEN
*. Inactive part
            ZERO = 0.0D0
            CALL SETVEC(SCR1,ZERO,NPX*NPPX)
            TWO = 2.0D0
            CALL SETDIA(SCR1,TWO,NPX,0)
          ELSE
            CALL GETD1(SCR1,IPXSM,IPXTP,IPPXSM,IPPXTP,1)
          END IF
          DO IP = 1, NP
          DO IQ = 1, NQ
            IF(IPQ.EQ.1) THEN
              IPX = IP
              IQX = IQ
            ELSE
              IPX = IQ
              IQX = IP
            END IF
            DO IPP = 1, NPP
            DO IQP = 1, NQP
              IF(IPPQP.EQ.1) THEN
                IPPX = IPP
                IQPX = IQP
              ELSE
                IPPX = IQP
                IQPX = IPP
              END IF
              IPQPQ = ((IQP-1)*NPP+IPP-1)*NQ*NP + (IQ-1)*NP + IP
              E2BLK(IPQRS) = E2BLK(IPQRS)     
     &       +2.0D0*XPQ*XRS*SCR1((IRR-1)*NPP+IPP)*SCR2((ISS-1)*NQQ+IQQ)
C    &       +2.0D0*XPQ*XRS*SCR1((IPPX-1)*NPX+IPX)*SCR2((IQX-1)*NQX+IQX)
* h-terms: 2D(pp')h(qq')
            END DO
            END DO
*           ^ End of loop over r,s
          END DO
          END DO
*         ^ End of loop over p,q
        END IF
*       ^ End if h-terms contributed
*
* Y-terms
*
*. Contributes only if px and ppx are occupied 
        IF(IPXTP.LE.NGAS.AND.IPPXTP.LE.NGAS) THEN
*. Construct Y-matrix
          CALL GET_YBLK_NORT(SCR1,NPX,IPXTP,IPXSM,NQX,IQXTP,IQXSM,
     &                       NPPX,IPPXTP,IPPXSM,NQPX,IQPXTP,IQPXSM,
     &                       FI,FA,H,
     &                       SCR2,SCR3,SCR4)
          DO IP = 1, NP
          DO IQ = 1, NQ
            IF(IPQ.EQ.1) THEN
              IPX = IP
              IQX = IQ
            ELSE
              IPX = IQ
              IQX = IP
            END IF
            DO IPP = 1, NPP
            DO IQP= 1, NQP
              IF(IPPPQ.EQ.1) THEN
                IPPX = IPP
                IQPX = IQP
              ELSE
                IPPX = IQP
                IQPX = IPP
              END IF
              IPQPQ = ((IQP-1)*NPP+IPP-1)*NQ*NP + (IQ-1)*NP + IP
              I2PQPQ = ((IQPX-1)*NPPX+I-1)*NPP*NQQ+(IQQ-1)*NPP+IPP
C?            WRITE(6,*) '  Updated with Y, IPQRS, I2PQRS, Y(I2PQRS)',
C?   &        IPQRS, I2PQRS, SCR1(I2PQRS)
              TWO = 2.0D0
              E2BLK(IPQRS)= E2BLK(IPQRS) + TWO*XPQ*XRS*SCR1(I2PQRS)
            END DO
            END DO
*           ^ End of loop over r,s
          END DO
          END DO
*         ^ End of loop over p,q
        END IF
*       ^ End if y-terms contributed
       END DO
*      ^ End of loop over IRS
      END DO
*     ^ End of loop over IPQ
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' Hessian block E2(p,q,r,s) '
        WRITE(6,*) ' =========================='
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' Sym  of p,q,r,s', IPSM,IQSM,IRSM,ISSM
        WRITE(6,'(A,4I3)') ' Type of p,q,r,s', IPTP,IQTP,IRTP,ISTP
        WRITE(6,*) ' Matrix in form E2(PQ,RS)'   
        WRITE(6,*)
        CALL WRTMAT(E2BLK,NP*NQ,NR*NS,NP*NQ,NR*NS)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E2BLK ')
*
      RETURN
      END
      SUBROUTINE GET_YBLK(Y,NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &                         NR,IRTP,IRSM,NS,ISTP,ISSM,
     &                         FI,FA,H,
     &                         SCR1,SCR2,SCR3)
*
* Obtain block of Y matrix defined as
*
* Y(pqrs) = Sum_(mn) [d(pmrn)+d(pmnr)](mq!sn) + d(mnpr)(mn!qs)
*
* Jeppe Olsen, Jan. 99
*              July 2010: Modified to treat inactive orbitals
*.             May 2011: Small changes in access of (og!og) integrals
*
* SCR1,SCR2, SCR3, should be able to hold largest 4index block
      IMPLICIT REAL*8(A-H,O-Z)
      REAL * 8 INPROD 
*
*. General input 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
CNW   DIMENSION FI(*), FA(*), H(*)
      integer FI, FA, H
*. Scratch
      DIMENSION SCR1(*),SCR2(*),SCR3(*)
*.Output
      DIMENSION Y(NP*NQ*NR*NS)
*. Y delivered as matrix Y(P,Q,R,S) = Y(PQ,RS)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ===================='
       WRITE(6,*) '  GET_YBLK in action '
       WRITE(6,*) ' ===================='
       WRITE(6,*)
       WRITE(6,'(A,4I3)') ' Sym  of p,q,r,s', IPSM,IQSM,IRSM,ISSM
       WRITE(6,'(A,4I3)') ' Type of p,q,r,s', IPTP,IQTP,IRTP,ISTP
       WRITE(6,*)
      END IF
    
*
      ZERO = 0.0D0
      ONE = 1.0D0
*
      CALL SETVEC(Y(1),ZERO,NP*NQ*NR*NS)
*
      IF(IPTP.EQ.0.AND.IRTP.EQ.0) THEN
*
* ========================
* inactive-inactive block:
* ========================
*
*      Y(PQRS) = 8(PQ!RS) - 2(PS!RQ) - 2(PR!QS) 
*              + 2Delta(P,R)(FI(QS) + FA(QS) - H(QS))
*. 8(PQ!RS)
        CALL LUCIAGETINT(Y,IPTP,IPSM,IQTP,IQSM,IRTP,IRSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
        EIGHT = 8.0D0
        CALL SCALVE(Y,EIGHT,NP*NQ*NR*NS)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' In-In, Y after 8(PQ!RS) '
          CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
        END IF
*. -2(PS!RQ)
        CALL LUCIAGETINT(SCR3,IPTP,IPSM,ISTP,ISSM,IRTP,IRSM,IQTP,IQSM,
     &              0,0,0,1,ONE,ONE)
        DO IS = 1, NS
         DO IR = 1, NR
          DO IQ = 1, NQ
           IOFF_1QRS = (IS-1)*NP*NQ*NR  
     &               + (IR-1)*NP*NQ    
     &               + (IQ-1)*NP
     &               + 1
           IOFF_1SRQ = (IQ-1)*NP*NS*NR 
     &               + (IR-1)*NP*NS    
     &               + (IS-1)*NP 
     &               + 1
           TWOM = -2.0D0
           ONE = 1.0D0
           CALL VECSUM(Y(IOFF_1QRS),Y(IOFF_1QRS),SCR3(IOFF_1SRQ),
     &                 ONE,TWOM,NP)
          END DO
         END DO
        END DO
*       ^ End of loops over Q,R,S
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' In-In, Y after -2(PS!RQ) '
          CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
        END IF
*-2(PR!QS)
        CALL LUCIAGETINT(SCR3,IPTP,IPSM,IRTP,IRSM,IQTP,IQSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
        DO IS = 1, NS
         DO IR = 1, NR
          DO IQ = 1, NQ
           TWOM = -2.0D0
           ONE = 1.0D0
*
           IOFF_1QRS = (IS-1)*NP*NQ*NR  
     &               + (IR-1)*NP*NQ    
     &               + (IQ-1)*NP
     &               + 1
*
           IOFF_1RQS = (IS-1)*NP*NR*NQ  
     &               + (IQ-1)*NP*NR    
     &               + (IR-1)*NP
     &               + 1
           CALL VECSUM(Y(IOFF_1QRS),Y(IOFF_1QRS),SCR3(IOFF_1RQS),
     &                 ONE,TWOM,NP)
          END DO
         END DO
        END DO
*       ^ End of loops over Q,R,S
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' In-In, Y after 2-el terms '
          CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
        END IF
*.+ 2Delta(P,R)(FI(QS) + FA(QS) - H(QS))
        IF(IPSM.EQ.IRSM.AND.IPTP.EQ.IRTP) THEN
*. Offset to sym 1-block with sym IQSM = (ISSM in this case)
        IOFF_H = 1
        DO ISM = 1, IQSM-1
         IOFF_H = IOFF_H + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
        END DO
       
        DO IQ = 1, NQ
         DO IS = 1, NS
*
           IQTO = IOBPTS_GN(IQTP,IQSM) -1 + IQ
           IQSO = IREOTS(IQTO)
           IQREL = IQSO - IBSO(IQSM) + 1
*
           ISTO = IOBPTS_GN(ISTP,ISSM) -1 + IS
           ISSO = IREOTS(ISTO)
           ISREL = ISSO - IBSO(ISSM) + 1
*
           IQS = 
     &     MAX(IQREL,ISREL)*(MAX(IQREL,ISREL)-1)/2 + MIN(IQREL,ISREL)
           IADR = IOFF_H-1+IQS
           call ga_get(FA,IADR,IADR,1,1,XVAL1,1)
           call ga_get(FI,IADR,IADR,1,1,XVAL2,1)
           call ga_get(H,IADR,IADR,1,1,XVAL3,1)
CNW        XFACTOR = 2.0D0*(FA(IADR)+FI(IADR)-H(IADR))
           XFACTOR = 2.0D0*(XVAL1+XVAL2-XVAL3)
*
           DO IP = 1, NP
            IPQPS = (IS-1)*NP*NQ*NP
     &            + (IP-1)*NQ*NP
     &            + (IQ-1)*NP
     &            + IP
            Y(IPQPS) = Y(IPQPS) + XFACTOR
          END DO
*         ^ End of loop over P
         END DO
        END DO
*       ^ End of loop over Q,S
       END IF
*      ^ End if P could be equal to R
      END IF
*     ^ End if P and R are inactive
      IF(0.LT.IPTP.AND.IPTP.LE.NGAS.AND.IRTP.EQ.0) THEN
*
* =================
*. Active-inactive
* =================
*
* sum(m:a) D(PM) (4(MQ!RS)-(MS!RQ)-(MR!QS))
* 
       IF(NTEST.GE.100) WRITE(6,*) ' Active-inactive block'
       IMSM = IPSM
       DO IMTP = 1, NGAS
*. (MQ!RS) in SCR1
        NM = NOBPTS(IMTP,IMSM)
        CALL lUCIAGETINT(SCR1,IMTP,IMSM,IQTP,IQSM,IRTP,IRSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. (MS!RQ) in SCR2
        CALL LUCIAGETINT(SCR2,IMTP,IMSM,ISTP,ISSM,IRTP,IRSM,IQTP,IQSM,
     &              0,0,0,1,ONE,ONE)
*. (MR!QS) in SCR3
        CALL LUCIAGETINT(SCR3,IMTP,IMSM,IRTP,IRSM,IQTP,IQSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. 4(MQ!RS)-(MS!RQ)-(MR!QS) in scr1
        DO IS = 1, NS
         DO IR = 1, NR
          DO IQ = 1, NQ
           IOFF_1QRS = (IS-1)*NM*NQ*NR
     &               + (IR-1)*NM*NQ
     &               + (IQ-1)*NM
     &               + 1
           IOFF_1SRQ = (IQ-1)*NM*NS*NR
     &               + (IR-1)*NM*NS
     &               + (IS-1)*NM
     &               + 1
           IOFF_1RQS = (IS-1)*NM*NR*NQ
     &               + (IQ-1)*NM*NR
     &               + (IR-1)*NM
     &               + 1
           FOUR =  4.0D0
           ONEM = -1.0D0
           CALL VECSUM(SCR1(IOFF_1QRS),SCR1(IOFF_1QRS),SCR2(IOFF_1SRQ),
     &                 FOUR,ONEM,NM)
           ONE = 1.0D0
           CALL VECSUM(SCR1(IOFF_1QRS),SCR1(IOFF_1QRS),SCR3(IOFF_1RQS),
     &                 ONE, ONEM, NM)
          END DO
         END DO
        END DO
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' (M!QRS) = 4(MQ!RS)-(MS!RQ)-(MR!QS)'
          CALL WRTMAT(SCR1,NM,NQ*NR*NS,NM,NQ*NR*NS)
        END IF
*. Obtain Density block D(P,M)
        ONE = 1.0D0
        CALL GETD1(SCR2,IPSM,IPTP,IMSM,IMTP,1)
* sum(m:a) D(PM) (4(MQ!RS)-(MS!RQ)-(MR!QS))
        CALL MATML7(Y,SCR2,SCR1,NP,NQ*NR*NS,NP,NM,NM,NQ*NR*NS,
     &              ONE,ONE,0)
       END DO
*      ^ End of loop over MTP
      END IF
*
      IF(IPTP.EQ.0.AND.0.LT.IRTP.AND.IRTP.LE.NGAS) THEN
*
* ===============
* Inactive-active 
* ===============
*
* 
* sum(m:a) D(RM) (4(MS!PQ)-(MQ!PS)-(MP!QS))
* 
       IF(NTEST.GE.100) WRITE(6,*) ' Inactive-active block'
       IMSM = IRSM
       DO IMTP = 1, NGAS
        NM = NOBPTS(IMTP,IMSM)
*. (MS!PQ) in SCR1
        CALL LUCIAGETINT(SCR1,IMTP,IMSM,ISTP,ISSM,IPTP,IPSM,IQTP,IQSM,
     &              0,0,0,1,ONE,ONE)
*. (MQ!PS) in SCR2
        CALL LUCIAGETINT(SCR2,IMTP,IMSM,IQTP,IQSM,IPTP,IPSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. (MP!QS) in SCR3
        CALL LUCIAGETINT(SCR3,IMTP,IMSM,IPTP,IPSM,IQTP,IQSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. 4(MS!PQ)-(MQ!PS)-(MP!QS) in scr1
        DO IQ = 1, NQ
         DO IP = 1, NP
          DO IS = 1, NS
           IOFF_1SPQ = (IQ-1)*NM*NS*NP
     &               + (IP-1)*NM*NS
     &               + (IS-1)*NM
     &               + 1
           IOFF_1QPS = (IS-1)*NM*NQ*NP
     &               + (IP-1)*NM*NQ
     &               + (IQ-1)*NM
     &               + 1
           IOFF_1PQS = (IS-1)*NM*NP*NQ
     &               + (IQ-1)*NM*NP
     &               + (IP-1)*NM
     &               + 1
           FOUR =  4.0D0
           ONEM = -1.0D0
C?         WRITE(6,*) ' IOFF_1SPQ, IOFF_1QPS, NM = ', 
C?   &                  IOFF_1SPQ, IOFF_1QPS, NM
           CALL VECSUM(SCR1(IOFF_1SPQ),SCR1(IOFF_1SPQ),SCR2(IOFF_1QPS),
     &                 FOUR,ONEM,NM)
           ONE = 1.0D0
           CALL VECSUM(SCR1(IOFF_1SPQ),SCR1(IOFF_1SPQ),SCR3(IOFF_1PQS),
     &                 ONE, ONEM, NM)
          END DO
         END DO
        END DO
*       ^ End of loop over QPS
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' (M!SPQ) =  4(MS!PQ)-(MQ!PS)-(MP!QS) '
          CALL WRTMAT(SCR1,NM,NS*NP*NQ,NM,NS*NP*NQ)
        END IF
*. Obtain Density block D(R,M)
        CALL GETD1(SCR2,IRSM,IRTP,IMSM,IMTP,1)
* sum(m:a) D(RM) (4(MS!PQ)-(MQ!PS)-(MP!QS))
        ONE = 1.0D0
        CALL MATML7(Y,SCR2,SCR1,NR,NS*NP*NQ,NR,NM,NM,NS*NP*NQ,
     &              ONE,ONE,0)
       END DO
*      ^ End of loop over MTP
*. We now have Y(RS,PQ), transpose
       CALL TRPMT3(Y,NR*NS,NP*NQ,SCR1)
       CALL COPVEC(SCR1,Y,NP*NQ*NR*NS)
      END IF
      IF(0.LT.IPTP.AND.IPTP.LE.NGAS.AND.0.LT.IRTP.AND.IRTP.LE.NGAS)THEN
*
* =============
* Active-active
* =============
*
       I_OLD_OR_NEW = 2
*. With I_OLD_OR_NEW = 2, integrals with two general indices are
*. accessed as (oo!gg) or (og!og)

       ONE = 1.0D0
*. Loop over symmetries and types of intermediate indeces M,N
       DO MSM = 1, NSMOB
*. symmetry of N defined by densities and integrals are assumed total sym.
         MPSM = MULTD2H(MSM,IPSM)
         MPRSM = MULTD2H(MPSM,IRSM)
         NSM = MULTD2H(1,MPRSM)
         DO MTP = 1, NGAS
           NM = NOBPTS(MTP,MSM)
           DO NTP = 1, NGAS
            NN = NOBPTS(NTP,NSM)
            IF(I_IAD(MTP).NE.3.AND.I_IAD(NTP).NE.3.AND.
     &         NM*NN.NE.0) THEN
*
* Y(pqrs) <= Sum_(mn) [d(mpnr)+d(mprn)](mq!sn) 
*
*. fetch blocks of integrals and integrals
C LUCIAGETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                  IXCHNG,IKSM,JLSM,ICOUL)
*. 
C     GETD2(RHO2B,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,ICOUL)
*. d(pmrn)
             CALL GETD2(SCR1(1),IPSM,IPTP,MSM,MTP,IRSM,IRTP,NSM,NTP,0)
*. d(pmnr)
             CALL GETD2(SCR2(1),IPSM,IPTP,MSM,MTP,NSM,NTP,IRSM,IRTP,0)
             IF(I_OLD_OR_NEW.EQ.1) THEN
*. (mq!sn)
               CALL LUCIAGETINT(SCR3(1),MTP,MSM,IQTP,IQSM,ISTP,ISSM,NTP,
     &                     NSM,
     &                     0,0,0,1,ONE,ONE)
             ELSE
*. (mq!ns)
               CALL LUCIAGETINT(SCR3(1),MTP,MSM,IQTP,IQSM,NTP,NSM,ISTP,
     &                     ISSM,
     &                     0,0,0,1,ONE,ONE)
             END IF
*
             DO IR = 1, NR
             DO IS = 1, NS
             DO N  = 1, NN
*. Address Y(1,1,R,S)
               IY11RS = ((IS-1)*NR+ IR-1)*NP*NQ+1
*. Address of d(11nr)
               ID11NR = ((IR-1)*NN + N-1)*NM*NP+1
*. Address of d(11rn)
               ID11RN = ((N-1)*NR + IR-1)*NM*NP+1
*. Address of (1 1 ! S N) or (1 1! N S)
               IF(I_OLD_OR_NEW.EQ.1) THEN
                 IG11 = ((N-1)*NS + IS -1 )*NM*NQ+1
               ELSE 
                 IG11 = ((IS-1)*NN + N - 1)*NM*NQ + 1
               END IF
*. And do the matrix mult ( it would of course be more efficient 
*  to first add the two density blocks, but that is not done p.t)
* sum(m) d(pmrn) (mq!sn)
               CALL MATML7(Y(IY11RS),SCR1(ID11RN),SCR3(IG11),
     &                     NP,NQ,NP,NM,NM,NQ,ONE,ONE,0)           
C?            WRITE(6,*) ' Updated Y(1),1', Y(1)
* sum(m) d(pmnr) (mq!sn)
               CALL MATML7(Y(IY11RS),SCR2(ID11NR),SCR3(IG11),
     &                     NP,NQ,NP,NM,NM,NQ,ONE,ONE,0)           
C?            WRITE(6,*) ' Updated Y(1),2', Y(1)
             END DO
             END DO
             END DO
*           ^ End of loops over r,s,n
*
* Y(pqrs) <= Sum_(mn) d(mnpr)(mn!qs)
*
             CALL GETD2(SCR1,MSM,MTP,NSM,NTP,IPSM,IPTP,IRSM,IRTP,0)
C?           WRITE(6,*) ' Block of D2'
C?           CALL WRTMAT(SCR1,NM*NN,NP*NR,NM*NN,NP*NR)
             CALL LUCIAGETINT(SCR3,MTP,MSM,NTP,NSM,IQTP,IQSM,ISTP,ISSM,
     &                   0,0,0,1,ONE,ONE)
C?           WRITE(6,*) ' Block of integrals'
C?           CALL WRTMAT(SCR3,NM*NN,NQ*NS,NM*NN,NQ*NS)
             DO IP = 1, NP
             DO IQ = 1, NQ
             DO IR = 1, NR
             DO IS = 1, NS
                IPQRS = ((IS-1)*NR + IR-1)*NP*NQ + (IQ-1)*NP + IP 
               ID11PR = ((IR-1)*NP + IP -1)* NM*NN + 1
               IG11QS = ((IS-1)*NQ+IQ-1)*NM*NN + 1
               Y(IPQRS) = Y(IPQRS) 
     &                  + INPROD(SCR1(ID11PR),SCR3(IG11QS),NM*NN)
C?            IF(IP.EQ.1. AND. IQ.EQ.1 .AND. IR.EQ.1 .AND .IS.EQ.1) 
C?   &           WRITE(6,*) ' Updated Y(1),3', Y(1)
             END DO
             END DO
             END DO
             END DO
*            ^ End of loops over p,q,r,s
            END IF
           END DO
         END DO
*       ^ End of loop  over types of m,n
       END DO
*      ^ End of loop over symmetry of m
*.+ D(P,R)(FI(QS) - H(QS))
       IF(IPSM.EQ.IRSM) THEN
*. Fetch block of D
        CALL GETD1(SCR1,IPSM,IPTP,IRSM,IRTP,1)
*. Offset to sym 1-block with sym IQSM = (ISSM in this case)
        IOFF_H = 1
        DO ISM = 1, IQSM-1
         IOFF_H = IOFF_H + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
        END DO
       
        DO IQ = 1, NQ
         DO IS = 1, NS
*
           IQTO = IOBPTS_GN(IQTP,IQSM) -1 + IQ
           IQSO = IREOTS(IQTO)
           IQREL = IQSO - IBSO(IQSM) + 1
*
           ISTO = IOBPTS_GN(ISTP,ISSM) -1 + IS
           ISSO = IREOTS(ISTO)
           ISREL = ISSO - IBSO(ISSM) + 1
*
           IQS = 
     &     MAX(IQREL,ISREL)*(MAX(IQREL,ISREL)-1)/2 + MIN(IQREL,ISREL)
           IADR = IOFF_H-1+IQS
           call ga_get(FI,IADR,IADR,1,1,XVAL1,1)
           call ga_get(H,IADR,IADR,1,1,XVAL2,1)
           XFACTOR = XVAL1 - XVAL2
CNW        XFACTOR = (FI(IADR)-H(IADR))
C          WRITE(6,*) ' XFACTOR = ', XFACTOR
*
           DO IR = 1, NR
           DO IP = 1, NP
            IPQRS = (IS-1)*NP*NQ*NR
     &            + (IR-1)*NP*NQ
     &            + (IQ-1)*NP
     &            + IP
            IPR = (IR-1)*NP + IP
            Y(IPQRS) = Y(IPQRS) + XFACTOR*SCR1(IPR)
          END DO
          END DO
*         ^ End of loops over P,R
         END DO
        END DO
*       ^ End of loop over Q,S
       END IF
*      ^ End P and R have same symmetry
      END IF
*     ^ End if P and R are active
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ===================='
       WRITE(6,*) ' Y block as Y(PQ,RS) '
       WRITE(6,*) ' ===================='
       WRITE(6,*)
       WRITE(6,'(A,4I3)') ' Sym  of p,q,r,s', IPSM,IQSM,IRSM,ISSM
       WRITE(6,'(A,4I3)') ' Type of p,q,r,s', IPTP,IQTP,IRTP,ISTP
       WRITE(6,*)
       CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
      END IF
*
      RETURN
      END
* 
      SUBROUTINE DIAG_ORBHES(DIAHES,FOCK,
     &     IOOEXC,NOOEXC,IOOSM,ITTACT,ISPC)
*
* Obtain diagonal part of orbital-orbital part of Hessian
* for non-redundant operators of symmetry IOOSM
*
*  from Jeppe's ORBHES routine: Andreas, summer 2004
*
* Argument list
*
* DIAHES : Space for diagonal of orbital-orbital hessian
* IOOEXC : Index matrix for indeces of orbital excitations
* NOOEXC : Total number of orbital excitations 
* IOOSM  : Required symmetry of orbital excitations
* ITTACT : Matrix giving nonredundant type-type excitations
*
*  ISPC  : 0 for closed shell cases
*          (be sure that I_UNRORB==0 on /OPER/)
*          1 get alpha-alpha contributions
*                   (please provide alpha-part of FOCK)
*          2 get beta-beta contributions (provide
*                   (please provide beta-part of FOCK)
*          3 get alpha-beta contributions
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
*. Input
      INTEGER IOOEXC(NTOOB,NTOOB)
      INTEGER ITTACT(0:NGAS+1,0:NGAS+1)
*. Output
      DIMENSION DIAHES(*)
*
      NTEST = 00
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'ORBHES')
*. Scratch : five two index blocks and the generalized Fock matrix 
*
      MXTSOB = 0
      DO ISM = 1, NSMOB
       DO IGAS = 1, NGAS
         MXTSOB = MAX(MXTSOB,NOBPTS(IGAS,ISM))
       END DO
      END DO
      LEN4 = MXTSOB * MXTSOB
      CALL MEMMAN(KLSCR1,MAX(NINT1,LEN4),'ADDL  ',2,'4IND_1')
      CALL MEMMAN(KLSCR2,LEN4,'ADDL  ',2,'4IND_2')
      CALL MEMMAN(KLSCR3,LEN4,'ADDL  ',2,'4IND_3')
      CALL MEMMAN(KLSCR4,LEN4,'ADDL  ',2,'4IND_4')
      CALL MEMMAN(KLSCR5,LEN4,'ADDL  ',2,'4IND_5')
*. Set up generalized fock matrix. It is assumed that the two-electron density
*  has been calculated 
*  the generalized fock matrix is already on FOCK
c      CALL FOCK_MAT(WORK(KLFG),2) 
*
      DO IPTP = 1, NGAS
      DO IQTP = 1, NGAS
        IPQIND = (IQTP-1)*NGAS+IPTP
        IF(ITTACT(IPTP,IQTP).EQ.1) THEN
          DO IPSM = 1, NSMOB
            IQSM = MULTD2H(IPSM,IOOSM)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Block of Hessian to be calculated : '
              WRITE(6,'(A,8I4)') 
     &            ' Type and sym for P,Q (IPTP IPSM ..) ',  
     &               IPTP, IPSM, IQTP, IQSM
            END IF
*. Obtain Hessian block E2(P,Q,R,S)
            NP = NOBPTS(IPTP,IPSM)
            NQ = NOBPTS(IQTP,IQSM)
            CALL GET_E2BLK_DIA(WORK(KLSCR5),NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &          FOCK,WORK(KINT1O),
     &          WORK(KLSCR1),WORK(KLSCR2),
     &          WORK(KLSCR3),WORK(KLSCR4),ISPC)

            IPOFF = IOBPTS(IPTP,IPSM)
            IQOFF = IOBPTS(IQTP,IQSM)
*     
            DO IP = 1, NP
              DO IQ = 1, NQ
                IPABS = IP + IPOFF-1
                IQABS = IQ + IQOFF-1
                IPQ = IOOEXC(IPABS,IQABS)
                I2PQ = (IP-1)*NQ + IQ
                DIAHES(IPQ) = WORK(KLSCR5-1+I2PQ)
              END DO
            END DO
*           ^ End of loop over p,q

          END DO
*         ^ End of loop over IPSM
        END IF
*       ^ End if pq op is nonredundant
      END DO
      END DO
*     ^ End of loop over ipsm, iqsm
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================================='
        WRITE(6,*) '  Diagonal of Orbital-Orbital Hessian '
        WRITE(6,*) ' ====================================='
        WRITE(6,*)
        CALL WRTMAT(DIAHES,NOOEXC,1,NOOEXC,1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ORBHES')
*
      RETURN
      END
      SUBROUTINE GET_E2BLK_DIA(E2BLK,NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &                         F,H,SCR1,SCR2,SCR3,SCR4,ISPC)
*
* Obtain diagonal block of orbital-orbital Hessian as
*
* E2(p,q,p,q) = 2 (1-P(p,q))
*                [2D(pp)h(qq)- 2F(pp)*delta(qq) + 2Y(pqpq)]
*
* Where 
*
* Y(pqpq) = Sum_(mn) [d(pmpn)+d(pmnp)](mq!qn) + d(mnpp)(mn!qq)
* 
* Andreas, adapted from Jeppes original
*
* F is the generalized Fock matrix and H is the one-electron Hamiltonian
*
*     ISPC: see driver routine DIAG_ORBHES
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
*. To access blocks of one-electron integrals
      INCLUDE 'intform.inc'
      INCLUDE 'oper.inc'
*. Output 
      DIMENSION E2BLK(NP*NQ)
      INCLUDE 'cintfo.inc'
*. Scratch : each SCR matrix should be able to hold largest 4 index matrix
*            of orbitals with given type-symmetry 
*
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Hessian block E2(pq),(pq) will be calculated for '
       WRITE(6,'(A,2I4)') ' Type and symmetry of P ' ,IPTP, IPSM
       WRITE(6,'(A,2I4)') ' Type and symmetry of Q ' ,IQTP, IQSM
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E2BLK ')
*
      ZERO = 0.0D0
      CALL SETVEC(E2BLK,ZERO,NP*NQ)
*
*. Loop over permutations P(pq),P(rs)
      DO IPQ = 1, 2
        IF(IPQ .EQ. 1 ) THEN
          NPP = NP
          IPPSM = IPSM
          IPPTP = IPTP
          NQQ = NQ
          IQQSM = IQSM
          IQQTP = IQTP 
          XPQ = 1.0D0
        ELSE
          NPP = NQ
          IPPSM = IQSM
          IPPTP = IQTP
          NQQ = NP
          IQQSM = IPSM
          IQQTP = IPTP
          XPQ = -1.0D0
        END IF
 
*
* F terms
*
*. Obtain F(p,r)
        IF (ISPC.LT.3) THEN

          IH1FORM = 2
          CALL COPVEC(WORK(KINT1),SCR1,NINT1) 
          LF = 0
          DO ISM = 1, NSMOB
            LF = LF + NTOOBS(ISM)**2
          END DO
c          IF (LF.GT.NINT1) STOP 'dimensions in GET_E2BLK_DIA'
          CALL COPVEC(F,WORK(KINT1),LF)
*
          ISPCAS = 1  ! alpha/beta-block was provided by user
*
          DO IP = 1, NP
            DO IQ = 1, NQ
              IF(IPQ.EQ.1) THEN
                IPP = IP
                IQQ = IQ
              ELSE
                IPP = IQ
                IQQ = IP
              END IF
              IPQPQ = (IP-1)*NQ+IQ  
              E2BLK(IPQPQ) = E2BLK(IPQPQ)
c     &           + XPQ*(2D0*SCR2((IPP-1)*NPP+IPP))
     &             -1D0*(2D0*GETH1E(IPP,IPPTP,IPPSM,IPP,IPPTP,IPPSM))
            END DO
          END DO
*         ^ End of loop  over p,q,r,s
*
          CALL COPVEC(SCR1,WORK(KINT1),NINT1)
          IH1FORM = 1
* h-terms
*
*. Nonvanishing if indeces PP and RR corresponds to occupied orbitals
*
          IF(I_IAD(IPPTP).NE.3) THEN
c          IF(.FALSE..AND.I_IAD(IPPTP).NE.3) THEN
            FAC = 2D0
            IF (ISPC.GT.0) FAC = 2D0
*. h(qq,qq)
c          CALL GETH1(SCR2,IQQSM,IQQTP,IQQSM,IQQTP)
*. D(pp,pp)
c          CALL GETD1(SCR1,IPPSM,IPPTP,IPPSM,IPPTP,1)
*
            ISPCAS = ISPC ! here, *we* have to take care
*
            DO IP = 1, NP
            DO IQ = 1, NQ
              IF(IPQ.EQ.1) THEN
                IPP = IP
                IQQ = IQ
              ELSE
                IPP = IQ
                IQQ = IP
              END IF
              IPQPQ = (IP-1)*NQ + IQ
              E2BLK(IPQPQ) = E2BLK(IPQPQ)     
     &             +FAC*GETD1E(IPP,IPPTP,IPPSM,IPP,IPPTP,IPPSM)
     &                     *GETH1E(IQQ,IQQTP,IQQSM,IQQ,IQQTP,IQQSM)

            END DO
            END DO
*           ^ End of loop over p,q
          END IF
*         ^ End if h-terms contributed
*
        END IF  ! ISPC.LT.3
*
* Y-terms
*
*. Contributes only if pp and rr are occupied 
        IF(I_IAD(IPPTP).NE.3) THEN
c        IF(.FALSE..AND.I_IAD(IPPTP).NE.3) THEN
          FAC = 2D0
          IF (ISPC.GT.0) FAC = 2D0  ! 2D0 is the factor
*. Construct Y-matrix
          CALL GET_YBLK_DIA(SCR1,NPP,IPPTP,IPPSM,NQQ,IQQTP,IQQSM,
     &                       SCR2,SCR3,SCR4,ISPC)
          DO IP = 1, NP
          DO IQ = 1, NQ
            IF(IPQ.EQ.1) THEN
              IPP = IP
              IQQ = IQ
            ELSE
              IPP = IQ
              IQQ = IP
            END IF
            IR = NR
            IS = NS
            IPQRS = (IP-1)*NQ + IQ
            I2PQRS = (IPP-1)*NQQ+IQQ
C?            WRITE(6,*) '  Updated with Y, IPQRS, I2PQRS, Y(I2PQRS)',
C?   &        IPQRS, I2PQRS, SCR1(I2PQRS)
            E2BLK(IPQRS)= E2BLK(IPQRS) + FAC*SCR1(I2PQRS)
          END DO
          END DO
*         ^ End of loop over p,q
        END IF
*       ^ End if y-terms contributed
      END DO
*     ^ End of loop over IPQ
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' Hessian block E2(p,q,p,q) '
        WRITE(6,*) ' =========================='
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' Sym  of p,q ', IPSM,IQSM
        WRITE(6,'(A,4I3)') ' Type of p,q ', IPTP,IQTP
        WRITE(6,*) ' Matrix in form E2(PQ,PQ)'   
        WRITE(6,*)
        CALL WRTMAT(E2BLK,NP,NQ,NP,NQ)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E2BLK ')
*
      RETURN
      END
      SUBROUTINE GET_YBLK_DIA(Y,NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &                         SCR1,SCR2,SCR3,ISPC)
*
* Obtain diagonal block of Y matrix defined as
*
* Y(pqpq) = Sum_(mn) [d(pmrn)+d(pmnr)](mq!sn) + d(mnpr)(mn!qs)
*
* actually, we return Y'(pqpq) = Y(pqpq) - Y(pqqp)
*
* AK based on Jeppe Olsens routine, June 2004
*
* ISCR1,ISCR2, ISCR3, should be able to hold largest 2index block
      IMPLICIT REAL*8(A-H,O-Z)
      REAL * 8 INPROD 
*
*. General input 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
c      INCLUDE 'oper.inc'
*. Scratch
      DIMENSION SCR1(*),SCR2(*),SCR3(*)
*.Output
      DIMENSION Y(NP*NQ)
*. Y delivered as matrix Y(P,Q) = Y(PQ,PQ)
*
      ZERO = 0.0D0
      CALL SETVEC(Y(1),ZERO,NP*NQ)
*
      ONE = 1.0D0
*. Loop over symmetries and types of intermediate indeces M,N
      DO MSM = 1, NSMOB
*. symmetry of N defined by densities and integrals are assumed total sym.
        MPSM = MULTD2H(MSM,IPSM)
        MPPSM = MSM
        NSM = MULTD2H(1,MPPSM)
        DO MTP = 1, NGAS
          NM = NOBPTS(MTP,MSM)
          DO NTP = 1, NGAS
            NN = NOBPTS(NTP,NSM)
*
* Y(pqpq) <=  Sum_(mn) [d(mpnp)+d(mppn)](mq!qn) 
* Y(pqpq) <= -Sum_(mn) [d(mpnq)+d(mpqn)](mq!pn) 
*
            ISPCAS = ISPC
            FAC = 1D0
            IF(ISPC.GT.0) FAC = 1.0D0
            DO IP = 1, NP
*. fetch blocks of integrals and integrals
*. d(pmpn)
              CALL GETD2_A(SCR1(1),IP,IPSM,IPTP,0,MSM,MTP,
     &                             IP,IPSM,IPTP,0,NSM,NTP,ISPCAS)
*. d(pmnp)
              CALL GETD2_A(SCR2(1),IP,IPSM,IPTP,0,  MSM, MTP,
     &                             0, NSM, NTP, IP,IPSM,IPTP,ISPCAS)
              CALL VECSUM(SCR1,SCR1,SCR2,1D0,1D0,NM*NN)

              DO IQ = 1, NQ
*. (mq!qn)
                CALL GETH2_A(SCR2(1),0,  MSM, MTP, IQ,IQSM,IQTP,
     &                               IQ,IQSM,IQTP,0,  NSM, NTP, ISPCAS)
                IYPQ = (IP-1)*NQ + IQ
                Y(IYPQ) = Y(IYPQ) + FAC*INPROD(SCR1,SCR2,NM*NN)
                
              END DO

            END DO
* and permutation term:
            ISPCAS = ISPC
            DO IP = 1, NP
              DO IQ = 1, NQ
*. d(pmqn)
                CALL GETD2_A(SCR1(1),IP,IPSM,IPTP,0,MSM,MTP,
     &                             IQ,IQSM,IQTP,0,NSM,NTP,ISPCAS)
*. d(pmnq)
                CALL GETD2_A(SCR2(1),IP,IPSM,IPTP,0,  MSM, MTP,
     &                             0, NSM, NTP, IQ,IQSM,IQTP,ISPCAS)
                CALL VECSUM(SCR1,SCR1,SCR2,1D0,1D0,NM*NN)
*. (mq!pn)
                CALL GETH2_A(SCR2(1),0,  MSM, MTP, IQ,IQSM,IQTP,
     &                               IP,IPSM,IPTP,0,  NSM, NTP, ISPCAS)
                IYPQ = (IP-1)*NQ + IQ
                Y(IYPQ) = Y(IYPQ) - FAC*INPROD(SCR1,SCR2,NM*NN)
                
              END DO

            END DO
*
            IF (ISPC.EQ.3) CYCLE ! no further contributions
*
* Y(pqpq) <=  Sum_(mn) d(mnpp)(mn!qq)
* Y(pqpq) <= -Sum_(mn) d(mnpq)(mn!qp)
*         
*
            ISP_STRT = ISPC

            IF (ISPC.EQ.0) ISP_STOP = 0 ! only one turn
            IF (ISPC.EQ.1) ISP_STOP = 3 ! fetch alpha-alpha and alpha-beta
            IF (ISPC.EQ.2) ISP_STOP = 4 ! fetch beta-beta and beta-alpha

            FAC = 1D0
            IF(ISPC.GT.0) FAC = 1.0D0
            DO ISPCAS = ISP_STRT, ISP_STOP, 2

              ! well, ok, aehem ...
              ! if fetch the integral and densities as (M,N,P,P),
              ! so we have to invert alpha/beta to beta/alpha and v.v.
              ISPCAS_ = ISPCAS
              IF (ISPCAS.EQ.3) ISPCAS_ = 4
              IF (ISPCAS.EQ.4) ISPCAS_ = 3

             DO IP = 1, NP
              CALL GETD2_A(SCR1,0,MSM,MTP,0,NSM,NTP,
     &             IP,IPSM,IPTP,IP,IPSM,IPTP,ISPCAS_)
              DO IQ = 1, NQ
                CALL GETH2_A(SCR3,0,MSM,MTP,0,NSM,NTP,
     &                       IQ,IQSM,IQTP,IQ,IQSM,IQTP,ISPCAS_)
                IYPQ = (IP-1)*NQ + IQ
                Y(IYPQ) = Y(IYPQ) 
     &                 + FAC*INPROD(SCR1,SCR3,NM*NN)
              END DO
             END DO
* and permutation term
             DO IP = 1, NP
              DO IQ = 1, NQ
                CALL GETD2_A(SCR1,0,MSM,MTP,0,NSM,NTP,
     &                     IP,IPSM,IPTP,IQ,IQSM,IQTP,ISPCAS_)
                CALL GETH2_A(SCR3,0,MSM,MTP,0,NSM,NTP,
     &                       IQ,IQSM,IQTP,IP,IPSM,IPTP,ISPCAS_)
                IYPQ = (IP-1)*NQ + IQ
                Y(IYPQ) = Y(IYPQ) 
     &                 - FAC*INPROD(SCR1,SCR3,NM*NN)
              END DO
             END DO
*            ^ End of loops over p,q
            END DO ! ISPCAS-loop
          END DO
        END DO
*      ^ End of loop  over types of m,n
      END DO
*     ^ End of loop over symmetry of m
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ========================'
       WRITE(6,*) ' Diag. Y block as Y(P,Q) '
       WRITE(6,*) ' ========================'
       WRITE(6,*)
       WRITE(6,*) ' IPSM, IPTP, IQSM,IQTP : '
       WRITE(6,'(8I5)') IPSM, IPTP, IQSM,IQTP
       WRITE(6,*)
       CALL WRTMAT(Y,NP,NQ,NP,NQ)
      END IF
*
      RETURN
      END
* 
      SUBROUTINE FOCK_MAT2(F,I12,IVSM)
*
* Construct Fock matrix for general one-electron operator with symmetry IVSM
*
* F(i,j) = SUM(K) V(i,K) * RHO1(j,K)
*
* Modified FOCK_MAT, January 1999
*
* Only I12 = 1 has been programmed 
c      IMPLICIT REAL*8(A-H,O-Z)
*. Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
*
      INCLUDE 'cintfo.inc'
*. Output
      DIMENSION F(*)
*
      NTEST = 0
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'FOO2  ')
*
      ONE = 1.0D0
      ZERO = 0.0D0
*. Largest set of orbitals with given symmetry and type
      MXTSOB = 0
      DO ISM = 1, NSMOB
      DO IGAS = 1, NGAS
        MXTSOB = MAX(MXTSOB,NOBPTS(IGAS,ISM))
      END DO
      END DO
C?    WRITE(6,*) 'MXTSOB = ', MXTSOB
*. Allocate scratch space for 2-electron integrals and 
*. two-electron densities
      MX4IBLK = MXTSOB ** 4
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT ')
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN ')
*. And a block of F
      MX2IBLK = MXTSOB ** 2
      CALL MEMMAN(KLFBLK,MX2IBLK,'ADDL  ',2,'KLFBL ')
*. 
      ONE = 1.0D0
      IJSM = IVSM
      IFOFF = 1
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IVSM)
        NIS = NOCOBS(ISM)
        NJS = NOCOBS(JSM)
*
        DO JGAS = 1, NGAS
          IF(JGAS.EQ.1) THEN
            JJ = 1
          ELSE 
            JJ = JJ + NOBPTS(JGAS-1,JSM)
          END IF
          NJ = NOBPTS(JGAS,IJSM)
          DO IGAS = 1, NGAS
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) 
     &        ' Fock matrix for ISM IGAS JGAS',ISM,IGAS,JGAS
            END IF
            NI = NOBPTS(IGAS,ISM)
            IF(IGAS.EQ.1) THEN
              II = 1
            ELSE 
              II = II + NOBPTS(IGAS-1,ISM)
            END IF
*
*  =======================
*. block F(ism,igas,jsm,jgas)
*  =======================
*
            CALL SETVEC(WORK(KLFBLK),ZERO,NI*NJ)
* 1 : One-electron part 
            KSM = JSM
            DO KGAS = 1, NGAS
              NK = NOBPTS(KGAS,KSM)
*. blocks of one-electron integrals and one-electron density
              CALL GETD1(WORK(KLDEN),JSM,JGAS,KSM,KGAS,1)
              CALL GETH1(WORK(KLINT),ISM,IGAS,KSM,KGAS)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) 
     &          ' 1-e ints for ISM IGAS KGAS ',ISM,IGAS,KGAS
                CALL WRTMAT(WORK(KLINT),NI,NK,NI,NK)
                WRITE(6,*) 
     &          ' 1-e densi for ISM JGAS KGAS ',ISM,JGAS,KGAS
                CALL WRTMAT(WORK(KLDEN),NJ,NK,NJ,NK)
              END IF
*. And then a matrix multiply( they are pretty much in fashion 
*. these days )
              CALL MATML7(WORK(KLFBLK),WORK(KLINT),WORK(KLDEN),
     &                    NI,NJ,NI,NK,NJ,NK,ONE,ONE,2)
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated block '
                 CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
               END IF
 
            END DO
*           ^ End of loop over KGAS
*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One-electron contributions'
              WRITE(6,*) ' =========================='
              CALL WRTMAT(WORK(KLFBLK),NI,NJ,NI,NJ)
            END IF
*
            IF(I12.EQ.2) THEN
*. 2 : Two-electron part
              WRITE(6,*) ' Only I12 = 1 allowed in FOCK_MAT2'
              STOP 'Only I12 = 1 allowed in FOCK_MAT2'
            END IF
*. Block has been constructed , transfer to -complete- 
*. symmetry blocked Fock matrix
            DO J = 1, NJ
              DO I = 1, NI
                F(IFOFF-1+(J+JJ-1-1)*NIS + I+II-1 ) = 
     &          WORK(KLFBLK-1+(J-1)*NI+I)
              END DO
            END DO
*           ^ End of loop over JI
          END DO
        END DO
*       ^ End of loop over IGAS,JGAS
*. Update pointer to start of symmetry block
        IFOFF = IFOFF+NOCOBS(ISM)*NOCOBS(JSM)
      END DO
*     ^ End of loop over symmetry ISM
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from FOO '
        WRITE(6,*) ' ================'
        CALL APRBLM2(F,NOCOBS,NOCOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'FOO   ')
      RETURN
      END
*
      SUBROUTINE WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
*
* Write one-electron integrals with symmetry IVSM
* ISYM = 1 => Only lower triangular matrix included
*
* Jeppe Olsen, Jan. 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER NRPSM(NSMOB),NCPSM(NSMOB)
      INCLUDE 'multd2h.inc'
*. Specific input
      DIMENSION H(*)
*
      IOFF = 1
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IHSM)
        NI = NRPSM(ISM)
        NJ = NCPSM(JSM)
        IF(ISYM.EQ.0.OR.ISM.GT.JSM) THEN
*. Complete block
          WRITE(6,*) ' Block with symmetry ISM, JSM ',ISM,JSM
          CALL WRTMAT(H(IOFF),NI,NJ,NI,NJ)
          IOFF = IOFF + NI*NJ
        ELSE IF (ISYM.EQ.1.AND.ISM.EQ.JSM) THEN
          CALL PRSYM(H(IOFF),NI)
          IOFF = IOFF + NI*(NI+1)/2
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE GET_FN(FN,DEN1N,DEN2N,MAXN,LFOCK)
*
* Construct n'th order contribution to Fock matrix                  
*
*. Contains two terms : 1 : Zero order Hamiltonian with N'th order densities
*                       2 : Perturbation with N-1'th order densities
*
*  Jeppe Olsen, Sping of 99
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
*. Input : one and two-body densities
      DIMENSION DEN1N(*),DEN2N(*)
*. Output 
      DIMENSION FN(*)
*
      NTEST = 000
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_FN')
      CALL MEMMAN(KLFSCR,LFOCK,'ADDL  ',2,'FSCR  ')
*
*. The Zero'th order operator is assumed in place in KFI
      WRITE(6,*) ' Operator in KFI '
      CALL APRBLM2(WORK(KFI),NTOOBS,NTOOBS,NSMOB,1)
*. Eliminate off diagonal elements   
      CALL ZERO_OFFDIAG_BLM(WORK(KFI),NSMOB,NTOOBS,1)
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB **2*(NTOOB**2+1)/2
*
      ONE = 1.0D0
      ONEM = -1.0D0
      DO K = 0, MAXN
        ZERO = 0.0D0
        CALL SETVEC(FN((K-0)*LFOCK+1),ZERO,LFOCK)
        IF(K.GT.0) THEN
*. Full Hamiltonian with K-1 order density
          CALL COPVEC(DEN1N(1+(K-1-0)*LRHO1),WORK(KRHO1),LRHO1)
          CALL COPVEC(DEN2N(1+(K-1-0)*LRHO2),WORK(KRHO2),LRHO2)
          CALL FOCK_MAT(WORK(KLFSCR),2)
          CALL COPVEC(WORK(KLFSCR),FN(1+(K-0)*LFOCK),LFOCK)
*. Subtract Zero order operator with K-1 order densities
          CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)  
          CALL FOCK_MAT(WORK(KLFSCR),1)
          CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)  
          CALL VECSUM(FN(1+(K-0)*LFOCK),FN(1+(K-0)*LFOCK),
     &                WORK(KLFSCR),ONE,ONEM,LFOCK)
        END IF
*. Zero order hamiltonian with k'th order densities
        CALL COPVEC(DEN1N(1+(K-0)*LRHO1),WORK(KRHO1),LRHO1)
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)  
        CALL FOCK_MAT(WORK(KLFSCR),1)
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)  
        CALL VECSUM(FN(1+(K-0)*LFOCK),FN(1+(K-0)*LFOCK),
     &              WORK(KLFSCR),ONE,ONE,LFOCK)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 'Correction to Fock matrix of order =',K
          CALL APRBLM2(FN(1+(K-0)*LFOCK),NOCOBS,NOCOBS,NSMOB,0)
        END IF
*
      END DO
*
      IF(NTEST.GE.100) THEN
* Accumulate corrections to Fock matrix
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLFSCR),ZERO,LFOCK)
        ONE = 1.0D0
        DO K = 0, MAXN
          CALL VECSUM(WORK(KLFSCR),WORK(KLFSCR),
     &         FN(1+(K-0)*LFOCK),ONE,ONE,LFOCK)
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' ============== '
        WRITE(6,*) ' sum(k) Fock(k) '
        WRITE(6,*) ' ============== '
        WRITE(6,*)
        CALL APRBLM2(WORK(KLFSCR),NOCOBS,NOCOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_FN')
*
      RETURN
      END
      FUNCTION IS_SX_OP_REDUNDANT(IGAS_C,IGAS_A,NOCCLS,IOCCLS,NGAS)
*
* A single excitation from space IGAS_A to space IGAS_C 
* is considered. Check if this operator is redundant 
* for a CI space defined by the NOCCLS occupation classes IOCCLS
* pointer IPCISPC
*
* Jeppe Olsen, Aug. 31, 2003
*
* Output 
*  IS_SX_OP_REDUNDANT = 1 => excitation is redundant
*  IS_SX_OP_REDUNDANT = 0 => excitation is notredundant
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'

*. General input
      INTEGER IOCCLS(NGAS,*)
*. Local scratch
      INTEGER JOC(MXPNGAS)
*
      IRED = 1
      DO I = 1, NOCCLS
        CALL ICOPVE(IOCCLS(1,I),JOC,NGAS)
        JOC(IGAS_C) = JOC(IGAS_C) + 1
        JOC(IGAS_A) = JOC(IGAS_A) - 1
*. Vanishing ?
        IF(JOC(IGAS_A).LT.0.OR.JOC(IGAS_C).GT.2*NOBPT(IGAS_C)) THEN
          IVANISH = 1
        ELSE
          IVANISH = 0
        END IF
        INCLUDED = 1
        IF(IVANISH.EQ.0) THEN
* See if list is in block
          INCLUDED = IS_IVEC_IN_LIST(JOC,NGAS,IOCCLS,NOCCLS)
        END IF
        IF(INCLUDED.EQ.0) IRED = 0
      END DO
      IS_SX_OP_REDUNDANT = IRED
*
      RETURN
      END

      FUNCTION IS_IVEC_IN_LIST(IVEC,NELMNT,LIST,LLIST)
* An array of NELMNT integers are given in IVEC
* Is this array included in the vectors given in LIST
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IVEC(NELMNT)
*. General input
      INTEGER LIST(NELMNT,LLIST)
*
      INCLUDED = 0
      DO JLIST = 1, LLIST
        IDENTICAL = 1
        DO JELMNT = 1, NELMNT
          IF(IVEC(JELMNT).NE.LIST(JELMNT,JLIST)) IDENTICAL = 0
        END DO
        IF(IDENTICAL.EQ.1) INCLUDED = 1
      END DO
*
      IS_IVEC_IN_LIST = INCLUDED
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector being tested '
        CALL IWRTMA(IVEC,1,NELMNT,1,NELMNT)
        IF(INCLUDED.EQ.1) THEN
          WRITE(6,*) ' Vector is in list '
        ELSE
          WRITE(6,*) ' Vector is not in list '
        END IF
      END IF
*
      RETURN
      END

      SUBROUTINE FGEN_RELAX(FGEN,F0,DREL,IHPVGAS_AB,LDIM)
*
*     update general F-matrix with relaxation contributions
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'

      LEN=0
      DO ISM = 1,NSMOB
        LEN=LEN+NTOOBS(ISM)*NTOOBS(ISM)
      END DO
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'FGRELX')
      CALL MEMMAN(KLBUFF,LEN,'ADDL  ',2,'BUFFER')
      CALL MEMMAN(KLBUFF2,LEN,'ADDL  ',2,'BUFF2')
      CALL FGEN_RELAX_S(FGEN,F0,DREL,
     &     IHPVGAS_AB,LDIM,WORK(KLBUFF),IBSO,NSMOB,
     &     NTOOBS,NACOB,NTOOB,IREOST,ITPFSO)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'FGRELX')

      RETURN
      END

      SUBROUTINE FGEN_RELAX_S(FGEN,F0,DREL,IHPVGAS_AB,LDIM,BUFF,
     &            IBSO,NSMOB,NTOOBS,NACOB,NTOOB,IREOST,ITPFSO)

      INCLUDE 'implicit.inc'
      DIMENSION FGEN(*),DREL(*),F0(*),BUFF(*)
      INTEGER IBSO(*),NTOOBS(*),IREOST(*),ITPFSO(*),IHPVGAS_AB(LDIM,*)

      NTEST = 0
      IF(NTEST.NE.0) THEN
*
       WRITE(6,*) ' Initial FGEN in symmetry blocked form '
       WRITE(6,*) ' ======================================'
       WRITE(6,*) 
       ISYM = 0
       CALL APRBLM2(FGEN,NTOOBS,NTOOBS,NSMOB,ISYM)

       WRITE(6,*) ' Initial F0 in symmetry blocked form '
       WRITE(6,*) ' ======================================'
       WRITE(6,*) 
       ISYM = 1
       CALL APRBLM2(F0,NTOOBS,NTOOBS,NSMOB,ISYM)

       WRITE(6,*) ' Initial DREL in symmetry blocked form '
       WRITE(6,*) ' ======================================'
       WRITE(6,*) 
       ISYM = 0
       CALL APRBLM2(DREL,NTOOBS,NTOOBS,NSMOB,ISYM)
      END IF
      
*
*.  Assume spatial symmetric fock matrix
*
      IRT = 1

      IF (IRT==0) THEN

      IJSM = 1

      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
        IF(JSM.EQ.0) CYCLE
        BUFF(IJ+1:IJ+NTOOBS(ISM)**2) = 0D0
        DO I = IBSO(ISM),IBSO(ISM) + NTOOBS(ISM)-1
          DO J = IBSO(JSM),IBSO(JSM) + NTOOBS(JSM)-1
            IJ = IJ+1
            IP = IREOST(I)
            JP = IREOST(J)
            KL = 0
            DO KSM = 1, NSMOB
              CALL SYMCOM(2,6,KSM,LSM,IJSM)
              IF (LSM.EQ.0) CYCLE
              DO K = IBSO(KSM),IBSO(KSM) + NTOOBS(KSM)-1
                DO L = IBSO(LSM),IBSO(LSM) + NTOOBS(LSM)-1
                  KL = KL+1
                  KP = IREOST(K)
                  LP = IREOST(L)
                  BUFF(IJ) = BUFF(IJ)
     &                 + DREL(KL)*(GTIJKL(IP,JP,KP,LP)
     &                      -0.5D0*GTIJKL(IP,LP,KP,JP))
                END DO
              END DO
            END DO
            
            IF (ISM.NE.JSM) STOP 'not prepared for that'
            IOFF = 0
            IOFF2 = 0
            DO KSM = 1, ISM-1
              IOFF = IOFF + NTOOBS(KSM)**2
              IOFF2 = IOFF2 + NTOOBS(KSM)*(NTOOBS(KSM)+1)/2
            END DO
            ILEN = NTOOBS(KSM)
            KSM = ISM
            IORB = I - IBSO(ISM) + 1
            JORB = J - IBSO(JSM) + 1
            ! loop over k
            DO K = IBSO(KSM),IBSO(KSM) + NTOOBS(KSM)-1
              KP = IREOST(K)
              KORB = K - IBSO(KSM) + 1
            ! F(ij) += F0(ik)drel(kj)
            !        + Buf(ik)d0(kj)
              FGEN(IJ) = FGEN(IJ) +
     &             F0(IOFF2+(MAX(IORB,KORB)*(MAX(IORB,KORB)-1)/2)+
     &                       MIN(IORB,KORB))*
     &             DREL(IOFF+(JORB-1)*ILEN+KORB)
c     &             BUFF(IOFF+(KORB-1)*ILEN+IORB)*
c     &             D0(KP,JP)
            END DO
            JGAS=ITPFSO(J)
            FAC = 0D0
            IF (IHPVGAS_AB(JGAS,1).EQ.1) FAC = 1.0D0
            IF (IHPVGAS_AB(JGAS,2).EQ.1) FAC = FAC+1.0D0
            FGEN(IJ) = FGEN(IJ) + FAC*BUFF(IJ)

          END DO  ! J
        END DO  ! I

      END DO ! symmetry blocks

      ELSE ! IRT

      IJSM = 1

      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
        IF(JSM.EQ.0) CYCLE
        DO I = IBSO(ISM),IBSO(ISM) + NTOOBS(ISM)-1
          DO J = IBSO(JSM),IBSO(JSM) + NTOOBS(JSM)-1
            IJ = IJ+1
            IP = IREOST(I)
            JP = IREOST(J)

            JGAS=ITPFSO(J)
            FAC = 0D0
            IF (IHPVGAS_AB(JGAS,1).EQ.1) FAC = 1.0D0
            IF (IHPVGAS_AB(JGAS,2).EQ.1) FAC = FAC+1.0D0

            IF (FAC.NE.0D0) THEN

              KL = 0
              DO KSM = 1, NSMOB
                CALL SYMCOM(2,6,KSM,LSM,IJSM)
                IF (LSM.EQ.0) CYCLE
                DO K = IBSO(KSM),IBSO(KSM) + NTOOBS(KSM)-1
                  DO L = IBSO(LSM),IBSO(LSM) + NTOOBS(LSM)-1
                    KL = KL+1
                    KP = IREOST(K)
                    LP = IREOST(L)
                    FGEN(IJ) = FGEN(IJ)
     &                   + FAC*DREL(KL)*(GTIJKL(IP,JP,KP,LP)
     &                                    -0.5D0*GTIJKL(IP,LP,KP,JP))
                  END DO
                END DO
              END DO
            END IF
            
            IF (ISM.NE.JSM) STOP 'not prepared for that'
            IOFF = 0
            IOFF2 = 0
            DO KSM = 1, ISM-1
              IOFF = IOFF + NTOOBS(KSM)*NTOOBS(KSM)
              IOFF2 = IOFF2 + NTOOBS(KSM)*(NTOOBS(KSM)+1)/2
            END DO
            ILEN = NTOOBS(KSM)
            KSM = ISM
            IORB = I - IBSO(ISM) + 1
            JORB = J - IBSO(JSM) + 1
            ! loop over k
            DO K = IBSO(KSM),IBSO(KSM) + NTOOBS(KSM)-1
              KP = IREOST(K)
              KORB = K - IBSO(KSM) + 1
            ! F(ij) += F0(ik)drel(kj)
              FGEN(IJ) = FGEN(IJ) +
     &             F0(IOFF2+(MAX(IORB,KORB)*(MAX(IORB,KORB)-1)/2)+
     &                       MIN(IORB,KORB))*
     &             DREL(IOFF+(JORB-1)*ILEN+KORB)
            END DO

          END DO  ! J
        END DO  ! I

      END DO ! symmetry blocks

      END IF ! IRT
*
      IF(NTEST.NE.0) THEN
*
       WRITE(6,*) ' Updated FGEN in symmetry blocked form '
       WRITE(6,*) ' ======================================'
       WRITE(6,*) 
       ISYM = 0
       CALL APRBLM2(FGEN,NTOOBS,NTOOBS,NSMOB,ISYM)
      END IF
* 
      RETURN
      END
      SUBROUTINE NONRED_OO_EXC(NOOEXC,IOOEXC,IOOEXCC,IPHSM,ITTACT,IFLAG)
*
* Obtain number of non-redundant orbital excitations
* (a+i aj and a+j ai counts as one excitation )
*
* Construct matrix IOOEXC giving index of non-redundant single excitation
* of symmetry IPHSM 
*
* If IFLAG .EQ. 1, only the number of nonvanishing excitations are 
* obtained
*
* IPHMAT(IORB,JORB) = 0 => excitation between Iorb and Jorb is redundant
* IPHMAT(IORB,JORB) = IJ => a+IORB A JORB is redundant and has number 
*                           ABS(IJ)
*
* Obtain compact list of orbital excitation operators IOOEXCC
*
* Orbital numbers refers to type ordered numbers 
*
* IJ > 0 => Excitation a+Iorb aJorb is excitation operator
*    < 0 => Excitation a+Iorb aJorb is deexcication operator
*
* Jeppe Olsen, Jan. 99
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INTEGER ITTACT(0:NGAS+1,0:NGAS+1) 
*. Output
      INTEGER IOOEXC(NTOOB,NTOOB), IOOEXCC(2,NTOOB*NTOOB)
*
      IF(IFLAG.NE.1) THEN
        IZERO = 0
        CALL ISETVC(IOOEXC,IZERO,NTOOB ** 2 )
      END IF
*
      IJEXC = 0
*. Loop over I orbitals, in type-ordered form
      DO IORB = 1, NTOOB
*. Symmetry and type of IORB
        ISM = ISMFTO(IORB)
        ITP = ITPFTO(IORB)
        DO JORB = 1, IORB - 1
          JSM = ISMFTO(JORB)
          JTP = ITPFTO(JORB)
*. Right symmetry ?
          IJSM = MULTD2H(ISM,JSM)
          IF( IJSM.EQ.IPHSM .AND. ITTACT(ITP,JTP).EQ.1) THEN            
*. Correct symmetry and nonredundant
            IJEXC = IJEXC + 1 
            IF(IFLAG.NE.1) THEN
              IOOEXC(IORB,JORB) = IJEXC
              IOOEXC(JORB,IORB) = -IJEXC
*
              IOOEXCC(1,IJEXC) = IORB
              IOOEXCC(2,IJEXC) = JORB
            END IF
          END IF
        END DO
      END DO
      NOOEXC = IJEXC
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
     &  ' Number of non-redundant orbital excitations', NOOEXC
      END IF
*
      IF(IFLAG.NE.1.AND. NTEST.GE.100) THEN
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' List of non-redundant orbital excitations '
          CALL IWRTMA(IOOEXC,NTOOB,NTOOB,NTOOB,NTOOB)
        END IF
*
        WRITE(6,*) ' Orbital excitations in compact form '
        WRITE(6,*) ' ===================================='
        WRITE(6,*)
        DO IEXC = 1, NOOEXC
          WRITE(6,'(2I5)') IOOEXCC(1,IEXC),IOOEXCC(2,IEXC)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE E1_FROM_F_NORT(E1,F1,F2,IOPSM,IOOEXC,IOOEXCC,
     &           NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)                 
*
* Obtain gradient E1 for non-orthogonal expansion from Fock matrices as
*
*    E1(ij) = 2F1(ij)-2F2(ji)
*
* Jeppe Olsen, June 2012
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER IOOEXC(NTOOB,NTOOB), IOOEXCC(2,NOOEXC)
      INTEGER NTOOBS(NSMOB),IBSO(NSMOB), IREOST(*)
      INCLUDE 'multd2h.inc'
*. Specific input
      DIMENSION F1(*),F2(*)
*. Output
      DIMENSION E1(*)
*
      NTEST = 00
*
      ZERO = 0.0D0
      CALL SETVEC(E1,ZERO,NOOEXC)
*
      IOFF = 1
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IOPSM)
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)
        IB = IBSO(ISM)
        JB = IBSO(JSM)
        IF(NTEST.GE.1000)
     &  WRITE(6,*) ' ISM,JSM,NI,NJ,IB,JB ', ISM,JSM,NI,NJ,IB,JB
*
        DO IORB = 1, NI
          DO JORB = 1, NJ
            IORBABS = IB + IORB - 1
            JORBABS = JB + JORB - 1
*. And in Type-order as expected in IOOEXC
            IORBABS_T = IREOST(IORBABS)
            JORBABS_T = IREOST(JORBABS)
            IJEXC = IOOEXC(IORBABS_T,JORBABS_T)
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' IORB,JORB,IJEXC', IORB,JORB,IJEXC
            IF(IJEXC.GT.0) THEN
              E1(IJEXC) = E1(IJEXC) + 2.0D0*F1(IOFF-1+(JORB-1)*NI+IORB)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' Updated element E1(IJEXC) ', E1(IJEXC)
                WRITE(6,*) ' Element in F1: ', IOFF-1+(JORB-1)*NI+IORB
              END IF
            ELSE IF(IJEXC.LT.0) THEN
              E1(-IJEXC) = E1(-IJEXC) -2.0D0*F2(IOFF-1+(JORB-1)*NI+IORB)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' Updated element E1(-IJEXC) ', E1(-IJEXC)
              END IF
            END IF
          END DO
        END DO
        IOFF = IOFF + NI*NJ
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' E1 vector '
        WRITE(6,*) ' ========='
C       CALL WRTMAT(E1,1,NOOEXC,1,NOOEXC)
C            WRT_EXCVEC(VEC,IOOEXCC,NOOEXC)
        CALL WRT_EXCVEC(E1,IOOEXCC,NOOEXC)
      END IF
*
      RETURN
      END
      SUBROUTINE PRINT_ORBEXC_LIST(IOOEXC,NOOEXC_A,NOOEXC_S)
*
* Print list of excitations, defined by IOOEXC. Symmetric
* rotations may be included - assumed at the end of IOOEXC
*
*. Jeppe Olsen, June 2012
*
      INCLUDE 'implicit.inc'
*
      INTEGER IOOEXC(2,*)
*  
C?    WRITE(6,*)  ' TEST: NOOEXC_A,NOOEXC_S = ', 
C?   &                    NOOEXC_A,NOOEXC_S
      IEXC = 0
      DO IAS = 1,2 
        IF(IAS.EQ.1) THEN
         NEXC_L = NOOEXC_A
        ELSE
         NEXC_L = NOOEXC_S
        END IF
        IF(IAS.EQ.1.AND.NEXC_L.NE.0) THEN  
          WRITE(6,*) ' Antisymmetric excitations: '
          WRITE(6,*)
          WRITE(6,*) ' Excitation Index 1   Index 2 '
          WRITE(6,*) ' ============================='
        END IF
        IF(IAS.EQ.2.AND. NEXC_L.NE.0) THEN  
          WRITE(6,*) ' Symmetric excitations: '
          WRITE(6,*)
          WRITE(6,*) ' Excitation Index 1   Index 2 '
          WRITE(6,*) ' ============================='
        END IF
        DO IIEXC = 1, NEXC_L
          IEXC = IEXC + 1
          WRITE(6,'(3X,I5,5X,I3,5X,I3)')
     &    IEXC, IOOEXC(1,IEXC),IOOEXC(2,IEXC)
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE ORBHES_NORT(OOHES,IOOEXC,NOOEXC_A,NOOEXC_S,
     &           IOOSM,ITTACT)
*
* Obtain orbital part of orbital-orbital part of Hessian
* for non-redundant operators of symmetry IOOSM
*
* It is assumed that the Fock matrix, inactive Fock matrix, 
* active matrices have been calculated and are stored in 
* KF1, KF2, KFI/KINT, KFA  has been calculated and is 
* stored in KF whereas the raw one-electron integrals are 
* in KINTO
*
* Nonorthogonal version
*
* At the moment, only the antisymmetric component is calculated here
*
* Jeppe Olsen, Nonorthogonal version, June 2012 
*              Finished, June 2013
*
*
* Argument list
*
* OOHES : Space for orbital-orbital hessian, triangular packed
* IOOEXC : Index matrix for indeces of orbital excitations
* NOOEXC : Total number of orbital excitations 
* IOOSM  : Required symmetry of orbital excitations
* ITTACT : Matrix giving nonredundant type-type excitations
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
*. Input
      INTEGER IOOEXC(NTOOB,NTOOB)
      INTEGER ITTACT(0:NGAS+1,0:NGAS+1)
*. Output
      DIMENSION OOHES(*)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' -----------------------------'
        WRITE(6,*) ' Information from ORBHES_NORT '
        WRITE(6,*) ' -----------------------------'
        WRITE(6,*)
      END IF
*
      IDUM = 0
      CALL LUCIAQENTER('ORBHE')
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'ORBHES')
*. Scratch : five four index blocks and the generalized Fock matrix 
*
      MXTSOB = 0
      DO ISM = 1, NSMOB
       DO IGAS = 0, NGAS+1
         MXTSOB = MAX(MXTSOB,NOBPTS_GN(IGAS,ISM))
       END DO
      END DO
      LEN4 = MXTSOB * MXTSOB * MXTSOB * MXTSOB
      CALL MEMMAN(KLSCR1,LEN4,'ADDL  ',2,'4IND_1')
      CALL MEMMAN(KLSCR2,LEN4,'ADDL  ',2,'4IND_2')
      CALL MEMMAN(KLSCR3,LEN4,'ADDL  ',2,'4IND_3')
      CALL MEMMAN(KLSCR4,LEN4,'ADDL  ',2,'4IND_4')
      CALL MEMMAN(KLSCR5,LEN4,'ADDL  ',2,'4IND_5')
*
      DO IPTP = 0, NGAS+1
      DO IQTP = 0, NGAS+1
        IPQIND = IQTP*(NGAS+2)+IPTP
C?      WRITE(6,'(A,3I4)') ' IPTP, IQTP, ITTACT( = ',
C?   &                       IPTP, IQTP, ITTACT(IPTP,IQTP)
        IF(ITTACT(IPTP,IQTP).EQ.1) THEN
          DO IPSM = 1, NSMOB
            IQSM = MULTD2H(IPSM,IOOSM)
            DO IPPTP = 0, NGAS+1
            DO IQPTP = 0, NGAS+1

              IPPPQIND = IQPTP*(NGAS+2)+IPPTP
              DO IPPSM = 1, NSMOB
C?            WRITE(6,'(A,3I4)') ' IPPTP, IQPTP, ITTACT( = ',
C?   &                             IPPTP, IQPTP, ITTACT(IPPTP,IQPTP)
              IF(ITTACT(IPPTP,IQPTP).EQ.1.AND.
     &           (IPQIND.GT.IPPPQIND.OR.
     &            (IPQIND.EQ.IPPPQIND.AND.IPSM.GE.IPPSM)) )THEN
                IQPSM = MULTD2H(IPPSM,IOOSM)
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Block of Hessian to be calculated : '
                  WRITE(6,'(A,8I4)') 
     &            ' Type and sym for P,Q,P'',Q'' (IPTP IPSM ..) ',  
     &              IPTP, IPSM, IQTP, IQSM, IPPTP, IPPSM, IPQTP, IPQSM   
                END IF
*. Obtain Hessian block E2(P,Q,P',Q')
                NP = NOBPTS_GN(IPTP,IPSM)
                NQ = NOBPTS_GN(IQTP,IQSM)
                NPP = NOBPTS_GN(IPPTP,IPPSM)
                NQP = NOBPTS_GN(IQPTP,IQPSM)
                IF(NP*NQ*NPP*NQP.NE.0) THEN
                  CALL GET_E2BLK_NORT(WORK(KLSCR5),
     &            NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &            NPP,IPPTP,IPPSM,NPQ,IPQTP,IPQSM,
     &            WORK(KF),WORK(KFI),WORK(KFA),WORK(KINT1O),
     &            WORK(KLSCR1),WORK(KLSCR2),WORK(KLSCR3),WORK(KLSCR4))
                END IF
*. Scatter out
                IPOFF = IOBPTS_GN(IPTP,IPSM)
                IQOFF = IOBPTS_GN(IQTP,IQSM)
                IPPOFF = IOBPTS_GN(IPPTP,IPPSM)
                IQPOFF = IOBPTS_GN(IQPTP,IQPSM)
C?              WRITE(6,'(A,4I3)') ' IPOFF, IQOFF, IPPOFF, IQPOFF',
C?   &                               IPOFF, IQOFF, IPPOFF, IQPOFF
*
                DO IP = 1, NP
                DO IQ = 1, NQ
                DO IPP = 1, NPP
                DO IQP = 1, NQP
                  IPABS = IP + IPOFF-1
                  IQABS = IQ + IQOFF-1
                  IPPABS = IPP + IPPOFF-1
                  IQPABS = IQP + IQPOFF-1
                  IPQ = IOOEXC(IPABS,IQABS)
                  IPPQP = IOOEXC(IPPABS,IQPABS)
                  IF(IPQ.GT.0.AND.IPPQP.GT.0) THEN
C?                WRITE(6,'(A,4I3)') ' IPABS, IQABS, IPPABS, IQPABS ',
C?   &            IPABS, IQABS, IPPABS, IQPABS
C?                WRITE(6,*) ' IPQ, IPPQP = ', IPQ, IPPQP
                  IPQPQ = MAX(IPQ,IPPQP)*(MAX(IPQ,IPPQP)-1)/2+
     &                    MIN(IPQ,IPPQP)
                  I2PQPQ = ((IQP-1)*NPP+IPP-1)*NP*NQ+(IQ-1)*NP + IP
C?                WRITE(6,*) ' IPQPQ, I2PQPQ', IPQPQ,I2PQPQ
                  OOHES(IPQPQ) = WORK(KLSCR5-1+I2PQPQ)
                  END IF
                END DO
                END DO
                END DO
                END DO
*               ^ End of loop over p,q,r,s
              END IF
*              ^ End if block should be calculated
              END DO
*             ^ End of loop over irsm
            END DO
            END DO
*           ^ End of loop over irtp, istp
          END DO
*         ^ End of loop over IPSM
        END IF
*       ^ End if pq op is nonredundant
      END DO
      END DO
*     ^ End of loop over ipsm, iqsm
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================'
        WRITE(6,*) ' Orbital-Orbital Hessian '
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        CALL PRSYM(OOHES,NOOEXC_A)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ORBHES')
      CALL LUCIAQEXIT('ORBHE')
*
      RETURN
      END
      SUBROUTINE GET_YBLK_NORT(Y,NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &                         NPP,IPPTP,IPPSM,NQP,IQPTP,IQPSM,
     &                         FI,FA,H,
     &                         SCR1,SCR2,SCR3)
*
* Obtain block of non-orthogonal Y matrix defined as
*
* Y(pqp'q') = Sum_(rs:ina+act) 
*    d(qsrq') (ps!rp')
* +  d(qq'rs) (pp'!rs)
*+1/2d(qrq's) (pr!p's)
*+1/2d(rqsq') (rp!sp')
*
* Jeppe Olsen, June 2012
*
* SCR1,SCR2, SCR3, should be able to hold largest 4index block
      IMPLICIT REAL*8(A-H,O-Z)
      REAL * 8 INPROD 
*
*. General input 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      DIMENSION FI(*), FA(*), H(*)
*. Scratch
      DIMENSION SCR1(*),SCR2(*),SCR3(*)
*.Output
      DIMENSION Y(NP*NQ*NPP*NQP)
*. Y delivered as matrix Y(P,Q,P',Q') = Y(PQ,P'Q')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ===================='
       WRITE(6,*) '  GET_YBLK in action '
       WRITE(6,*) ' ===================='
       WRITE(6,*)
       WRITE(6,'(A,4I3)') ' Sym  of p,q,p'',q'' ', IPSM,IQSM,IPPSM,IQPSM
       WRITE(6,'(A,4I3)') ' Type of p,q,r,s', IPTP,IQTP,IPPTP,IQPTP
       WRITE(6,*)
      END IF
    
*
      ZERO = 0.0D0
      ONE = 1.0D0
*
      CALL SETVEC(Y(1),ZERO,NP*NQ*NPP*NQP)
*
      IF(IQTP.EQ.0.AND.IQPTP.EQ.0) THEN
*
* ========================
* inactive-inactive block:
* ========================
*
*. (standard/2 :  Y(PQRS) = 8(PQ!P'Q') - 2(PQ'!PQ) - 2(PP'!QQ') 
*      Y(PQP'Q') = 8(PQ!Q'P') 
*                - 4(PP'!Q'Q)
*                + 4(PQ! P'Q')
*                - 2(PQ'! P'Q)
*                + 4(QP!Q'P')
*                - 2(Q'P!QP') 
*                + 4Delta(QQ')(FI(PP') + FA(PP') - H(PP'))
*. 8(PQ!Q'P')
        CALL LUCIAGETINT(Y,IPTP,IPSM,IQTP,IQSM,IQPTP,IQPSM,IPPTP,IPPSM,
     &              0,0,0,1,ONE,ONE)
        EIGHT = 8.0D0
        CALL SCALVE(Y,EIGHT,NP*NQ*NPP*NQP)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' In-In, Y after 8(PQ!Q''P'') '
          CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
        END IF
*. -4(PP'!Q'Q)
        CALL LUCIAGETINT(SCR3,IPTP,IPSM,IPPTP,IPPSM,IQPTP,IQPSM,IQTP,
     &              IQSM,
     &              0,0,0,1,ONE,ONE)
        DO IQ = 1, NS
         DO IQP = 1, NQP
          DO IPP = 1, NPP
           IOFF_1QPPQP = (IQP-1)*NP*NQ*NPP
     &               + (IPP-1)*NP*NQ    
     &               + (IQ-1)*NP
     &               + 1
           IOFF_1PPQPQ = (IQ-1)*NP*NPP*NQP
     &               + (IQP-1)*NP*NPP
     &               + (IPP-1)*NP 
     &               + 1
           TWOM = -2.0D0
           ONE = 1.0D0
           CALL VECSUM(Y(IOFF_1QRS),Y(IOFF_1QRS),SCR3(IOFF_1SRQ),
     &                 ONE,TWOM,NP)
          END DO
         END DO
        END DO
*       ^ End of loops over Q,R,S
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' In-In, Y after  -4(PP'!Q'Q) '
          CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
        END IF
*+ 4(PQ! P'Q')
* HER ER JEG
        CALL LUCIAGETINT(SCR3,IPTP,IPSM,IQTP,IQSM,IPTP,IQPSM,IQTP,IQSM,
     &              0,0,0,1,ONE,ONE)
        DO IQ = 1, NS
         DO IQP = 1, NQP
          DO IPP = 1, NPP
           IOFF_1QPPQP = (IQP-1)*NP*NQ*NPP
     &               + (IPP-1)*NP*NQ    
     &               + (IQ-1)*NP
     &               + 1
           IOFF_1PPQPQ = (IQ-1)*NP*NPP*NQP
     &               + (IQP-1)*NP*NPP
     &               + (IPP-1)*NP 
     &               + 1
           TWOM = -2.0D0
           ONE = 1.0D0
           CALL VECSUM(Y(IOFF_1QRS),Y(IOFF_1QRS),SCR3(IOFF_1SRQ),
     &                 ONE,TWOM,NP)
          END DO
         END DO
        END DO
*       ^ End of loops over Q,R,S
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' In-In, Y after -2(PS!RQ) '
          CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
        END IF
*-2(PR!QS)
        CALL LUCIAGETINT(SCR3,IPTP,IPSM,IRTP,IRSM,IQTP,IQSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
        DO IS = 1, NS
         DO IR = 1, NR
          DO IQ = 1, NQ
           TWOM = -2.0D0
           ONE = 1.0D0
*
           IOFF_1QRS = (IS-1)*NP*NQ*NR  
     &               + (IR-1)*NP*NQ    
     &               + (IQ-1)*NP
     &               + 1
*
           IOFF_1RQS = (IS-1)*NP*NR*NQ  
     &               + (IQ-1)*NP*NR    
     &               + (IR-1)*NP
     &               + 1
           CALL VECSUM(Y(IOFF_1QRS),Y(IOFF_1QRS),SCR3(IOFF_1RQS),
     &                 ONE,TWOM,NP)
          END DO
         END DO
        END DO
*       ^ End of loops over Q,R,S
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' In-In, Y after 2-el terms '
          CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
        END IF

*.+ 2Delta(P,R)(FI(QS) + FA(QS) - H(QS))
        IF(IPSM.EQ.IRSM.AND.IPTP.EQ.IRTP) THEN
*. Offset to sym 1-block with sym IQSM = (ISSM in this case)
        IOFF_H = 1
        DO ISM = 1, IQSM-1
         IOFF_H = IOFF_H + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
        END DO
       
        DO IQ = 1, NQ
         DO IS = 1, NS
*
           IQTO = IOBPTS_GN(IQTP,IQSM) -1 + IQ
           IQSO = IREOTS(IQTO)
           IQREL = IQSO - IBSO(IQSM) + 1
*
           ISTO = IOBPTS_GN(ISTP,ISSM) -1 + IS
           ISSO = IREOTS(ISTO)
           ISREL = ISSO - IBSO(ISSM) + 1
*
           IQS = 
     &     MAX(IQREL,ISREL)*(MAX(IQREL,ISREL)-1)/2 + MIN(IQREL,ISREL)
           IADR = IOFF_H-1+IQS
           XFACTOR = 2.0D0*(FA(IADR)+FI(IADR)-H(IADR))
*
           DO IP = 1, NP
            IPQPS = (IS-1)*NP*NQ*NP
     &            + (IP-1)*NQ*NP
     &            + (IQ-1)*NP
     &            + IP
            Y(IPQPS) = Y(IPQPS) + XFACTOR
          END DO
*         ^ End of loop over P
         END DO
        END DO
*       ^ End of loop over Q,S
       END IF
*      ^ End if P could be equal to R
      END IF
*     ^ End if P and R are inactive
      IF(0.LT.IPTP.AND.IPTP.LE.NGAS.AND.IRTP.EQ.0) THEN
*
* =================
*. Active-inactive
* =================
*
* sum(m:a) D(PM) (4(MQ!RS)-(MS!RQ)-(MR!QS))
* 
       IF(NTEST.GE.100) WRITE(6,*) ' Active-inactive block'
       IMSM = IPSM
       DO IMTP = 1, NGAS
*. (MQ!RS) in SCR1
        NM = NOBPTS(IMTP,IMSM)
        CALL LUCIAGETINT(SCR1,IMTP,IMSM,IQTP,IQSM,IRTP,IRSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. (MS!RQ) in SCR2
        CALL LUCIAGETINT(SCR2,IMTP,IMSM,ISTP,ISSM,IRTP,IRSM,IQTP,IQSM,
     &              0,0,0,1,ONE,ONE)
*. (MR!QS) in SCR3
        CALL LUCIAGETINT(SCR3,IMTP,IMSM,IRTP,IRSM,IQTP,IQSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. 4(MQ!RS)-(MS!RQ)-(MR!QS) in scr1
        DO IS = 1, NS
         DO IR = 1, NR
          DO IQ = 1, NQ
           IOFF_1QRS = (IS-1)*NM*NQ*NR
     &               + (IR-1)*NM*NQ
     &               + (IQ-1)*NM
     &               + 1
           IOFF_1SRQ = (IQ-1)*NM*NS*NR
     &               + (IR-1)*NM*NS
     &               + (IS-1)*NM
     &               + 1
           IOFF_1RQS = (IS-1)*NM*NR*NQ
     &               + (IQ-1)*NM*NR
     &               + (IR-1)*NM
     &               + 1
           FOUR =  4.0D0
           ONEM = -1.0D0
           CALL VECSUM(SCR1(IOFF_1QRS),SCR1(IOFF_1QRS),SCR2(IOFF_1SRQ),
     &                 FOUR,ONEM,NM)
           ONE = 1.0D0
           CALL VECSUM(SCR1(IOFF_1QRS),SCR1(IOFF_1QRS),SCR3(IOFF_1RQS),
     &                 ONE, ONEM, NM)
          END DO
         END DO
        END DO
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' (M!QRS) = 4(MQ!RS)-(MS!RQ)-(MR!QS)'
          CALL WRTMAT(SCR1,NM,NQ*NR*NS,NM,NQ*NR*NS)
        END IF
*. Obtain Density block D(P,M)
        ONE = 1.0D0
        CALL GETD1(SCR2,IPSM,IPTP,IMSM,IMTP,1)
* sum(m:a) D(PM) (4(MQ!RS)-(MS!RQ)-(MR!QS))
        CALL MATML7(Y,SCR2,SCR1,NP,NQ*NR*NS,NP,NM,NM,NQ*NR*NS,
     &              ONE,ONE,0)
       END DO
*      ^ End of loop over MTP
      END IF
*
      IF(IPTP.EQ.0.AND.0.LT.IRTP.AND.IRTP.LE.NGAS) THEN
*
* ===============
* Inactive-active 
* ===============
*
* 
* sum(m:a) D(RM) (4(MS!PQ)-(MQ!PS)-(MP!QS))
* 
       IF(NTEST.GE.100) WRITE(6,*) ' Inactive-active block'
       IMSM = IRSM
       DO IMTP = 1, NGAS
        NM = NOBPTS(IMTP,IMSM)
*. (MS!PQ) in SCR1
        CALL LUCIAGETINT(SCR1,IMTP,IMSM,ISTP,ISSM,IPTP,IPSM,IQTP,IQSM,
     &              0,0,0,1,ONE,ONE)
*. (MQ!PS) in SCR2
        CALL LUCIAGETINT(SCR2,IMTP,IMSM,IQTP,IQSM,IPTP,IPSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. (MP!QS) in SCR3
        CALL LUCIAGETINT(SCR3,IMTP,IMSM,IPTP,IPSM,IQTP,IQSM,ISTP,ISSM,
     &              0,0,0,1,ONE,ONE)
*. 4(MS!PQ)-(MQ!PS)-(MP!QS) in scr1
        DO IQ = 1, NQ
         DO IP = 1, NP
          DO IS = 1, NS
           IOFF_1SPQ = (IQ-1)*NM*NS*NP
     &               + (IP-1)*NM*NS
     &               + (IS-1)*NM
     &               + 1
           IOFF_1QPS = (IS-1)*NM*NQ*NP
     &               + (IP-1)*NM*NQ
     &               + (IQ-1)*NM
     &               + 1
           IOFF_1PQS = (IS-1)*NM*NP*NQ
     &               + (IQ-1)*NM*NP
     &               + (IP-1)*NM
     &               + 1
           FOUR =  4.0D0
           ONEM = -1.0D0
C?         WRITE(6,*) ' IOFF_1SPQ, IOFF_1QPS, NM = ', 
C?   &                  IOFF_1SPQ, IOFF_1QPS, NM
           CALL VECSUM(SCR1(IOFF_1SPQ),SCR1(IOFF_1SPQ),SCR2(IOFF_1QPS),
     &                 FOUR,ONEM,NM)
           ONE = 1.0D0
           CALL VECSUM(SCR1(IOFF_1SPQ),SCR1(IOFF_1SPQ),SCR3(IOFF_1PQS),
     &                 ONE, ONEM, NM)
          END DO
         END DO
        END DO
*       ^ End of loop over QPS
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' (M!SPQ) =  4(MS!PQ)-(MQ!PS)-(MP!QS) '
          CALL WRTMAT(SCR1,NM,NS*NP*NQ,NM,NS*NP*NQ)
        END IF
*. Obtain Density block D(R,M)
        CALL GETD1(SCR2,IRSM,IRTP,IMSM,IMTP,1)
* sum(m:a) D(RM) (4(MS!PQ)-(MQ!PS)-(MP!QS))
        ONE = 1.0D0
        CALL MATML7(Y,SCR2,SCR1,NR,NS*NP*NQ,NR,NM,NM,NS*NP*NQ,
     &              ONE,ONE,0)
       END DO
*      ^ End of loop over MTP
*. We now have Y(RS,PQ), transpose
       CALL TRPMT3(Y,NR*NS,NP*NQ,SCR1)
       CALL COPVEC(SCR1,Y,NP*NQ*NR*NS)
      END IF
      IF(0.LT.IPTP.AND.IPTP.LE.NGAS.AND.0.LT.IRTP.AND.IRTP.LE.NGAS)THEN
*
* =============
* Active-active
* =============
*
       I_OLD_OR_NEW = 2
*. With I_OLD_OR_NEW = 2, integrals with two general indices are
*. accessed as (oo!gg) or (og!og)

       ONE = 1.0D0
*. Loop over symmetries and types of intermediate indeces M,N
       DO MSM = 1, NSMOB
*. symmetry of N defined by densities and integrals are assumed total sym.
         MPSM = MULTD2H(MSM,IPSM)
         MPRSM = MULTD2H(MPSM,IRSM)
         NSM = MULTD2H(1,MPRSM)
         DO MTP = 1, NGAS
           NM = NOBPTS(MTP,MSM)
           DO NTP = 1, NGAS
            NN = NOBPTS(NTP,NSM)
            IF(I_IAD(MTP).NE.3.AND.I_IAD(NTP).NE.3.AND.
     &         NM*NN.NE.0) THEN
*
* Y(pqrs) <= Sum_(mn) [d(mpnr)+d(mprn)](mq!sn) 
*
*. fetch blocks of integrals and integrals
C LUCIAGETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                  IXCHNG,IKSM,JLSM,ICOUL)
*. 
C     GETD2(RHO2B,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,ICOUL)
*. d(pmrn)
             CALL GETD2(SCR1(1),IPSM,IPTP,MSM,MTP,IRSM,IRTP,NSM,NTP,0)
*. d(pmnr)
             CALL GETD2(SCR2(1),IPSM,IPTP,MSM,MTP,NSM,NTP,IRSM,IRTP,0)
             IF(I_OLD_OR_NEW.EQ.1) THEN
*. (mq!sn)
               CALL LUCIAGETINT(SCR3(1),MTP,MSM,IQTP,IQSM,ISTP,ISSM,NTP,
     &                     NSM,
     &                     0,0,0,1,ONE,ONE)
             ELSE
*. (mq!ns)
               CALL lUCIAGETINT(SCR3(1),MTP,MSM,IQTP,IQSM,NTP,NSM,ISTP,
     &                     ISSM,
     &                     0,0,0,1,ONE,ONE)
             END IF
*
             DO IR = 1, NR
             DO IS = 1, NS
             DO N  = 1, NN
*. Address Y(1,1,R,S)
               IY11RS = ((IS-1)*NR+ IR-1)*NP*NQ+1
*. Address of d(11nr)
               ID11NR = ((IR-1)*NN + N-1)*NM*NP+1
*. Address of d(11rn)
               ID11RN = ((N-1)*NR + IR-1)*NM*NP+1
*. Address of (1 1 ! S N) or (1 1! N S)
               IF(I_OLD_OR_NEW.EQ.1) THEN
                 IG11 = ((N-1)*NS + IS -1 )*NM*NQ+1
               ELSE 
                 IG11 = ((IS-1)*NN + N - 1)*NM*NQ + 1
               END IF
*. And do the matrix mult ( it would of course be more efficient 
*  to first add the two density blocks, but that is not done p.t)
* sum(m) d(pmrn) (mq!sn)
               CALL MATML7(Y(IY11RS),SCR1(ID11RN),SCR3(IG11),
     &                     NP,NQ,NP,NM,NM,NQ,ONE,ONE,0)           
C?            WRITE(6,*) ' Updated Y(1),1', Y(1)
* sum(m) d(pmnr) (mq!sn)
               CALL MATML7(Y(IY11RS),SCR2(ID11NR),SCR3(IG11),
     &                     NP,NQ,NP,NM,NM,NQ,ONE,ONE,0)           
C?            WRITE(6,*) ' Updated Y(1),2', Y(1)
             END DO
             END DO
             END DO
*           ^ End of loops over r,s,n
*
* Y(pqrs) <= Sum_(mn) d(mnpr)(mn!qs)
*
             CALL GETD2(SCR1,MSM,MTP,NSM,NTP,IPSM,IPTP,IRSM,IRTP,0)
C?           WRITE(6,*) ' Block of D2'
C?           CALL WRTMAT(SCR1,NM*NN,NP*NR,NM*NN,NP*NR)
             CALL lUCIAGETINT(SCR3,MTP,MSM,NTP,NSM,IQTP,IQSM,ISTP,ISSM,
     &                   0,0,0,1,ONE,ONE)
C?           WRITE(6,*) ' Block of integrals'
C?           CALL WRTMAT(SCR3,NM*NN,NQ*NS,NM*NN,NQ*NS)
             DO IP = 1, NP
             DO IQ = 1, NQ
             DO IR = 1, NR
             DO IS = 1, NS
                IPQRS = ((IS-1)*NR + IR-1)*NP*NQ + (IQ-1)*NP + IP 
               ID11PR = ((IR-1)*NP + IP -1)* NM*NN + 1
               IG11QS = ((IS-1)*NQ+IQ-1)*NM*NN + 1
               Y(IPQRS) = Y(IPQRS) 
     &                  + INPROD(SCR1(ID11PR),SCR3(IG11QS),NM*NN)
C?            IF(IP.EQ.1. AND. IQ.EQ.1 .AND. IR.EQ.1 .AND .IS.EQ.1) 
C?   &           WRITE(6,*) ' Updated Y(1),3', Y(1)
             END DO
             END DO
             END DO
             END DO
*            ^ End of loops over p,q,r,s
            END IF
           END DO
         END DO
*       ^ End of loop  over types of m,n
       END DO
*      ^ End of loop over symmetry of m
*.+ D(P,R)(FI(QS) - H(QS))
       IF(IPSM.EQ.IRSM) THEN
*. Fetch block of D
        CALL GETD1(SCR1,IPSM,IPTP,IRSM,IRTP,1)
*. Offset to sym 1-block with sym IQSM = (ISSM in this case)
        IOFF_H = 1
        DO ISM = 1, IQSM-1
         IOFF_H = IOFF_H + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
        END DO
       
        DO IQ = 1, NQ
         DO IS = 1, NS
*
           IQTO = IOBPTS_GN(IQTP,IQSM) -1 + IQ
           IQSO = IREOTS(IQTO)
           IQREL = IQSO - IBSO(IQSM) + 1
*
           ISTO = IOBPTS_GN(ISTP,ISSM) -1 + IS
           ISSO = IREOTS(ISTO)
           ISREL = ISSO - IBSO(ISSM) + 1
*
           IQS = 
     &     MAX(IQREL,ISREL)*(MAX(IQREL,ISREL)-1)/2 + MIN(IQREL,ISREL)
           IADR = IOFF_H-1+IQS
           XFACTOR = (FI(IADR)-H(IADR))
C          WRITE(6,*) ' XFACTOR = ', XFACTOR
*
           DO IR = 1, NR
           DO IP = 1, NP
            IPQRS = (IS-1)*NP*NQ*NR
     &            + (IR-1)*NP*NQ
     &            + (IQ-1)*NP
     &            + IP
            IPR = (IR-1)*NP + IP
            Y(IPQRS) = Y(IPQRS) + XFACTOR*SCR1(IPR)
          END DO
          END DO
*         ^ End of loops over P,R
         END DO
        END DO
*       ^ End of loop over Q,S
       END IF
*      ^ End P and R have same symmetry
      END IF
*     ^ End if P and R are active
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ===================='
       WRITE(6,*) ' Y block as Y(PQ,RS) '
       WRITE(6,*) ' ===================='
       WRITE(6,*)
       WRITE(6,'(A,4I3)') ' Sym  of p,q,r,s', IPSM,IQSM,IRSM,ISSM
       WRITE(6,'(A,4I3)') ' Type of p,q,r,s', IPTP,IQTP,IRTP,ISTP
       WRITE(6,*)
       CALL WRTMAT(Y,NP*NQ,NR*NS,NP*NQ,NR*NS)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_E2BLK(E2BLK,NP,IPTP,IPSM,NQ,IQTP,IQSM,
     &                         NR,IRTP,IRSM,NS,ISTP,ISSM,
     &                         F, FI, FA, H, 
     &                         SCR1,SCR2,SCR3,SCR4)
*
* Obtain block of orbital-orbital Hessian as
*
* E2(p,q,r,s) = (1-P(p,q))(1-P(r,s)
*                [2D(pr)h(qs)-(F(pr)+F(rp))*delta(qs) + 2Y(pqrs)]
*
* Where 
*
* Y(pqrs) = Sum_(mn:ina+act) [d(pmrn)+d(pmnr)](mq!sn) + d(mnpr)(mn!qs)
*
* Modifications arises from the reduction of Y to inactive/active
* orbitals
* 
* Jeppe Olsen, Jan. 99
*              July 2010: Updated for inactive..
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
*. To access blocks of one-electron integrals
      INCLUDE 'intform.inc'
CNW   DIMENSION F(*), FI(*), FA(*), H(*)
      integer F, FI, FA, H
*. Output 
      DIMENSION E2BLK(NP*NQ*NR*NS)
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Scratch : each SCR matrix should be able to hold largest 4 index matrix
*            of orbitals with given type-symmetry 
*
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*)
*
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' Hessian block E2(pq),(rs) will be calculated for '
       WRITE(6,'(A,2I4)') ' Type and symmetry of P ' ,IPTP, IPSM
       WRITE(6,'(A,2I4)') ' Type and symmetry of Q ' ,IQTP, IQSM
       WRITE(6,'(A,2I4)') ' Type and symmetry of R ' ,IRTP, IRSM
       WRITE(6,'(A,2I4)') ' Type and symmetry of S ' ,ISTP, ISSM
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E2BLK ')
*
      ZERO = 0.0D0
      CALL SETVEC(E2BLK,ZERO,NP*NQ*NR*NS)
*
*. Loop over permutations P(pq),P(rs)
      DO IPQ = 1, 2
       IF(IPQ .EQ. 1 ) THEN
         NPP = NP
         IPPSM = IPSM
         IPPTP = IPTP
         NQQ = NQ
         IQQSM = IQSM
         IQQTP = IQTP 
         XPQ = 1.0D0
       ELSE
         NPP = NQ
         IPPSM = IQSM
         IPPTP = IQTP
         NQQ = NP
         IQQSM = IPSM
         IQQTP = IPTP
         XPQ = -1.0D0
       END IF
       DO IRS = 1, 2
        IF(IRS.EQ.1) THEN
          NRR = NR
          IRRSM = IRSM
          IRRTP = IRTP
          NSS = NS  
          ISSSM = ISSM
          ISSTP = ISTP
          XRS = 1.0D0
        ELSE
          NRR = NS
          IRRSM = ISSM
          IRRTP = ISTP
          NSS = NR  
          ISSSM = IRSM
          ISSTP = IRTP
          XRS = -1.0D0
        END IF
*
* F terms
*
        IF(IQQSM.EQ.ISSSM.AND.IQQTP.EQ.ISSTP) THEN
*. Obtain F(p,r), F(r,p) blocks 
          IH1FORM = 2
          LF = 0
          DO ISM = 1, NSMOB
            LF = LF + NTOOBS(ISM)**2
          END DO
CNW       CALL COPVEC(F,WORK(KINT1),LF)
          call ga_copy_patch('n',F,1,LF,1,1,KINT1,1,LF,1,1)
*. F(pr)
          CALL GETH1(SCR2,IPPSM,IPPTP,IRRSM,IRRTP)
*. F(rp)
          CALL GETH1(SCR3,IRRSM,IRRTP,IPPSM,IPPTP)
*. Cleanup
          IH1FORM = 1
*
          DO IP = 1, NP
          DO IR = 1, NR
          DO IQ = 1, NQ
          DO IS = 1, NS
            IF(IPQ.EQ.1) THEN
              IPP = IP
              IQQ = IQ
            ELSE
              IPP = IQ
              IQQ = IP
            END IF
            IF(IRS.EQ.1) THEN
              IRR = IR 
              ISS = IS
            ELSE
              IRR = IS
              ISS = IR
            END IF
            IF(IQQ.EQ.ISS) THEN
              IPQRS = ((IS-1)*NR+IR-1)*NQ*NP + (IQ-1)*NP+IP  
              E2BLK(IPQRS) = E2BLK(IPQRS)
     &      - XPQ*XRS*(SCR3((IPP-1)*NRR+IRR)+SCR2((IRR-1)*NPP+IPP))
CERR &      - XPQ*XRS*(SCR2((IPP-1)*NRR+IRR)+SCR3((IRR-1)*NPP+IPP))
            END IF
          END DO
          END DO
          END DO
          END DO
*         ^ End of loop  over p,q,r,s
        END IF
*       ^ End if F- terms are active
*
* h-terms
*
*. Nonvanishing if indeces PP and RR corresponds to occupied orbitals
*
        IF(IPPSM.EQ.IRRSM.AND.
     &     ((1.LE.IPPTP.AND.IPPTP.LE.NGAS.AND.
     &       1.LE.IRRTP.AND.IRRTP.LE.NGAS).OR.
     &        IPPTP.EQ.0.AND.IRRTP.EQ.0       ) ) THEN 
*.h(qq,ss)
CNW       CALL COPVEC(H,WORK(KINT1),NINT1)
          call ga_copy_patch('n',H,1,NINT1,1,1,KINT1,1,NINT1,1,1)
          CALL GETH1(SCR2,IQQSM,IQQTP,ISSSM,ISSTP)
*. D(pp,rr)
          IF(IPPTP.EQ.0.AND.IRRTP.EQ.0) THEN
*. Inactive part
            ZERO = 0.0D0
            CALL SETVEC(SCR1,ZERO,NPP*NRR)
            TWO = 2.0D0
            CALL SETDIA(SCR1,TWO,NPP,0)
          ELSE
            CALL GETD1(SCR1,IPPSM,IPPTP,IRRSM,IRRTP,1)
          END IF
          DO IP = 1, NP
          DO IQ = 1, NQ
            IF(IPQ.EQ.1) THEN
              IPP = IP
              IQQ = IQ
            ELSE
              IPP = IQ
              IQQ = IP
            END IF
            DO IR = 1, NR
            DO IS = 1, NS
              IF(IRS.EQ.1) THEN
                IRR = IR
                ISS = IS
              ELSE
                IRR = IS
                ISS = IR
              END IF
              IPQRS = ((IS-1)*NR+IR-1)*NQ*NP + (IQ-1)*NP + IP
              E2BLK(IPQRS) = E2BLK(IPQRS)     
     &       +2.0D0*XPQ*XRS*SCR1((IRR-1)*NPP+IPP)*SCR2((ISS-1)*NQQ+IQQ)
            END DO
            END DO
*           ^ End of loop over r,s
          END DO
          END DO
*         ^ End of loop over p,q
        END IF
*       ^ End if h-terms contributed
*
* Y-terms
*
*. Contributes only if pp and rr are occupied 
        IF(IPPTP.LE.NGAS.AND.IRRTP.LE.NGAS) THEN
*. Construct Y-matrix
          CALL GET_YBLK(SCR1,NPP,IPPTP,IPPSM,NQQ,IQQTP,IQQSM,
     &                       NRR,IRRTP,IRRSM,NSS,ISSTP,ISSSM,
     &                       FI,FA,H,
     &                       SCR2,SCR3,SCR4)
          DO IP = 1, NP
          DO IQ = 1, NQ
            IF(IPQ.EQ.1) THEN
              IPP = IP
              IQQ = IQ
            ELSE
              IPP = IQ
              IQQ = IP
            END IF
            DO IR = 1, NR
            DO IS = 1, NS
              IF(IRS.EQ.1) THEN
                IRR = IR
                ISS = IS
              ELSE
                IRR = IS
                ISS = IR
              END IF
              IPQRS = ((IS-1)*NR+IR-1)*NQ*NP + (IQ-1)*NP + IP
              I2PQRS = ((ISS-1)*NRR+IRR-1)*NPP*NQQ+(IQQ-1)*NPP+IPP
C?            WRITE(6,*) '  Updated with Y, IPQRS, I2PQRS, Y(I2PQRS)',
C?   &        IPQRS, I2PQRS, SCR1(I2PQRS)
              TWO = 2.0D0
              E2BLK(IPQRS)= E2BLK(IPQRS) + TWO*XPQ*XRS*SCR1(I2PQRS)
            END DO
            END DO
*           ^ End of loop over r,s
          END DO
          END DO
*         ^ End of loop over p,q
        END IF
*       ^ End if y-terms contributed
       END DO
*      ^ End of loop over IRS
      END DO
*     ^ End of loop over IPQ
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' Hessian block E2(p,q,r,s) '
        WRITE(6,*) ' =========================='
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' Sym  of p,q,r,s', IPSM,IQSM,IRSM,ISSM
        WRITE(6,'(A,4I3)') ' Type of p,q,r,s', IPTP,IQTP,IRTP,ISTP
        WRITE(6,*) ' Matrix in form E2(PQ,RS)'   
        WRITE(6,*)
        CALL WRTMAT(E2BLK,NP*NQ,NR*NS,NP*NQ,NR*NS)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E2BLK ')
*
      RETURN
      END
      SUBROUTINE REO_4INDMAT(AIN,AUT,IAC,N1IN,N2IN,N3IN,N4IN,
     &                       IN1, IN2, IN3, IN4,
     &                       N1UT,N2UT,N3UT,N4UT)
*
* A 4 dimensional matrix (AIN(I1,I2,I3,I4) is given 
* with dimensions N1IN, N2IN, N3IN, N4IN, respectively
*
* Obtain matrix where the indeces have been reordered,so
* IN1 is old index of new index1..
*
* IAC = 1: Add reordered elements to AUT
* IAC = 2: Copy reordered elements to AUT
*
* Jeppe Olsen; June3, 2013; Lugano, trying to get the orbital E2
*
      INCLUDE 'implicit.inc'
* Input 
      DIMENSION AIN(N1IN,N2IN,N3IN,N4IN)
*. Output
      DIMENSION AUT(*)
*. Local scratch for easing reorder
*
      INTEGER IO_TO_N(4), IN_TO_O(4), IARR(4)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from REO_4INDMAT '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
      END IF
*. New to Old
      IN_TO_O(1) = IN1
      IN_TO_O(2) = IN2
      IN_TO_O(3) = IN3
      IN_TO_O(4) = IN4
*. Old to New
      IO_TO_N(IN_TO_O(1)) = 1
      IO_TO_N(IN_TO_O(2)) = 2
      IO_TO_N(IN_TO_O(3)) = 3
      IO_TO_N(IN_TO_O(4)) = 4
*             
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' IN_TO_O and IO_TO_N: '
        CALL IWRTMA3(IN_TO_O,1,4,1,4)
        WRITE(6,*)
        CALL IWRTMA3(IO_TO_N,1,4,1,4)
      END IF
*
* Dimensions 
*
      IARR(1) = N1IN
      IARR(2) = N2IN
      IARR(3) = N3IN
      IARR(4) = N4IN
      N1UT = IARR(IN_TO_O(1))
      N2UT = IARR(IN_TO_O(2))
      N3UT = IARR(IN_TO_O(3))
      N4UT = IARR(IN_TO_O(4))
*
*. Loop over elements of output matrix
*
      IADDR_UT = 0
      DO I4UT = 1, N4UT
       DO I3UT = 1, N3UT
        DO I2UT = 1, N2UT
         DO I1UT = 1, N1UT
           IARR(1) = I1UT
           IARR(2) = I2UT
           IARR(3) = I3UT
           IARR(4) = I4UT
           I1IN = IARR(IO_TO_N(1))
           I2IN = IARR(IO_TO_N(2))
           I3IN = IARR(IO_TO_N(3))
           I4IN = IARR(IO_TO_N(4))
           IADDR_UT = IADDR_UT + 1
           IF(IAC.EQ.1) THEN
             AUT(IADDR_UT) = 
     &       AUT(IADDR_UT) + AIN(I1N, I2IN, I3IN, I4IN)
           ELSE IF (IAC.EQ.2) THEN
             AUT(IADDR_UT) = AIN(I1N, I2IN, I3IN, I4IN)
           END IF
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 
     &  ' Input and output 4 index matrices as MAT(12,34)'
        WRITE(6,*)
        N12 = N1IN*N2IN
        N34 = N3IN*N4IN
        CALL WRTMAT(AIN, N12, N34,N12, N34)
        WRITE(6,*)
        N12 = N1UT*N2UT
        N34 = N3UT*N4UT
        CALL WRTMAT(AUT, N12, N34,N12, N34)
      END IF
*
      RETURN
      END
      SUBROUTINE TRNS_IND_4INDMAT(AIN,AUT,N1,N2,N3,N4,IINDT,C)
*
* Transform one index of a 4-dimensional matrix with matric C
*
* A  transformation of an index by a matrix is defined as 
* new(i) = sum_i' old(i') c_i'i
*
*. Jeppe Olsen; June 3, 2013; Lugonao for the non-ort Hessian
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION AIN(N1*N2*N3*N4), C(*)
*. Output
      DIMENSION AUT(N1*N2*N3*N4)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TRNS_IND_4INDMAT '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Index to be transformed: IINDT'
      END IF
*. To eliminate compiler warning  
      NC = 1
*   
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      IF(IINDT.EQ.1) THEN
* AUT(i,j,k,l) = sum_i' AIN(i',j,k,l) C_i'i
        NRUT = N1
        NCUT = N2*N3*N4
        CALL MATML7(AUT,C,AIN,
     &       NRUT,NCUT,NRUT,NRUT,NRUT,NCUT,
     &       FACTORC, FACTORAB, 1)
      ELSE IF (IINDT.EQ.2) THEN
* AUT(i,j,k,l) = sum_j' AIN(i,j',k,l) C_j'j
        NR = N1
        NC = N2
        DO IKL = 1, N3*N4
          IB_KL = (IKL-1)*N1*N2 + 1
          CALL MATML7(AUT(IB_KL),AIN(IB_KL),C,
     &                NR,NC,NR,NC,NC,NC,
     &                FACTORC, FACTORAB,0)
        END DO
      ELSE IF (IINDT.EQ.3) THEN
* AUT(i,j,k,l) = sum_k' AIN(i,j,k',l) C_k'k
        NR = N1*N2
        NC = N3
        DO IL = 1, N4
          IB_L = (IL-1)*N1*N2*N3 + 1
          CALL MATML7(AUT(IB_L),AIN(IB_L),C,
     &                NR,NC,NR,NC,NC,NC,
     &                FACTORC, FACTORAB,0)
        END DO
      ELSE IF (IINDT.EQ.4) THEN
* AUT(i,j,k,l) = sum_l' AIN(i,j,k,l') C_l'l
        NR = N1*N2*N3
        NC = N3
        CALL MATML7(AUT,AIN,C,
     &              NR,NC,NR,NC,NC,NC,
     &                FACTORC, FACTORAB,0)
      ELSE
        WRITE(6,*) ' Illegal value of IINDT = ', IINDT 
        STOP '  Illegal value of IINDT'
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 
     &  ' Input as (12, 34)  and transformation matrix'
        N12 = N1*N2
        N34 = N3*N4
        CALL WRTMAT(AIN,N12,N34,N12,N34)
        WRITE(6,*)
        CALL WRTMAT(C,NC,NC,NC,NC)
        WRITE(6,*) ' Output matrix as (12,34) '
        CALL WRTMAT(AUT,N12,N34,N12,N34)
      END IF
*
      RETURN
      END
      SUBROUTINE OUT_PROD_2MAT(AB,A,B,NAB1,NAB2,NAB3,NAB4,
     &           NA1,NA2,NB1,NB2,
     &           I1MAT,I2MAT,I3MAT,I4MAT,I1IND,I2IND,I3IND,I4IND)
*
* Obtain outer product of two two-dimensional matrices and
* store
*
*. Jeppe Olsen; June 3, 2013; Lugano for E2 non-ort 
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(NA1,NA2), B(NB1,NB2)
*. Output
      DIMENSION AB(NAB1*NAB2*NAB3*NAB4)
*. Scratch to ease indexing
      DIMENSION IMAT_FOR_4IND(4),I4IND_FOR_2IND(2,2), IARR(4)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from OUT_PROD_2MAT '
       WRITE(6,*) ' ========================'
       WRITE(6,*) 
      END IF
*
      IMAT_FOR_4IND(1) = I1MAT
      IMAT_FOR_4IND(2) = I2MAT
      IMAT_FOR_4IND(3) = I3MAT
      IMAT_FOR_4IND(4) = I4MAT
*
      I4IND_FOR_2IND(I1IND,I1MAT) = 1
      I4IND_FOR_2IND(I2IND,I2MAT) = 2
      I4IND_FOR_2IND(I3IND,I3MAT) = 3
      I4IND_FOR_2IND(I4IND,I4MAT) = 4
*
      I4ELMNT = 0
      DO I4 = 1, NAB4
       DO I3 = 1, NAB3
        DO I2 = 1, NAB2
         DO I1 = 1, NAB1
           I4ELMNT = I4ELMNT + 1
           IARR(1) = I1
           IARR(2) = I2
           IARR(3) = I3
           IARR(4) = I4
           IA1 = IARR(I4IND_FOR_2IND(1,1))
           IA2 = IARR(I4IND_FOR_2IND(1,2))
           IB1 = IARR(I4IND_FOR_2IND(2,1))
           IB2 = IARR(I4IND_FOR_2IND(2,2))
           AB(IELMNT) = A(IA1,IA2)*B(IB1,IB2)
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output AB matrix as AB(12,34) '
        N12 = NAB1*NAB2
        N34 = NAB3*NAB4
        CALL WRTMAT(AB,N12,N34,N12,N34)
      END IF
*
      RETURN
      END
c $Id$
