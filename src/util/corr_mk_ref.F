C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     Corr_Mk_Ref -- Create reference for correlated calculation
C
C REVISION
C     $Id: corr_mk_ref.F,v 1.2 1998-05-11 03:09:01 d3e129 Exp $
C
C SYNOPSIS
      Subroutine Corr_Mk_Ref(RTDB, Mod, Geom, NAO, RefSpin, RefNMO,
     $     G_RefVec, RChunk, CChunk, MxSpin, H_RefOcc, H_RefEval,
     $     NFrzOcc, NFrzVir, KeepFrz)
C
C ARGUMENTS
      Implicit NONE
      Integer RTDB              ! [in]  Database handle
      Character*(*) Mod         ! [in]  Module name
      Integer Geom              ! [in]  Geometry handle
      Integer NAO               ! [out] Number of AO basis functions
      Integer RefSpin           ! [out] Number of spin cases
      Integer MxSpin            ! [in]  Maximum number of spin cases expected
      Integer RefNMO(MxSpin)    ! [out] Number of MOs in each spin case
      Integer G_RefVec(MxSpin)  ! [out] GA handles of reference eigenvectors
      Integer RChunk            ! [in]  Desired row chunksize for vector GA 
                                !       (0 = NAO)
      Integer CChunk            ! [in]  Desired col chunksize for vector GA 
                                !       (0 = NMO)
      Integer H_RefOcc          ! [out] MA handle of reference occupation 
                                !       vector (heap)
      Integer H_RefEVal         ! [out] MA handle of reference eigenvalue 
                                !       vector (heap)
      Integer NFrzOcc           ! [out] Number of frozen occupied orbitals
      Integer NFrzVir           ! [out] Number of frozen virtual orbitals
      Logical KeepFrz           ! [in]  Keep frozen orbitals in reference 
                                !       returned
C
C INCLUDE FILES
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "util.fh"
#include "stdio.fh"
C
C EXTERNAL ROUTINES
      Logical  MOVecs_Read_Header, MOVecs_Read
      External MOVecs_Read_Header, MOVecs_Read
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Character*(*) VecDef
      Parameter( VecDef = 'movecs')
C
      Integer MyMxSpin
      Parameter (MyMxSpin = 2) 
C 
C LOCAL VARIABLES
C
      Integer MLen, MyRChunk, MyCChunk, Spin, DType, Lo, Hi, Me, I
      Character*(80) VecFil, VecTitle, VecBas, VecType
      Character*(80) Key, Key2, Type, Date
      Integer VecSpin, VecNMO(MyMxSpin), G_Vec(MyMxSpin)
      Integer SpinMap(MyMxSpin)
      Double Precision TotOcc(MyMxSpin)
      Integer H_Occ, I_Occ, H_EVal, I_EVal
c$$$      Integer H_Occ0, I_Occ0, H_EVal0, I_EVal0
      Integer I_RefOcc, I_RefEVal
      Logical Status, PrInfo
C
      MLen = Inp_StrLen(Mod)
      PrInfo = Util_Print('ref info', Print_Low)
      Me = GA_NodeID()
C
C     Simple sanity checks
C
      If ( MxSpin .gt. MyMxSpin) Call ErrQuit(
     $     'corr_mk_ref: caller expects too many spin cases', MxSpin)
C
C     Precedence for source of vectors:
C     1) VECTORS INPUT for this module
C     2) VECTORS OUTPUT for SCF module
C     3) <job_prefix>.movecs
C
      Write( Key, '(a,'':'',a)' ) Mod(1:MLen), 'input vectors'
      If ( RTDB_CGet(RTDB, Key, 1, VecFil) ) then
      ElseIf ( RTDB_CGet(RTDB, 'scf:output vectors', 1, VecFil) ) then
      Else
         Call Util_File_Name(VecDef, .False., .False., VecFil)
      EndIf
C
C     Get basic dimension info from VecFil
C
      If ( .NOT.  MOVecs_Read_Header(VecFil, VecTitle, VecBas, 
     $     VecType, NAO, VecSpin, VecNMO, MyMxSpin) )
     $   Call ErrQuit('corr_mk_ref: can''t read vectors (header)', 0)
C
C     Allocate space for vectors, eigenvalues, and occupation
C
      Do Spin = 1, VecSpin
         MyRChunk = RChunk
         MyCChunk = CChunk
         If (MyRChunk .eq. 0) MyRChunk = NAO
         If (MyCChunk .eq. 0) MyCChunk = VecNMO(Spin)
C
         If ( .NOT. GA_Create(MT_Dbl, NAO, VecNMO(Spin),
     $        'SCF Vectors', MyRChunk, MyCChunk, G_Vec(Spin) ) )
     $        Call ErrQuit(
     $        'corr_mk_ref: failed to create GA for vectors',
     $        NAO*VecNMO(Spin) )
      EndDo
C
C     For eigenvalues and occupation vector, each spin requires NAO
C     elements (not VecNMO(Spin)), and the entire vector (both spins)
C     must be contiguous (for MOVecs_Swap).
C     
      If ( .NOT. MA_Push_Get(MT_Dbl, NAO*VecSpin,
     $     'SCF Occupation', H_Occ, I_Occ) ) Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for occupation',
     $     NAO*VecSpin)
      If ( .NOT. MA_Push_Get(MT_Dbl, NAO*VecSpin,
     $     'SCF Eigenvalues', H_EVal, I_EVal) ) Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for eigenvalues',
     $     NAO*VecSpin)
c$$$C
c$$$C     Another set so we can preserve the originals for a while
c$$$C
c$$$      If ( .NOT. MA_Push_Get(MT_Dbl, NAO*VecSpin,
c$$$     $     'SCF Occupation', H_Occ0, I_Occ0) ) Call ErrQuit(
c$$$     $     'corr_mk_ref: failed to allocate for occupation',
c$$$     $     NAO*VecSpin)
c$$$      If ( .NOT. MA_Push_Get(MT_Dbl, NAO*VecSpin,
c$$$     $     'SCF Eigenvalues', H_EVal0, I_EVal0) ) Call ErrQuit(
c$$$     $     'corr_mk_ref: failed to allocate for eigenvalues',
c$$$     $     NAO*VecSpin)
C
C     Read in everything
C
      Do Spin = 1, VecSpin
         If ( .NOT. MOVecs_Read(VecFil, Spin,
     $        Dbl_MB(I_Occ + NAO*(Spin-1)),
     $        Dbl_MB(I_EVal+ NAO*(Spin-1)), G_Vec(Spin) ) )
     $        Call ErrQuit(
     $        'corr_mk_ref: can''t read vectors', 0)
      EndDo
c$$$C
c$$$C     Preserve original eigenvalues & vectors
c$$$C
c$$$      Call DCopy(NAO*VecSpin, Dbl_MB(I_Occ),  1, Dbl_MB(I_Occ0),  1)
c$$$      Call DCopy(NAO*VecSpin, Dbl_MB(I_EVal), 1, Dbl_MB(I_EVal0), 1)
C
      Call GA_Sync
      Call Inp_UCase(VecType)   ! Presumptive, but reasonable treatment
      If (PrInfo .AND. Me .eq. 0) then
         Write (LuOut, '(1X, A, 1X, A)')
     $        VecType(1:Inp_StrLen(VecType)), 'Reference Information:'
         Write (LuOut, 9000) 'File:', VecFil(1:Inp_StrLen(VecFil))
         Write (LuOut, 9000) 'Title:', VecTitle(1:Inp_StrLen(VecTitle))
C
         If ( VecSpin .eq. 1) then
            Write (LuOut, 9020) 'Size (AOxMO):', NAO, VecNMO(1)
         ElseIf ( VecSpin .eq. 2) then
            Write (LuOut, 9025) 'Size (AOxMO):',
     $           NAO, (VecNMO(i), i=1, VecSpin)
         Else                   ! Just in case
            Write (LuOut, 9026) 'Size (AOxMO):',
     $           NAO, (VecNMO(i), i=1, VecSpin)
         EndIf
         Do Spin = 1, VecSpin
            TotOcc(Spin) = 0.0d0
            Do I = 1, VecNMO(Spin)
               TotOcc(Spin) = TotOcc(Spin) +
     $              Dbl_MB(I_Occ + NAO*(Spin-1) + I-1) 
            EndDo
         EndDo
         If ( VecSpin .eq. 1) then
            Write (LuOut, 9030) 'Total Occ:', TotOcc(1)
         ElseIf (VecSpin .eq. 2) then
            Write (LuOut, 9035) 'Total Occ:',
     $           (TotOcc(i), i = 1, VecSpin)
         Else
            Write (LuOut, 9036) 'Total Occ:',
     $           (TotOcc(i), i = 1, VecSpin)
         EndIf
      EndIf
 9000 Format(4X, A13, 1X, '''', A, '''')
 9010 Format(4X, A13, 10(1X, I7))
 9020 Format(4X, A13, I7, ' x', I7)
 9025 Format(4X, A13, I7, ' x ', '(', I7, ',', I7, ')')
 9026 Format(4X, A13, I7, ' x ', 8(I7, ','))
 9030 Format(4X, A13, 12x,  F7.2)
 9035 Format(4X, A13, 14x, F7.2, ',', F7.2)
 9036 Format(4X, A13, 14x, 8(F7.2, ','))
C
      Call GA_Sync
C
C     Process VECTORS SWAP directive.  Doing SWAP before FREEZE
C     allows arbitrary orbitals to be frozen by swapping them to 
C     the bottom or top of the vector.
C
      If ( VecSpin .eq. 2) then
         Type = 'UHF'
      Else
         Type = 'RHF'
      EndIf
      Call MOVecs_Swap(RTDB, Mod, Type, G_Vec, Dbl_MB(I_Occ),
     $     Dbl_MB(I_EVal) )
C
C     Now handle FREEZE input directive.  Simply freeze the lowest or
C     highest orbitals.  Use SWAP to change around orbitals and you
C     have complete generality.
C
C     First for frozen core
C
      Write (Key, '(a,'':'',a,'':'',a)') Mod(1:MLen), 'frozen core',
     $     'freeze by atoms'
      Write (Key2, '(a,'':'',a)') Mod(1:MLen), 'frozen core'
      If ( RTDB_Get(RTDB, Key(1:Inp_StrLen(Key)), MT_Log, 1,
     $     Status) ) Then
         If (.NOT. Geom_Num_Core(Geom, NFrzOcc) ) Call ErrQuit(
     $        'corr_mk_ref: problem with geom_num_core', 0)
      ElseIf ( RTDB_Get(RTDB, Key2(1:Inp_StrLen(Key2)), MT_Int, 1,
     $        NFrzOcc) ) Then
      Else
         NFrzOcc = 0
      EndIf
C
C     Then for frozen virtuals
C
      Write (Key, '(a,'':'',a,'':'',a)') Mod(1:MLen), 'frozen virtual',
     $     'freeze by atoms'
      Write (Key2, '(a,'':'',a)') Mod(1:MLen), 'frozen virtual'
      If ( RTDB_Get(RTDB, Key(1:Inp_StrLen(Key)), MT_Log, 1,
     $     Status) ) Then
         If (.NOT. Geom_Num_Core(Geom, NFrzVir) ) Call ErrQuit(
     $        'corr_mk_ref: problem with geom_num_core', 0)
      ElseIf ( RTDB_Get(RTDB, Key2(1:Inp_StrLen(Key2)), MT_Int, 1,
     $        NFrzVir) ) Then
      Else
         NFrzVir = 0
      EndIf
C
C     Now handle spin case mapping if required
C
      Write( Key, '(a,'':'',a)' ) Mod, 'reference spin mapping'
      If ( RTDB_Get_Info( RTDB, Key(1:Inp_StrLen(Key)), DType, RefSpin,
     $     Date) ) then 
         If ( RefSpin .gt. MyMxSpin) Call ErrQuit(
     $      'corr_mk_ref: too many elements in spin mapping',
     $      RefSpin)
C
         If ( .NOT. RTDB_Get( RTDB, Key(1:Inp_StrLen(Key)),
     $      MT_Int, RefSpin, SpinMap) ) Call ErrQuit(
     $      'corr_mk_ref: problems reading runtime database', 0)
C
      Else                      ! NO RTDB ENTRY -- defaults
         RefSpin = VecSpin
         Do Spin = 1, RefSpin
            SpinMap(Spin) = Spin ! 1-1 mapping
         EndDo
      EndIf
C
C     Sanity check number of spin cases.  Don't do this earlier to allow
C     for using the spin mapping to give choose a subset of spin cases
C     to give to caller.  For example, feed on spin case of a UHF to
C     caller expecting RHF reference.  Potentially useful for debugging
C     and other unusual things.
C
      If ( RefSpin .gt. MxSpin ) Call ErrQuit(
     $     'corr_mk_ref: too many spin cases for caller!', RefSpin)
C
C     Determine final sizes and allocate the output arrays
C
      Do Spin = 1, RefSpin
         If ( SpinMap(Spin) .lt. 1 .OR. SpinMap(Spin) .gt. VecSpin)
     $        Call ErrQuit('corr_mk_ref: bad spin map entry',
     $        SpinMap(Spin))
C
         RefNMO(Spin) = VecNMO( SpinMap(Spin) )
         If ( .NOT. KeepFrz )
     $        RefNMO(Spin) = RefNMO(Spin) - NFrzOcc - NFrzVir
C
         MyRChunk = RChunk
         MyCChunk = CChunk
         If (MyRChunk .eq. 0) MyRChunk = NAO
         If (MyCChunk .eq. 0) MyCChunk = RefNMO(Spin)
C
         If ( .NOT. GA_Create(MT_Dbl, NAO, RefNMO(Spin),
     $        'Reference Vectors', MyRChunk, MyCChunk, G_RefVec(Spin)))
     $        Call ErrQuit(
     $        'corr_mk_ref: failed to create GA for vectors',
     $        NAO*RefNMO(Spin) )
      EndDo
C     
      If ( .NOT. MA_Alloc_Get(MT_Dbl, NAO*RefSpin,
     $     'Reference Occupation', H_RefOcc, I_RefOcc) ) Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for occupation',
     $     NAO*RefSpin)
      If ( .NOT. MA_Alloc_Get(MT_Dbl, NAO*RefSpin,
     $     'Reference Eigenvalues', H_RefEVal, I_RefEVal) )
     $     Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for eigenvalues',
     $     NAO*RefSpin)
C
C     Now copy everything into its final destination
C
      Do Spin = 1, RefSpin
         Lo = 1
         Hi = VecNMO(SpinMap(Spin))
         If ( .NOT. KeepFrz ) then
            Lo = Lo + NFrzOcc
            Hi = Hi - NFrzVir
         EndIf
         Call GA_Copy_Patch('No transpose',
     $        G_Vec(SpinMap(Spin)), 1, NAO, Lo, Hi,
     $        G_RefVec(Spin),       1, NAO, 1,  RefNMO(Spin) )
C
         Call dcopy(RefNMO(Spin),
     $        Dbl_MB(I_Occ    + NAO*(SpinMap(Spin)-1) + Lo - 1), 1,
     $        Dbl_MB(I_RefOcc + NAO*(Spin-1)),                   1)
C
         Call dcopy(RefNMO(Spin),
     $        Dbl_MB(I_EVal    + NAO*(SpinMap(Spin)-1) + Lo - 1), 1,
     $        Dbl_MB(I_RefEVal + NAO*(Spin-1)),                   1)
      EndDo
C
c$$$      Write (LuOut, '(/(16F5.2))')
c$$$     $     (Dbl_MB(I_Occ + i), i=0, (VecSpin*NAO)-1)
c$$$      Write (LuOut, '(/(16F5.2))')
c$$$     $     (Dbl_MB(I_RefOcc + i), i=0, (RefSpin*NAO)-1)
C
C     Free up memory
C
c$$$      If ( .NOT. MA_Pop_Stack(H_Eval0) ) Call ErrQuit(
c$$$     $     'corr_mk_ref: failed to pop stack', 0)
c$$$      If ( .NOT. MA_Pop_Stack(H_Occ0) ) Call ErrQuit(
c$$$     $     'corr_mk_ref: failed to pop stack', 0)
C
      If ( .NOT. MA_Pop_Stack(H_Eval) ) Call ErrQuit(
     $     'corr_mk_ref: failed to pop stack', 0)
      If ( .NOT. MA_Pop_Stack(H_Occ) ) Call ErrQuit(
     $     'corr_mk_ref: failed to pop stack', 0)
      Do Spin = VecSpin, 1, -1
         If ( .NOT. GA_Destroy(G_Vec(Spin) ) ) Call ErrQuit(
     $        'corr_mk_ref: failed to destroy ga', 0)
      EndDo
C
      Return
      End
