C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     Corr_Mk_Ref -- Create reference for correlated calculation
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function Corr_Mk_Ref(RTDB, Mod, Geom, Basis, DoSCF, 
     $     NAO, RefSpin, RefNMO, G_RefVec, RChunk, CChunk, MxSpin,
     $     H_RefOcc, H_RefEval, NFrzOcc, NFrzVir, KeepFrz)
C
C ARGUMENTS
      Implicit NONE
#include "errquit.fh"
      Integer RTDB              ! [in]  Database handle
      Character*(*) Mod         ! [in]  Module name
      Integer Geom              ! [in]  Geometry handle
      Integer Basis             ! [in]  Basis handle
      Logical DoSCF             ! [in]  Invoke SCF if needed
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
      logical dftmp2            ! use DFT reference orbitals and energy
C
C DESCRIPTION
C     Prepares reference wavefunction for post-SCF calculations.
C     Raw reference info is read from (in order of preference)
C     1) VECTORS INPUT for current module
C     2) VECTORS OUTPUT for SCF module
C     3) <job_prefix>.movecs
C     Swap and freeze directives are applied to produce the reference
C     returned to the caller.  Note that since swaps are processed 
C     before freezes, it is possible to freeze arbitrary orbitals.
C
C     If KeepFrz is .TRUE., frozen occupied and virtual orbitals are
C     eliminated from the occupation, eigenvalues, and eigenvectors
C     returned to the user.
C
C     NFrzOcc and NFrzVir always return the number of orbitals
C     which were frozen regardless of the value of KeepFrz.
C
C     If DoSCF is .TRUE. the SCF module will be invoked under these
C     conditions:
C     1) SCF reference does not exist
C     2) SCF reference exists, but does not satisfy the requirements
C        of movecs_converged (basis & geometry checksums, scftype, and
C        date of movecs file must match current RTDB and active
C        calculation)
C
C     Note that this mechanism is imperfect:
C     1) Input decks for which the SCF module inputs don't match those
C        for the correlated module (primarily w.r.t. movecs filenames)
C        will probably run the SCF then fail.
C     2) No check is make of the convergence threshold to insure it is
C        appropriate for a correlated calculation.  Some modules which
C        do NOT use corr_mk_ref check the convergence threshold, and may
C        even invoke the SCF if they consider the threshold inadequate.
C
C     When orbitals are frozen, a simple test is made to try to warn
C     users if they are splitting a degeneracy.  It is possible to
C     fool this check by swaping orbitals, but this is simple and
C     useful as a sanity check.
C
C RTDB ENTRIES USED
C <Mod>:input vectors, scf:output vectors
C     Used to figure out what movecs file to read
C
C Entries understood by movecs_swap
C     Reorder vectors prior to freezing
C
C <Mod>:frozen core, <Mod>:frozen virtual
C     Determines which orbitals to freeze
C
C <Mod>:reference spin mapping
C     Allows switching or cloning of references by spin case
C
C PRINT CONTROLS
C 'reference info' (low)
C     Print basic information about the reference: filename, title
C     string, dimensions, total occupation.
C 'raw reference'  (debug)
C     Print raw reference information: occupation, eigenvalues,
c     eigenvectors
C 'reference'      (debug)
C     Print final reference information: occupation, eigenvalues,
c     eigenvectors
C
C INCLUDE FILES
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "util.fh"
#include "stdio.fh"
#include "msgids.fh"
C
C EXTERNAL ROUTINES
      Logical  MOVecs_Read_Header, MOVecs_Read, MOVecs_Converged, SCF
      External MOVecs_Read_Header, MOVecs_Read, MOVecs_Converged, SCF
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C CONSTANTS
      Character*(*) VecDef
      Parameter( VecDef = 'movecs')
C
      Character*(*) Routine
      Parameter ( Routine = 'Corr_Mk_Ref: ')
C
      Integer MyMxSpin
      Parameter (MyMxSpin = 2)
C
C     Warn user when freeze directives split degeneracies
      Double Precision Thresh_Degen
      Parameter (Thresh_Degen = 1.0d-4)
C 
C LOCAL VARIABLES
C
      Integer MLen, MyRChunk, MyCChunk, Spin, DType, Lo, Hi, Me, I
      integer mitob1
      Character*(256) VecFil, VecTitle, VecBas, VecType
      Character*(80) Key, Key2, Type, Date
      Integer VecSpin, VecNMO(MyMxSpin), G_Vec(MyMxSpin), MxVecNMO,
     $     MxRefNMO
      Integer SpinMap(MyMxSpin)
      Double Precision TotOcc(MyMxSpin), Delta
      Integer H_Occ, I_Occ, H_EVal, I_EVal
c$$$      Integer H_Occ0, I_Occ0, H_EVal0, I_EVal0
      Integer I_RefOcc, I_RefEVal
      Logical Status, PrInfo, Header
C
C STATEMENT FUNCTIONS
      Character*1 SpinItoA
      SpinItoA(I) = Char( I + IChar('a') - 1 )
C
      MLen = Inp_StrLen(Mod)
      Me = GA_NodeID()
      PrInfo = Util_Print('reference info', Print_Low) .AND. Me .eq. 0
c
c     Are we doing a double-hybrid?
c
      if (.not. rtdb_get(rtdb,'dft:dftmp2',mt_log,1,dftmp2))
     $     dftmp2 = .false.
C
C     Simple sanity checks
C
      If ( MxSpin .gt. MyMxSpin) Call ErrQuit(
     $     'corr_mk_ref: caller expects too many spin cases', MxSpin,
     &       INPUT_ERR)
C
C     ****************************
C     * Figure out input vectors *
C     ****************************
C     Precedence for source of vectors:
C     1) VECTORS INPUT for this module
C     2) VECTORS OUTPUT for SCF module
C     3) <job_prefix>.movecs
C
      if (dftmp2) then
        if (.not. rtdb_cget(rtdb,'dft:output vectors',1,VecFil))
     $       call errquit('rimp2: no DFT MO vectors',0, RTDB_ERR)
      else
        Write( Key, '(a,'':'',a)' ) Mod(1:MLen), 'input vectors'
        If ( RTDB_CGet(RTDB, Key, 1, VecFil) ) then
        ElseIf ( RTDB_CGet(RTDB, 'scf:output vectors', 1, VecFil) ) then
        Else
           Call Util_File_Name(VecDef, .False., .False., VecFil)
        EndIf
      endif
C
C     ***************************************
C     * Try to insure vectors are converged *
C     ***************************************
C
      If (Me .eq. 0) Inquire( File=VecFil, Exist=Status)
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      Call GA_Brdcst(Msg_CMR_File, Status, MIToB1, 0)
C
      If ( Status )
     $     Status = MOVecs_Converged( RTDB, Basis, 'scf', VecFil)
C
      If ( .NOT. Status .AND. PrInfo ) Write (LUOut, 8910) Routine,
     $     VecFil(1:Inp_StrLen(VecFil) )
 8910 Format(1X, A, 'MO vectors ''', A, ''' don''t exist or aren''t ',
     $     'converged.')
C
      If ( DoSCF .AND. .NOT. Status ) then
         If ( PrInfo ) Write (LuOut, 8920) Routine
         If ( .NOT. SCF( RTDB ) ) Then
            Corr_Mk_Ref = .FALSE.
            if(Prinfo) Call Util_Flush(LuOut)
            Return
         EndIf
      EndIf
 8920 Format(1X, A, 'Invoking SCF to produce converged vectors.')
C
C     Get basic dimension info from VecFil
C
      If ( .NOT.  MOVecs_Read_Header(VecFil, VecTitle, VecBas, 
     $     VecType, NAO, VecSpin, VecNMO, MyMxSpin) )
     $   Call ErrQuit('corr_mk_ref: can''t read vectors (header)', 0,
     &       DISK_ERR)
C
C     Allocate space for vectors, eigenvalues, and occupation
C
      MxVecNMO = 0
      Do Spin = 1, VecSpin
         MxVecNMO = Max( MxVecNMO, VecNMO(Spin) )
C
         MyRChunk = RChunk
         MyCChunk = CChunk
         If (MyRChunk .eq. 0) MyRChunk = NAO
         If (MyCChunk .eq. 0) MyCChunk = VecNMO(Spin)
C
         Key = 'Raw SCF Vectors spin ' // SpinItoA(Spin)
         If ( .NOT. GA_Create(MT_Dbl, NAO, VecNMO(Spin),
     $        Key(:Inp_StrLen(Key)), MyRChunk, MyCChunk, G_Vec(Spin) ) )
     $        Call ErrQuit(
     $        'corr_mk_ref: failed to create GA for vectors',
     $        NAO*VecNMO(Spin) , GA_ERR)
      EndDo
C
C     For eigenvalues and occupation vector, each spin requires NAO
C     elements (not VecNMO(Spin)), and the entire vector (both spins)
C     must be contiguous (for MOVecs_Swap).
C     
      If ( .NOT. MA_Push_Get(MT_Dbl, NAO*VecSpin,
     $     'SCF Occupation', H_Occ, I_Occ) ) Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for occupation',
     $     NAO*VecSpin, MA_ERR)
      If ( .NOT. MA_Push_Get(MT_Dbl, NAO*VecSpin,
     $     'SCF Eigenvalues', H_EVal, I_EVal) ) Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for eigenvalues',
     $     NAO*VecSpin, MA_ERR)
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
     $        'corr_mk_ref: can''t read vectors', 0, DISK_ERR)
      EndDo
c$$$C
c$$$C     Preserve original eigenvalues & vectors
c$$$C
c$$$      Call DCopy(NAO*VecSpin, Dbl_MB(I_Occ),  1, Dbl_MB(I_Occ0),  1)
c$$$      Call DCopy(NAO*VecSpin, Dbl_MB(I_EVal), 1, Dbl_MB(I_EVal0), 1)
C
      Call GA_Sync
      Call Inp_UCase(VecType)   ! Presumptive, but reasonable treatment
      If (PrInfo ) then
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
C     Before we mess up the raw reference, print it for the user
C
      If ( Util_Print('raw reference', Print_Debug) ) then
         If ( Me .eq. 0) then
            Call DGEWR(MxVecNMO, VecSpin, Dbl_MB(I_Occ), NAO, 0, 0,
     $           LuOut, 'Raw Occupation Vectors', 80, 'F6.3')
            Call DGEWR(MxVecNMO, VecSpin, Dbl_MB(I_Eval), NAO, 0, 0,
     $           LuOut, 'Raw Eigenvalues', 80, 'F12.6')
         EndIf
         Do Spin = 1, VecSpin
            Call GA_Print( g_Vec(Spin) )
         EndDo
      EndIf
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
C      Write (Key, '(a,'':'',a,'':'',a)') Mod(1:MLen), 'frozen core',
C     $     'freeze by atoms'
C      Write (Key2, '(a,'':'',a)') Mod(1:MLen), 'frozen core'
C      If ( RTDB_Get(RTDB, Key(1:Inp_StrLen(Key)), MT_Log, 1,
C     $     Status) ) Then
C         If (.NOT. Geom_Num_Core(rtdb,Geom, mod,NFrzOcc) ) 
C     .       Call ErrQuit(
C     $        'corr_mk_ref: problem with geom_num_core', 0)
C      ElseIf ( RTDB_Get(RTDB, Key2(1:Inp_StrLen(Key2)), MT_Int, 1,
C     $        NFrzOcc) ) Then
C      Else
C         NFrzOcc = 0
C      EndIf

       if(.not. geom_freeze(rtdb, geom, Mod, NFrzOcc)) then
          call errquit('corr_mk_ref: problem with geom_freeze',0,
     &       GEOM_ERR)
       endif

C
C     Then for frozen virtuals
C
      Write (Key, '(a,'':'',a,'':'',a)') Mod(1:MLen), 'frozen virtual',
     $     'freeze by atoms'
      Write (Key2, '(a,'':'',a)') Mod(1:MLen), 'frozen virtual'
      If ( RTDB_Get(RTDB, Key(1:Inp_StrLen(Key)), MT_Log, 1,
     $     Status) ) Then
         If (.NOT. Geom_Num_Core(rtdb,Geom,mod,NFrzVir) ) 
     C        Call ErrQuit(
     $        'corr_mk_ref: problem with geom_num_core', 0, GEOM_ERR)
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
     $      RefSpin, INPUT_ERR)
C
         If ( .NOT. RTDB_Get( RTDB, Key(1:Inp_StrLen(Key)),
     $      MT_Int, RefSpin, SpinMap) ) Call ErrQuit(
     $      'corr_mk_ref: problems reading runtime database', 0,
     &       RTDB_ERR)
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
     $     'corr_mk_ref: too many spin cases for caller!', RefSpin,
     &       INPUT_ERR)
C
C     Determine final sizes and allocate the output arrays
C
      MxRefNMO = 0
      Do Spin = 1, RefSpin
         If ( SpinMap(Spin) .lt. 1 .OR. SpinMap(Spin) .gt. VecSpin)
     $        Call ErrQuit('corr_mk_ref: bad spin map entry',
     $        SpinMap(Spin), INPUT_ERR)
C
         RefNMO(Spin) = VecNMO( SpinMap(Spin) )
         If ( .NOT. KeepFrz )
     $        RefNMO(Spin) = RefNMO(Spin) - NFrzOcc - NFrzVir
         MxRefNMO = Max( MxRefNMO, RefNMO(Spin) )
C
         MyRChunk = RChunk
         MyCChunk = CChunk
         If (MyRChunk .eq. 0) MyRChunk = NAO
         If (MyCChunk .eq. 0) MyCChunk = RefNMO(Spin)
C
         Key = 'Reference Vectors spin ' // SpinItoa(Spin)
         If ( .NOT. GA_Create(MT_Dbl, NAO, RefNMO(Spin),
     $        Key(:Inp_StrLen(Key)), MyRChunk, MyCChunk,
     $        G_RefVec(Spin))) Call ErrQuit(
     $        'corr_mk_ref: failed to create GA for vectors',
     $        NAO*RefNMO(Spin) , GA_ERR)
      EndDo
C     
      If ( .NOT. MA_Alloc_Get(MT_Dbl, NAO*RefSpin,
     $     'Reference Occupation', H_RefOcc, I_RefOcc) ) Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for occupation',
     $     NAO*RefSpin, MA_ERR)
      If ( .NOT. MA_Alloc_Get(MT_Dbl, NAO*RefSpin,
     $     'Reference Eigenvalues', H_RefEVal, I_RefEVal) )
     $     Call ErrQuit(
     $     'corr_mk_ref: failed to allocate for eigenvalues',
     $     NAO*RefSpin, MA_ERR)
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
C     Check for broken degeneracies in frozen/active split, which can
C     cause denominators to blow up. This is a simplistic check, but 
C     it should catch most cases.
C
      Header = .FALSE.
      Do Spin = 1, RefSpin
         If ( NFrzOcc .gt. 0 ) then
            Delta = Abs ( Dbl_MB(I_EVal+NAO*(SpinMap(Spin)-1) + NFrzOcc)
     $           - Dbl_MB(I_EVal+NAO*(SpinMap(Spin)-1) + NFrzOcc - 1) )
            If ( Delta .lt. Thresh_Degen ) then
               If ( .NOT. Header ) Write (LuOut, 8950)
               Header = .TRUE.
               Write (LuOut, 8960) 'occupied', Spin, Delta
            EndIf
         EndIF
         If ( NFrzVir .gt. 0 ) then
            Delta = Abs( Dbl_MB(I_EVal+NAO*(SpinMap(Spin)-1) + NFrzVir)
     $           - Dbl_MB(I_EVal+NAO*(SpinMap(Spin)-1) + NFrzVir + 1) )
            If ( Delta .lt. Thresh_Degen ) then
               If ( .NOT. Header ) Write (LuOut, 8950)
               Header = .TRUE.
               Write (LuOut, 8960) 'virtual ', Spin, Delta
            EndIf
         EndIF
      EndDo
 8950 Format(/1X, '>>>>> WARNING!  FREEZE directive appears to spit ',
     $     'degeneracy    <<<<<' )

 8960 Format(1X, '>>>>> WARNING!  in ', A, ' space, spin ', I2,
     $     ' (eval diff ', 1P E7.1, ') <<<<<')
C
C     Free up memory
C
c$$$      If ( .NOT. MA_Pop_Stack(H_Eval0) ) Call ErrQuit(
c$$$     $     'corr_mk_ref: failed to pop stack', 0)
c$$$      If ( .NOT. MA_Pop_Stack(H_Occ0) ) Call ErrQuit(
c$$$     $     'corr_mk_ref: failed to pop stack', 0)
C
      If ( .NOT. MA_Pop_Stack(H_Eval) ) Call ErrQuit(
     $     'corr_mk_ref: failed to pop stack', 0, MA_ERR)
      If ( .NOT. MA_Pop_Stack(H_Occ) ) Call ErrQuit(
     $     'corr_mk_ref: failed to pop stack', 0, MA_ERR)
      Do Spin = VecSpin, 1, -1
         If ( .NOT. GA_Destroy(G_Vec(Spin) ) ) Call ErrQuit(
     $        'corr_mk_ref: failed to destroy ga', 0, GA_ERR)
      EndDo
C
C     Print results for the user
C
      If ( Util_Print('reference', Print_Debug) ) then
         If ( Me .eq. 0) then
            Call DGEWR(MxRefNMO, RefSpin, Dbl_MB(I_RefOcc), NAO, 0, 0,
     $           LuOut, 'Reference Occupation Vectors', 80, 'F6.3')
            Call DGEWR(MxRefNMO, RefSpin, Dbl_MB(I_RefEval), NAO, 0, 0,
     $           LuOut, 'Reference Eigenvalues', 80, 'F12.6')
         EndIf
         Do Spin = 1, RefSpin
            Call GA_Print( g_RefVec(Spin) )
         EndDo
      EndIf
C
      Corr_Mk_Ref = .TRUE.
      if(me.eq.0) Call Util_Flush(LuOut)
      Return
      End
