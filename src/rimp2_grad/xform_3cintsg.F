C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     Xform_3cintsG -- Produce transformed three-center integrals
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Logical Function Xform_3cintsG(RTDB, Geom, Type3ci, AOBas, FitBas,
     $   ALo, AHi, ILo, IHi, G_C, XFFit, G_V,
     $   Dest, G_Out, IntFil, Neglect)
      Implicit NONE
#include "errquit.fh"
      Integer RTDB, AOBas, FitBas
      Integer IntFil
      Integer ALo, AHi, ILo, IHi
      Integer G_Out, G_C, G_V, Geom
      Character*1 Type3ci, XFFit, Dest
      Double precision Neglect
C
C ARGUMENTS
C     Type3ci  Type of 3-center integrals to be used ('E'RI or 'O'verlap) [IN]
C     AOBas    Handle to basis for two of the indices ("AO basis") [IN]
C     FitBas   Handle to basis for the third index ("fitting basis") [IN]
C     ALo      Lower limit of target index range for first AO basis index [IN]
C     AHi      Upper limit of target index range for first AO basis index [IN]
C     ILo      Lower limit of target index range for 2nd AO basis index [IN]
C     IHi      Lower limit of target index range for 2nd AO basis index [IN]
C     G_C      GA handle for transformation matrix for AO basis indices [IN]
C     XFFit    'Y'es to transform fitting basis index, 'N'o to skip it [IN]
C     G_V      GA handle for transformation matrix for fitting basis index [IN]
C     Dest     Send transformed integrals to 'D'isk or 'M'emory [IN]
C     G_Out    GA handle for output buffer of integrals [IN]
C     IntFil   LUN for file to hold transformed integrals [IN]
C     Neglect  Neglect raw integrals smaller than this value [IN]
C
C RETURN VALUE
C     .TRUE.   Successful completion
C     .FALSE.  Errors in routine or dependents
C
C DESCRIPTION
C     Produce transformed three-center integrals either in disk or in memory.
C
C     The global arrays G_Out, G_C, G_V must be created and filled
C     as appropriate before entry.
C
C NOTES
C     Since the Schwarz routines aren't really general enough to
c     accomodate 3 center integrals yet, the NEGLECT threshold passed in
c     may have been manipulated to account for the size of the (s|s)
c     integrals that form the second term of the limit.  Neglect may in
c     fact be different on each node, since that's the obvious way to
c     get the best possible screening in the current circumstances.  But
c     only node 0 reports Neglect to the user, so it may be a little
c     misleading.
C
C LIMITATIONS
C     Currently DEST must be 'D'isk.
C     G_Out must be the same size
C
C MEMORY USE
C     Not here, but our dependents do use some.  I'd rather not talk about
C     it at the moment.
C
C INCLUDE FILES
#include "numerical_constants.fh"
#include "stdio.fh"
#include "bas.fh"
#include "global.fh"
#include "rimp2g_timers.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "dra.fh"
C
C EXTERNAL ROUTINES
      Logical LSame     ! From BLAS-[23]
      External LSame
      Logical  XF3cI_Step12, XF3cI_Step3, Put_3cInts
      External XF3cI_Step12, XF3cI_Step3, Put_3cInts
      Integer LNBlnk
      External LNBlnk
      Logical GA_Create_One_by_Shells
      External GA_Create_One_by_Shells
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Logical Status, TransFit, Disk
      Integer Dim1, Dim2, GAType, NFit, G_Int
      Integer IBlkSiz, OutBlkSiz, NA, NI, NJ, NK, JLo, JHi, KLo, KHi
      Character*256 String1
      Integer Req
C
C     Since we'll synchronize pretty soon on the GA creations,
C     let's do it here to make the timings look right
C
      Call GA_Sync
      Call PStat_On( PS_XForm)
C
C     ***********************
C     * Sanity check inputs *
C     ***********************
C     Type3ci is handled by XF3cI_Step12, and it takes little time to
C     get to that point, so leave that check out here.
C
      Disk = LSame( Dest, 'D')
      If ( .NOT. (Disk .or. LSame( Dest, 'M'))) then
         Write (LuErr, *) 'Xform_3cintsG: Can only place result on ',
     $      'Disk or in Memory.'
         Xform_3cintsG = .False.
         Return
      EndIf
      TransFit = LSame( XFFit, 'Y')
      If ( .NOT. (TransFit .or. LSame( XFFit, 'N'))) then
         Write (LuErr, *) 'Xform_3cintsG: Transform fitting basis ',
     $      'is a Y/N question.'
         Xform_3cintsG = .False.
         Return
      EndIf
C
C     Basis sets
C
      Status = .TRUE.
      Status = Status .AND.
     $   Bas_Check_Handle( AOBas, 'Xform_3cintsG: AO basis')
      Status = Status .AND.
     $   Bas_Check_Handle( FitBas, 'Xform_3cintsG: Fitting basis')
C
C     Message already printed by the checking routine
C
      If ( .NOT. Status ) then
         Xform_3cintsG = .False.
         Return
      EndIf
C
C     Global arrays which should be already created
C
      Call GA_Check_Handle( G_C, 'Xform_3cintsG: AO transformation')
      If ( TransFit )
     $   Call GA_Check_Handle( G_V, 'Xform_3cintsG: Fit transformation')
      If ( TransFit .or. .not. Disk)
     $   Call GA_Check_Handle( G_Out, 'Xform_3cintsG: output array')
C
C     Now sanity check the index ranges provided.
C
      Status = .TRUE.
      Status = Status .AND. ALo .gt. 0 .AND. AHi .gt. 0
     $   .AND. AHi .ge. ALo
      Status = Status .AND. ILo .gt. 0 .AND. IHi .gt. 0
     $   .AND. IHi .ge. ILo
C
      If ( .NOT. Status ) then
         Write (LuErr, 9000) ALo, AHi, ILo, IHi
         Xform_3cintsG = .False.
         Return
      EndIf
 9000 Format (1X, 'Xform_3cintsG: Bad range for transformation: (',
     $   I5, ':', I5, ', ', I5, ':', I5, ')')
C
C     Compute some basic quantities
C
      NA = AHi - Alo + 1
      NI = IHi - Ilo + 1
C
      Status = Bas_NumBF( FitBas, NFit)
      If ( .NOT. Status ) then
         Write (LuErr, *) 'Xform_3cintsG: Unable to get basis info.'
         Xform_3cintsG = .FALSE.
         Return
      EndIf
C
C     ******************************
C     * Determine number of passes *
C     ******************************
C     Use the sizes of the GAs passed in to determine block sizes.
C     This may have to change somewhat when we want to put results
C     in GAs instead of to disk.
C
      Call GA_Inquire( G_Out, GAType, Dim1, Dim2)
      IBlkSiz = Dim1 / NA
C
C     Also verify the size
C
      If ( Dim1 .lt. NA .OR. Dim2 .lt. NFit) then
         Write (LuErr, 9010) Dim1, Dim2, NA, NFit
         Xform_3cintsG = .FALSE.
         Return
      EndIf
C
      If ( TransFit ) then
         Call GA_Inquire( G_Out, GAType, Dim1, Dim2)
         OutBlkSiz = Dim1 / NA
C
C        Also verify the size
C
         If ( Dim1 .lt. NA .OR. Dim2 .lt. NFit) then
            Write (LuErr, 9010) Dim1, Dim2, NA, NFit
            Xform_3cintsG = .FALSE.
            Return
         EndIf
      Else
         OutBlkSiz = IBlkSiz
      EndIf
C
 9010 Format( 1X, 'Xform_3cintsG: GA matrix dimensions ',
     $   'inconsistent with request:'
     $   /1X, I5, 'x', I5, ' vs. ', I5, 'x', I5, '.')
C
C     Provide some info to the user
C
      If ( Util_Print('information', Print_Low)
     $   .AND. GA_NodeID() .eq. 0) then
C
c$$$         Write (LuOut, 9110) IBlkSiz, (NI + IBlkSiz - 1) / IBlkSiz
      EndIf
C
 9110 Format(1X, 'Transformation block size', I5, ', requiring ',
     $   I5, ' integral passes.')
C
C     ****************************************************
C     * Cover the requested range of I in several passes *
C     ****************************************************
C     The block size is determined by the global array G_Out.
C     The current range of I is [JLo,JHi].
C
      Do JLo = ILo, IHi, IBlkSiz
         JHi = Min( JLo+IBlkSiz-1, IHi)
         NJ = JHi - JLo + 1
C
C        Create the GA needed for the first 2/3 transformation
C
         If ( .NOT. GA_Create_One_By_Shells( MT_Dbl, NA * NJ,
     $      FitBas, 'xf3ci: 2/3 integral matrix', G_Int) )
     $        Call ErrQuit(
     $        'Xform_3cintsG: Unable to allocate temporary GA', 0,
     &       GA_ERR)
C
C        Generate blocks of (ai|s) (or (ais) ) integrals for the
C        specified range of I.
C
         Status = XF3cI_Step12( RTDB, Geom, Type3ci, ALo, AHi,
     $      JLo, JHi, AOBas, G_C, FitBas, G_Int, Neglect)
         If ( .NOT. Status ) then
            Write (LuErr, 9020) JLo, JHi, Ilo, Ihi
            Xform_3cintsG = .FALSE.
            Return
         EndIf
 9020    Format( 1x, 'Xform_3cintsG: Error occurred in integral block ',
     $      I5, ':', I5, ' of requested ', I5, ':', I5, '.')
C
         If ( Util_Print('2/3 ints', Print_Debug) ) Call GA_Print(G_Int)
C
C        Transform fitting basis & write out the result
C
         If ( TransFit ) then
C
C           For this transformation, we need another global array in 
C           which to deposit the results.  It may not be the same size
C           as G_Int -- for space efficiency when the results ultimately
C           go to disk, for example, G_Out should hold only a single I.
C
            Do KLo = JLo, JHi, OutBlkSiz
               KHi = Min( KLo+OutBlkSiz-1, JHi)
               NK = KHi - KLo + 1
C
C              Do the final 1/3 transformation.  Note that here
C              we specify the region to transform relative to the
C              existing global array, not relative to absolute
C              indices of the occupied & virtual orbitals.
C
               Call GA_Zero(G_Out)
               Status = XF3cI_Step3( G_Int, NA, NJ, 1, NA,
     $            KLo-JLo+1, KHi-JLo+1, 1, NFit,
     $            G_V, 1, NFit, G_Out, NA, NK, 1, 1, 1)
               If ( .NOT. Status ) then
                  Write (LuErr, 9030) KLo, KHi, JLo, JHi, ILo, IHi
                  Xform_3cintsG = .FALSE.
                  Return
               EndIf
C
C              Now we're back to actual indices of the occupied
C              and virtual orbitals.
C
               If ( Disk ) then
                  Call PStat_On( PS_Wr_3ci)
c$$$                  Status = Put_3cInts(IntFil, G_Out, NA, KLo-ILo+1,
c$$$     $                 KHi-ILo+1)
c$$$                  If ( .NOT. Status ) then
c$$$                     Write (LuErr, 9030) KLo, KHi, JLo, JHi, ILo, IHi
c$$$                     Xform_3cintsG = .FALSE.
c$$$                     Return
c$$$                  EndIf
                  If ( DRA_Write_Section(.FALSE., G_Out, 1,
     $                 NA*(KHi-KLo+1), 1, NFit, IntFil,
     $                 1+NA*(KLo-ILo), NA*(KHi-ILo+1), 1, NFit, Req)
     $                 .ne. 0) Call ErrQuit(
     $                 'Put_3cInts: Error from DRA_Write_Section',
     $                 KLo-ILo+1, DISK_ERR)
C     
                  If ( DRA_Wait(Req) .ne. 0) Call ErrQuit(
     $                 'Put_3cInts: Error from DRA_Wait', KLo-ILo+1,
     &       DISK_ERR)
                  Call PStat_Off( PS_Wr_3ci)
               EndIf
C
               If ( Util_Print('mo ints', Print_Debug) )
     $            Call GA_Print(G_Out)
            EndDo               ! KLo
C
C           If we are skipping the final 1/3-transformation, we
C           can just write directly from G_Int in one go.
C
         Else
            If ( Disk ) then
               Call PStat_On( PS_Wr_3ci)
c$$$               Status = Put_3cInts(IntFil, G_Int, NA, JLo-ILo+1,
c$$$     $              JHi-ILo+1)
c$$$               If ( .NOT. Status ) then
c$$$                  Write (LuErr, 9030) JLo, JHi, JLo, JHi, ILo, IHi
c$$$                  Xform_3cintsG = .FALSE.
c$$$                  Return
c$$$               EndIf
                  If ( DRA_Write_Section(.FALSE., G_Int, 1,
     $                 NA*(JHi-JLo+1), 1, NFit, IntFil,
     $                 1+NA*(JLo-ILo), NA*(JHi-ILo+1), 1, NFit, Req)
     $                 .ne. 0) Call ErrQuit(
     $                 'Put_3cInts: Error from DRA_Write_Section',
     $                 JLo-ILo+1, DISK_ERR)
C     
                  If ( DRA_Wait(Req) .ne. 0) Call ErrQuit(
     $                 'Put_3cInts: Error from DRA_Wait', JLo-ILo+1,
     &       DISK_ERR)
               Call PStat_Off( PS_Wr_3ci)
            EndIf
C
            If ( Util_Print('mo ints', Print_Debug) ) then
               If ( TransFit ) then
                  Call GA_Print(G_Out)
               Else
                  Call GA_Print(G_Int)
               EndIf
            EndIf
C
         EndIf
 9030    Format( 1X, 'Xform_3cintsG: Error occurred in output block ',
     $      I5, ':', I5, ' of integral block ', I5, ':', I5,
     $      ' of requested ', I5, ':', I5, '.')
C
C        Free 2/3 integral array to allow it to be resized next pass
C
         If ( .NOT. GA_Destroy( G_Int) ) Call ErrQuit(
     $        'Xform_3cintsG: Unable to deallocate temporary GA', 0,
     &       GA_ERR)
C
C        Print out a progress report.  Eventually want to make this
C        dependent on user desires and time required by each pass.
C
c$$$         If ( Util_Print('progress reports', Print_Default) ) then
c$$$            Write (String1, 9200) JLo, JHi
c$$$ 9200       Format('Xform_3cintsG: Completed ', I5, ':', I5, ' pass.')
c$$$            Call Progress_Report(LuOut, String1(:LNBlnk(String1)) )
c$$$         EndIf
      EndDo                     ! JLo
C
      Call GA_Sync
C
      Xform_3cintsG = .TRUE.
C
      Call PStat_Off( PS_XForm)
C
      Return
      End
