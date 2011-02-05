C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     RIMP2_IJ_Int -- Get integrals for IJ block
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RIMP2_IJ_Int(IJSame, D_I, ILo, IHi, NA,
     $                        D_J, JLo, JHi, NB, FitBas, NFit,
     $     G_I, ReqI, G_J, ReqJ)
      Implicit NONE
#include "errquit.fh"
      Logical IJSame![in]
      Integer D_I   ![in]
      Integer ILo   ![in]
      Integer IHi   ![in]
      Integer NA    ![in]
      Integer D_J   ![in]
      Integer JLo   ![in]
      Integer JHi   ![in]
      Integer NB    ![in]
      Integer FitBas![in]
      Integer NFit  ![in]
      Integer G_I   ![out]
      Integer ReqI  ![out]
      Integer G_J   ![out]
      Integer ReqJ  ![out]
C
C DESCRIPTION
C
C WARNING
C     This routine does NOT guarantee that the DRA requests to fill
C     the integral GAs have completed at the time it returns to the
C     caller!  It is the caller's responsibility to DRA_Wait() on 
C     these requests before G_I and G_J must actually be used.
C
C EXTERNAL ROUTINES
      Logical  GA_Create_One_By_Shells
      External GA_Create_One_By_Shells
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
#include "util.fh"
#include "rimp2g_timers.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer NI, NJ
C
      NI = IHi - ILo + 1
      NJ = JHi - JLo + 1
C     
C     Create GAs for 3c2e integrals
C     
      If ( .NOT. GA_Create_One_By_Shells( MT_Dbl, NI*NA, FitBas,
     $     '3c2e ints I', G_I) ) Call ErrQuit(
     $     'RIMP2_IJ_INT: can''t allocate integral GA', NI*NA*NFit,
     &       GA_ERR)
C     
      If ( .NOT. GA_Create_One_By_Shells( MT_Dbl, NJ*NB, FitBas,
     $     '3c2e ints J', G_J) ) Call ErrQuit(
     $     'RIMP2_IJ_INT: can''t allocate integral GA', NJ*NB*NFit,
     &       GA_ERR)
C     
C     Get integrals (ia|s') for this I
C     
      Call PStat_On( PS_Rd_3ci)
      If ( DRA_Read_Section(.FALSE.,
     $     G_I, 1,            NA*NI,  1, NFit,
     $     D_I, 1+NA*(ILo-1), NA*IHi, 1, NFit, ReqI) .ne. 0)
     $     Call ErrQuit( 'RIMP2_IJ_Int: Error reading I integrals', ILo,
     &       INT_ERR)
      If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
     $     'RIMP2_IJ_Int: Error from DRA_Wait on I integrals', 1,
     &       INT_ERR)
      Call PStat_Off( PS_Rd_3ci)
C     
C     Find the integrals (jb|s') for this J.  Use G_I instead of
C     disk if J is a subset.  We could go further and make use of
C     any overlap between I and J, but that might require two read
C     requests, which slightly messes up the idea of allowing
C     this routine to return with ONE request id for each array.
C     Not insurmountable, but not something I want to deal with now.
C     
      If ( IJSame .AND. JLo .ge. ILo .AND. JHi .le. IHi) then
C
C        Must be sure G_I is full before we copy
C
c$$$         Call PStat_On( PS_Rd_3ci)
c$$$         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_IJ_Int: Error from DRA_Wait', ILo)
c$$$         Call PStat_Off( PS_Rd_3ci)
         Call GA_Sync
C
         Call GA_Copy_Patch('N',
     $        G_I, (JLo-ILo)*NB+1, (JHi-ILo+1)*NB, 1, NFit,
     $        G_J, 1,              NJ*NB,          1, NFit )
      Else
         Call PStat_On( PS_Rd_3ci)
         If ( DRA_Read_Section(.FALSE.,
     $        G_J, 1,            NB*NJ,  1, NFit,
     $        D_J, 1+NB*(JLo-1), NB*JHi, 1, NFit, ReqJ) .ne. 0)
     $        Call ErrQuit( 'RIMP2_IJ_Int: Error reading J integrals',
     $        JLo, DISK_ERR)
         If ( DRA_Wait(ReqJ) .ne. 0) Call ErrQuit(
     $        'RIMP2_IJ_Int: Error from DRA_Wait on J integrals', 2,
     &       CALC_ERR)
         Call PStat_Off( PS_Rd_3ci)
      EndIf
C     
C     Must wait until GAs are full before we can print
C
      If ( Util_Print('3c ints', Print_Debug) ) then
c$$$         Call PStat_On( PS_Rd_3ci)
c$$$         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_IJ_Int: Error from DRA_Wait', ILo)
c$$$         Call PStat_Off( PS_Rd_3ci)
         Call GA_Print( G_I )
C
c$$$         If ( DRA_Wait(ReqI) .ne. 0) Call ErrQuit(
c$$$     $        'RIMP2_IJ_Int: Error from DRA_Wait', ILo)
c$$$         Call PStat_Off( PS_Rd_3ci)
         Call GA_Print( G_J )
      EndIf
C
      Return
      End
