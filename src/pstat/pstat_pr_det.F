      Subroutine PStat_Print_Detail(handles,n_handles)
c $Id$
      Implicit NONE
#include "errquit.fh"
C
#include "pstatP.fh"
#include "global.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "stdio.fh"
C
      integer n_handles
      integer handles(n_handles)
      integer current_handle
c
C
      Integer I, Off, N, H_Collect, I_Collect
C
      Call GA_Sync()
C
      If ( GA_NodeID() .eq. 0) then
         Call Banner(luout, 'Detailed Performance Information', '-',
     $   .TRUE., .TRUE., .FALSE.)
         Write (luout, 9000)
      EndIf
 9000 Format(1X, 'Name', 18X, 'Node', 6X, 'Calls', 1X,
     $     'CPU Time (s)', 2X, 'Wall Time (s)')
C
C     Create a local array to collect the info into
C
      If (.NOT. MA_Push_Get( MT_Dbl, 3*N_Handles*GA_NNodes(),
     $     'pstat temporary', H_Collect, I_Collect) ) Call ErrQuit(
     $     'Pstat_Print_Detail: Unable to allocate local memory', 0,
     &       MA_ERR)
C
      Call DFill(3*N_Handles*GA_NNodes(), 0.0d0, Dbl_MB(I_Collect), 1)
C
C     Fill with our local data
C
      Off = 3 * GA_NodeID()
C
      Do I = 1, n_handles
         current_handle = handles(i)
         if (pstat_active(current_handle)) then
C     
            Dbl_MB( I_Collect+Off)  =
     $           Dble( PStat_Nr_Calls( current_handle ) )
            Dbl_MB( I_Collect+Off+1) = PStat_Times( current_handle , 1)
            Dbl_MB( I_Collect+Off+2) = PStat_Times( current_handle , 2)
C
            Off = Off + 3 * GA_NNodes()
          endif
      EndDo
C
      Call ga_DGOp(Msg_PStat_Avg, Dbl_MB(I_Collect),
     $     3*N_Handles*GA_NNodes(), '+')
C
      If ( GA_NodeID() .eq. 0) then
         Off = 0
         Do I = 1, n_handles
            current_handle = handles(i)
            if (pstat_active(current_handle)) then
C
               Do N = 0, GA_NNodes()-1
                  Write (luout, 9010)
     &                 PStat_Name(current_handle), N,
     $                 Int( Dbl_MB( I_Collect + Off )),
     $                 Dbl_MB( I_Collect + Off + 1 ),
     $                 Dbl_MB( I_Collect + Off + 2 )
                  Off = Off + 3
               EndDo
            EndIF
         EndDo
         call util_flush(luout)
      EndIf
C
      If ( .NOT. MA_Pop_Stack( H_Collect) ) Call ErrQuit(
     $   'PStat_Print_Detail: unable to free memory', 0, MA_ERR)
C
 9010 Format(1X, A20, 1X, I5, 1X, I10, 1X, F13.3, 1X, F13.3)
C
      Return
      End
