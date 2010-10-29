      Subroutine PStat_Print_Handles(handles,n_handles)
c $Id$
      Implicit NONE
C
#include "pstatP.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
C
      integer n_handles
      integer handles(n_handles)
      integer current_handle
c
      Integer I, IMAM(3), icalls
      Double Precision DMAM(3, 2)
C
      Call GA_Sync()
C
      If ( GA_NodeID() .eq. 0) then
         Call Banner(luout, 'Performance information', '-',
     $   .TRUE., .TRUE., .FALSE.)
         Write (luout, 9000)
      EndIf
 9000 Format(1X, 21X,
     $   '             Nr. of calls          ',1X,
     $   '      CPU time (s)  ',1X, '           Wall time (s)  ',/
     $   1X, 21X,
     $   '       -------------------------      ---------------------',
     $     '      ---------------------'/
     $   1X, 'Name', 17X, '       Min', 1X, '       Avg', 1X,
     $                     '       Max',1X,
     $  '     Min',1X, '     Avg',1X, '     Max',1X,
     $  '     Min',1X, '     Avg',1X, '     Max')
C
      Do I = 1, n_handles
         current_handle = handles(i)
         if (pstat_active(current_handle)) then
C     
C     Determine the min/avg/max for the number of calls
C     
            iMAM(1) = PStat_Nr_Calls( current_handle )
            iMAM(2) = PStat_Nr_Calls( current_handle )
            iMAM(3) = PStat_Nr_Calls( current_handle )
C     
            Call ga_IGOp(Msg_PStat_Min, iMAM(1), 1, 'min')
            Call ga_IGOp(Msg_PStat_Avg, iMAM(2), 1, '+')
            iMAM(2) = iMAM(2) / GA_NNodes()
            Call ga_IGOp(Msg_PStat_Max, iMAM(3), 1, 'max')
C     
C     Determine min/avg/max for the times
C     
            dMAM(1, 1) = PStat_Times( current_handle , 1)
            dMAM(2, 1) = PStat_Times( current_handle , 1)
            dMAM(3, 1) = PStat_Times( current_handle , 1)
            dMAM(1, 2) = PStat_Times( current_handle , 2)
            dMAM(2, 2) = PStat_Times( current_handle , 2)
            dMAM(3, 2) = PStat_Times( current_handle , 2)
C     
            Call ga_DGOp(Msg_PStat_Min, dMAM(1, 1), 1, 'min')
            Call ga_DGOp(Msg_PStat_Avg, dMAM(2, 1), 1, '+')
            dMAM(2, 1) = dMAM(2, 1) / GA_NNodes()
            Call ga_DGOp(Msg_PStat_Max, dMAM(3, 1), 1, 'max')
C     
            Call ga_DGOp(Msg_PStat_Min, dMAM(1, 2), 1, 'min')
            Call ga_DGOp(Msg_PStat_Avg, dMAM(2, 2), 1, '+')
            dMAM(2, 2) = dMAM(2, 2) / GA_NNodes()
            Call ga_DGOp(Msg_PStat_Max, dMAM(3, 2), 1, 'max')
C     
            If ( GA_NodeID() .eq. 0) then
              icalls = iMAM(1) + iMAM(2) + iMAM(3)
              if (icalls.gt.0)
     &            Write (luout, 9010)
     &            PStat_Name(current_handle), iMAM, dMAM
            endif
          endif
      EndDo
c
      if (ga_nodeid() .eq. 0) call util_flush(luout)
C
 9010 Format(1X, A20, 1X, 3(I10,1X), 1P, 3(G8.2, 1X), 2X, 3(G8.2, 1X))
C
      Return
      End
