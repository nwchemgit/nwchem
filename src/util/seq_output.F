C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     begin_seq_output -- Mark the start of a section of sequentialized 
C     output
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Begin_Seq_Output
      Implicit NONE
#include "errquit.fh"
C
C DESCRIPTION
C     Indicates the start of a section of sequentialized output.
C     Should be followed by any number of Write_Seq calls terminated
C     by an End_Seq_Output.  All output will be sent to node 0 and
C     written there in order my node number.
C
C     All nodes must participate in a sequential output section.
C
C INCLUDE FILES
#include "msgids.fh"
#include "global.fh"
#include "mafdecls.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer GoAhead, Len, Node, Me, TCG_Master,mitob1
C
C     Translate from GA_NodeID to TCGMSG NODEID
C     Only need Node 0 at the moment
C
      Call GA_List_NodeID(TCG_Master, 1)
C
C     Expect node 0 to tell us to go ahead with our output.  If we
C     are node 0, this is a no-op.
C
      Me = ga_nodeid()
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      If ( Me .ne. TCG_Master) then
         Call RCV(Msg_SeqOut_Ctrl, GoAhead, MItoB1, Len,
     $      TCG_Master, Node, 1)
         If ( GoAHead .ne. Me .OR. Node .ne. TCG_Master)
     $      Call ErrQuit(
     $        'Begin_Seq_Output: Incorrect signal from node 0: ',
     $        GoAhead, CALC_ERR)
      EndIf
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     end_seq_output -- Mark the end of a section of sequentialized 
C     output
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine End_Seq_Output
      Implicit NONE
#include "errquit.fh"
C
C DESCRIPTION
C     Indicates the end of a section of sequentialized output.
C
C     All nodes must participate in a sequential output section.
C
C BUGS
C     Because we have to declare a fixed length string, it is possible
C     for some transmissions to be truncated.  In practice, however,
C     we choose something rather longer than typical line lengths and
C     it should not be a serious problem.
C
C INCLUDE FILES
#include "msgids.fh"
#include "global.fh"
#include "mafdecls.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C PARAMETERS
C
C     Text longer than MaxLen will be truncated at this end.
C
      Integer MaxLen
      Parameter (MaxLen = 256)
C
C LOCAL VARIABLES
      Integer Len, Me, LUN, Node, FrNode, I, TCG_Master, H_Map, I_Map
      integer mitob1
      Character*(MaxLen) Text
C
C     Only GA nodes can participate, but we need TCGMSG node ids.
C
      if ( .NOT. MA_Push_Get(MT_Int, GA_NNodes(), 'Node list',
     $   H_Map, I_Map) )
     $   Call ErrQuit('End_Seq_Output: Insufficient memory',
     $   GA_NNodes() , MA_ERR)
      Call GA_List_NodeID(Int_MB(I_Map), GA_NNodes())
      TCG_Master = Int_MB(I_Map)
C
      Me = ga_NodeID()
C
C     If we are not node 0, simply say we are done
C
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      If ( Me .ne. TCG_Master) then
         Call SND(Msg_SeqOut_Ctrl, -Me, MItoB1, TCG_Master, 1)
         If ( .NOT. MA_Pop_Stack(H_Map)) Call ErrQuit(
     $      'End_Seq_Output: Unable to free local memory', 0, MA_ERR)
         Return
      EndIf
C
C     If we are node 0, we have more responsibilities.  At this point
C     we have printed all of our stuff.  Now we need to listen to each
C     of the other nodes in turn.
C
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      Do I = 1, GA_NNodes() - 1
         Node = Int_MB(I_Map + I)
C
C        Tell them to go ahead
C
         Call SND(Msg_SeqOut_Ctrl, Node, MItoB1, Node, 1)
C
C        Listen for their output until we get -Node back for the
C        control signal.
C
 1000    Call RCV(Msg_SeqOut_Ctrl, LUN, MItoB1, Len, Node, FrNode, 1)
C
         If ( LUN .ge. 0) then
            Call FCRCV(Msg_SeqOut_Text, Text, Len, Node, FrNode, 1)
            Write (LUN, '(A)') Text(:Min(Len,MaxLen))
            Goto 1000
         ElseIf ( LUN .ne. -Node) then
            Call ErrQuit('End_Seq_Output: Incorrect signal from node',
     $         Node, UNKNOWN_ERR)
         EndIf
      EndDo
C
      If ( .NOT. MA_Pop_Stack(H_Map)) Call ErrQuit(
     $   'End_Seq_Output: Unable to free local memory', 0, MA_ERR)
C
      Return
      End
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     write_seq -- Write sequentialized output
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine Write_Seq(LUN, Text)
      Implicit NONE
#include "errquit.fh"
      Integer LUN
      Character*(*) Text
C
C DESCRIPTION
C     Write formatted output on the chosen logical unit on node 0.
C     If invoked on nodes other than 0 sends data to node 0.  Must
C     be surrounded by Begin_Seq_Output/End_Seq_Output calls.
C
C     Observe that LUN is the unit on node 0, not the invoking node!
C
C INCLUDE FILES
#include "msgids.fh"
#include "global.fh"      
#include "mafdecls.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer TCG_Master,mitob1
C
C     Safety check the node number because it is also a signal
C     to node 0.
C
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      If ( LUN .lt. 0) Call ErrQuit(
     $   'Write_Seq: Negative unit numbers are invalid', LUN, INPUT_ERR)
C
C     Translate from GA_NodeID to TCGMSG NODEID
C     Only need Node 0 at the moment
C
      Call GA_List_NodeID(TCG_Master, 1)
C
C     Of course node 0 prints what it has, everyone else sends to node 0.
C
      If (ga_NodeID() .eq. TCG_Master) then
         Write (LUN, '(A)') Text
      Else
         Call SND(Msg_SeqOut_Ctrl, LUN, MItoB1, TCG_Master, 1)
         Call FCSND(Msg_SeqOut_Text, Text, TCG_Master, 1)
      EndIf
C
      Return
      End
