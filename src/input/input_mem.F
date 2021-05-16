      subroutine input_mem_size(stack_size, heap_size, global_size,
     $   overify, ohardfail)
C$Id$
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer stack_size        ! [output]
      integer heap_size         ! [output]
      integer global_size       ! [output]
      logical overify           ! [output]
      logical ohardfail         ! [output]
C
      integer Total_Size,mitob1
      Character*1 Total_Src, Heap_Src, Stack_Src, Global_Src
      integer util_mynodeid
      external util_mynodeid
c
c     Scan the input file for the memory directive and call
c     memory_input() to parse the line without storing the data
c     in the database (which may not be available for use)
c
c     Return also flag to enable/disable automatic memory verification.
c
c     Search for memory keyword from beginning of input file
c
      overify    = .true.       ! Verify unless told otherwise!
      ohardfail  = .false.      ! Don't hard fail
C
      if (util_mynodeid() .eq. 0) then ! MPI not GA due to init. order
         rewind LuIn
         call inp_init(LuIn,LuOut)
         if (inp_search(.false., 'memory', 1)) then
            call memory_input(-1, .false., stack_size, heap_size, 
     $           global_size, Total_Size, overify, ohardfail)
         Else
            Total_Src = 'N'
            Heap_Src = 'N'
            Stack_Src = 'N'
            Global_Src = 'N'
            Call Memory_Defaults( Total_Size, Total_Src, Heap_Size,
     $         Heap_Src, Stack_Size, Stack_Src, Global_Size,
     $         Global_Src )
         EndIf
c 2GB check on 32bit archs 
#ifndef EXT_INT
      if((stack_size+heap_size)*8d0.ge.2d0**31) call errquit
     c     ( ' sum of stack and heap must be < 2 GB on 32bit archs',
     c     0    ,MA_ERR)
      if(global_size*8d0.ge.2d0**31) call errquit
     c     ( ' global memory must be < 2 GB on 32bit archs',0,MA_ERR)
#endif      
         rewind LuIn
         call inp_init(LuIn,LuOut)
      endif
c
c     Has to be TCGMSG broadcast since GA is not yet initialized
c
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      call brdcst(Msg_HeapMem, heap_size,  mitob1, 0)
      call brdcst(Msg_StkMem, stack_size,  mitob1, 0)
      call brdcst(Msg_GblMem, global_size, mitob1, 0)
      call brdcst(Msg_StkMem,    overify,  mitob1, 0)
c
      end
