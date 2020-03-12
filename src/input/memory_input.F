      subroutine memory_input(rtdb, ostore, stack_size, heap_size, 
     $     global_size, total_size, overify, ohardfail)
C$Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
      integer rtdb              ! [input]
      logical ostore            ! [input]
      integer stack_size        ! [output]
      integer heap_size         ! [output]
      integer global_size       ! [output]
      integer total_size        ! [output]
      logical overify           ! [output]
      logical ohardfail         ! [output]
c
      character*8 user_units, ctest
      double precision user_stack, user_heap, user_global, user_total
      character*1 stack_Src, heap_Src, global_Src, total_Src
      integer mem(3)            ! Used to stuff data into the rtdb
      integer nunits
      integer ind
      integer meg,gig
      character*255 name        ! Used to store info in the rtdb
      logical memory_units
      external memory_units
      parameter (nunits = 9)
      parameter (meg = 1024*1024)
      parameter (gig = 1024*1024*1024)
      integer conv_units(nunits)
      character*8 units(nunits)
      data units /'double','real','integer','byte','kb',
     &            'mb','mw','gb','gw'/
      data conv_units /1,1,1,1,1024,
     &     meg,meg,gig,gig/
c
c     Parse memory directive
c
c     memory [[total] <total_size>] [stack <stack_size>] 
c            [heap <heap_size>] [global <global_size>]
c            [<units>] [verify|noverify] [hardfail|nohardfail]
c
c     input units default to doubles
c
c     Return required sizes in units of double precision words
c     or default if not specified in the input
c
c     If (ostore) store the values in the rtdb using the current context
c     (so that in the near future each module can control its own memory)
c
      user_stack  =  0.0
      user_heap   =  0.0
      user_global =  0.0
      user_total  =  0.0
      user_units  = 'double'
      stack_Src    = 'N'
      heap_Src     = 'N'
      global_Src   = 'N'
      total_Src    = 'N'
c
c     units conversion
c
      conv_units(3)= ma_sizeof(mt_int, conv_units(3),   mt_dbl)
      conv_units(4)= ma_sizeof(mt_byte, conv_units(4),   mt_dbl)
      conv_units(5)= ma_sizeof(mt_byte, conv_units(5),   mt_dbl)
      conv_units(6)= ma_sizeof(mt_byte, conv_units(6),   mt_dbl)
      conv_units(8)= ma_sizeof(mt_byte, conv_units(8),   mt_dbl)
c
c     Check input line is vaguely sensible
c
      call inp_set_field(0)
      if (.not. inp_a(ctest)) goto 1000
      if (.not. inp_compare(.false., ctest, 'memory')) goto 1000
      if (inp_n_field() .lt. 2) goto 1000
c
c     Parse the input line
c
 10   if (inp_a_trunc(ctest)) then
         if (inp_compare(.false.,'verify',ctest)) then
            overify = .true.
         else if (inp_compare(.false.,'noverify',ctest)) then
            overify = .false.
         else if (inp_compare(.false.,'trace',ctest)) then
            call ma_trace(1)
         else if (inp_compare(.false.,'hardfail',ctest)) then
            ohardfail = .true.
         else if (inp_compare(.false.,'nohardfail',ctest)) then
            ohardfail = .false.
         else if (inp_compare(.false.,'heap', ctest)) then
            if (.not. memory_units(heap_size,
     n           units,nunits,conv_units)) goto 1000
            heap_Src = 'U'
         else if (inp_compare(.false.,'stack', ctest)) then
            if (.not. memory_units(stack_size,
     n           units,nunits,conv_units)) goto 1000
            Stack_Src = 'U'
         else if (inp_compare(.false.,'global', ctest)) then
            if (.not. memory_units(global_size,
     n           units,nunits,conv_units)) goto 1000
            Global_Src = 'U'
         else if (inp_compare(.false.,'total', ctest)) then
            if (.not. memory_units(total_size,
     n           units,nunits,conv_units)) goto 1000
            Total_Src = 'U'
         else
            call inp_prev_field
            if (.not. memory_units(total_size,
     n           units,nunits,conv_units)) goto 1000
            Total_Src = 'U'
         endif
         goto 10
      endif
C
C     User inputs must be non-negative
C
      If ( User_Stack .lt. 0 .OR. User_Heap .lt. 0
     $   .OR. User_Global .lt. 0 .OR. User_Total .lt. 0) Call ErrQuit(
     $   'Memory_Input: Memory limits must be non-negative', 0,
     &       INPUT_ERR)
c$$$c
c$$$c     debug
c$$$c
c$$$      write(0,*) ' total ', user_total, ' heap ', user_heap, ' stack ',
c$$$     $     user_stack, ' global ', user_global, ' ',
c$$$     $     user_units, ' verify ', overify
C
C     Take what the user has entered and turn it into a complete
C     memory specification.  
C
      Call Memory_Defaults(Total_Size, Total_Src, Heap_Size, Heap_Src,
     $   Stack_Size, Stack_Src, Global_Size, Global_Src)
C
c$$$      write(0,*) ' total ', total_size, ' heap ', heap_size, ' stack ',
c$$$     $     stack_size, ' global ', global_size, ' doubles ',
c$$$     $     ' verify ', overify
c
c     Store data to the rtdb
c
      if(global_size.eq.0) call errquit(
     E     'memory_input: no global memory defined ',0,MEM_ERR)
      if (ostore) then
         mem(1) = heap_size
         mem(2) = stack_size
         mem(3) = global_size
         if (.not. context_prefix('memory',name))
     $        call errquit('memory_input: context buffer too small',0,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb, name, mt_int, 3, mem))
     $        call errquit('memory_input: rtdb_put failed', 0,
     &       RTDB_ERR)
      endif
      return
c
 1000 call errquit('input_mem_size: memory [[total] <total_size>] '//
     $   '[heap <heap_size>] [stack <stack_size>] '//
     $   '[global <global_size>] [<units>] [verify|noverify] '//
     $   '[hardfail|nohardfail]', 0, INPUT_ERR)
c
      end
      logical function memory_units(memory_size,
     n     units,nunits,conv_units)
      implicit none
#include "errquit.fh"
#include "inp.fh"
c
      integer memory_size ! [out]
      character*8 units(*) ! [in]
      integer nunits,conv_units(*) ! [in]
c
      double precision user_memory
      character*8 ctest
      integer ind
c      
      memory_units=.false.
      if (.not. inp_f(user_memory)) goto 1000
      if (.not.inp_a(ctest)) goto 1000
      if (inp_match(nunits, .false., ctest, units, ind)) then
         memory_size  = int(user_memory)*conv_units(ind)
         memory_units=.true.
      else
         call errquit('inp_match: units missing ',0,INPUT_ERR)
      endif
 1000 continue
      return
      end
