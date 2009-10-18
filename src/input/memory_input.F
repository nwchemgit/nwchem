      subroutine memory_input(rtdb, ostore, stack_size, heap_size, 
     $     global_size, total_size, overify, ohardfail)
C$Id: memory_input.F,v 1.9 2003-10-17 22:55:18 carlfahl Exp $
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
      parameter (nunits = 9)
      parameter (meg = 1024*1024)
      parameter (gig = 1024*1024*1024)
      character*8 units(nunits)
      data units /'double','real','integer','byte','kb',
     &            'mb','mw','gb','gw'/
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
            if (.not. inp_f(user_heap)) goto 1000
            heap_Src = 'U'
         else if (inp_compare(.false.,'stack', ctest)) then
            if (.not. inp_f(user_stack)) goto 1000
            Stack_Src = 'U'
         else if (inp_compare(.false.,'global', ctest)) then
            if (.not. inp_f(user_global)) goto 1000
            Global_Src = 'U'
         else if (inp_compare(.false.,'total', ctest)) then
            if (.not. inp_f(user_total)) goto 1000
            Total_Src = 'U'
         else if (inp_match(nunits, .false., ctest, units, ind)) then
            user_units = units(ind)
         else
            call inp_prev_field
            if (.not. inp_f(user_total)) goto 1000
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
c
c     Convert user units to doubles/reals
c
      if (user_units.eq.units(1) .or. user_units.eq.units(2)) then
         heap_size   = user_heap
         stack_size  = user_stack
         global_size = user_global
         total_size  = user_total
      else if (user_units .eq. units(3)) then
         heap_size   = ma_sizeof(mt_int, int(user_heap),   mt_dbl)
         stack_size  = ma_sizeof(mt_int, int(user_stack),  mt_dbl)
         global_size = ma_sizeof(mt_int, int(user_global), mt_dbl)
         total_size  = ma_sizeof(mt_int, int(user_total),  mt_dbl)
      else if (user_units .eq. units(4)) then
         heap_size   = ma_sizeof(mt_byte, int(user_heap),   mt_dbl)
         stack_size  = ma_sizeof(mt_byte, int(user_stack),  mt_dbl)
         global_size = ma_sizeof(mt_byte, int(user_global), mt_dbl)
         total_size  = ma_sizeof(mt_byte, int(user_total),  mt_dbl)
      else if (user_units .eq. units(5)) then
         heap_size   = ma_sizeof(mt_byte, 1024*int(user_heap),   mt_dbl)
         stack_size  = ma_sizeof(mt_byte, 1024*int(user_stack),  mt_dbl)
         global_size = ma_sizeof(mt_byte, 1024*int(user_global), mt_dbl)
         total_size  = ma_sizeof(mt_byte, 1024*int(user_total),  mt_dbl)
      else if (user_units .eq. units(6)) then
         heap_size   = ma_sizeof(mt_byte, meg*int(user_heap),   mt_dbl)
         stack_size  = ma_sizeof(mt_byte, meg*int(user_stack),  mt_dbl)
         global_size = ma_sizeof(mt_byte, meg*int(user_global), mt_dbl)
         total_size  = ma_sizeof(mt_byte, meg*int(user_total),  mt_dbl)
      else if (user_units .eq. units(7)) then
         heap_size   = meg*user_heap
         stack_size  = meg*user_stack
         global_size = meg*user_global
         total_size  = meg*user_total
      else if (user_units .eq. units(8)) then
         heap_size   = ma_sizeof(mt_byte,gig*int(user_heap),   mt_dbl)
         stack_size  = ma_sizeof(mt_byte,gig*int(user_stack),  mt_dbl)
         global_size = ma_sizeof(mt_byte,gig*int(user_global), mt_dbl)
         total_size  = ma_sizeof(mt_byte,gig*int(user_total),  mt_dbl)
      else if (user_units .eq. units(9)) then
         heap_size   = gig*user_heap
         stack_size  = gig*user_stack
         global_size = gig*user_global
         total_size  = gig*user_total
      endif
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
