      subroutine memory_input(rtdb, ostore, stack_size, heap_size, 
     $     overify)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"
      integer rtdb              ! [input]
      logical ostore            ! [input]
      integer stack_size        ! [output]
      integer heap_size         ! [output]
      logical overify           ! [output]
c
      character*8 user_units, ctest
      double precision user_stack, user_heap
      integer mem(2)            ! Used to stuff data into the rtdb
      integer nunits
      integer ind
      integer default
      character*255 name        ! Used to store info in the rtdb
      parameter (nunits = 7)
      character*8 units(nunits)
      data units /'double', 'real', 'integer', 'byte', 'kb', 'mb', 'mw'/
#ifdef DELTA
      data default /917504/
#else
      data default /2097152/
#endif
c
c     Parse memory directive
c
c     memory [[heap] <heap_size>] [stack <stack_size>] [<units>] [verify|noverify]
c
c     input units default to doubles
c
c     Return required sizes in units of double precision words
c     or default if not specified in the input
c
c     Delta default = 7  Mbytes =   917,504 doubles
c     Other default = 16 Mbytes = 2,097,152 doubles
c
c     If (ostore) store the values in the rtdb using the current context
c     (so that in the near future each module can control its own memory)
c
      user_stack = 0.0
      user_heap  = 0.0
      user_units = 'double'
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
         else if (inp_compare(.false.,'heap', ctest)) then
            if (.not. inp_f(user_heap)) goto 1000
         else if (inp_compare(.false.,'stack', ctest)) then
            if (.not. inp_f(user_stack)) goto 1000
         else if (inp_match(nunits, .false., ctest, units, ind)) then
            user_units = units(ind)
         else
            call inp_prev_field
            if (.not. inp_f(user_heap)) goto 1000
         endif
         goto 10
      endif
c
c     debug
c
*      write(6,*) ' heap ', user_heap, ' stack ', user_stack,' ',
*     $     user_units, ', verify ', overify
c
c     Convert user units to doubles/reals
c
      if (user_units.eq.units(1) .or. user_units.eq.units(2)) then
         heap_size = user_heap
         stack_size= user_stack
      else if (user_units .eq. units(3)) then
         heap_size = ma_sizeof(mt_int, int(user_heap), mt_dbl)
         stack_size= ma_sizeof(mt_int, int(user_stack), mt_dbl)
      else if (user_units .eq. units(4)) then
         heap_size = ma_sizeof(mt_byte, int(user_heap), mt_dbl)
         stack_size= ma_sizeof(mt_byte, int(user_stack), mt_dbl)
      else if (user_units .eq. units(5)) then
         heap_size = ma_sizeof(mt_byte, 1024*int(user_heap), mt_dbl)
         stack_size= ma_sizeof(mt_byte, 1024*int(user_stack), mt_dbl)
      else if (user_units .eq. units(6)) then
         heap_size = ma_sizeof(mt_byte,1024*1024*int(user_heap), mt_dbl)
         stack_size= ma_sizeof(mt_byte,1024*1024*int(user_stack),mt_dbl)
      else if (user_units .eq. units(7)) then
         heap_size = 1024*1024*user_heap
         stack_size= 1024*1024*user_stack
      endif
c
c     Make sure that defaults OK even if units were changed
c
      if (heap_size .le. 0) heap_size = -1
      if (stack_size .le. 0) stack_size = -1
      if (heap_size.eq.-1 .and. stack_size.eq.-1) then
         heap_size = default / 2
         stack_size= heap_size
      endif
#ifdef DELTA
c
c     The RTDB will not function if it cannot malloc some memory
c     and 7MB is the upper limit of available user memory
c
      if ((heap_size+stack_size).gt.default) call errquit
     $     ('memory_input: delta memory maximum is 7 Mb', 
     $     heap_size+stack_size)
#endif
c
c     Store data to the rtdb
c
      if (ostore) then
         mem(1) = heap_size
         mem(2) = stack_size
         if (.not. context_prefix('memory',name))
     $        call errquit('memory_input: context buffer too small',0)
         if (.not. rtdb_put(rtdb, name, mt_int, 2, mem))
     $        call errquit('memory_input: rtdb_put failed', 0)
      endif
      return
c
 1000 call errquit('input_mem_size: memory [[heap] <heap_size>] '//
     $     '[stack <stack_size>] [<units>] [verify|noverify]', 0)
c
      end
