      subroutine input_set(rtdb)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input]
      character*255 field       ! input buffer
      character*255 name        ! name of variable being set
      character*255 type        ! type of data
      integer ma_type
      integer ntypes            ! no. of known data types
      parameter (ntypes = 3)
      character*8 types(ntypes) ! known data types
      integer ma_types(3)
      logical status            ! accumulates return codes
      integer ind               ! used to match types
      integer i, nvalues, max_values
      parameter (max_values = 10)
      double precision dvals(max_values)
      character*16 cvals(max_values)
      integer ivals(max_values)
      data types /'string', 'real', 'integer'/
c
c     Set entries in the rtdb from the input file
c
c     set <name> [<type = string>] <list of values>
c
      ma_types(1) = MT_BYTE
      ma_types(2) = MT_DBL
      ma_types(3) = MT_INT
c
      field   = ' '
      type    = types(1)
      ma_type = ma_types(1)
c
c     Crude checking
c
      call inp_set_field(0)
      if (inp_n_field() .lt. 3) goto 1000
      if (.not. inp_a(field)) call errquit
     $     ('input_set: very first input failed!', 0)
      if (.not. inp_compare(.false.,'set',field)) call errquit
     $     ('input_set: not a set directive?', 0)
      if (.not. inp_a(name)) call errquit
     $     ('input_set: reading name failed', 0)
c
      if (inp_n_field() .gt. 3) then
c
c     A type may be specified ... see if we can match something
c
         if (.not. inp_a(type)) call errquit
     $        ('input_set: reading type failed', 0)
         if (inp_match(ntypes, .false., type, types, ind)) then 
            type = types(ind)
         else
            type = types(1)     ! No match ... back up
            call inp_prev_field()
         endif
      endif
      if (inp_match(ntypes, .false., type, types, ind)) then 
         type = types(ind)
         ma_type = ma_types(ind)
      else
         call errquit('input_set: very wierd error matching type', 0)
      endif
c
c     How many values are provided ?  For character variables
c     we will always have statically allocatate buffers, but
c     for integer/real we could make some with ma?
c
      nvalues = inp_n_field() - inp_cur_field()
      if (nvalues .gt. max_values) call errquit
     $     ('input_set: too many values', nvalues)
c
      do i = 1, nvalues
         if (ma_type .eq. mt_int) then
            status = inp_i(ivals(i))
         else if (ma_type .eq. mt_dbl) then
            status = inp_f(dvals(i))
         else 
            status = inp_a(cvals(i))
         endif
         if (.not. status) call errquit
     $        ('input_set: failed reading value ', i)
      enddo
c
      if (ma_type .eq. mt_int) then
         status = rtdb_put(rtdb, name, mt_int, nvalues, ivals)
      else if (ma_type .eq. mt_dbl) then
         status = rtdb_put(rtdb, name, mt_dbl, nvalues, dvals)
      else
         status = rtdb_cput(rtdb, name, nvalues, cvals)
      endif
c
      if (.not. status) call errquit
     $     ('input_set: failed writing to rtdb', 0)
c
      return
c
 1000 call errquit
     $     ('input_set: set <name> [<type = string>] <list of values>',
     $     0)
c
      end
