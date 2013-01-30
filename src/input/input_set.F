      subroutine input_set(rtdb)
C$Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb                   !< [Input] The RTDB handle
      character*(rtdb_max_key) field ! input buffer
      character*(rtdb_max_key) name  ! name of variable being set
      integer ma_type
      integer ntypes                 ! no. of known data types
      parameter (ntypes = 5)
      character*8 type, types(ntypes) ! known data types
      integer ma_types(5)
      logical status            ! accumulates return codes
      integer ind               ! used to match types
      integer itest
      double precision dtest
      integer i, nvalues, max_values, H_Values, I_Values, Space_Avail
      logical input_is_ilist
      parameter (max_values = 100)
      character*(rtdb_max_file) cvals(max_values) ! may contain file names
      character*8 test
      data types /'string', 'real', 'integer', 'logical', 'double'/
c
c     Set entries in the rtdb from the input file
c
c     set <name> [<type = automatic>] <list of values>
c
      ma_types(1) = MT_BYTE
      ma_types(2) = MT_DBL
      ma_types(3) = MT_INT
      ma_types(4) = MT_LOG
      ma_types(5) = MT_DBL
c
      field   = ' '
      type    = ' '             ! This setting tested below
      ma_type = -1
c
c     Crude checking
c
      call inp_set_field(0)
      if (inp_n_field() .lt. 3) goto 1000
      if (.not. inp_a(field)) call errquit
     $     ('input_set: very first input failed!', 0, INPUT_ERR)
      if (.not. inp_compare(.false.,'set',field)) call errquit
     $     ('input_set: not a set directive?', 0, INPUT_ERR)
      if (.not. inp_a(name)) call errquit
     $     ('input_set: reading name failed', 0, INPUT_ERR)
c
      if (inp_n_field() .gt. 3) then
c
c     A type may be specified ... see if we can match something
c
         if (.not. inp_a(field)) call errquit
     $        ('input_set: reading type failed', 0, INPUT_ERR)
         if (inp_match(ntypes, .false., field, types, ind)) then 
            type = types(ind)
         else
            call inp_prev_field() ! No match, backup
         endif
      endif
      if (type .eq. ' ') then
c
c     No type was specified.  Attempt to infer the type
c     from reading the first datum.
c
c     Examine in sequence integer, real, logical, integerlist, string.
c
         if (inp_i(itest)) then
            type = types(3)
         else if (inp_f(dtest)) then
            type = types(5)
         else if (inp_a(field)) then
            if (inp_compare(.false.,field,'true') .or.
     $           inp_compare(.false.,field,'.true.') .or.
     $           inp_compare(.false.,field,'t') .or.
     $           inp_compare(.false.,field,'false') .or.
     $           inp_compare(.false.,field,'.false.') .or.
     $           inp_compare(.false.,field,'f')) then
               type = types(4)
            else if (input_is_ilist(field)) then
               type = types(3)
            else
               type = types(1)
            endif
         else
            goto 1000
         endif
         call inp_prev_field()
      endif
c            
      if (inp_match(ntypes, .false., type, types, ind)) then 
         type = types(ind)
         ma_type = ma_types(ind)
      else
         call errquit('input_set: very wierd error matching type', 0,
     &       INPUT_ERR)
      endif
C
C     How much data is there to read? (For integers there may be more,
C     since each field may actually be a range.
C
      nvalues = inp_n_field() - inp_cur_field()
*      write(LuOut,*) name(1:inp_strlen(name)),' is ', type,' ',nvalues
C
C     Read the data according to the specified type.  For integers, we
C     allocate an obnoxiously large buffer in the hopes of ensuring
C     sufficient space to expand any ranges.  For reals and logicals
C     we can just dynamically allocate NValues.  Characters must be
C     handled statically.
C     
      if ( ma_type .eq. mt_int) then
c
         space_avail = min( 100000, ma_inquire_avail(ma_type) )
         if ( .not. ma_push_get( ma_type, space_avail, 'input scratch',
     $      h_values, i_values) ) call errquit(
     $      'input_set: unable to allocate input scratch space',
     $      space_avail, MA_ERR)
c
         status = inp_ilist(space_avail, int_mb(i_values), nvalues)
         if ( .not. status) call errquit(
     $      'input_set: failed reading integer list ', 0, INPUT_ERR)
c
      elseif (ma_type .eq. mt_dbl .or. ma_type .eq. mt_log) then
         if ( .not. ma_push_get( ma_type, nvalues, 'input scratch',
     $      h_values, i_values) ) call errquit(
     $      'input_set: unable to allocate input scratch space',
     $      nvalues, MA_ERR)
c
         do i = 1, nvalues
            if (ma_type .eq. mt_dbl) then
c              real
               status = inp_f( dbl_MB(i_Values + i - 1) )
            else if (ma_type .eq. mt_log) then
c              logical
               status = inp_a(test)
               if (inp_compare(.false.,test,'true') .or.
     $            inp_compare(.false.,test,'.true.') .or.
     $            inp_compare(.false.,test,'t')) then
                  Log_MB( I_Values + i - 1) = .true.
               else if (inp_compare(.false.,test,'false') .or.
     $               inp_compare(.false.,test,'.false.') .or.
     $               inp_compare(.false.,test,'f')) then
                  Log_MB( I_Values + i - 1) = .false.
               else
                  call errquit('input_set: error parsing logical', i,
     &       INPUT_ERR)
               endif
            endif
            if (.not. status) call errquit
     $         ('input_set: failed reading value ', i, INPUT_ERR)
         enddo
c
      else                      ! character strings
         if ( nvalues .gt. max_values) call errquit(
     $      'input_set: too many values', nvalues, INPUT_ERR)
c
         do i = 1, nvalues
            status = inp_a(cvals(i))
            if (.not. status) call errquit
     $         ('input_set: failed reading value ', i, INPUT_ERR)
         enddo
      endif
c
      if (ma_type .eq. mt_int) then
         status = rtdb_put(rtdb, name, mt_int, nvalues,
     $      int_mb(i_values))
      else if (ma_type .eq. mt_dbl) then
         status = rtdb_put(rtdb, name, mt_dbl, nvalues,
     $      dbl_mb(i_values))
      else if (ma_type .eq. mt_log) then
         status = rtdb_put(rtdb, name, mt_log, nvalues,
     $      log_mb(i_values))
      else
         status = rtdb_cput(rtdb, name, nvalues, cvals)
      endif
c
      if (.not. status) call errquit
     $     ('input_set: failed writing to rtdb', 0, RTDB_ERR)
c
      if ( ma_type .eq. mt_int .or. ma_type .eq. mt_dbl
     $   .or. ma_type .eq. mt_log) then
         if ( .not. ma_pop_stack(h_values) ) call errquit(
     $      'input_set: unable to free scratch space', h_values, MA_ERR)
      endif
C
      return
c
 1000 call errquit
     $     ('input_set: set <name> [<type = automatic>] <list>', 0,
     &       INPUT_ERR)
c
      end
      subroutine input_unset(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "stdio.fh"
      integer rtdb
c     
      character*255 test, name
      integer ltest
      logical status
c     
      if (.not. inp_a(test)) call errquit
     $     ('input_unset: first inp_a failed?',0, INPUT_ERR)
      if (.not. inp_compare(.false.,'unset',test))
     $     call errquit('input_unset: not unset directive',0, INPUT_ERR)
c     
      if (.not. inp_a(test)) call errquit
     $     ('input_unset: name missing',0, INPUT_ERR)
c     
c     Check for wildcard
c     
      ltest = inp_strlen(test)
      if (test(ltest:ltest) .eq. '*') then
	 ltest = ltest - 1
	 status = rtdb_first(rtdb, name)
 10      if (status) then
	    if (inp_compare(.true.,test(1:ltest),name(1:ltest))) then
               write(LuOut,*) ' unset: ', name(1:inp_strlen(name))
	       if (.not. rtdb_delete(rtdb,name)) call errquit
     $              ('input_unset: failed deleting existing entry',0,
     &       RTDB_ERR)
	    endif
            status = rtdb_next(rtdb, name)
            goto 10
         endif
      else if (.not. rtdb_delete(rtdb,test)) then
         write(LuOut,*) ' unset: warning: ',test(1:inp_strlen(test)),
     $        ' is not in the database '
      endif
c     
      end
      logical function input_is_ilist(a)
      implicit none
#include "inp.fh"
      character*(*) a
c
      integer lena, i, j
      character*1 ok(11)
      data ok/'0','1','2','3','4','5','6','7','8','9',':'/
c
      lena = inp_strlen(a)
      do i = 1, lena
         do j = 1, 11
            if (a(i:i).eq.ok(j)) goto 10
         enddo
         input_is_ilist = .false.
         return
 10      continue
      enddo
c
      input_is_ilist = .true.
c
      end
         
         
