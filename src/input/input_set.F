      subroutine input_set(rtdb)
C$Id: input_set.F,v 1.10 1996-07-18 07:17:00 d3g681 Exp $
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
      parameter (ntypes = 4)
      character*8 types(ntypes) ! known data types
      integer ma_types(4)
      logical status            ! accumulates return codes
      integer ind               ! used to match types
      integer i, nvalues, max_values, H_Values, I_Values, Space_Avail
      parameter (max_values = 100)
      character*32     cvals(max_values)
      character*8 test
      data types /'string', 'real', 'integer', 'logical'/
c
c     Set entries in the rtdb from the input file
c
c     set <name> [<type = string>] <list of values>
c
      ma_types(1) = MT_BYTE
      ma_types(2) = MT_DBL
      ma_types(3) = MT_INT
      ma_types(4) = MT_LOG
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
C
C     How much data is there to read? (For integers there may be more,
C     since each field may actually be a range.
C
      nvalues = inp_n_field() - inp_cur_field()
C
C     Read the data according to the specified type.  For integers, we
C     allocate an obnoxiously large buffer in the hopes of insuring
C     sufficient space to expand any ranges.  For reals and logicals
C     we can just dynamically allocate NValues.  Characters must be
C     handled statically.
C     
      If ( Ma_Type .eq. MT_Int) then
C
         Space_Avail = Min( 100000, MA_Inquire_Avail(ma_type) )
         If ( .NOT. MA_Push_Get( ma_type, Space_Avail, 'Input scratch',
     $      H_Values, I_Values) ) Call ErrQuit(
     $      'input_set: unable to allocate input scratch space',
     $      Space_avail)
C
         Status = Inp_IList(Space_Avail, Int_MB(I_Values), NValues)
         If ( .NOT. Status) Call ErrQuit(
     $      'input_set: failed reading integer list ', 0)
C
      ElseIf (MA_Type .eq. MT_Dbl .OR. MA_Type .eq. MT_Log) then
         If ( .NOT. MA_Push_Get( MA_Type, NValues, 'Input scratch',
     $      H_Values, I_Values) ) Call ErrQuit(
     $      'input_set: unable to allocate input scratch space',
     $      NValues)
C
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
                  call errquit('input_set: error parsing logical', i)
               endif
            endif
            if (.not. status) call errquit
     $         ('input_set: failed reading value ', i)
         enddo
C
      Else                      ! Character strings
         If ( NValues .gt. Max_Values) Call ErrQuit(
     $      'input_set: too many values', NValues)
C
         do i = 1, nvalues
            status = inp_a(cvals(i))
            if (.not. status) call errquit
     $         ('input_set: failed reading value ', i)
         enddo
      EndIf
c
      if (ma_type .eq. mt_int) then
         status = rtdb_put(rtdb, name, mt_int, nvalues,
     $      Int_MB(I_Values))
      else if (ma_type .eq. mt_dbl) then
         status = rtdb_put(rtdb, name, mt_dbl, nvalues,
     $      Dbl_MB(I_Values))
      else if (ma_type .eq. mt_log) then
         status = rtdb_put(rtdb, name, mt_log, nvalues,
     $      Log_MB(I_Values))
      else
         status = rtdb_cput(rtdb, name, nvalues, cvals)
      endif
c
      if (.not. status) call errquit
     $     ('input_set: failed writing to rtdb', 0)
c
      If ( Ma_type .eq. MT_Int .OR. Ma_Type .eq. MT_Dbl
     $   .OR. Ma_Type .eq. MT_Log) then
         if ( .NOT. MA_Pop_Stack(H_Values) ) Call ErrQuit(
     $      'input_set: unable to free scratch space', H_Values)
      EndIf
C
      return
c
 1000 call errquit
     $     ('input_set: set <name> [<type = string>] <list of values>',
     $     0)
c
      end
      subroutine input_unset(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
      integer rtdb
c     
      character*255 test, name
      integer ltest
      logical status
c     
      if (.not. inp_a(test)) call errquit
     $     ('input_unset: first inp_a failed?',0)
      if (.not. inp_compare(.false.,'unset',test))
     $     call errquit('input_unset: not unset directive',0)
c     
      if (.not. inp_a(test)) call errquit
     $     ('input_unset: name missing',0)
c     
c     Check for wildcard
c     
      ltest = inp_strlen(test)
      if (test(ltest:ltest) .eq. '*') then
	 ltest = ltest - 1
	 status = rtdb_first(rtdb, name)
 10      if (status) then
	    if (inp_compare(.true.,test(1:ltest),name(1:ltest))) then
               write(6,*) ' unset: ', name(1:inp_strlen(name))
	       if (.not. rtdb_delete(rtdb,name)) call errquit
     $              ('input_unset: failed deleting existing entry',0)
	    endif
            status = rtdb_next(rtdb, name)
            goto 10
         endif
      else if (.not. rtdb_delete(rtdb,test)) then
         write(6,*) ' unset: warning: ',test(1:inp_strlen(test)),
     $        ' is not in the database '
      endif
c     
      end
