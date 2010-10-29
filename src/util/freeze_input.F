c $Id$
      subroutine freeze_input(rtdb, module)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer rtdb
      character*(*) module      ! [input]
c
c     .  [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c     .          [virtual <integer nfzv default 0>]]
c
c     Freeze directive has already been detected.  Just process
c     this line
c     External Routines

      character*80 key
      character*20 test, space
      integer mlen
      integer n

      character*16 tag
      character*20 symbol
      character*20 element
      integer atn
      integer tagfrz(60)
      integer numcore(60)
      integer num
      integer num_elem             !loop counter
      logical find_atn
      logical is_atom

      num_elem = 1
      is_atom = .false.
c
      mlen = inp_strlen(module)
      space = 'frozen core'
c
c     Remove previous entries so empty directive forces freeze nothing.
c
      key = ' '
      write(key,'(a,'':'',a)') module(1:mlen),'frozen core'
      call rtdb_clean(rtdb,key)
      key = ' '
      write(key,'(a,'':'',a,'':'')') module(1:mlen),'frozen virtual'
      call rtdb_clean(rtdb,key)
c
 10   if (.not. inp_a(test)) return
c
      if (inp_compare(.false.,test,'core')) then
         space = 'frozen core'
      else if (inp_compare(.false.,test,'virtual')) then
         space = 'frozen virtual'
      else if (inp_compare(.false.,test,'atomic')) then
         key = ' '
         write(key,'(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms'
         if (.not. rtdb_put(rtdb, key, mt_log, 1, .true.))
     $        call errquit('freeze_input: failed writing rtbd',0,0)
       
 20      if(inp_a(tag)) then
           find_atn = geom_tag_to_element(tag, symbol, element, atn)
           if( find_atn ) then
             tagfrz(num_elem) = atn

             if(inp_i(num)) then
               numcore(num_elem) = num
               num_elem = num_elem + 1
               is_atom = .true.
               goto 20
            end if
          end if
         endif

         if(is_atom) then
           num_elem = num_elem - 1 
           space = 'number tag elements'
           key = ' '
           write(key, '(a,'':'',a,'':'',a)') module(1:mlen),
     $       space(1:inp_strlen(space)), 'freeze by atoms'
           call rtdb_clean(rtdb, key)
           if(.not. rtdb_put(rtdb, key, mt_int, 1, num_elem)) then
             call errquit('freeze_input: failed writing to rtdb', 0,0)
           endif
 
           space = 'tag freeze'
           key = ' '
           write(key,'(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms'
           call rtdb_clean(rtdb, key)
           if(.not. rtdb_put(rtdb, key, mt_int, num_elem, tagfrz)) then
             call errquit('freeze_input: failed writing to rtdb', 0,0)
           endif
           
           space = 'num freeze'
           key = ' '
           write(key,'(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms'
           call rtdb_clean(rtdb, key)
           if(.not. rtdb_put(rtdb, key, mt_int, num_elem, numcore)) then
             call errquit('freeze_input: failed writing to rtdb', 0,0)
           endif
         endif

         space = 'frozen core'         !restore default

      else 
         call inp_prev_field()
	 if (.not. inp_i(n)) call errquit
     $      ('freeze_input: looking for nfzc/v',0,0)
         key = ' '
         write(key, '(a,'':'',a)') module(1:mlen),
     $           space(1:inp_strlen(space))
         if (.not. rtdb_put(rtdb, key, mt_int, 1, n))
     $        call errquit('freeze_input: failed writing rtdb',0,0)
         space = 'frozen core'         ! restore default
      endif
c      
      goto 10
c
      end
      



