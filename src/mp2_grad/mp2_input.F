      subroutine mp2_input(rtdb)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
c
c     MP2
c     .  [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c     .          [virtual <integer nfzv default 0>]]
c     .  [AOTOL2E <real aotol2e default 1d-9>]
c     .  [MOTOL2E <real motol2e default 1d-9>]
c     .  [PRINT]
c     .  [NOPRINT]
c     .  [VECTORS <string filename default scf-output-vectors> \
c     .    [swap [(alpha||beta)] <integer pair-list>] ]
c     .  [RIAPPROX <string riapprox default V>]
c     .  [FILE3C <string filename default $file_prefix$.mo3cint">]
c     END
c
      character*20 test, riapprox
      character*255 rifile
      double precision aotol2e, motol2e
c
 10   if (.not. inp_read()) 
     $     call errquit('mp2_input: failed reading input',0)
      if (.not. inp_a(test)) 
     $     call errquit('mp2_input: failed reading keyword',0)
c
      if (inp_compare(.false.,test,'freeze')) then
         call freeze_input(rtdb,'mp2')
      else if (inp_compare(.false.,test,'aotol2e')) then
         if (.not. inp_f(aotol2e)) call errquit
     $        ('mp2_input: failed reading aotol2e',0)
         if (.not. rtdb_put(rtdb, 'mp2:aotol2e', mt_dbl, 1, aotol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'motol2e')) then
         if (.not. inp_f(aotol2e)) call errquit
     $        ('mp2_input: failed reading motol2e',0)
         if (.not. rtdb_put(rtdb, 'mp2:motol2e', mt_dbl, 1, motol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'print') .or.
     $        inp_compare(.false.,test,'noprint')) then
         call util_print_input(rtdb,'mp2')
      else if (inp_compare(.false.,test,'vectors')) then
         call vectors_input(rtdb,'mp2')
      else if (inp_compare(.false.,test,'riapprox')) then
         if (.not. inp_a(riapprox)) call errquit
     $        ('mp2_input: failed reading ri approximation',0)
         if (.not. rtdb_cput(rtdb, 'mp2:ri approximation', 1, riapprox))
     $        call errquit('mp2_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'file3c')) then
         if (.not. inp_a(rifile)) call errquit
     $        ('mp2_input: failed reading ri 3c file',0)
         if (.not. rtdb_cput(rtdb,'mp2:mo 3-center integral file', 1, 
     $        rifile)) 
     $        call errquit('mp2_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'end')) then
         goto 20
      else
         call errquit('mp2_input: unknown directive',0)
      endif
      goto 10
c
 20   return
      end
      subroutine freeze_input(rtdb, module)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      character*(*) module      ! [input]
c
c     .  [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c     .          [virtual <integer nfzv default 0>]]
c
c     Freeze directive has already been detected.  Just process
c     this line
c
      character*80 key
      character*20 test, space
      integer mlen
      integer maxlist
      parameter (maxlist=1024)
      integer list(maxlist)
      integer n, i
      logical status
c
      mlen = inp_strlen(module)
      space = 'core'
c
 10   if (.not. inp_a(test)) return
c
      if (inp_compare(.false.,test,'core')) then
         space = 'core'
      else if (inp_compare(.false.,test,'virtual')) then
         space = 'virtual'
      else if (inp_compare(.false.,test,'atomic')) then
         key = ' '
         write(key,'(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms'
         write(6,*) ' THE KEY IS ', key
         if (.not. rtdb_put(rtdb, key, mt_log, 1, .true.))
     $        call errquit('freeze_input: failed writing rtbd',0)
         space = 'core'         ! restore default
      else 
         call inp_prev_field()
         status = inp_ilist(maxlist, list, n)
         if (n .eq. 0) then
            call errquit('freeze_input: failed reading list',0)
         else if (n.le.maxlist) then
            key = ' '
            write(key, '(a,'':'',a)') module(1:mlen),
     $           space(1:inp_strlen(space))
            write(6,*) ' KEY ', key, (list(i),i=1,n)
            if (.not. rtdb_put(rtdb, key, mt_int, n, list))
     $           call errquit('freeze_input: failed writing rtdb',0)
         else
            call errquit('freeze_input: too many frozen',n)
         endif
      endif
c      
      goto 10
c
      end
      
