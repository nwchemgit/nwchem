      subroutine mp2_input(rtdb)
*
* $Id: mp2_input.F,v 1.5 1998-08-21 00:41:29 d3g681 Exp $
*
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
c     .  [TIGHT]
c     .  [PRINT]
c     .  [NOPRINT]
c     .  [VECTORS <string filename default scf-output-vectors> \
c     .    [swap [(alpha||beta)] <integer pair-list>] ]
c     .  [RIAPPROX <string riapprox default V>]
c     .  [FILE3C <string filename default "$file_prefix$.mo3cint">]
c     END
c
      character*20 test, riapprox
      character*255 rifile
      double precision aotol2e, motol2e, scfthresh, cphfthresh, backtol
      logical status
c
 10   if (.not. inp_read()) 
     $     call errquit('mp2_input: failed reading input',0)
      if (.not. inp_a(test)) 
     $     call errquit('mp2_input: failed reading keyword',0)
c
      if (inp_compare(.false.,test,'freeze')) then
         call freeze_input(rtdb,'mp2')
      else if (inp_compare(.false.,test,'tight')) then
         aotol2e = 1e-11
         motol2e = 1e-11
         backtol = 1e-11
         if (.not. rtdb_put(rtdb, 'mp2:aotol2e', mt_dbl, 1, aotol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0)
         if (.not. rtdb_put(rtdb, 'mp2:motol2e', mt_dbl, 1, motol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0)
         if (.not. rtdb_put(rtdb, 'mp2:backtol', mt_dbl, 1, backtol))
     $        call errquit('mp2_input: failed writing to rtdb',0)
         if (.not. rtdb_get(rtdb,'scf:thresh',mt_dbl,1,scfthresh))
     $        scfthresh = 1.0d-4
         if (scfthresh .gt. 1.0d-8) then
            scfthresh = 1.0d-8
            if (.not. rtdb_put(rtdb,'scf:thresh',mt_dbl,1,scfthresh))
     $           call errquit('mp2_input: failed setting scf:thresh',0)
            status = rtdb_delete(rtdb,'scf:converged')
         endif
         if (.not. rtdb_get(rtdb,'cphf:thresh',mt_dbl,1,cphfthresh))
     $        cphfthresh = 1.0d-4
         if (cphfthresh .gt. 1.0d-6) then
            cphfthresh = 1.0d-6
            if (.not. rtdb_put(rtdb,'cphf:thresh',mt_dbl,1,cphfthresh))
     $           call errquit('mp2_input: failed setting cphf:thresh',0)
         endif
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
      integer n
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
     $        call errquit('freeze_input: failed writing rtbd',0)
         space = 'frozen core'         ! restore default
      else 
         call inp_prev_field()
	 if (.not. inp_i(n)) call errquit
     $      ('freeze_input: looking for nfzc/v',0)
         key = ' '
         write(key, '(a,'':'',a)') module(1:mlen),
     $           space(1:inp_strlen(space))
         if (.not. rtdb_put(rtdb, key, mt_int, 1, n))
     $        call errquit('freeze_input: failed writing rtdb',0)
         space = 'frozen core'         ! restore default
      endif
c      
      goto 10
c
      end
      
