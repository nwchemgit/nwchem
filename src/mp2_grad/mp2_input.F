      subroutine mp2_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
c
c     MP2
c     .  [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c     .          [virtual <integer nfzv default 0>]]
c     .  [AOTOL2E <real aotol2e default 1d-9>]
c     .  [MOTOL2E <real motol2e default 1d-9>] ... DEFUNCT
c     .  [TIGHT]
c     .  [PRINT]
c     .  [NOPRINT]
c     .  [VECTORS <string filename default scf-output-vectors> \
c     .    [swap [(alpha||beta)] <integer pair-list>] ]
c     .  [RIAPPROX <string riapprox default V>]
c     .  [FILE3C <string filename default "$file_prefix$.mo3cint">]
c     .  [SCRATCHDISK <integer>]
c     .  [FSS <real same spin scaling factor for SCS-MP2 default 1/3>]
c     .  [FOS <real opposite spin scaling factor for SCS-MP2 default 6/5>]
c     END
c
      character*20 test, riapprox
      character*255 rifile
      double precision aotol2e, motol2e, scfthresh, cphfthresh, backtol
      double precision fss,fos
      logical status
      integer scrdisk
c
 10   if (.not. inp_read()) 
     $     call errquit('mp2_input: failed reading input',0, INPUT_ERR)
      if (.not. inp_a(test)) 
     $     call errquit('mp2_input: failed reading keyword',0,
     &       INPUT_ERR)
c
      if (inp_compare(.false.,test,'freeze')) then
         call freeze_input(rtdb,'mp2')
      else if (inp_compare(.false.,test,'tight')) then
         aotol2e = 1e-11
         motol2e = 1e-11
         backtol = 1e-11
         if (.not. rtdb_put(rtdb, 'mp2:aotol2e', mt_dbl, 1, aotol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'mp2:aotol2e fock', mt_dbl, 1, 
     $        motol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'mp2:backtol', mt_dbl, 1, backtol))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
         if (.not. rtdb_get(rtdb,'scf:thresh',mt_dbl,1,scfthresh))
     $        scfthresh = 1.0d-4
         if (scfthresh .gt. 1.0d-8) then
            scfthresh = 1.0d-8
            if (.not. rtdb_put(rtdb,'scf:thresh',mt_dbl,1,scfthresh))
     $           call errquit('mp2_input: failed setting scf:thresh',0,
     &       RTDB_ERR)
            status = rtdb_delete(rtdb,'scf:converged')
         endif
         if (.not. rtdb_get(rtdb,'cphf:thresh',mt_dbl,1,cphfthresh))
     $        cphfthresh = 1.0d-4
         if (cphfthresh .gt. 1.0d-6) then
            cphfthresh = 1.0d-6
            if (.not. rtdb_put(rtdb,'cphf:thresh',mt_dbl,1,cphfthresh))
     $           call errquit('mp2_input: failed setting cphf:thresh',0,
     &       RTDB_ERR)
         endif
      else if (inp_compare(.false.,test,'aotol2e')) then
         if (.not. inp_f(aotol2e)) call errquit
     $        ('mp2_input: failed reading aotol2e',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'mp2:aotol2e', mt_dbl, 1, aotol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'mp2:aotol2e fock',mt_dbl,1,aotol2e))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (inp_compare(.false.,test,'motol2e')) then
         write(6,*) ' mp2_input: "MOTOL2E" is no longer used.',
     $        ' Use "TIGHT" instead.'
      else if (inp_compare(.false.,test,'print') .or.
     $        inp_compare(.false.,test,'noprint')) then
         call util_print_input(rtdb,'mp2')
      else if (inp_compare(.false.,test,'vectors')) then
         call vectors_input(rtdb,'mp2')
      else if (inp_compare(.false.,test,'riapprox')) then
         if (.not. inp_a(riapprox)) call errquit
     $        ('mp2_input: failed reading ri approximation',0,
     &       INPUT_ERR)
         if (.not. rtdb_cput(rtdb, 'mp2:ri approximation', 1, riapprox))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (inp_compare(.false.,test,'file3c')) then
         if (.not. inp_a(rifile)) call errquit
     $        ('mp2_input: failed reading ri 3c file',0, INPUT_ERR)
         if (.not. rtdb_cput(rtdb,'mp2:mo 3-center integral file', 1, 
     $        rifile)) 
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (inp_compare(.false.,test,'scratchdisk')) then
         if (.not. inp_i(scrdisk)) call errquit
     $        ('mp2_input: failed reading scrdisk',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'mp2:scratchdisk', mt_dbl, 1,
     $        dble(scrdisk))) 
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (inp_compare(.false.,test,'fss')) then
         if (.not. inp_f(fss)) call errquit
     $        ('mp2_input: failed reading fss',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (inp_compare(.false.,test,'fos')) then
         if (.not. inp_f(fos)) call errquit
     $        ('mp2_input: failed reading fos',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (inp_compare(.false.,test,'scs')) then
         if (.not. rtdb_put(rtdb, 'mp2:scs', mt_log, 1, .true.))
     $        call errquit('mp2_input: failed writing to rtdb',0,
     &       RTDB_ERR)
      else if (inp_compare(.false.,test,'end')) then
         goto 20
      else
         call errquit('mp2_input: unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
 20   return
      end
