       subroutine tce_mrcc_input(rtdb)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "tce_mrcc.fh"
#include "tce_mrcc_diis.fh"
      integer rtdb
      character*1025 test
c      integer nref 
c     integer targetroot
      integer i,j
      character*8 t1
      character*3 t2
      character*4 ds
      integer nullc,icheckpoint
      logical nodezero
ckbn-7
      logical se4t
      logical se4t_offdiagonal
      logical no_aposteriori
      logical uncoup_pt3
      logical uncoup_se4t
      logical uncoup_not1t2

      integer mrccdebuglevel
      double precision dvalue
      integer iactel,iactorb,itotel
      logical lcas
c      integer iignore,isubsize
      integer isubsize
      integer istartmk
c
c ----------
c Initialization
c ----------
c
      nodezero = (ga_nodeid().eq.0)

      targetroot = 0
      iroottooverlap = 0
      nref = 0
         i = 0
        t1 = 'bwcc:ref'
      nullc = ICHAR('0')
      if(nodezero)write(LuOut,"('Entering MRCC input module',/)")
         j = 0
         lcas = .false.
c
c ----------
c Read input
c ----------
c
 10   if (.not. inp_read())
     1  call errquit('tce_mrcc_input: failed reading input',0,
     2  RTDB_ERR)
      if (.not. inp_a(test))
     1  call errquit('tce_mrcc_input: failed reading keyword',0,
     2  RTDB_ERR)
c
c ------------------------------------
c Read references and other parameters
c ------------------------------------
c
      if (inp_compare(.false.,test,'end')) then
          goto 20
      endif
      if (inp_compare(.false.,test,'nref')) then
        if (.not. inp_i(nref))
     1    call errquit('tce_mrcc_input: failed reading value',0,
     2    RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'root')) then
        if (.not. inp_i(targetroot))
     1    call errquit('tce_mrcc_input: failed reading value',1,
     2    RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'fermivacuum')) then
       if (.not.rtdb_put(rtdb,'mrcc:usescffermiv',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'save2e')) then
       if (.not.rtdb_put(rtdb,'mrcc:save2e',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'read2e')) then
       if (.not.rtdb_put(rtdb,'mrcc:read2e',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'fullheff')) then
       if (.not.rtdb_put(rtdb,'mrcc:fullheff',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'forcedegen')) then
       if (.not.rtdb_put(rtdb,'mrcc:forcedegen',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'improvetiling')) then
       if (.not.rtdb_put(rtdb,'mrcc:improvetiling',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'subgroups')) then
       if (.not.rtdb_put(rtdb,'mrcc:subgroups',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
        if (inp_i(isubsize)) then
        if (.not.rtdb_put(rtdb,'mrcc:subsize',mt_int,1,isubsize))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
        endif
          goto 10
      endif
      if (inp_compare(.false.,test,'startmk')) then
        if (.not. inp_i(istartmk))
     1    call errquit('tce_mrcc_input: failed reading value',1,
     2    RTDB_ERR)
       if (.not.rtdb_put(rtdb,'mrcc:istartmk',mt_int,1,istartmk))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'zignore')) then
        if (.not. inp_i(iignore))
     1    call errquit('tce_mrcc_input: failed reading value',1,
     2    RTDB_ERR)
       if (.not.rtdb_put(rtdb,'mrcc:zignore',mt_int,1,iignore))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
cxxx
      if (inp_compare(.false.,test,'mrccdebug')) then
          if (.not.rtdb_put(rtdb,'mrcc:debug',mt_log,1,.true.))
     +     call errquit('tce_mrcc_input: failed writing to mrccdebug',
     +       0,RTDB_ERR)
          goto 10
      endif

ckbn-7
      if (inp_compare(.false.,test,'se4t')) then
       if (.not.rtdb_put(rtdb,'mrcc:se4t',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
ckbn-9
      if (inp_compare(.false.,test,'no_aposteriori')) then
       if (.not.rtdb_put(rtdb,'mrcc:no_aposteriori',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif

ckbn-9
      if (inp_compare(.false.,test,'se4t_offdiagonal')) then
       if (.not.rtdb_put(rtdb,'mrcc:se4t_offdiagonal',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
ckbn-9
      if (inp_compare(.false.,test,'uncoup_pt3')) then
       if (.not.rtdb_put(rtdb,'mrcc:uncoup_pt3',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
ckbn-9
      if (inp_compare(.false.,test,'uncoup_se4t')) then
       if (.not.rtdb_put(rtdb,'mrcc:uncoup_se4t',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
ckbn-9
      if (inp_compare(.false.,test,'uncoup_not1t2')) then
       if (.not.rtdb_put(rtdb,'mrcc:uncoup_not1t2',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'usesamefock_it')) then
       if (.not.rtdb_put(rtdb,'mrcc:usesamefock_it',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'usesamefock_nonit')) then
       if (.not.rtdb_put(rtdb,'mrcc:usesamefock_nonit',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif

ckbn roottooverlap
      if (inp_compare(.false.,test,'roottooverlap')) then
        if (.not. inp_i(iroottooverlap))
     +    call errquit('tce_mrcc_input: failed reading value',1,
     +    RTDB_ERR)
c       if (.not.rtdb_put(rtdb,'mrcc:roottooverlap',mt_int,1,
c     +    iroottooverlap))
c     +   call errquit('tce_mrcc_input: failed writing to rtdb',0,
c     +   RTDB_ERR)
          goto 10
      endif
      write(*,*) "roottooverlap", iroottooverlap

ckbn ignore complex root has to be used with caution 
      if (inp_compare(.false.,test,'ignorecomplex')) then
       if (.not.rtdb_put(rtdb,'mrcc:ignorecomplex',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif

ckbn  savet when converged
      if (inp_compare(.false.,test,'savet')) then
       if (.not.rtdb_put(rtdb,'mrcc:savet',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif


ckbn zora
      if (inp_compare(.false.,test,'zorafromdft')) then
      if (.not.rtdb_put(rtdb,'mrcc:zorafromdft',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      



cxxx
      if (inp_compare(.false.,test,'diistype')) then
        if (.not. inp_i(idiis))
     1    call errquit('tce_mrcc_input: failed reading value',1,
     2    RTDB_ERR)
       if (.not.rtdb_put(rtdb,'mrcc:diistype',mt_int,1,idiis))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'bwcorr')) then
       if (.not.rtdb_put(rtdb,'mrcc:bwcorr',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'readt')) then
       if (.not.rtdb_put(rtdb,'mrcc:readt',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'useeaf2el')) then
       if (.not.rtdb_put(rtdb,'mrcc:useeaf2el',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          goto 10
      endif
      if (inp_compare(.false.,test,'checkpoints')) then
       if (.not.rtdb_put(rtdb,'mrcc:checkpoints',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
        if (inp_i(icheckpoint)) then
        if (.not.rtdb_put(rtdb,'mrcc:icheckpoint',mt_int,1,icheckpoint))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
        endif
          goto 10
      endif
      if (inp_compare(.false.,test,'rootmucilage')) then
          j = 0
777       if(inp_f(dvalue)) then
          j = j + 1
          write(ds,"(I3.3)")j
          write(6,*)ds,dvalue
          drootmuc(j) = dvalue
       if (.not.rtdb_put(rtdb,'mrcc:rootmuc'//ds,mt_dbl,1,dvalue))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
       else
         goto 10
       endif
         goto 777
      endif
c      if (inp_compare(.false.,test,'mrccdebug')) then
c        if (.not. inp_i(mrccdebuglevel))
c     1    mrccdebuglevel=1
c          if (.not.rtdb_put(rtdb,'mrcc:debug',mt_log,1,.true.))
c     1    call ERRQUIT('tce_mrcc_input',3,RTDB_ERR)
c          goto 10
c      endif
      if (inp_compare(.false.,test,'cas')) then
        if (.not. inp_i(iactel))
     1    call ERRQUIT('tce_mrcc_input',3,RTDB_ERR)
        if (.not. inp_i(iactorb))
     1    call ERRQUIT('tce_mrcc_input',3,RTDB_ERR)
       if (.not.rtdb_put(rtdb,'mrcc:cas',mt_log,1,.true.))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
       if (.not.rtdb_put(rtdb,'mrcc:iactel',mt_int,1,iactel))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
       if (.not.rtdb_put(rtdb,'mrcc:iactorb',mt_int,1,iactorb))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
          lcas = .true.
          goto 10
      endif
      i = i + 1
      if(nodezero)write(LuOut,"(T3,'Reference ',I3,1X,':',A50)")i,test
      write(t2,"(I3.3)")i
      if(lcas) call
     1ERRQUIT('CAS and ref input cannot be defined together',4,RTDB_ERR)
      if (.not.rtdb_cput(rtdb,t1//t2,1,test))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
      goto 10
c
c ------
c Return
c ------
c
 20   continue 

      if (.not.rtdb_put(rtdb,'bwcc:nref',mt_int,1,nref))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)

ckbn test not root and roottooverlap are "on" together
      write(*,*) "roottooverlap", iroottooverlap
      if(targetroot.eq.0  .and. iroottooverlap.eq.0) then
        write(*,'(A)') "Use root or roottooverlap option"
        call errquit('Rerun with root or roottooverlap ',INPUT_ERR,30)
      endif      

       if(iroottooverlap .gt.0) then
       if (.not.rtdb_put(rtdb,'mrcc:roottooverlap',mt_int,1,
     +    iroottooverlap))
     +   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     +   RTDB_ERR)
       endif

      if(targetroot.gt.0) then
       if (.not.rtdb_put(rtdb,'bwcc:targetroot',mt_int,1,targetroot))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
      endif
      if(j.gt.0) then
      if (.not.rtdb_put(rtdb,'mrcc:rootmuc',mt_int,1,j))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)
      endif
       return
       end


