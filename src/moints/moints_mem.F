c
c Return the maximum number of operator matrices
c that will fit in memory. (lots of assumptions)
c
c
c
       integer function moints_numoper( basis, algo, nbf, nocc,
     $                                  nvir, blen )
C$Id: moints_mem.F,v 1.2 1995-02-02 23:18:48 d3g681 Exp $
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
c
c
       integer MSG_MIN_LOCAL_MEM
       parameter(MSG_MIN_LOCAL_MEM=1971)
c
       integer basis, nbf, nocc, nvir, blen
       character*(*) algo
       integer max2e, mem2, maxbfsh, bsize, mem_req
       integer mem_avail, mem_min_local, noper_node
       integer mem_int, mem_mo, mem_blk, mem_buf1, mem_buf2
c
c
c
       call int_mem_2e4c(max2e, mem2)
       if (.not.bas_nbf_cn_max(basis,maxbfsh))
     $     call errquit('moints_numoper: cannot get basis info',0)
c
c Get minimum local memory that any
c
       mem_avail = ma_inquire_avail(MT_DBL)
       mem_min_local = mem_avail
       call ga_sync()
       call ga_igop(MSG_MIN_LOCAL_MEM,mem_min_local,1,'min')
c
c
c
       bsize = max(blen,maxbfsh)
       mem_int = max2e + mem2
       mem_mo = nbf*nbf
       mem_blk = maxbfsh*maxbfsh*bsize*bsize
       mem_blk = mem_blk + max(mem_blk,(nocc*nocc))
       mem_buf1 = max((bsize*maxbfsh*maxbfsh*nocc),(maxbfsh*nbf))
       if (algo.eq.'sixfold') then
         mem_buf2 = maxbfsh*maxbfsh*nbf*nocc
       else
         mem_buf2 = maxbfsh*maxbfsh*(nbf*nocc + nvir)
       endif
       mem_req = mem_int + mem_mo + mem_blk + mem_buf1 + mem_buf2
c
c 10% safety margin
c
       mem_req = mem_req*1.1
c
c
c
       if (mem_req.gt.mem_min_local)
     $   call errquit(
     $   'not enough local memory to do transformation',mem_req)
       noper_node = (mem_min_local-mem_req)/(nbf*nbf)
       moints_numoper = ga_nnodes()*noper_node

#ifdef MOINTS_STATS
       if (ga_nodeid().eq.0) write(6,933) mem_int,
     $                                    mem_mo,
     $                                    mem_blk,
     $                                    mem_buf1,
     $                                    mem_buf2,
     $                                    mem_req,
     $                                    mem_min_local,
     $                                    noper_node,
     $                                    moints_numoper
 933   format(/,24x,'Memory Summary',
     $        /,24x,14('-'),
     $        /,10x,'AO integral:',19x,i10,
     $        /,10x,'MO Coefficients:',15x,i10,
     $        /,10x,'Block:',25x,i10,
     $        /,10x,'Buffer 1:',22x,i10,
     $        /,10x,'Buffer 2:',22x,i10,
     $        /,10x,'Local memory required:',9x,i10,
     $        /,10x,'Local memory available:',8x,i10,
     $        /,10x,'Operator matrices per node:',9x,i5,
     $        /,10x,'Max number of operators:',12x,i5)
#endif
c
c
c
       return
       end


