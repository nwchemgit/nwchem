c
c Check if the calculation will fit in memory
c
       logical function ccsd_memory_fit( basis, algo, 
     $            nbf, nocc, nvir, blen, mxvec, do_trip)
C$Id: ccsd_mem.F,v 2.1 2003-11-15 00:31:26 bert Exp $
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "ccsd_len.fh"
c
c
       integer MSG_MIN_LOCAL_MEM
       parameter(MSG_MIN_LOCAL_MEM=1971)
c
       integer basis, nbf, nocc, nvir, blen
       character*(*) algo
       logical do_trip
       integer max2e, mem2, maxbfsh, bsize, ma_req, nsh
       integer mem_avail, mem_min_local, noper_node, noper
       integer mem_int, mem_mo, mem_blk, mem_buf1, mem_buf2
       integer gmem_max, noccc, nvirc
       integer ma_req_ccsd, ga_req_ccsd
       integer ma_req_trip, ga_req_trip
       integer ga_req, mxvec
       integer moints_trp_mem
       external moints_trp_mem
c
       call int_mem_2e4c(max2e, mem2)
       if (.not.bas_nbf_cn_max(basis,maxbfsh))
     $     call errquit('ccsd_numoper: cannot get basis info',0,
     &       BASIS_ERR)
c
c MA memory for integral transformation
c
       bsize = max(blen,maxbfsh)
       mem_int = max2e + mem2
       mem_mo =  2*nbf*nbf + 10*nbf
       mem_blk = maxbfsh*maxbfsh*bsize*bsize
       mem_blk = mem_blk + max(mem_blk,(nocc*nocc))
       mem_buf1 = max((bsize*maxbfsh*maxbfsh*nocc),(maxbfsh*nbf))
       if (algo.eq.'sixfold') then
         mem_buf2 = maxbfsh*maxbfsh*nbf*nocc
       else
         mem_buf2 = maxbfsh*maxbfsh*(nbf+nvir+nocc)*nocc
       endif
       ma_req = mem_int + mem_mo + mem_blk + mem_buf1 + mem_buf2
c
c GA memory for integral transformation              
c
       noper = (nocc*(nocc+1))/2
       ga_req = 2*noper*nbf*nbf
c
c MA memory for CCSD calculation
c
       ma_req_ccsd = 3*lnov+2*lnoo+2*lnvv+4*nvir*nvir
       ma_req_ccsd = ma_req_ccsd + mxvec*(mxvec+1)/2 
     &                + (mxvec+1)*(mxvec+2)
       if (.not.bas_numcont(basis,nsh))
     $     call errquit('ccsd_numoper: cannot get basis info',0,
     &       BASIS_ERR)
       ma_req_ccsd = ma_req_ccsd + nsh*(nsh+1)
       ma_req=max(ma_req,ma_req_ccsd)
c
c GA memory for CCSD calculation
c
       ga_req_ccsd = 2*nvir*nvir*nocc*nocc
       ga_req_ccsd = ga_req_ccsd + nocc*(nocc+1)*nvir*(nvir+1)*mxvec
       ga_req=max(ga_req,ga_req_ccsd)

       if (do_trip) then 
c
c MA memory for triples calculation
c
       ma_req_trip = lnov + 2*max(2*lnoov+lnooo,3*lnoov)
       ma_req_trip = ma_req_trip + 14*lnvv + 9*lnov + 2*nvir           
       ma_req_trip = ma_req_trip + moints_trp_mem(basis,1,blen)
       ma_req=max(ma_req,ma_req_trip)
c
c GA memory for triples calculation
c
       ga_req_trip = max(2*lnoov+lnooo,3*lnoov)
       ga_req_trip = ga_req_trip + 2*(lnovv+2*lnoov)*nocc
       ga_req=max(ga_req,ga_req_trip)
c
       endif
c
c Get minimum amount of memory available on nodes
c
       mem_avail = ma_inquire_avail(MT_DBL)
       mem_min_local = mem_avail
       call ga_sync()
       call ga_igop(MSG_MIN_LOCAL_MEM,mem_min_local,1,'min')
       if (ga_uses_ma()) then
         gmem_max = mem_min_local-ma_req
       else
         gmem_max = ga_memory_avail()
       endif
c
c Check if the memory allocation will fit
c
       if (ma_req.gt.mem_min_local) then
          if(ga_nodeid().eq.0) then
             write(6,*) 'The CCSD calculation needs at least 
     $                  ',ma_req*8,'of local MA memory !'
             write(6,*) 'The CCSD calculation needs at least 
     $                  ',ga_req*8,'of global GA memory !'
          endif
          call errquit
     $     ('ccsd_mem: not MA local memory to do transformation',
     $      ma_req*8, MEM_ERR)
       endif
c
       if (ga_req.gt.gmem_max) then
          if(ga_nodeid().eq.0) then
             write(6,*) 'The CCSD calculation needs at least 
     $                  ',ga_req*8,'of global GA memory !'
             write(6,*) 'The CCSD calculation needs at least 
     $                  ',ma_req*8,'of local MA memory !'
          endif
          call errquit
     $     ('ccsd_mem: not GA local memory to do transformation',
     $      ga_req*8, GA_ERR)
       endif
c
       ccsd_memory_fit = .true.
       return
       end
