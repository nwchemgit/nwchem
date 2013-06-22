        subroutine tce_mrcc_gacp(g_a,g_b,pga,size)
        implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_mrcc_diis.fh"
#include "tce_main.fh"

        integer g_a,g_b
        integer pga
        integer size
        integer l_buff,k_buff
        integer chunksize
        integer i,ntask,j
        integer istart,iend
        integer count,next,nprocs
        integer nxtasksub
        external nxtasksub

        nprocs = GA_PGROUP_NNODES(pga)
        count = 0
        next = nxtasksub(nprocs,1,pga)

        chunksize = pgbsize*pgbsize*pgbsize
        ntask = size/chunksize

        if(ntask.eq.0)ntask=1
        if((ntask*chunksize).lt.size)ntask=ntask+1

c        if(ga_pgroup_nodeid(pga).eq.0)
c     1 write(LuOut,"('CHUNKSIZE:',I10,' SIZE:',I12,' NTASK:',I10)")
c     2 chunksize,size,ntask

        do i=1,ntask

        if(count.eq.next) then


        if (.not.ma_push_get(mt_dbl,chunksize,'buff',l_buff,
     1  k_buff)) call errquit('tce_mrcc_gacp: MA problem',1,
     2  MA_ERR)

        do j=1,chunksize
          dbl_mb(k_buff+j-1)=0.0d0
        enddo

        istart = (i-1)*chunksize+1
        iend = i*chunksize

        if(iend.gt.size)iend=size

        call ga_get(g_a,istart,iend,1,1,dbl_mb(k_buff),1)
        call ga_put(g_b,istart,iend,1,1,dbl_mb(k_buff),1)

        if (.not.ma_pop_stack(l_buff))
     1  call errquit('tce_mrcc_gacp: MA problem',2,MA_ERR)

        next = nxtasksub(nprocs,1,pga)

        endif

        count = count + 1

        enddo

        next = nxtasksub(-nprocs,1,pga)
        call ga_pgroup_sync(pga)

        return
        end

        subroutine tce_mrcc_daxp(n,coefficients,g_a,g_b,size,iref)
        implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_mrcc_diis.fh"
#include "tce_main.fh"

        integer g_a(maxdiism,maxref),g_b
        integer pga
        integer size
        integer l_buff,k_buff,l_buff2,k_buff2
        integer chunksize
        integer i,ntask,j,k
        integer istart,iend
        integer count,next,nprocs
        integer nxtasksub
        external nxtasksub
        integer n,iref
        double precision coefficients(*)

        pga = int_mb(k_innodes+ga_nnodes()+ga_nodeid()) 

        nprocs = GA_PGROUP_NNODES(pga)
        count = 0
        next = nxtasksub(nprocs,1,pga)

        chunksize = pgbsize*pgbsize*pgbsize
        ntask = size/chunksize

        if(ntask.eq.0)ntask=1
        if((ntask*chunksize).lt.size)ntask=ntask+1

        do i=1,ntask

        if(count.eq.next) then

        if (.not.ma_push_get(mt_dbl,chunksize,'buff',l_buff,
     1  k_buff)) call errquit('tce_mrcc_gacp: MA problem',1,
     2  MA_ERR)

        if (.not.ma_push_get(mt_dbl,chunksize,'buff2',l_buff2,
     1  k_buff2)) call errquit('tce_mrcc_gacp: MA problem',3,
     2  MA_ERR)

        do j=1,chunksize
          dbl_mb(k_buff+j-1)=0.0d0
          dbl_mb(k_buff2+j-1)=0.0d0
        enddo

        istart = (i-1)*chunksize+1
        iend = i*chunksize

        if(iend.gt.size)iend=size

        do k=1,n

        call ga_get(g_a(k,iref),istart,iend,1,1,dbl_mb(k_buff),1)

        do j=1,chunksize
          dbl_mb(k_buff2+j-1) = dbl_mb(k_buff2+j-1) +
     1 coefficients(k) * dbl_mb(k_buff+j-1)
        enddo
        enddo

        call ga_put(g_b,istart,iend,1,1,dbl_mb(k_buff2),1)

        if (.not.ma_pop_stack(l_buff2))
     1  call errquit('tce_mrcc_gacp: MA problem',4,MA_ERR)

        if (.not.ma_pop_stack(l_buff))
     1  call errquit('tce_mrcc_gacp: MA problem',2,MA_ERR)

        next = nxtasksub(nprocs,1,pga)

        endif

        count = count + 1

        enddo

        next = nxtasksub(-nprocs,1,pga)
        call ga_pgroup_sync(pga)

        return
        end



c $Id$
