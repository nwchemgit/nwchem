      subroutine eaf_w_mrcc(fh,d_v2,size_xx)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
#include "eaf.fh"
#include "tce_mrcc.fh"
       integer fh,d_v2,size_xx
       integer i,j,k,l
       integer me,total_np
       integer big_chunk 
       integer offset,size_node
       integer total_proc
       integer chunk_t,current_chunk
       integer num_task
       integer l_aux,k_aux
       integer l_s,k_s
       integer rec_size
       integer length_r
       double precision xx
       logical nodezero
       integer nxtask
       integer next
       integer nprocs
       integer count
       external nxtask 
c global structure
       nodezero = (ga_nodeid().eq.0)
       me = ga_nodeid()
       total_proc = ga_nnodes()
       big_chunk = size_xx/total_proc
c local  structure
       chunk_t=(rec_mem*1024)/8
c
       do i = 1,total_proc  ! --- main do loop

        if(me.eq.i-1) then  ! --- main if statement
         offset=(i-1)*big_chunk
         if(i.ne.total_proc) then
          size_node=big_chunk
         else
          size_node=size_xx-(i-1)*big_chunk
         end if 
c local
         num_task=size_node/chunk_t  
         if(num_task*chunk_t.lt.size_node) num_task=num_task+1
         do j=1,num_task ! local do loop
c
          if(j.lt.num_task) then
           current_chunk=chunk_t
          else
           current_chunk=size_node-(j-1)*chunk_t
          end if
c
          if (.not.ma_push_get(mt_dbl,current_chunk,'l_aux',
     1        l_aux,k_aux))
     1    call errquit('eaf write:eaf l_aux 1',0,MA_ERR)
          call dfill(current_chunk, 0.0d0, dbl_mb(k_aux), 1)
          call ga_get(d_v2,offset+(j-1)*chunk_t+1,
     1         offset+(j-1)*chunk_t+current_chunk,1,1,dbl_mb(k_aux),1)
c  write to eaf
         if (eaf_write(fh,dble(bytes*((j-1)*chunk_t+offset)),
     1       dbl_mb(k_aux),current_chunk*bytes)
     1       .ne.0) 
     1       call errquit('eaf write: write failed',1,DISK_ERR)
c
          if (.not.ma_pop_stack(l_aux))
     1     call errquit('eaf write:eaf l_aux 2',0,MA_ERR)

         enddo ! --- end of local do loop
        end if ! --- end of if
       end do  ! --- end of main do loop
c
       return
       end
c
c
c
      subroutine eaf_r_mrcc(fh,d_v2,size_xx)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
#include "eaf.fh"
#include "tce_mrcc.fh"
       integer fh,d_v2,size_xx
       integer i,j,k,l
       integer me,total_np
       integer big_chunk 
       integer offset,size_node
       integer total_proc
       integer chunk_t,current_chunk
       integer num_task
       integer l_aux,k_aux
       integer l_s,k_s
       integer rec_size
       integer length_r
       double precision xx
       logical nodezero
       integer nxtask
       integer next
       integer nprocs
       integer count
       external nxtask 
c global structure
       nodezero = (ga_nodeid().eq.0)
       me = ga_nodeid()
       total_proc = ga_nnodes()
       big_chunk = size_xx/total_proc
c local  structure
       chunk_t=(rec_mem*1024)/8
c
       do i = 1,total_proc  ! --- main do loop

        if(me.eq.i-1) then  ! --- main if statement
         offset=(i-1)*big_chunk
         if(i.ne.total_proc) then
          size_node=big_chunk
         else
          size_node=size_xx-(i-1)*big_chunk
         end if 
c local
         num_task=size_node/chunk_t  
         if(num_task*chunk_t.lt.size_node) num_task=num_task+1
         do j=1,num_task ! local do loop
c
          if(j.lt.num_task) then
           current_chunk=chunk_t
          else
           current_chunk=size_node-(j-1)*chunk_t
          end if
c
          if (.not.ma_push_get(mt_dbl,current_chunk,'l_aux',
     1        l_aux,k_aux))
     1    call errquit('eaf write:eaf l_aux 1',0,MA_ERR)
          call dfill(current_chunk, 0.0d0, dbl_mb(k_aux), 1)
c  read from eaf
         if (eaf_read(fh,dble(bytes*((j-1)*chunk_t+offset)),
     1       dbl_mb(k_aux),current_chunk*bytes)
     1       .ne.0) 
     1       call errquit('eaf read: read failed',1,DISK_ERR)
c
          call ga_put(d_v2,offset+(j-1)*chunk_t+1,
     1         offset+(j-1)*chunk_t+current_chunk,1,1,dbl_mb(k_aux),1)
c
          if (.not.ma_pop_stack(l_aux))
     1     call errquit('eaf write:eaf l_aux 2',0,MA_ERR)

         enddo ! --- end of local do loop
        end if ! --- end of if
       end do  ! --- end of main do loop
c
       return
       end
c $Id$
