      subroutine write_tensor(filename,d_xx,size_xx,unitn)
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
       integer i,j,k,l
       integer d_xx
       integer size_xx
       integer unitn
       integer num_task
       integer chunk_t,current_chunk
       integer l_aux,k_aux
       integer l_s,k_s
       integer rec_size
       integer length_r
       double precision xx
       character*255 filename
       logical nodezero
       integer nxtask
       integer next
       integer nprocs
       integer count
       external nxtask 
       nodezero=(ga_nodeid().eq.0)
c
       rec_size=(rec_mem*1024)/8 
       chunk_t=rec_size
       num_task=size_xx/chunk_t
       if(num_task*chunk_t.lt.size_xx) num_task=num_task+1
c --
       if(nodezero) then
        write(LuOut,*)'WRITE TENSOR'
        write(LuOut,10)filename
        write(LuOut,9) unitn
        write(LuOut,12)size_xx
        write(LuOut,13)rec_mem
        write(LuOut,11)rec_size
        write(LuOut,14)num_task
        call util_flush(6)
       end if
    9  format(1x,' unit nr:  ',2x,i5)
   10  format(1x,' filename: ',A80)
   11  format(1x,' rec_size: ',2x,i12) 
   12  format(1x,' file size: ',2x,i12)
   13  format(1x,' rec_mem (KB): ',2x,i6)
   14  format(1x,' number of tasks: ',2x,i10)
c
c
       if(nodezero) then ! -------------------------------
c open file
        length_r=8*rec_size
        open(unit=unitn,file=filename, status='unknown',
     &  access='direct',form='unformatted',recl=length_r,err=1000,
     *  iostat=k)
c
        do i=1,num_task
c ---
         if (.not.ma_push_get(mt_dbl,chunk_t,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('write tensor:1MA l_aux',0,MA_ERR)
         call dfill(chunk_t, 0.0d0, dbl_mb(k_aux), 1)
c
         if(i.lt.num_task) then  ! ----
          call ga_get(d_xx,(i-1)*chunk_t+1,(i-1)*chunk_t+chunk_t,
     &                1,1,dbl_mb(k_aux),1)
c *** debug ***
c          write(6,*)'write reg chunk: ',i
c          write(6,*)'chunk_t: ',chunk_t
c          do j=1,chunk_t
c           write(6,100) dbl_mb(k_aux+j-1)
c          enddo
c          call util_flush(6)
c 100      format(f12.5)
c *************
         else                    ! ---
          current_chunk=size_xx-(i-1)*chunk_t
          if (.not.ma_push_get(mt_dbl,current_chunk,'l_s',
     1        l_s,k_s))
     1        call errquit('write tensor:1MA l_s',0,MA_ERR)
          call dfill(current_chunk, 0.0d0, dbl_mb(k_s), 1)
          call ga_get(d_xx,(i-1)*chunk_t+1,(i-1)*chunk_t+current_chunk,
     &                1,1,dbl_mb(k_s),1)
          do j=1,current_chunk
           dbl_mb(k_aux+j-1)=dbl_mb(k_s+j-1)
          enddo
c *** debug ***
c          write(6,*)'write last chunk: ',i
c          write(6,*)'current_chunk: ',current_chunk
c          do j=1,current_chunk
c           write(6,101) dbl_mb(k_aux+j-1)
c          enddo
c          call util_flush(6)
c 101      format(f12.5)
c *************
          if (.not.ma_pop_stack(l_s))
     1     call errquit('write tensor: 2MA l_s',0,MA_ERR)
         end if                  ! --
c
c         write(unitn,rec=i,err=150) dbl_mb(k_aux)
          write(unitn,rec=i,err=150) dbl_mb(k_aux:k_aux+chunk_t-1)
c
          if (.not.ma_pop_stack(l_aux))
     1     call errquit('write tensor: 2MA l_aux',0,MA_ERR)
c
       enddo
       go to 200
c
 1000  write(6,*)'failed to open the file: ',filename
       call util_flush(6)
       go to 200
  150  write(6,*)'write problem IOSTAT: ',k
       call util_flush(6)
       go to 200
c close file
 200   continue
       close(unitn)
c
       end if ! ----------------------------------
c
       call ga_sync()
c
       return
       end
c
c
c
c
c
c
c
c
      subroutine read_tensor(filename,d_xx,size_xx,unitn)
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
       integer i,j,k,l
       integer d_xx
       integer size_xx
       integer unitn
       integer num_task
       integer chunk_t,current_chunk
       integer l_aux,k_aux
       integer l_s,k_s
       integer rec_size
       integer length_r
       double precision xx
       character*255 filename
       logical nodezero
       integer nxtask
       integer next
       integer nprocs
       integer count
       external nxtask 
       nodezero=(ga_nodeid().eq.0)
c
c
       rec_size=(rec_mem*1024)/8 
       chunk_t=rec_size
       num_task=size_xx/chunk_t
       if(num_task*chunk_t.lt.size_xx) num_task=num_task+1
c --
      if(nodezero) then
        write(6,*)'READ TENSOR'
        write(6,10)filename
        write(6,9) unitn
        write(6,12)size_xx
        write(6,13)rec_mem
        write(6,11)rec_size
        write(6,14)num_task
        call util_flush(6)
       end if
    9  format(1x,' unit nr:  ',2x,i5)
   10  format(1x,' filename: ',A80)
   11  format(1x,' rec_size: ',2x,i12)
   12  format(1x,' file size: ',2x,i12)
   13  format(1x,' rec_mem (KB): ',2x,i6)
   14  format(1x,' number of tasks: ',2x,i10)
c
c
       if(nodezero) then ! -------------------------------
c open file
        length_r=8*rec_size
        open(unit=unitn,file=filename, status='unknown',
     &  access='direct',form='unformatted',recl=length_r,err=1000,
     &  iostat=k)
c
        do i=1,num_task
c ---
         if (.not.ma_push_get(mt_dbl,chunk_t,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('read tnsor:1MA l_aux',0,MA_ERR)
         call dfill(chunk_t, 0.0d0, dbl_mb(k_aux), 1)
c
c         read(unitn,rec=i,err=150) dbl_mb(k_aux)
         read(unitn,rec=i,err=150) dbl_mb(k_aux:k_aux+chunk_t-1)
c
         if(i.lt.num_task) then ! -------
c *** debug ***
c          write(6,*)'read reg chunk: ',i
c          write(6,*)'chunk_t: ',chunk_t
c          do j=1,chunk_t
c           write(6,100) dbl_mb(k_aux+j-1)
c          enddo
c          call util_flush(6)
c 100      format(f12.5)
c *************
          call ga_put(d_xx,(i-1)*chunk_t+1,(i-1)*chunk_t+chunk_t,
     &                1,1,dbl_mb(k_aux),1)
         else                   ! -------
          current_chunk=size_xx-(i-1)*chunk_t
          if (.not.ma_push_get(mt_dbl,current_chunk,'l_s',
     1        l_s,k_s))
     1        call errquit('read tensor:1MA l_s',0,MA_ERR)
          call dfill(current_chunk, 0.0d0, dbl_mb(k_s), 1)
          do j=1,current_chunk
           dbl_mb(k_s+j-1)=dbl_mb(k_aux+j-1)
          enddo
c *** debug ***
c          write(6,*)'read last chunk: ',i
c          write(6,*)'current_chunk: ',current_chunk
c          do j=1,current_chunk
c           write(6,101) dbl_mb(k_aux+j-1)
c          enddo
c          call util_flush(6)
c 101      format(f12.5)
c *************
          call ga_put(d_xx,(i-1)*chunk_t+1,(i-1)*chunk_t+current_chunk,
     &                1,1,dbl_mb(k_s),1)
          if (.not.ma_pop_stack(l_s))
     1     call errquit('read tensor: 2MA l_s',0,MA_ERR)
         end if                 ! -------
c
c
          if (.not.ma_pop_stack(l_aux))
     1     call errquit('write tensor: 2MA l_aux',0,MA_ERR)
c
       enddo
       go to 200
c
 1000  write(6,*)'failed to open the file: ',filename
       call util_flush(6)
       go to 200
  150  write(6,*)'write problem IOSTAT: ',k
       call util_flush(6)
       go to 200
c close file
 200   continue
c close file
       close(unitn)
c
       end if ! ----------------------------------
c
       call ga_sync()
c
       return
       end
c
c
c
c
c $Id$
