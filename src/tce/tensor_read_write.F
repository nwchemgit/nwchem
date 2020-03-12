#define USE_EAF 1
      subroutine write_tensor(filename,d_xx,size_xx,unitn)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "inp.fh"
#include "eaf.fh"
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
       character*255 errmsg
       integer lgth
       double precision offset
       integer max_ioprocs
       integer size_xx_perproc,ga_offset
       max_ioprocs=min(max_ioprocs_in,ga_nnodes())
       nodezero=(ga_nodeid().lt.max_ioprocs)
c
       rec_size=(rec_mem*1024)/8 
       size_xx_perproc=size_xx/max_ioprocs
       ga_offset=ga_nodeid()*size_xx_perproc
c     last io proc does the remainder
       if(max_ioprocs*size_xx_perproc.lt.size_xx) then
          if(ga_nodeid().eq.max_ioprocs-1) then
             size_xx_perproc=size_xx-(max_ioprocs-1)*size_xx_perproc
          endif
       endif
       chunk_t=rec_size
       num_task=size_xx_perproc/chunk_t
       if(num_task*chunk_t.lt.size_xx_perproc) num_task=num_task+1
c --
       if(nodezero) then
         write(6,'(I3,A,I19,A,I19,A,I4)') 
     G ga_nodeid(),' ga offset ',ga_offset,
     C      ' size_xx_perproc ',size_xx_perproc,'mx ',max_ioprocs
         if(ga_nodeid().eq.0) then
            write(LuOut,*)'WRITE TENSOR'
            write(LuOut,10)filename
            write(LuOut,9) unitn
            write(LuOut,12)size_xx_perproc
            write(LuOut,13)rec_mem
            write(LuOut,11)rec_size
            write(LuOut,14)num_task
         endif
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
#ifdef USE_EAF
         offset=0d0
c delete existing files to avoid appending
        k=eaf_delete(filename(1:inp_strlen(filename)))
c
        k=eaf_open(filename(1:inp_strlen(filename)),
     ,        eaf_rw, unitn)
        if (k .ne. 0) then
           call eaf_errmsg(k,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('write_tensor: open failed for proc',
     G          ga_nodeid(), DISK_ERR)
        endif
c     write header
          k = eaf_write(unitn, offset, max_ioprocs,
     M       ma_sizeof(MT_INT,1,MT_BYTE))
          offset=offset+ma_sizeof(MT_INT,1,MT_BYTE)
          k = eaf_write(unitn, offset, size_xx_perproc,
     M       ma_sizeof(MT_INT,1,MT_BYTE))
          offset=offset+ma_sizeof(MT_INT,1,MT_BYTE)

#else
        open(unit=unitn,file=filename, status='unknown',
     &  access='direct',form='unformatted',recl=length_r,err=1000,
     *  iostat=k)
#endif
c
         if (.not.ma_push_get(mt_dbl,chunk_t,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('write tensor:1MA l_aux',0,MA_ERR)
        do i=1,num_task
c ---
c
         if(i.eq.num_task) then  ! ----
            current_chunk=size_xx_perproc-(i-1)*chunk_t
         else
            current_chunk=chunk_t
         endif
          call ga_get(d_xx,ga_offset+(i-1)*chunk_t+1,
     g        ga_offset+(i-1)*chunk_t+current_chunk,
     &                1,1,dbl_mb(k_aux),current_chunk)

c
c         write(unitn,rec=i,err=150) dbl_mb(k_aux)
#ifdef USE_EAF
          lgth=ma_sizeof(MT_DBL,1,MT_BYTE)*chunk_t
         if(i.eq.num_task) lgth=(size_xx_perproc-(i-1)*chunk_t)*
     *        ma_sizeof(MT_DBL,1,MT_BYTE)
          k = eaf_write(unitn, offset, dbl_mb(k_aux),lgth)
        if (k .ne. 0) then
           call eaf_errmsg(k,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('write_tensor: write failed for proc',
     G          ga_nodeid(), DISK_ERR)
        endif
          offset=offset+lgth
#else
          write(unitn,rec=i,err=150) dbl_mb(k_aux:k_aux+chunk_t-1)
#endif
c
       enddo
c
          if (.not.ma_pop_stack(l_aux))
     1     call errquit('write tensor: 2MA l_aux',0,MA_ERR)
       go to 200
c
 1000  write(6,*)'failed to open the file: ',filename
       call util_flush(6)
       go to 200
  150  write(6,*)'write problem IOSTAT: ',k
       call util_flush(6)
       go to 200
c close file
  200  continue
#ifdef USE_EAF
       k=eaf_close(unitn)
        if (k .ne. 0) then
           call eaf_errmsg(k,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
        endif
c        call util_fadvise(filename(1:inp_strlen(filename)))
#else
       close(unitn)
#endif
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
#include "inp.fh"
#include "eaf.fh"
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
       integer rec_size,rec_read
       integer length_r
       double precision xx
       character*255 filename
       logical nodezero
       integer nxtask
       integer next
       integer nprocs
       integer count
       integer nbhandle
       character*255 errmsg
       integer lgth
       double precision offset
       external nxtask 
       integer ga_offset
       integer max_ioprocs,size_xx_perproc
       integer max_ioprocs_r,size_xx_perproc_r
       max_ioprocs=min(max_ioprocs_in,ga_nnodes())
       nodezero=(ga_nodeid().lt.max_ioprocs)
c
c
       rec_size=(rec_mem*1024)/8 
       chunk_t=rec_size
       size_xx_perproc=size_xx/max_ioprocs
       ga_offset=ga_nodeid()*size_xx_perproc
c     last io proc does the remainder
       if(max_ioprocs*size_xx_perproc.lt.size_xx) then
          if(ga_nodeid().eq.max_ioprocs-1) then
             size_xx_perproc=size_xx-(max_ioprocs-1)*size_xx_perproc
          endif
       endif
       chunk_t=rec_size
       num_task=size_xx_perproc/chunk_t
       if(num_task*chunk_t.lt.size_xx_perproc) num_task=num_task+1
c --
      if(nodezero) then
         write(6,'(I3,A,I19,A,I19)') 
     G        ga_nodeid(),' ga offset ',ga_offset,
     C      ' size_xx_perproc ',size_xx_perproc
         if(ga_nodeid().eq.0) then
            write(6,*)'READ TENSOR'
            write(6,10)filename
            write(6,9) unitn
            write(6,12)size_xx_perproc
            write(6,13)rec_mem
            write(6,11)rec_size
            write(6,14)num_task
         endif
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
#ifdef USE_EAF
         offset=0d0
        k=eaf_open(filename(1:inp_strlen(filename)),
     ,        eaf_r, unitn)
        if (k .ne. 0) then
           call eaf_errmsg(k,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('read_tensor: open failed for proc',ga_nodeid(),
     D          DISK_ERR)
        endif
        call util_fadvise_noreuse(filename(1:inp_strlen(filename)))
c     read header
          k = eaf_read(unitn, offset, max_ioprocs_r,
     M       ma_sizeof(MT_INT,1,MT_BYTE))
          offset=offset+ma_sizeof(MT_INT,1,MT_BYTE)
          if(max_ioprocs_r.ne.max_ioprocs) call errquit(
     R         'read_tensor: wrong no. ioprocs ',
     D         max_ioprocs_r,DISK_ERR)
          k = eaf_read(unitn, offset, size_xx_perproc_r,
     M       ma_sizeof(MT_INT,1,MT_BYTE))
          if(size_xx_perproc_r.ne.size_xx_perproc) call errquit(
     R         'read_tensor: size_xx_perproc ',
     D         ga_nodeid(),DISK_ERR)
          offset=offset+ma_sizeof(MT_INT,1,MT_BYTE)
#else
        open(unit=unitn,file=filename, status='unknown',
     &  access='direct',form='unformatted',recl=length_r,err=1000,
     &  iostat=k)
#endif
c
         if (.not.ma_push_get(mt_dbl,chunk_t,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('read tnsor:1MA l_aux',0,MA_ERR)
        do i=1,num_task
         call dfill(chunk_t, 0.0d0, dbl_mb(k_aux), 1)
#ifdef USE_EAF
         lgth=ma_sizeof(MT_DBL,1,MT_BYTE)*chunk_t
         if(i.eq.num_task) lgth=(size_xx_perproc-(i-1)*chunk_t)*
     *        ma_sizeof(MT_DBL,1,MT_BYTE)
          k = eaf_read(unitn, offset, dbl_mb(k_aux),lgth)
        if (k .ne. 0) then
           call eaf_errmsg(k,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('read_tensor: read failed for proc',ga_nodeid(),
     D          DISK_ERR)
        endif
          offset=offset+lgth
#else
         read(unitn,rec=i,err=150) dbl_mb(k_aux:k_aux+chunk_t-1)
#endif
         rec_read=chunk_t
         if(i.eq.num_task) rec_read=size_xx_perproc-(i-1)*chunk_t
c
          call ga_put(d_xx,
     g        ga_offset+(i-1)*chunk_t+1,
     g        ga_offset+(i-1)*chunk_t+rec_read,
     &                1,1,dbl_mb(k_aux),rec_read)
c
c
       enddo
          if (.not.ma_chop_stack(l_aux))
     1     call errquit('write tensor: 2MA l_aux',0,MA_ERR)
       go to 200
c
 1000  write(6,*)'failed to open the file: ',filename
       call util_flush(6)
       go to 200
  150  write(6,*)'read problem IOSTAT: ',k
       call util_flush(6)
       go to 200
c close file
  200  continue
c close file
#ifdef USE_EAF
       k=eaf_close(unitn)
        if (k .ne. 0) then
           call eaf_errmsg(k,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
        endif
        call util_fadvise_dontneed(filename(1:inp_strlen(filename)))
#else
       close(unitn)

#endif
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
