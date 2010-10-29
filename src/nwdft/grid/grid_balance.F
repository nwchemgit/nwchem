      logical function grid_file_balance()
*
* $Id$
*  
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "eaf.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
c     
c     Attempt to load balance the grid pts file
c
c     Eliminate intermediate syncrhonizations, attempt to overlap
c     big transfers, fuzzy termination.
c     
cedo#if !(defined(SP) || defined(SP1) || defined(SOLARIS))
cedo      return
cedo#else 
      integer nrec, me, nproc, him, msglen, from, record
      integer tmp, i, j, p, q, test, nexch, nmoved, npass
      integer his_count, my_count, total, shift
      double precision newfilelen, used
c
      integer maxproc
      parameter (maxproc = 2048)
      integer info(2,0:maxproc-1) ! 1=node#, 2=nrec
      logical grid_packed_buf_read, grid_packed_buf_write
     $     , grid_file_rewind
      external grid_packed_buf_read, grid_packed_buf_write
     $     , grid_file_rewind
      logical oprint, odebug
      intrinsic mod
c     
      npass = 0                 ! Used for paranoid infinite loop avoidance
      me = ga_nodeid()          ! Not true on workstation clusters
      nproc = ga_nnodes()
      grid_file_balance=.true.
      if (nproc.eq.1) return
c
      call ga_sync()            ! So can use message passing
      used = util_wallsec()
      odebug = util_print('balance debug',print_never) 
      oprint = util_print('balance', print_default) .and.  me.eq.0
      if (odebug) 
     $     write(LuOut,*) me, ' initial n_re_in_file ', n_rec_in_file
c
      nexch = 0                 ! Just for printing statistics
      nmoved= 0
      grid_file_balance=.false.
c     
c     Fill in info
c     
      call ifill(2*nproc, 0, info, 1)
      info(1,me) = me
      info(2,me) = n_rec_in_file
      call ga_igop(33, info, 2*nproc, '+')
      total = 0
      do i = 0, nproc-1
         total = total + info(2,i)
      enddo
      nrec = (total - 1) / nproc + 1
      test = max(2,nint(0.01*nrec))
c
c     Adjust info so the record count is difference from the mean
c
      do i = 0, nproc-1
         info(2,i) = info(2,i) - nrec
      enddo
c
c     The logic with shifts is to ensure that have mod(total,nproc)
c     processes with nrec, with rest with nrec-1 (rather than
c     nproc-1 with nrec and one with nrec-nproc-1)
c
      do i = 0, nproc-1
         shift = 0
         if (i.ge.mod(total,nproc)) shift = 1
         info(2,i) = info(2,i) + shift
      enddo
c
 10   continue                  ! While loop
      npass = npass + 1
c
c     Sort nodes by #records suplus/deficit
c
      do i = 0, nproc-1
         do j = 0, i-1
            if (info(2,i).gt.info(2,j)) then
               tmp = info(1,i)
               info(1,i) = info(1,j)
               info(1,j) = tmp
               tmp = info(2,i)
               info(2,i) = info(2,j)
               info(2,j) = tmp
            endif
         enddo
      enddo
c
c     Find my place in the list and thus my counterpart (note i and j
c     used from here to 100)
c
      do i = 0, nproc-1
         if (info(1,i) .eq. me) goto 11
      enddo
      call errquit('grid_file_balance: algorithm failure?',0, CALC_ERR)
 11   j = nproc-1-i
c
      if (odebug .and. me.eq.0) then
         write(LuOut,1) (info(1,tmp), info(2,tmp), tmp=0,nproc-1)
 1       format(400(6(4x,i4,i4:)/))
         call util_flush(LuOut)
      endif
c
c     If discrepancy is small enough then quit
c
      if (info(2,0).le.test .and. info(2,nproc-1).ge.-test) goto 100
      if (npass .gt. ((nproc-1)/2+1)) goto 100
      grid_file_balance=.true.
c
c     Pair processor with largest surplus with one with largest deficit
c
      him = info(1,j)
      his_count = info(2,j)
      my_count = info(2,i)
c
c     The case of nproc odd is automatically handled by the logic.
c
      if (my_count.lt.0 .and. his_count.gt.0) then
c     
c     He will send to me
c     
         my_count = min(abs(my_count), his_count)
         if (odebug) then
            write(LuOut,*) me, ' receiving ', my_count, ' from ', him
            call util_flush(LuOut)
         endif
         do record = 1, my_count
            call rcv(2,buf,8*rec_len_in_dbl,msglen,him,from,1)
            if (.not. grid_packed_buf_write(.false.))
     $           call errquit('grid_file_balance:write failed',0,
     &       DISK_ERR)
         enddo
      else if (my_count.gt.0 .and. his_count.lt.0) then
c     
c     I'll send to him
c     
         his_count = min(my_count, abs(his_count))
         cur_rec = n_rec_in_file - his_count + 1
         do record = 1, his_count
            if (.not. grid_packed_buf_read())
     $         call errquit('grid_file_balance: read failed',1,
     &       DISK_ERR)
            call snd(2, buf, 8*rec_len_in_dbl, him, 1)
         enddo
         n_rec_in_file = n_rec_in_file - his_count
      endif
c     
c     Adjust everyone's counters by repeating above logic on all nodes
c
      do p = 0, (nproc-2)/2
         q = nproc-1-p
         if (info(2,p).le.0 .and. info(2,q).ge.0) then
            tmp = min(abs(info(2,p)),info(2,q))
            info(2,p) = info(2,p) + tmp
            info(2,q) = info(2,q) - tmp
            nexch  = nexch + 1
            nmoved = nmoved + tmp
         else if (info(2,p).ge.0 .and. info(2,q).le.0) then
            tmp = min(abs(info(2,q)),info(2,p))
            info(2,p) = info(2,p) - tmp
            info(2,q) = info(2,q) + tmp
            nexch  = nexch + 1
            nmoved = nmoved + tmp
         endif
      enddo
c            
      goto 10
c     
 100  if (n_rec_in_file .gt. max_cache_rec) then
         newfilelen = max(0.0d0,
     $        8.0d0*rec_len_in_dbl*(n_rec_in_file-max_cache_rec))
         if (odebug) then
            write(LuOut,*) me,' truncating at ', n_rec_in_file, 
     $                     newfilelen
            call util_flush(LuOut)
         endif
         if (eaf_truncate(fd, newfilelen) .ne. 0) call errquit
     $        ('grid_file_balance: file truncate failed', 3, DISK_ERR)
         cur_rec = n_rec_in_file
      endif
c
      if (oprint) then
         used = util_wallsec() - used
         write(LuOut,21) nexch, nmoved, used
 21      format(/'File balance: exchanges=',i6,'  moved=',i6,
     $        '  time=',f6.1/)
         call util_flush(LuOut)
      endif
c
      if (.not. grid_file_rewind()) call errquit
     $     ('grid_file_balance: failed rewind integral file', 0,
     &       DISK_ERR)
      if (odebug) 
     $     write(LuOut,*) me, ' final n_re_in_file ', n_rec_in_file
      call ga_sync()
c
cedo#endif
      end
