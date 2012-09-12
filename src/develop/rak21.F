*
* $Id$
*
C> \ingroup task
C> @{
      logical function task_dddd(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
c
      integer numthr
      parameter (numthr = 40)
      integer nvcount,nval(0:numthr), nvallast(0:numthr)
      common /nvalstuff/nvcount,nval,nvallast
c
      integer k_buf, h_buf, k_scr, h_scr
      integer k_buf2, h_buf2, k_diff, h_diff
      logical status
      integer maxg, maxs
      integer ip, ia, ii, ia2
      double precision xyz(3,4), prims(20,4), coefs(4)
      double precision accy(10), norm
      double precision ddot
      external ddot
c      
      call dfill(12,1.0d00,xyz,1)
      xyz(3,2) = -1.72d00*1.889725989d00
      xyz(3,3) =  1.72d00*1.889725989d00
      xyz(3,4) =  1.72d00*1.889725989d00
      call dfill(4,1.0d00,coefs,1)
      call dfill(20,1.0d00,prims(1,1),1)
      call dfill(20,0.5d00,prims(1,3),1)
      call dfill(20,0.5d00,prims(1,4),1)
      prims(1,2) = 100.d00
      do ip = 2,20
        prims(ip,2) = prims((ip-1),2) + 10.0d00
      enddo
*use reported problem
      prims(1,1) = 0.374506d00
      prims(1,2) = 11720.0d00
      prims(1,3) = 0.332
      prims(1,4) = 0.332
      prims(2,2) = 5000.0d00
      maxg = 1296
      maxs = 100 000
      status = ma_alloc_get(mt_dbl,maxg,'int buffer' ,h_buf,k_buf)
      status = status.and.
     &    ma_alloc_get(mt_dbl,maxg,'int buffer2' ,h_buf2,k_buf2)
      status = status.and.
     &    ma_alloc_get(mt_dbl,maxg,'int diff buffer' ,h_diff,k_diff)
      status = status.and.
     &    ma_alloc_get(mt_dbl,maxs,'int scratch',h_scr,k_scr)
      if (.not.status) call errquit('dddd: ma falied ',911, MA_ERR)
      call dfill(maxg,0.0d00,dbl_mb(k_buf),1)
      call dfill(maxg,0.0d00,dbl_mb(k_buf2),1)
      call dfill(maxs,0.0d00,dbl_mb(k_scr),1)
c
      call ifill((numthr+1),-565,nvallast,1)
c
      call defNxyz(3)
      call int_acc_std
      call int_acc_get(accy(1))
*      accy(1) = accy(1)*100.0d00
      do ia = 2,10
        accy(ia) = accy(ia-1)/10.0d00
      enddo
      do ip = 1,2 ! 0
        write(6,10000)
        do ia = 1,10
          nvcount=ia
          call int_acc_set(accy(ia))
          call dfill(maxg,0.0d00,dbl_mb(k_buf),1)
          call dfill(maxs,0.0d00,dbl_mb(k_scr),1)
*rak:          call hf2(
*rak:     &        xyz(1,1),prims(ip,1),coefs(1),1,1,2,
*rak:     &        xyz(1,2),prims(ip,2),coefs(2),1,1,2,
*rak:     &        xyz(1,3),prims(ip,3),coefs(3),1,1,2,
*rak:     &        xyz(1,4),prims(ip,4),coefs(4),1,1,2,
*rak:     &        dbl_mb(k_buf),maxg,.false.,.false.,.false.,
*rak:     &        .false.,dbl_mb(k_scr),maxs)
c      
          call dddd_cnt_zero(ip,maxg,dbl_mb(k_buf),1.0d-14,
     &        prims(ip,2),accy(ia))
          call dddd_histogram(ip,maxg,dbl_mb(k_buf),1.0d-14,
     &        prims(ip,2),accy(ia),(nvcount.eq.1))
          do ii = 0,numthr
            nvallast(ii)=nval(ii)
          enddo
          do ia2 = 1,10
            
            call int_acc_set(accy(ia2))
            call dfill(maxg,0.0d00,dbl_mb(k_buf2),1)
            call dfill(maxs,0.0d00,dbl_mb(k_scr),1)
*rak:            call hf2(
*rak:     &          xyz(1,1),prims(ip,1),coefs(1),1,1,2,
*rak:     &          xyz(1,2),prims(ip,2),coefs(2),1,1,2,
*rak:     &          xyz(1,3),prims(ip,3),coefs(3),1,1,2,
*rak:     &          xyz(1,4),prims(ip,4),coefs(4),1,1,2,
*rak:     &          dbl_mb(k_buf2),maxg,.false.,.false.,.false.,
*rak:     &          .false.,dbl_mb(k_scr),maxs)
            call dcopy(maxg,dbl_mb(k_buf),1,dbl_mb(k_diff),1)
            call daxpy(maxg,-1.0d00,dbl_mb(k_buf2),1,dbl_mb(k_diff),1)
            norm = ddot(maxg,dbl_mb(k_diff),1,dbl_mb(k_diff),1)
            write(6,10001)ia2,accy(ia),accy(ia2),prims(ip,2),norm
          enddo
        enddo
      enddo
c
      status = ma_free_heap(h_buf)
      status = status .and. ma_free_heap(h_buf2)
      status = status .and. ma_free_heap(h_diff)
      status = status .and. ma_free_heap(h_scr)
      task_dddd = status
10000 format(1x,40('=-'))
10001 format(1x,i3,1x,'<acc1=',1pd12.5,'> <acc2=',1pd12.5,'> ',
     &    '<bexp=',1pd12.5,'> norm:',1pd20.10)
      end
C> @}
      subroutine dddd_histogram(cnt,n,buf,thresh,bexp,accy,firstpass)
      implicit none
      integer n,cnt
      double precision buf(n),thresh,bexp,accy
      integer numthr
      parameter (numthr = 40)
      integer nvcount,nval(0:numthr), nvallast(0:numthr)
      common /nvalstuff/nvcount,nval,nvallast
      integer i,it
      double precision thr(numthr), value, valmax, valmin, range
      logical firstpass
      logical found_slot

c
      valmax = 1.0d-20
      valmin = 1.0d20
      do i = 1,n
        value = abs(buf(i))
        valmax = max(valmax,value)
        valmin = min(valmin,value)
      enddo
      valmin = max(valmin,1.0d-07)
      valmax = min(valmax,0.1d00)
      if (valmin.gt.valmax) then
        range = valmin
        valmin = valmax
        valmax = range
      endif
*      write(6,*)' valmax :',valmax
*      write(6,*)' valmin :',valmin
c
      call dfill(numthr,0.0d00,thr,1)
      range = valmax - valmin
*      write(6,*)' range 1:',range
      range = range/20.0d00
*      write(6,*)' range 2:',range
*      thr(1) = valmax
      thr(1) = 0.01d00
      do i = 2,numthr,2
        thr(i)  = thr(i-1)/2.0d00
        thr(i+1)= thr(i)/5.0d00
      enddo
c
      call ifill((numthr+1),0,nval,1)
c
      do i = 1,n
        value = abs(buf(i))
        found_slot = .false.
        do it = 1,numthr
          if (.not.found_slot.and.value.gt.thr(it)) then
            nval((it-1)) = nval((it-1)) + 1
            found_slot = .true.
          endif
        enddo
        if (.not.found_slot) nval(numthr) = nval(numthr) + 1
      enddo
      call dddd_print_nvalue(nval,thr,numthr,n,nvallast,firstpass)
      end
      subroutine dddd_print_nvalue(nval,thr,numthr,n,nvallast,firstpass)
      implicit none
      integer n
      integer numthr
      integer nval(0:numthr)
      integer nvallast(0:numthr)
      double precision thr(numthr)
      integer it, sum, diff
      logical firstpass
c
      diff = nvallast(0) - nval(0)
      if ((.not.firstpass).and.(diff.ne.0)) then
        write(6,10000)thr(1),nval(0),nvallast(0),diff
      elseif (firstpass) then
        diff = 0
        nvallast(0) = 0
        write(6,10000)thr(1),nval(0),nvallast(0),diff
      endif
      do it = 2,numthr
        diff = nvallast(it-1) - nval(it-1)
        if ((.not.firstpass).and.(diff.ne.0)) then
          write(6,10001)
     &        thr(it-1),thr(it),nval(it-1),nvallast(it-1),diff
        elseif (firstpass) then
          diff = 0
          nvallast(it-1) = 0
          write(6,10001)
     &        thr(it-1),thr(it),nval(it-1),nvallast(it-1),diff
        endif
      enddo
      diff = nvallast(numthr) - nval(numthr)
      if ((.not.firstpass).and.(diff.ne.0)) then
        write(6,10002)thr(numthr),nval(numthr),nvallast(numthr),diff
      elseif (firstpass) then
        diff = 0
        nvallast(numthr) = 0
        write(6,10002)thr(numthr),nval(numthr),nvallast(numthr),diff
      endif
10000 format(1x,20x,     '   x > ',1pd20.10, ' :',3(i10))
10001 format(1x,1pd20.10,' > x > ',1pd20.10, ' :',3(i10))
10002 format(1x,20x,     '   x < ',1pd20.10, ' :',3(i10))
c
      sum = 0
      do it = 0,numthr
        sum = sum + nval(it)
      enddo
      if ((sum-n).ne.0) write(6,*)' sum = ',sum, ' n = ',n
c      
      end
      subroutine dddd_cnt_zero(cnt,n,buf,thresh,bexp,accy)
      implicit none
      integer n,cnt
      double precision buf(n),thresh,bexp,accy
      integer i,nz
      nz = 0
      do i=1,n
        if (abs(buf(i)).lt.thresh)
     &      nz = nz + 1
      enddo
      write(6,*)' '
      write(6,10000)cnt,bexp,accy
      write(6,*)'dddd_cnt_zero: total    :',cnt,':',n
      write(6,*)'dddd_cnt_zero:  zero    :',cnt,':',nz
      write(6,*)'dddd_cnt_zero: non-zero :',cnt,':',(n-nz)
10000 format(1x,25('-+'),i3,2(1pd20.6))
      end
