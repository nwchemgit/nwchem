      logical function task_dddd(rtdb)
      implicit none
#include "mafdecls.fh"
      integer rtdb
c
      integer k_buf, h_buf, k_scr, h_scr
      logical status
      integer maxg, maxs
      integer ip, ia
      double precision xyz(3,4), prims(20,4), coefs(4)
      double precision accy(5)
c      
      call dfill(3,0.0d00,xyz(1,1),1)
      call dfill(3,0.5d00,xyz(1,2),1)
      call dfill(3,1.0d00,xyz(1,3),1)
      call dfill(3,1.5d00,xyz(1,4),1)
      call dfill(12,1.0d00,xyz,1)
      call dfill(4,1.0d00,coefs,1)
      call dfill(20,1.0d00,prims(1,1),1)
      call dfill(20,0.5d00,prims(1,3),1)
      call dfill(20,0.5d00,prims(1,4),1)
      prims(1,2) = 10.d00
      do ip = 2,20
        prims(ip,2) = prims((ip-1),2) + 10.0d00
      enddo
      maxg = 1296
      maxs = 100 000
      status = ma_alloc_get(mt_dbl,maxg,'int buffer' ,h_buf,k_buf)
      status = status.and.
     &    ma_alloc_get(mt_dbl,maxs,'int scratch',h_scr,k_scr)
      if (.not.status) call errquit('dddd: ma falied ',911)
      call dfill(maxg,0.0d00,dbl_mb(k_buf),1)
      call dfill(maxs,0.0d00,dbl_mb(k_scr),1)
c
      call defNxyz(3)
      call int_acc_std
      call int_acc_get(accy(1))
      accy(1) = accy(1)*100.0d00
      do ia = 2,5
        accy(ia) = accy(ia-1)/100.0d00
      enddo
      do ip = 1,20
        write(6,10000)
        do ia = 1,5
          call int_acc_set(accy(ia))
          call hf2(
     &        xyz(1,1),prims(ip,1),coefs(1),1,1,2,
     &        xyz(1,2),prims(ip,2),coefs(2),1,1,2,
     &        xyz(1,3),prims(ip,3),coefs(3),1,1,2,
     &        xyz(1,4),prims(ip,4),coefs(4),1,1,2,
     &        dbl_mb(k_buf),maxg,.false.,.false.,.false.,
     &        .false.,dbl_mb(k_scr),maxs)
c      
          call dddd_cnt_zero(ip,maxg,dbl_mb(k_buf),1.0d-14,
     &        prims(ip,2),accy(ia))
        enddo
      enddo
c
      status = ma_free_heap(h_buf)
      status = status .and. ma_free_heap(h_scr)
      task_dddd = status
10000 format(1x,40('=-'))
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
