      logical function task_dddd(rtdb)
      implicit none
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
      if (.not.status) call errquit('dddd: ma falied ',911)
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
      logical function task_accy(rtdb)
      implicit none
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "global.fh"
      logical int_normalize
      external int_normalize
      integer rtdb
      integer geom, basis
      integer szbuf, maxg, maxscr, lmax, nshell
      integer size_grad, size_scr, size_labs
      integer h_grad_1, h_grad_2, h_labs, h_scr
      integer k_grad_1, k_grad_2, k_labs, k_scr
      logical status
c
      if (.not.geom_create(geom,'geometry'))
     &    call errquit('raktest_accy: geom_create failed?',911)
c
      if (.not. geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('raktest_accy: geom_rtdb_load -ref failed',911)
      
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &      ('bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &      ('bas_rtdb_load failed',911)
c
      if (ga_nodeid().eq.0) then
        if (.not.geom_print(geom)) stop ' print error'
        if (.not.bas_print(basis)) stop ' print error'
        if (.not.gbs_map_print(basis)) stop ' gbs_map_print 2?'
      endif
      if (.not.int_normalize(basis)) stop ' norm error 1'
      if (.not.bas_numcont(basis,nshell)) stop 'bas_numcont fe'
      if (.not.bas_high_angular(basis,Lmax)) stop 'bas_ha fe'
c
      szbuf = (Lmax+1)*(Lmax+2)/2
      szbuf = szbuf**4
c
      call intd_init(rtdb,1,basis)
      call intb_mem_2e4c(maxg,maxscr)
      maxscr = maxscr + szbuf
      maxscr = maxscr + maxscr/5 + 1
      maxscr = max(51000,maxscr)
c
      size_grad = 12*szbuf ! 2 grad buffers
      size_scr  = maxscr
      size_labs = 4*szbuf ! for labels
c
      if (.not.ma_push_get(mt_dbl,size_grad,' grad buf 1 ',
     &    h_grad_1, k_grad_1)) stop ' push failed grad 1'
      if (.not.ma_push_get(mt_dbl,size_grad,' grad buf 2 ',
     &    h_grad_2, k_grad_2)) stop ' push failed grad 2'
      if (.not.ma_push_get(mt_dbl,size_scr,' scratch ',
     &    h_scr, k_scr)) stop ' push failed scratch '
      if (.not.ma_push_get(mt_int,size_labs,' labels ',
     &    h_labs, k_labs)) stop ' push failed labels'
      call task_accy_a(rtdb,geom,basis,
     &    size_grad,size_scr,size_labs,nshell,
     &    dbl_mb(k_grad_1),dbl_mb(k_grad_2),
     &    dbl_mb(k_scr),int_mb(k_labs))
c
      
c
      status = .true.
      status = status.and.ma_pop_stack(h_labs)
      status = status.and.ma_pop_stack(h_scr)
      status = status.and.ma_pop_stack(h_grad_1)
      status = status.and.ma_pop_stack(h_grad_2)
      task_accy = status
c
      end
      subroutine task_accy_a(rtdb,geom,basis,
     &    lgrad,lscr,llab,nshell,
     &    diff, grad1, grad2, scr, labs)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
c::functions
      logical intb_init4c
      logical intb_2e4c
      external intb_init4c
      external intb_2e4c
      integer rtdb,geom,basis
      integer lgrad, lscr, llab
      integer nshell
      double precision grad1(12,*)
      double precision grad2(12,*)
      double precision scr(lscr)
      double precision diff(12,*)
      integer labs(4,*)
c
      double precision q4
      double precision dummy
      integer me, nproc
      integer task_count
      integer i, j, k, l
      integer ish, jsh, ksh, lsh
      integer tish, tjsh, tksh, tlsh
      integer ilo, ihi, nbfi
      integer jlo, jhi, nbfj
      integer klo, khi, nbfk 
      integer llo, lhi, nbfl
      integer nint
c
      q4 = 1.0d00
      dummy = 1.0d00
      me    = ga_nodeid()
      nproc = ga_nnodes()
      task_count = (me-1)
      do ish = 1,nshell
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) stop 'cn2bfr i'        
        nbfi = ihi - ilo + 1
        do jsh = i,nshell
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) stop 'cn2bfr j'        
          nbfj = jhi - jlo + 1
          do ksh = 1,nshell
            if (.not.bas_cn2bfr(basis,ksh,klo,khi)) stop 'cn2bfr k'        
            nbfk = khi - klo + 1
            do lsh = 1,nshell
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi)) stop 'cn2bfr l'        
              nbfl = lhi - llo + 1
              nint = nbfi*nbfj*nbfk*nbfl
              task_count = task_count + 1
              if (mod(task_count,nproc).eq.0) then
                tish = ish
                tjsh = jsh
                tksh = ksh
                tlsh = lsh
                call ifill((4*nint),0,labs,1)
                call dfill(nint,0.0d00,grad1,1)
                if (.not.intb_init4c(
     &              basis,tish,tjsh,
     &              basis,tksh,tlsh,1,
     &              q4,.false.,lscr,scr,nint,dummy)) stop 'error'
00100           continue
              endif
            enddo
          enddo
        enddo
      enddo
c
      end
