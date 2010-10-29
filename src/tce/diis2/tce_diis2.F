      subroutine tce_diis2(transpose,iter,t1,t2,t3,t4,
     1  d_r1,d_t1,k_t1_offset,size_t1,
     2  d_r2,d_t2,k_t2_offset,size_t2,
     3  d_r3,d_t3,k_t3_offset,size_t3,
     4  d_r4,d_t4,k_t4_offset,size_t4,omega)
c
c $Id$
c
c     Copies residual and t1-t4 files in a round robin fashion
c     and also performs DIIS if necessary
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
      logical nodezero         ! True if node 0
      logical transpose
      logical t1, t2, t3, t4
      integer d_r1
      integer d_t1
      integer k_t1_offset
      integer d_r1_diis
      integer d_t1_diis
      integer size_t1
      integer d_r2
      integer d_t2
      integer k_t2_offset
      integer d_r2_diis
      integer d_t2_diis
      integer size_t2
      integer d_r3
      integer d_t3
      integer k_t3_offset
      integer d_r3_diis
      integer d_t3_diis
      integer size_t3
      integer d_r4
      integer d_t4
      integer k_t4_offset
      integer d_r4_diis
      integer d_t4_diis
      integer size_t4
      double precision a(maxdiis+1,maxdiis+1)
      double precision b(maxdiis+1)
      integer iwork(maxdiis+1)
      integer info
      character*255 filename
c      integer order
      integer i,j
c      integer tce_diis_switch
c      external tce_diis_switch
      double precision ddotfile,ddotfile_1,ddotfile_2
      external ddotfile,ddotfile_1,ddotfile_2
      character*7 r1suffixes(maxdiis)
      data r1suffixes/'rr1.d1','rr1.d2','rr1.d3','rr1.d4',
     1   'rr1.d5','rr1.d6','rr1.d7','rr1.d8','rr1.d9',
     2   'rr1.d10','rr1.d11','rr1.d12','rr1.d13',
     3   'rr1.d14','rr1.d15','rr1.d16','rr1.d17',
     4   'rr1.d18','rr1.d19','rr1.d20'/
      character*7 t1suffixes(maxdiis)
      data t1suffixes/'tr1.d1','tr1.d2','tr1.d3','tr1.d4',
     1   'tr1.d5','tr1.d6','tr1.d7','tr1.d8','tr1.d9',
     2   'tr1.d10','tr1.d11','tr1.d12','tr1.d13',
     3   'tr1.d14','tr1.d15','tr1.d16','tr1.d17',
     4   'tr1.d18','tr1.d19','tr1.d20'/
      character*7 r2suffixes(maxdiis)
      data r2suffixes/'rr2.d1','rr2.d2','rr2.d3','rr2.d4',
     1   'rr2.d5','rr2.d6','rr2.d7','rr2.d8','rr2.d9',
     2   'rr2.d10','rr2.d11','rr2.d12','rr2.d13',
     3   'rr2.d14','rr2.d15','rr2.d16','rr2.d17',
     4   'rr2.d18','rr2.d19','rr2.d20'/
      character*7 t2suffixes(maxdiis)
      data t2suffixes/'tr2.d1','tr2.d2','tr2.d3','tr2.d4',
     1   'tr2.d5','tr2.d6','tr2.d7','tr2.d8','tr2.d9',
     2   'tr2.d10','tr2.d11','tr2.d12','tr2.d13',
     3   'tr2.d14','tr2.d15','tr2.d16','tr2.d17',
     4   'tr2.d18','tr2.d19','tr2.d20'/
      character*7 r3suffixes(maxdiis)
      data r3suffixes/'rr3.d1','rr3.d2','rr3.d3','rr3.d4',
     1   'rr3.d5','rr3.d6','rr3.d7','rr3.d8','rr3.d9',
     2   'rr3.d10','rr3.d11','rr3.d12','rr3.d13',
     3   'rr3.d14','rr3.d15','rr3.d16','rr3.d17',
     4   'rr3.d18','rr3.d19','rr3.d20'/
      character*7 t3suffixes(maxdiis)
      data t3suffixes/'tr3.d1','tr3.d2','tr3.d3','tr3.d4',
     1   'tr3.d5','tr3.d6','tr3.d7','tr3.d8','tr3.d9',
     2   'tr3.d10','tr3.d11','tr3.d12','tr3.d13',
     3   'tr3.d14','tr3.d15','tr3.d16','tr3.d17',
     4   'tr3.d18','tr3.d19','tr3.d20'/
      character*7 r4suffixes(maxdiis)
      data r4suffixes/'rr4.d1','rr4.d2','rr4.d3','rr4.d4',
     1   'rr4.d5','rr4.d6','rr4.d7','rr4.d8','rr4.d9',
     2   'rr4.d10','rr4.d11','rr4.d12','rr4.d13',
     3   'rr4.d14','rr4.d15','rr4.d16','rr4.d17',
     4   'rr4.d18','rr4.d19','rr4.d20'/
      character*7 t4suffixes(maxdiis)
      data t4suffixes/'tr4.d1','tr4.d2','tr4.d3','tr4.d4',
     1   'tr4.d5','tr4.d6','tr4.d7','tr4.d8','tr4.d9',
     2   'tr4.d10','tr4.d11','tr4.d12','tr4.d13',
     3   'tr4.d14','tr4.d15','tr4.d16','tr4.d17',
     4   'tr4.d18','tr4.d19','tr4.d20'/
      double precision omega,shift
c
c     =====================
c     Jacobi if DIIS is off
c     =====================
c
      nodezero = (ga_nodeid() .eq. 0)
      if (nodezero.and.util_print('time',print_debug)) then
        write(LuOut,8000) 'createfile   ', cpusecs(1), cpusecs(51)
        write(LuOut,8000) 'deletefile   ', cpusecs(2), cpusecs(52)
        write(LuOut,8000) 'put_block    ', cpusecs(3), cpusecs(53)
        write(LuOut,8000) 'add_block    ', cpusecs(4), cpusecs(54)
        write(LuOut,8000) 'get_block    ', cpusecs(5), cpusecs(55)
        write(LuOut,8000) 'reconcilefile', cpusecs(6), cpusecs(56)
      endif
 8000 format(1x,'cpu & wall time in seconds for ',A,2f20.1)
      if (diis .le. 0) then
        if (transpose) then
          shift = -1.0d0*zlshift2(2)
          if (t1) call tce_jacobi_lr1(d_r1,d_t1,k_t1_offset,omega,shift)
          shift = -2.0d0*zlshift2(2)
          if (t2) call tce_jacobi_lr2(d_r2,d_t2,k_t2_offset,omega,shift)
          shift = -3.0d0*zlshift2(2)
          if (t3) call tce_jacobi_lr3(d_r3,d_t3,k_t3_offset,omega,shift)
          shift = -4.0d0*zlshift2(2)
          if (t4) call tce_jacobi_lr4(d_r4,d_t4,k_t4_offset,omega,shift)
        else
          shift = -1.0d0*zlshift2(1)
          if (t1) call tce_jacobi_tr1(d_r1,d_t1,k_t1_offset,omega,shift)
          shift = -2.0d0*zlshift2(1)
          if (t2) call tce_jacobi_tr2(d_r2,d_t2,k_t2_offset,omega,shift)
          shift = -3.0d0*zlshift2(1)
          if (t3) call tce_jacobi_tr3(d_r3,d_t3,k_t3_offset,omega,shift)
          shift = -4.0d0*zlshift2(1)
          if (t4) call tce_jacobi_tr4(d_r4,d_t4,k_t4_offset,omega,shift)
        endif
        return
      endif
c
c
c
c ATTENTION !!! ORTHOGONALIZATION !!! (for diis > 1) ---
c     write(6,*)'before tce_ortho iter=',iter
c     call util_flush(6) 
c      call tce_r_ortho(d_r1,d_r2,k_t1_offset,k_t2_offset,
c     &                 size_t1,size_t2,iter)
c      write(6,*)'after tce_ortho iter=',iter
c      call util_flush(6)
c ----------------------------------------
c
c     ================
c     Copy the r files
c     ================
c
      if (t1) then
      r1_pointer = r1_pointer + 1
      if (r1_pointer .gt. diis2) r1_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*) 'diis_exist',(r1_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.r1_diis_exist(r1_pointer)) then
        call tce_filename(r1suffixes(r1_pointer),filename)
        call createfile(filename,d_r1_diis,size_t1)         
        call gatoeaf(d_r1_diis)         
        call gatodra(d_r1_diis)
        r1_diis(r1_pointer) = d_r1_diis
        r1_diis_exist(r1_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(r1_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_r1_diis = r1_diis(r1_pointer)
      endif
c ---
c        call tce_r1_divide(d_r1,k_t1_offset)
c ---
      call copyfile(d_r1,d_r1_diis,size_t1)
      endif

      if (t2) then
      r2_pointer = r2_pointer + 1
      if (r2_pointer .gt. diis2) r2_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r2_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.r2_diis_exist(r2_pointer)) then
        call tce_filename(r2suffixes(r2_pointer),filename)
        call createfile(filename,d_r2_diis,size_t2)         
        call gatoeaf(d_r2_diis)         
        call gatodra(d_r2_diis)
        r2_diis(r2_pointer) = d_r2_diis
        r2_diis_exist(r2_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(r2_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_r2_diis = r2_diis(r2_pointer)
      endif
c ---
c        call tce_r2_divide(d_r2,k_t2_offset)
c ---
      call copyfile(d_r2,d_r2_diis,size_t2)
      endif

      if (t3) then
      r3_pointer = r3_pointer + 1
      if (r3_pointer .gt. diis2) r3_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r3_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.r3_diis_exist(r3_pointer)) then
        call tce_filename(r3suffixes(r3_pointer),filename)
        call createfile(filename,d_r3_diis,size_t3)         
        call gatoeaf(d_r3_diis)         
        call gatodra(d_r3_diis)
        r3_diis(r3_pointer) = d_r3_diis
        r3_diis_exist(r3_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(r3_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_r3_diis = r3_diis(r3_pointer)
      endif
      call copyfile(d_r3,d_r3_diis,size_t3)
      endif

      if (t4) then
      r4_pointer = r4_pointer + 1
      if (r4_pointer .gt. diis2) r4_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r4_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.r4_diis_exist(r4_pointer)) then
        call tce_filename(r4suffixes(r4_pointer),filename)
        call createfile(filename,d_r4_diis,size_t4)         
        call gatoeaf(d_r4_diis)         
        call gatodra(d_r4_diis)
        r4_diis(r4_pointer) = d_r4_diis
        r4_diis_exist(r4_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(r4_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_r4_diis = r4_diis(r4_pointer)
      endif
      call copyfile(d_r4,d_r4_diis,size_t4)
      endif
c
c     ================
c     Copy the t files
c     ================
c
      if (t1) then
      t1_pointer = t1_pointer + 1
      if (t1_pointer .gt. diis2) t1_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t1_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.t1_diis_exist(t1_pointer)) then
        call tce_filename(t1suffixes(t1_pointer),filename)
        call createfile(filename,d_t1_diis,size_t1)         
        call gatoeaf(d_t1_diis)         
        call gatodra(d_t1_diis)         
        t1_diis(t1_pointer) = d_t1_diis
        t1_diis_exist(t1_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(t1_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_t1_diis = t1_diis(t1_pointer)
      endif
      call copyfile(d_t1,d_t1_diis,size_t1)
      endif

      if (t2) then
      t2_pointer = t2_pointer + 1
      if (t2_pointer .gt. diis2) t2_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t2_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.t2_diis_exist(t2_pointer)) then
        call tce_filename(t2suffixes(t2_pointer),filename)
        call createfile(filename,d_t2_diis,size_t2)         
        call gatoeaf(d_t2_diis)         
        call gatodra(d_t2_diis)         
        t2_diis(t2_pointer) = d_t2_diis
        t2_diis_exist(t2_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(t2_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_t2_diis = t2_diis(t2_pointer)
      endif
      call copyfile(d_t2,d_t2_diis,size_t2)
      endif

      if (t3) then
      t3_pointer = t3_pointer + 1
      if (t3_pointer .gt. diis2) t3_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t3_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.t3_diis_exist(t3_pointer)) then
        call tce_filename(t3suffixes(t3_pointer),filename)
        call createfile(filename,d_t3_diis,size_t3)         
        call gatoeaf(d_t3_diis)         
        call gatodra(d_t3_diis)         
        t3_diis(t3_pointer) = d_t3_diis
        t3_diis_exist(t3_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(t3_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_t3_diis = t3_diis(t3_pointer)
      endif
      call copyfile(d_t3,d_t3_diis,size_t3)
      endif

      if (t4) then
      t4_pointer = t4_pointer + 1
      if (t4_pointer .gt. diis2) t4_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t4_diis_exist(i),i=1,diis2)
        call util_flush(LuOut)
      endif
      if (.not.t4_diis_exist(t4_pointer)) then
        call tce_filename(t4suffixes(t4_pointer),filename)
        call createfile(filename,d_t4_diis,size_t4)         
        call gatoeaf(d_t4_diis)         
        call gatodra(d_t4_diis)         
        t4_diis(t4_pointer) = d_t4_diis
        t4_diis_exist(t4_pointer) = .true.
        if (nodezero.and.util_print('diis',print_debug)) then
          write(LuOut,9010) 'diis',filename(1:20)
          write(LuOut,9020) 'diis',(t4_diis(i),i=1,diis2)
          call util_flush(LuOut)
        endif
      else
        d_t4_diis = t4_diis(t4_pointer)
      endif
      call copyfile(d_t4,d_t4_diis,size_t4)
      endif
c
c     ==============
c     DIIS or Jacobi
c     ==============
c
ccc      order = tce_diis_switch(iter)
      
      IF (r1_pointer .lt. diis2) THEN
        if (transpose) then
          shift = -1.0d0*zlshift2(2)
          if (t1) call tce_jacobi_lr1(d_r1,d_t1,k_t1_offset,omega,shift)
          shift = -2.0d0*zlshift2(2)
          if (t2) call tce_jacobi_lr2(d_r2,d_t2,k_t2_offset,omega,shift)
          shift = -3.0d0*zlshift2(2)
          if (t3) call tce_jacobi_lr3(d_r3,d_t3,k_t3_offset,omega,shift)
          shift = 0.0d0
          if (t4) call tce_jacobi_lr4(d_r4,d_t4,k_t4_offset,omega,shift)
        else
          shift = -1.0d0*zlshift2(1)
          if (t1) call tce_jacobi_tr1(d_r1,d_t1,k_t1_offset,omega,shift)
          shift = -2.0d0*zlshift2(1)
          if (t2) call tce_jacobi_tr2(d_r2,d_t2,k_t2_offset,omega,shift)
          shift = -3.0d0*zlshift2(1)
          if (t3) call tce_jacobi_tr3(d_r3,d_t3,k_t3_offset,omega,shift)
          shift = 0.0d0
          if (t4) call tce_jacobi_tr4(d_r4,d_t4,k_t4_offset,omega,shift)
        endif
      END IF
c
c     Make Matrix A (at the end of DIIS cycle)
c
      IF(r1_pointer.eq.diis2) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!
      if(nodezero) then
      write(6,*)'MICROCYCLE DIIS UPDATE:',
     &          iter,r1_pointer
      end if 
c
        do i = 1, diis2
          do j = i, diis2
            a(i,j) = 0.0d0
            if (t1) a(i,j) = a(i,j)
     1      + ddotfile_1(r1_diis(i),r1_diis(j),size_t1)
            if (t2) a(i,j) = a(i,j)
     1      + ddotfile_2(r2_diis(i),r2_diis(j),size_t2)
            if (t3) a(i,j) = a(i,j)
     1      + ddotfile(r3_diis(i),r3_diis(j),size_t3)
            if (t4) a(i,j) = a(i,j)
     1      + ddotfile(r4_diis(i),r4_diis(j),size_t4)
            a(j,i) = a(i,j)
          enddo
        enddo
        do i = 1,diis2
          a(i,diis2+1) = -1.0d0
          a(diis2+1,i) = -1.0d0
        enddo
        a(diis2+1,diis2+1) = 0.0d0
        if (nodezero.and.util_print('diis',print_debug)) then
          call ma_print(a,maxdiis+1,maxdiis+1,'A matrix')
          call util_flush(LuOut)
        endif
c
c     Make Vector B
c
        do i = 1,diis2
          b(i) = 0.0d0
        enddo
        b(diis2+1) = -1.0d0
        if (nodezero.and.util_print('diis',print_debug)) then
          call ma_print(b,maxdiis+1,1,'B vector')
          call util_flush(LuOut)
        endif
c
c     Solve AX = B
c
        call dgesv(diis2+1,1,a,maxdiis+1,iwork,b,maxdiis+1,info)
        if (info .ne. 0)
     1    call errquit('tce_diis: LU decomposition failed',info,
     2    UNKNOWN_ERR)
        if (nodezero.and.util_print('diis',print_debug)) then
          call ma_print(b,maxdiis+1,1,'X matrix')
          call util_flush(LuOut)
        endif
c
c     Make New t's
c
        if (t1) call daxfile(diis2,b,t1_diis,d_t1,size_t1)
        if (t2) call daxfile(diis2,b,t2_diis,d_t2,size_t2)
        if (t3) call daxfile(diis2,b,t3_diis,d_t3,size_t3)
        if (t4) call daxfile(diis2,b,t4_diis,d_t4,size_t4)
        call ga_sync()
      END IF
 9010 format(1x,A,' file name    = ',A)
 9020 format(1x,A,' file handles = ',20i5)
      return
      end
