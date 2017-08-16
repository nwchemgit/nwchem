      subroutine tce_diis(transpose,iter,t1,t2,t3,t4,
     1  d_r1,d_t1,k_t1_offset,size_t1,
     2  d_r2,d_t2,k_t2_offset,size_t2,
     3  d_r3,d_t3,k_t3_offset,size_t3,
     4  d_r4,d_t4,k_t4_offset,size_t4)
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
      logical nodezero
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
      integer i,j
      double precision ddotfile,ddotfile_1,ddotfile_2
      external ddotfile,ddotfile_1,ddotfile_2
      character*7 r1suffixes(maxdiis)
      data r1suffixes/'r1.d1','r1.d2','r1.d3','r1.d4',
     1   'r1.d5','r1.d6','r1.d7','r1.d8','r1.d9',
     2   'r1.d10','r1.d11','r1.d12','r1.d13',
     3   'r1.d14','r1.d15','r1.d16','r1.d17',
     4   'r1.d18','r1.d19','r1.d20'/
      character*7 t1suffixes(maxdiis)
      data t1suffixes/'t1.d1','t1.d2','t1.d3','t1.d4',
     1   't1.d5','t1.d6','t1.d7','t1.d8','t1.d9',
     2   't1.d10','t1.d11','t1.d12','t1.d13',
     3   't1.d14','t1.d15','t1.d16','t1.d17',
     4   't1.d18','t1.d19','t1.d20'/
      character*7 r2suffixes(maxdiis)
      data r2suffixes/'r2.d1','r2.d2','r2.d3','r2.d4',
     1   'r2.d5','r2.d6','r2.d7','r2.d8','r2.d9',
     2   'r2.d10','r2.d11','r2.d12','r2.d13',
     3   'r2.d14','r2.d15','r2.d16','r2.d17',
     4   'r2.d18','r2.d19','r2.d20'/
      character*7 t2suffixes(maxdiis)
      data t2suffixes/'t2.d1','t2.d2','t2.d3','t2.d4',
     1   't2.d5','t2.d6','t2.d7','t2.d8','t2.d9',
     2   't2.d10','t2.d11','t2.d12','t2.d13',
     3   't2.d14','t2.d15','t2.d16','t2.d17',
     4   't2.d18','t2.d19','t2.d20'/
      character*7 r3suffixes(maxdiis)
      data r3suffixes/'r3.d1','r3.d2','r3.d3','r3.d4',
     1   'r3.d5','r3.d6','r3.d7','r3.d8','r3.d9',
     2   'r3.d10','r3.d11','r3.d12','r3.d13',
     3   'r3.d14','r3.d15','r3.d16','r3.d17',
     4   'r3.d18','r3.d19','r3.d20'/
      character*7 t3suffixes(maxdiis)
      data t3suffixes/'t3.d1','t3.d2','t3.d3','t3.d4',
     1   't3.d5','t3.d6','t3.d7','t3.d8','t3.d9',
     2   't3.d10','t3.d11','t3.d12','t3.d13',
     3   't3.d14','t3.d15','t3.d16','t3.d17',
     4   't3.d18','t3.d19','t3.d20'/
      character*7 r4suffixes(maxdiis)
      data r4suffixes/'r4.d1','r4.d2','r4.d3','r4.d4',
     1   'r4.d5','r4.d6','r4.d7','r4.d8','r4.d9',
     2   'r4.d10','r4.d11','r4.d12','r4.d13',
     3   'r4.d14','r4.d15','r4.d16','r4.d17',
     4   'r4.d18','r4.d19','r4.d20'/
      character*7 t4suffixes(maxdiis)
      data t4suffixes/'t4.d1','t4.d2','t4.d3','t4.d4',
     1   't4.d5','t4.d6','t4.d7','t4.d8','t4.d9',
     2   't4.d10','t4.d11','t4.d12','t4.d13',
     3   't4.d14','t4.d15','t4.d16','t4.d17',
     4   't4.d18','t4.d19','t4.d20'/
      double precision shift
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
          shift = -1.0d0*zlshiftl
          if (t1) call tce_jacobi_l1(d_r1,d_t1,k_t1_offset,shift)
          shift = -2.0d0*zlshiftl
          if (t2) call tce_jacobi_l2(d_r2,d_t2,k_t2_offset,shift)
          shift = -3.0d0*zlshiftl
          if (t3) call tce_jacobi_l3(d_r3,d_t3,k_t3_offset,shift)
          shift = 0.0d0
          if (t4) call tce_jacobi_l4(d_r4,d_t4,k_t4_offset,shift)
        else
          shift = -1.0d0*zlshift
          if (t1) call tce_jacobi_t1(d_r1,d_t1,k_t1_offset,shift)
          shift = -2.0d0*zlshift
          if (t2) call tce_jacobi_t2(d_r2,d_t2,k_t2_offset,shift)
          shift = -3.0d0*zlshift
          if (t3) then
           if(activecalc) then
             call tce_jacobi_t3a(d_r3,d_t3,k_t3_offset,shift)
           else
             call tce_jacobi_t3(d_r3,d_t3,k_t3_offset,shift)
           end if
          end if
          shift = 0.0d0
          if (t4) call tce_jacobi_t4(d_r4,d_t4,k_t4_offset,shift)
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
      if (r1_pointer .gt. diis) r1_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*) 'diis_exist',(r1_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(r1_diis(i),i=1,diis)
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
      if (r2_pointer .gt. diis) r2_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r2_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(r2_diis(i),i=1,diis)
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
      if (r3_pointer .gt. diis) r3_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r3_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(r3_diis(i),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_r3_diis = r3_diis(r3_pointer)
      endif
      call copyfile(d_r3,d_r3_diis,size_t3)
      endif

      if (t4) then
      r4_pointer = r4_pointer + 1
      if (r4_pointer .gt. diis) r4_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(r4_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(r4_diis(i),i=1,diis)
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
      if (t1_pointer .gt. diis) t1_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t1_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(t1_diis(i),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_t1_diis = t1_diis(t1_pointer)
      endif
      call copyfile(d_t1,d_t1_diis,size_t1)
      endif

      if (t2) then
      t2_pointer = t2_pointer + 1
      if (t2_pointer .gt. diis) t2_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t2_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(t2_diis(i),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_t2_diis = t2_diis(t2_pointer)
      endif
      call copyfile(d_t2,d_t2_diis,size_t2)
      endif

      if (t3) then
      t3_pointer = t3_pointer + 1
      if (t3_pointer .gt. diis) t3_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t3_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(t3_diis(i),i=1,diis)
          call util_flush(LuOut)
        endif
      else
        d_t3_diis = t3_diis(t3_pointer)
      endif
      call copyfile(d_t3,d_t3_diis,size_t3)
      endif

      if (t4) then
      t4_pointer = t4_pointer + 1
      if (t4_pointer .gt. diis) t4_pointer = 1
      if (nodezero.and.util_print('diis',print_debug)) then
        write(LuOut,*)
        write(LuOut,*) 'diis_exist',(t4_diis_exist(i),i=1,diis)
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
          write(LuOut,9020) 'diis',(t4_diis(i),i=1,diis)
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
      IF (r2_pointer .lt. diis) THEN
        if (transpose) then
          shift = -1.0d0*zlshiftl
          if (t1) call tce_jacobi_l1(d_r1,d_t1,k_t1_offset,shift)
          shift = -2.0d0*zlshiftl
          if (t2) call tce_jacobi_l2(d_r2,d_t2,k_t2_offset,shift)
          shift = -3.0d0*zlshiftl
          if (t3) call tce_jacobi_l3(d_r3,d_t3,k_t3_offset,shift)
          shift = 0.0d0
          if (t4) call tce_jacobi_l4(d_r4,d_t4,k_t4_offset,shift)
        else
          shift = -1.0d0*zlshift
          if (t1) call tce_jacobi_t1(d_r1,d_t1,k_t1_offset,shift)
          shift = -2.0d0*zlshift
          if (t2) call tce_jacobi_t2(d_r2,d_t2,k_t2_offset,shift)
          shift = -3.0d0*zlshift
          if (t3) then
           if(activecalc) then
            call tce_jacobi_t3a(d_r3,d_t3,k_t3_offset,shift)
           else
            call tce_jacobi_t3(d_r3,d_t3,k_t3_offset,shift)
           end if
          end if
          shift = 0.0d0
          if (t4) call tce_jacobi_t4(d_r4,d_t4,k_t4_offset,shift)
        endif
      END IF
c
c     Make Matrix A (at the end of DIIS cycle)
c
      IF(r1_pointer.eq.diis) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!
      if(nodezero) then
      write(6,*)'MICROCYCLE DIIS UPDATE:',
     &          iter,r1_pointer
      end if 
c
        do i = 1, diis
          do j = i, diis
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
        do i = 1,diis
          a(i,diis+1) = -1.0d0
          a(diis+1,i) = -1.0d0
        enddo
        a(diis+1,diis+1) = 0.0d0
        if (nodezero.and.util_print('diis',print_debug)) then
          call ma_print(a,maxdiis+1,maxdiis+1,'A matrix')
          call util_flush(LuOut)
        endif
c
c     Make Vector B
c
        do i = 1,diis
          b(i) = 0.0d0
        enddo
        b(diis+1) = -1.0d0
        if (nodezero.and.util_print('diis',print_debug)) then
          call ma_print(b,maxdiis+1,1,'B vector')
          call util_flush(LuOut)
        endif
c
c     Solve AX = B
c
        call dgesv(diis+1,1,a,maxdiis+1,iwork,b,maxdiis+1,info)
        if (info .ne. 0)
     1    call errquit('tce_diis: LU decomposition failed',info,
     2    UNKNOWN_ERR)
        if (nodezero.and.util_print('diis',print_debug)) then
          call ma_print(b,maxdiis+1,1,'X matrix')
          call util_flush(LuOut)
        endif
c
c     Make New T amplitudes
c
        if (t1) call daxfile(diis,b,t1_diis,d_t1,size_t1)
        if (t2) call daxfile(diis,b,t2_diis,d_t2,size_t2)
        if (t3) call daxfile(diis,b,t3_diis,d_t3,size_t3)
        if (t4) call daxfile(diis,b,t4_diis,d_t4,size_t4)
        call ga_sync()
      END IF
 9010 format(1x,A,' file name    = ',A)
 9020 format(1x,A,' file handles = ',20i5)
      return
      end
