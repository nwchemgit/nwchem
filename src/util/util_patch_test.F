c $Id$
c
c
      subroutine util_patch_test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      logical status
      integer me
c
c***  Initialize a message passing library
c
#ifdef MPI
c     integer ierr
c     call mpi_init(ierr)
#else
c     call pbeginf
#endif
c
c     call ga_initialize()
c     if(ga_nodeid().eq.0)then
c        print *,' GA initialized'
c        call ffflush(6)
c     endif
c
c     status = ma_init(MT_DBL, 50000, 90000/ga_nnodes())
c     if (.not. status)call ga_error( 'ma_init failed', -1)
      me = ga_nodeid()
      if (me.eq.0) then
         write(6,'(A)') ' Checking doubles  '
         write(6,*)
         call ffflush(6)
      endif
      call dpatch_test()
c
      if (me.eq.0) then
         write(6,*)
         write(6,'(A)') ' Checking double complexes  '
         write(6,*)
         call ffflush(6)
      endif
      call zpatch_test()
c 
      if (me.eq.0) then
         write(6,*)
         write(6,'(A)') ' Checking single precisions '
         write(6,*)
         call ffflush(6)
      endif
      call spatch_test()
c
c     call ga_terminate()
c
#ifdef MPI
c     call mpi_finalize(ierr)
#else
c     call pend()
#endif
      end



      subroutine dpatch_test()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 128)
      parameter (m = n*n)
      double precision a(n,n), b(n,n), c(n,n), buf(m), val
      double precision alpha, beta
      double precision crap
      integer nproc, me 
      integer i, j, ailo, ajlo, bilo, bjlo, base
      integer aihi, ajhi, bihi, bjhi
      integer g_a, g_b, g_c 
      integer rows, cols, loop
      logical status, dist_same 
      integer iran
      external iran
      dist_same = .false.
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) = i+j 
         enddo
      enddo
c
c***  Create a global array
c
      status = ga_create(MT_DBL, n, n, 'a', 0, 0, g_a)
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_create ........................ OK'
        call ffflush(6)
      endif
c
      if(dist_same) then
         status = ga_duplicate(g_a, g_b, 'a_duplicated')
         if(.not.ga_compare_distr(g_a, g_b))
     $           call ga_error("g_b distribution different",0) 
         status = ga_duplicate(g_a, g_c, 'a_duplicated_again')
         if(.not.ga_compare_distr(g_a, g_c))
     $           call ga_error("g_c distribution different",0) 
      else
         status = ga_create(MT_DBL, n, n, 'b', 0, n, g_b)
         if (.not. status) call ga_error('ga_create failed:b',0) 
         status = ga_create(MT_DBL, n, n, 'c', n, 0, g_c)
         if (.not. status) call ga_error('ga_create failed:c',0) 
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_fill_patch ... '
c        call ffflush(6)
c     endif
c
      val = 1d0
      call ga_fill_patch(g_a, 2,n/2, 2,n, val)
*     call ga_print(g_a, 1)
      do j = 2+me, n, nproc 
         call ga_get(g_a, 1,n/2, j,j, buf,n/2) 
         do i = 2, n/2
            if(buf(i) .ne.val ) then
               print *,me, ' error ',i,j, buf(i),val 
               call ga_error('exiting ...',0)
            endif
         enddo
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_fill_patch .................... OK'
        call ffflush(6)
      endif
c  
      crap = util_drand(12345) ! Same seed for each process
      do j = 1+me, n, nproc 
         call ga_put(g_a,1,n,j,j,a(1,j),n) 
      enddo
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_copy_patch ... '
c        call ffflush(6)
c     endif
c
      do loop =1, 10
              ailo = iran(n/2)
              ajlo = iran(n/2)
              aihi = min(n, -1+ailo+n/2)
              ajhi = min(n, -1+ajlo+n/4)
              rows = aihi -ailo+1
              cols = ajhi -ajlo +1
c
              bilo = iran(n/3)
              bjlo = iran(n/3)
              bihi = bilo + rows -1
              bjhi = bjlo + cols -1
c             if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c
c               call ffflush(6)
c             endif

c
c--------------------------------------------------------------
c             call ga_print(g_a,1)
c             call ga_get(g_a,ailo,aihi, ajlo, ajhi, buf, 1)
c             call ga_put(g_b,bilo,bihi, bjlo, bjhi, buf, 1)
c             call ga_print(g_b,1)
c             call ga_zero(g_b)
c--------------------------------------------------------------
              call ga_copy_patch('n', g_a, ailo, aihi, ajlo, ajhi,
     &                                g_b, bilo, bihi, bjlo, bjhi)
c             call ga_print(g_a,1)
c             call ga_print(g_b,1)
              call ga_get(g_b,bilo,bihi, bjlo, bjhi, buf, rows)
              base = 0
              do j = ajlo, ajhi
                 if(Mod(j,nproc).eq.me) then
                    do i = ailo, aihi
                       base = base+1
                       if(buf(base) .ne. a(i,j)) then
                          print *,me, ' error ',i,j, buf(base), a(i,j)
                          call ga_error('exiting ...',0)
                       endif
                    enddo
                 else
                    base = base + rows
                 endif
              enddo
      enddo
c
      ailo = iran(n/2)
      ajlo = iran(n/2)
      bilo = iran(n/2)
      bjlo = iran(n/2)
c
c     if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c     endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1) 
      base = 0
      do j = ajlo, ajlo+n/3 
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2 
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  without transpose:   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_copy_patch .................... OK'
        call ffflush(6)
      endif
c 
c     if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,2h]~,5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c     endif
      call ga_copy_patch('t', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do i = ailo, ailo+n/2
         if(Mod(j,nproc).eq.me) then
         do j = ajlo, ajlo+n/3
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/3+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  transposed:   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_copy_patch (transpose) ........ OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_scale_patch ... '
c        call ffflush(6)
c     endif
      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/2, bjlo, bjlo+n/3)
      val = 1.d0
      call ga_scale_patch(g_b, bilo,bilo+n/2,bjlo, bjlo+n/3,val)
      call ga_get(g_b,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)*val) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_scale_patch ................... OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_add_patch ... '
c        call ffflush(6)
c     endif
      alpha = .1d0
      beta = .2d0
      call ga_zero(g_c)
      call ga_add_patch(alpha, g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                   beta, g_b, bilo,bilo+n/2, bjlo, bjlo+n/3,
     $                         g_c, bilo,bilo+n/2, bjlo, bjlo+n/3)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      val = val*beta + alpha
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(ABS(buf(base)- a(i,j)*val).gt.1d-5) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_add_patch ..................... OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_ddot_patch ... '
c        call ffflush(6)
c     endif
      alpha= ga_ddot_patch(g_a,'n', ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                     g_c,'n', bilo,bilo+n/2, bjlo, bjlo+n/3)
      beta = 0d0
      do j = ajlo, ajlo+n/3
            do i = ailo, ailo+n/2
               beta = beta + a(i,j)**2
            enddo
      enddo

      if(ABS(beta*val- alpha).gt.1d-6*alpha) then
             print *,me, ' error ', beta*val, alpha
             call ga_error('exiting ...',0)
      endif
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_ddot_patch .................... OK'
        call ffflush(6)
      endif
c
c......................................................
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_matmul_patch ... '
c        call ffflush(6)
c     endif
      do j = 1+me, n, nproc
         call ga_put(g_b,1,n,j,j,b(1,j),n)
      enddo
      call ga_sync()
      call ga_matmul_patch('n','n', 1d0, 0d0, 
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call DGEMM('n','n',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
     $            b(bilo,bjlo),n, 0d0, c, n)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  a*b: OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch .................. OK'
        call ffflush(6)
      endif
c
      call ga_sync()
      call ga_matmul_patch('t','n', 1d0, 0d0,
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call DGEMM('t','n',n/2+1,n/2+1,n/3+1,1d0,a(ajlo,ailo), n,
     $            b(bilo,bjlo),n, 0d0, c, n)
*     call ga_print(g_a,1) 
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  trans(a)*b: OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch (transposed A) ... OK'
        call ffflush(6)
      endif
c
      call ga_sync()
      call ga_matmul_patch('n','t', 1d0, 0d0,
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call DGEMM('n','t',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
     $            b(bjlo,bilo),n, 0d0, c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  a*trans(b): OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch (transposed B) ... OK'
        call ffflush(6)
      endif
c
      call ga_sync()
      call ga_matmul_patch('t','t', 1d0, 0d0,
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call DGEMM('t','t',n/2+1,n/2+1,n/3+1,1d0,a(ajlo,ailo), n,
     $            b(bjlo,bilo),n, 0d0, c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  trans(a)*trans(b): OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch (transposed A & B) OK'
        call ffflush(6)
      endif
      status = ga_destroy(g_a)
      status = status .and. ga_destroy(g_b)
      status = status .and. ga_destroy(g_c)
      if(.not. status) print *, 'ga_destroy failed'
      if (me.eq.0) then
        write(6,'(A)') ' ga_destory ....................... OK'
        call ffflush(6)
      endif
c
      end





      subroutine zpatch_test()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 128)
      parameter (m = n*n)
      double complex a(n,n), b(n,n), c(n,n), buf(m),  val
      double complex alpha, beta
      integer nproc, me 
      integer i, j, ailo, ajlo, bilo, bjlo, base
      integer aihi, ajhi, bihi, bjhi
      integer g_a, g_b, g_c 
      integer rows, cols, loop
      logical status, dist_same 
      integer iran
      external iran
      dist_same = .false.
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      do j = 1, n
         do i = 1, n
            a(i,j) = cmplx(dble(i-1), dble((j-1)*n))
            b(i,j) = cmplx(dble(i+j),1d0)
         enddo
      enddo
c
c***  Create a global array
c
      status = ga_create(MT_DCPL, n, n, 'a', 0, 0, g_a)
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_create ........................ OK'
        call ffflush(6)
      endif
c
      if(dist_same) then
         status = ga_duplicate(g_a, g_b, 'a_duplicated')
         if(.not.ga_compare_distr(g_a, g_b))
     $           call ga_error("g_b distribution different",0) 
         status = ga_duplicate(g_a, g_c, 'a_duplicated_again')
         if(.not.ga_compare_distr(g_a, g_c))
     $           call ga_error("g_c distribution different",0) 
      else
         status = ga_create(MT_DCPL, n, n, 'b', 0, n, g_b)
         if (.not. status) call ga_error('ga_create failed:b',0) 
         status = ga_create(MT_DCPL, n, n, 'c', n, 0, g_c)
         if (.not. status) call ga_error('ga_create failed:c',0) 
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_fill_patch ... '
c        call ffflush(6)
c     endif
c
      val = (1d0,-1d0)
      call ga_fill_patch(g_a, 2,n/2, 2,n, val)
*     call ga_print(g_a, 1)
      do j = 2+me, n, nproc 
         call ga_get(g_a, 1,n/2, j,j, buf,n/2) 
         do i = 2, n/2
            if(buf(i) .ne.val ) then
               print *,me, ' error ',i,j, buf(i),val 
               call ga_error('exiting ...',0)
            endif
         enddo
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_fill_patch .................... OK'
        call ffflush(6)
      endif
c  
      do j = 1+me, n, nproc 
         call ga_put(g_a,1,n,j,j,a(1,j),n) 
      enddo
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_copy_patch ... '
c        call ffflush(6)
c     endif
c
      do loop =1, 10
              ailo = iran(n/2)
              ajlo = iran(n/2)
              aihi = min(n, -1+ailo+n/2)
              ajhi = min(n, -1+ajlo+n/4)
              rows = aihi -ailo+1
              cols = ajhi -ajlo +1
c
              bilo = iran(n/3)
              bjlo = iran(n/3)
              bihi = bilo + rows -1
              bjhi = bjlo + cols -1
c             if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c
c               call ffflush(6)
c             endif

c
              call ga_copy_patch('n', g_a, ailo, aihi, ajlo, ajhi,
     &                                g_b, bilo, bihi, bjlo, bjhi)
c             call ga_print(g_a,1)
c             call ga_print(g_b,1)
              call ga_get(g_b,bilo,bihi, bjlo, bjhi, buf, rows)
              base = 0
              do j = ajlo, ajhi
                 if(Mod(j,nproc).eq.me) then
                    do i = ailo, aihi
                       base = base+1
                       if(buf(base) .ne. a(i,j)) then
                          print *,me, ' error ',i,j, buf(base), a(i,j)
                          call ga_error('exiting ...',0)
                       endif
                    enddo
                 else
                    base = base + rows
                 endif
              enddo
      enddo
c
      ailo = iran(n/2)
      ajlo = iran(n/2)
      bilo = iran(n/2)
      bjlo = iran(n/2)
c
c     if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c     endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1) 
      base = 0
      do j = ajlo, ajlo+n/3 
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2 
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  without transpose:   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_copy_patch .................... OK'
        call ffflush(6)
      endif
c 
c     if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,2h]~,5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c     endif
      call ga_copy_patch('t', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do i = ailo, ailo+n/2
         if(Mod(j,nproc).eq.me) then
         do j = ajlo, ajlo+n/3
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/3+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  transposed:   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_copy_patch (transpose) ........ OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_scale_patch ... '
c        call ffflush(6)
c     endif
      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/2, bjlo, bjlo+n/3)
      val = 1.d0
      call ga_scale_patch(g_b, bilo,bilo+n/2,bjlo, bjlo+n/3,val)
      call ga_get(g_b,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)*val) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_scale_patch ................... OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_add_patch ... '
c        call ffflush(6)
c     endif
      alpha = .1d0
      beta = .2d0
      call ga_zero(g_c)
      call ga_add_patch(alpha, g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                   beta, g_b, bilo,bilo+n/2, bjlo, bjlo+n/3,
     $                         g_c, bilo,bilo+n/2, bjlo, bjlo+n/3)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      val = val*beta + alpha
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(ABS(buf(base)- a(i,j)*val).gt.1d-5) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_add_patch ..................... OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_zdot_patch ... '
c        call ffflush(6)
c     endif
      alpha= ga_zdot_patch(g_a,'n', ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                     g_c,'n', bilo,bilo+n/2, bjlo, bjlo+n/3)
      beta = (0d0,0d0)
      do j = ajlo, ajlo+n/3
            do i = ailo, ailo+n/2
               beta = beta + a(i,j)*a(i,j)
            enddo
      enddo
      if(ABS(beta*val- alpha)/abs(alpha).gt.1d-6) then
             print *,me, ' error ', beta*val, alpha
             call ga_error('exiting ...',0)
      endif
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_zdot_patch .................... OK'
        call ffflush(6)
      endif
c
c......................................................
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_matmul_patch ... '
c        call ffflush(6)
c     endif
      do j = 1+me, n, nproc
         call ga_put(g_b,1,n,j,j,b(1,j),n)
      enddo
      call ga_sync()
      call ga_matmul_patch('n','n', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ZGEMM('n','n',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ailo,ajlo), n,
     $            b(bilo,bjlo),n, (0d0,0d0), c, n)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  a*b: OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch .................. OK'
        call ffflush(6)
      endif
c
      call ga_sync()
      call ga_matmul_patch('t','n', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ZGEMM('t','n',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ajlo,ailo), n,
     $            b(bilo,bjlo),n, (0d0,0d0), c, n)
*     call ga_print(g_a,1) 
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  trans(a)*b: OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch (transposed A) ... OK'
        call ffflush(6)
      endif
c
      call ga_sync()
      call ga_matmul_patch('n','t', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ZGEMM('n','t',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ailo,ajlo), n,
     $            b(bjlo,bilo),n, (0d0,0d0), c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  a*trans(b): OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch (transposed B) ... OK'
        call ffflush(6)
      endif
c
      call ga_sync()
      call ga_matmul_patch('t','t', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ZGEMM('t','t',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ajlo,ailo), n,
     $            b(bilo,bjlo),n, (0d0,0d0), c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'  trans(a)*trans(b): OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_matmul_patch (transposed A & B) OK'
        call ffflush(6)
      endif
      status = ga_destroy(g_a)
      status = status .and. ga_destroy(g_b)
      status = status .and. ga_destroy(g_c)
      if(.not. status) print *, 'ga_destroy failed'
      if (me.eq.0) then
        write(6,'(A)') ' ga_destroy ....................... OK'
        call ffflush(6)
      endif
c
      end



      subroutine spatch_test()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 128)
      parameter (m = n*n)
      real a(n,n), b(n,n), c(n,n), buf(m), val
      real alpha, beta
      integer nproc, me
      integer i, j, ailo, ajlo, bilo, bjlo, base
      integer aihi, ajhi, bihi, bjhi
      integer g_a, g_b, g_c
      integer rows, cols, loop
      logical status, dist_same
      integer iran
      external iran
      dist_same = .false.
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) = i+j
         enddo
      enddo
c
c***  Create a global array
c
      status = ga_create(MT_REAL, n, n, 'a', 0, 0, g_a)
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_create ........................ OK'
        call ffflush(6)
      endif
c
      if(dist_same) then
         status = ga_duplicate(g_a, g_b, 'a_duplicated')
         if(.not.ga_compare_distr(g_a, g_b))
     $           call ga_error("g_b distribution different",0)
         status = ga_duplicate(g_a, g_c, 'a_duplicated_again')
         if(.not.ga_compare_distr(g_a, g_c))
     $           call ga_error("g_c distribution different",0)
      else
         status = ga_create(MT_REAL, n, n, 'b', 0, n, g_b)
         if (.not. status) call ga_error('ga_create failed:b',0)
         status = ga_create(MT_REAL, n, n, 'c', n, 0, g_c)
         if (.not. status) call ga_error('ga_create failed:c',0)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_fill_patch ... '
c        call ffflush(6)
c     endif
c
      val = 1.0
      call ga_fill_patch(g_a, 2,n/2, 2,n, val)
*     call ga_print(g_a, 1)
      do j = 2+me, n, nproc
         call ga_get(g_a, 1,n/2, j,j, buf,n/2)
         do i = 2, n/2
            if(buf(i) .ne.val ) then
               print *,me, ' error ',i,j, buf(i),val
               call ga_error('exiting ...',0)
            endif
         enddo
      enddo
c
      call ga_sync()
c     if (me .eq. 0) then
c        write(6,*)'   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_fill_patch .................... OK'
        call ffflush(6)
      endif
c 
      do j = 1+me, n, nproc
         call ga_put(g_a,1,n,j,j,a(1,j),n)
      enddo
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_copy_patch ... '
c        call ffflush(6)
c     endif
c
      do loop =1, 10
              ailo = iran(n/2)
              ajlo = iran(n/2)
              aihi = min(n, -1+ailo+n/2)
              ajhi = min(n, -1+ajlo+n/4)
              rows = aihi -ailo+1
              cols = ajhi -ajlo +1
c
              bilo = iran(n/3)
              bjlo = iran(n/3)
              bihi = bilo + rows -1
              bjhi = bjlo + cols -1
c             if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi

c               call ffflush(6)
c             endif

c
              call ga_copy_patch('n', g_a, ailo, aihi, ajlo, ajhi,
     &                                g_b, bilo, bihi, bjlo, bjhi)
c             call ga_print(g_a,1)
c             call ga_print(g_b,1)
              call ga_get(g_b,bilo,bihi, bjlo, bjhi, buf, rows)
              base = 0
              do j = ajlo, ajhi
                 if(Mod(j,nproc).eq.me) then
                    do i = ailo, aihi
                       base = base+1
                       if(buf(base) .ne. a(i,j)) then
                          print *,me, ' error ',i,j, buf(base), a(i,j)
                          call ga_error('exiting ...',0)
                       endif
                    enddo
                 else
                    base = base + rows
                 endif
              enddo
      enddo
c
      ailo = iran(n/2)
      ajlo = iran(n/2)
      bilo = iran(n/2)
      bjlo = iran(n/2)
c
c     if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c     endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  without transpose:   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_copy_patch .................... OK'
        call ffflush(6)
      endif
c
c     if (me .eq. 0) then
c               write(6,'(2x,1h[,4i4,2h]~,5h --> ,1h[,4i4,1h])')
c    $                ailo,aihi,ajlo,ajhi,
c    $                bilo,bihi, bjlo, bjhi
c     endif
      call ga_copy_patch('t', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do i = ailo, ailo+n/2
         if(Mod(j,nproc).eq.me) then
         do j = ajlo, ajlo+n/3
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/3+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  transposed:   OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_copy_patch (transpose) ........ OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_scale_patch ... '
c        call ffflush(6)
c     endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/2, bjlo, bjlo+n/3)
      val = 1.0
      call ga_scale_patch(g_b, bilo,bilo+n/2,bjlo, bjlo+n/3,val)
      call ga_get(g_b,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)*val) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_scale_patch ................... OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_add_patch ... '
c        call ffflush(6)
c     endif
      alpha = 0.1
      beta = 0.2
      call ga_zero(g_c)

      call ga_add_patch(alpha, g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                   beta, g_b, bilo,bilo+n/2, bjlo, bjlo+n/3,
     $                         g_c, bilo,bilo+n/2, bjlo, bjlo+n/3)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      val = val*beta + alpha
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(ABS(buf(base)- a(i,j)*val).gt.1e-3) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_add_patch ..................... OK'
        call ffflush(6)
      endif
c
c***
c     if (me .eq. 0) then
c        print *, ' '
c        write(6,*)'> Checking ga_sdot_patch ... '
c        call ffflush(6)
c     endif
      alpha= ga_sdot_patch(g_a,'n', ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                     g_c,'n', bilo,bilo+n/2, bjlo, bjlo+n/3)
      beta = 0.0

      do j = ajlo, ajlo+n/3
            do i = ailo, ailo+n/2
               beta = beta + a(i,j)**2
            enddo
      enddo
      if(ABS(beta*val- alpha).gt.(1.0e-2*alpha)) then
             print *,me, ' error ', beta*val, alpha
             call ga_error('exiting ...',0)
      endif
      call ga_sync()
c
c     if (me .eq. 0) then
c        write(6,*)'  OK '
c        call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A)') ' ga_sdot_patch .................... OK'
        call ffflush(6)
      endif
c......................................................
ccc      if (me .eq. 0) then
ccc         print *, ' '
ccc         write(6,*)'> Checking ga_matmul_patch ... '
ccc         call ffflush(6)
ccc      endif
ccc      do j = 1+me, n, nproc
ccc         call ga_put(g_b,1,n,j,j,b(1,j),n)
ccc      enddo
ccc      call ga_sync()
ccc      call ga_matmul_patch('n','n', 1.0, 0.0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call DGEMM('n','n',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
ccc     $            b(bilo,bjlo),n, 0d0, c, n)
ccc*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  a*b: OK '
ccc         call ffflush(6)
ccc      endif
cccc
ccc      call ga_sync()
ccc      call ga_matmul_patch('t','n', 1d0, 0d0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call DGEMM('t','n',n/2+1,n/2+1,n/3+1,1.0,a(ajlo,ailo), n,
ccc     $            b(bilo,bjlo),n, 0.0, c, n)
ccc*     call ga_print(g_a,1)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  trans(a)*b: OK '
ccc         call ffflush(6)
ccc      endif
cccc
ccc      call ga_sync()
ccc      call ga_matmul_patch('n','t', 1.0, 0.0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call DGEMM('n','t',n/2+1,n/2+1,n/3+1,1.0,a(ailo,ajlo), n,
ccc     $            b(bjlo,bilo),n, 0.0, c, n)
ccc*     call ga_print(g_a,1)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  a*trans(b): OK '
ccc         call ffflush(6)
ccc      endif
cccc
ccc      call ga_sync()
ccc      call ga_matmul_patch('t','t', 1d0, 0d0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call DGEMM('t','t',n/2+1,n/2+1,n/3+1,1.0,a(ajlo,ailo), n,
ccc     $            b(bjlo,bilo),n, 0.0, c, n)
ccc*     call ga_print(g_a,1)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  trans(a)*trans(b): OK '
ccc         call ffflush(6)
ccc      endif
      status = ga_destroy(g_a)
      status = status .and. ga_destroy(g_b)
      status = status .and. ga_destroy(g_c)
      if(.not. status) print *, 'ga_destroy failed'
      if (me.eq.0) then
        write(6,'(A)') ' ga_destory ....................... OK'
        call ffflush(6)
      endif
      return
      end


