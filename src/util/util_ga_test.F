c $Id$
c vector boxes lack arithmetic precision 
#ifdef CRAY_YMP
# define THRESH 1d-10
# define THRESHF 1e-5
#elif defined(FUJITSU)
# define THRESH 1d-12
# define THRESHF 1e-5
#else
# define THRESH 1d-13
# define THRESHF 1e-5
#endif

#define MISMATCH(x,y) abs(x-y)/max(1d0,abs(x)).gt.THRESH
#define MISMATCHF(x,y) abs(x-y)/max(1.0,abs(x)).gt.THRESHF 


      subroutine util_ga_test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer heap, stack, fudge, ma_heap, me, nproc
      logical status
      parameter (heap=100*100*4, fudge=100, stack=100*100)
c     
c***  Intitialize a message passing library
c
#ifdef MPI
c      integer ierr
c      call mpi_init(ierr)
#else
c      call pbeginf
#endif
c
c***  Initialize GA
c
c     There are 2 choices: ga_initialize or ga_initialize_ltd.
c     In the first case, there is no explicit limit on memory usage.
c     In the second, user can set limit (per processor) in bytes.
c
c     call ga_initialize()
      nproc = ga_nnodes()
      me = ga_nodeid()
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      ma_heap = heap/nproc + fudge 
c     call GA_set_memory_limit(util_mdtob(ma_heap))
c
c     if(ga_nodeid().eq.0)then
c        print *,' GA initialized '
c        call ffflush(6)
c     endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
c     status = ma_init(MT_DCPL, stack, ma_heap)
c     if (.not. status) call ga_error('ma_init failed',-1) 
c
c     Uncomment the below line to register external memory allocator
c     for dynamic arrays inside GA routines.
c      call register_ext_memory()
c
c     if(me.eq.(nproc-1))then
c       print *, 'using ', nproc,' process(es) ', ga_cluster_nnodes(),
c    $           ' cluster nodes'
c       print *,'process ', me, ' is on node ',ga_cluster_nodeid(),
c    $          ' with ',  ga_cluster_nprocs(-1), ' processes' 
c       call ffflush(6)
c     endif
      if (me.eq.0) then
        write(6,'(A,I3)') ' Number of processes ..............',nproc
        write(6,'(A,I3)') ' Number of cluster nodes ..........',
     1                    ga_cluster_nnodes()
        call ffflush(6)
      endif
      call ga_sync()
c
c***  Check support for double precision arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,'(A)') ' Checking doubles  '
         write(6,*)
         call ffflush(6)
      endif

      call check_dbl()
c
c***  Check support for double precision complex arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,'(A)') ' Checking double complexes'
         write(6,*)
         call ffflush(6)
      endif

      call check_complex()
c
c***  Check support for integer arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,'(A)') ' Checking integers  '
         write(6,*)
         call ffflush(6)
      endif

      call check_int()
c
c
c***  Check support for single precision
c
      if (me.eq.0) then
         write(6,*)
         write(6,'(A)') ' Checking single precisions '
         write(6,*)
         call ffflush(6)
      endif
 
      call check_flt() 
c     
      if (me.eq.0) then
         write(6,*)
         write(6,'(A)')' Checking wrappers to MP collective operations'
         write(6,*)
         call ffflush(6)
      endif

      call check_wrappers
c
c***  Check if memory limits are enforced 
c
c     if(ga_memory_limited())
c    1  call check_mem
c
c     if(me.eq.0) call ga_print_stats()
c     if(me.eq.0) print *,' ' 
c     if(me.eq.0) print *,'All tests succesful ' 
c
c***  Tidy up the GA package
c
c     call ga_terminate()
c
c***  Tidy up after message-passing library
c
#ifdef MPI
c     call mpi_finalize(ierr)
#else
c     call pend()
#endif
c
      end


      subroutine check_dbl()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n,m
      parameter (n = 128)
      parameter (m = 2*n)
      double precision a(n,n), b(n,n), v(m),w(m)
      integer iv(m), jv(m)
      logical status
      integer g_a, g_b
      integer i,j, loop,nloop,maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, ij, inc, ii, jj
      parameter (maxloop = 100)
      integer maxproc
      parameter (maxproc = 128)
      double precision crap, sum1, sum2, x
      double precision nwords
      integer iran
      external iran
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
      nloop = Min(maxloop,n)
c     
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) =-1. 
         enddo
      enddo
*      write(6,*) ' correct '
*      call output(a, 1, n, 1, n, n, n, 1)
*      call ffflush(6)
c
c     Create a global array
c
*     print *,ga_nodeid(), ' creating array'
*     call ffflush(6)
c     call setdbg(1)
      status = ga_create(MT_DBL, n, n, 'a', 0, 0, g_a)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_create ........................ OK'
        else
        write(6,'(A)') ' ga_create ........................ Failed'
        stop
        endif
        call ffflush(6)
      endif
c
c     check if handle is valid. be quiet unless error
C
      if(.not.ga_valid_handle(g_a)) call ga_error("invalid handle",g_a)
c
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
      call ga_sync()
c     
c     Zero the array
c
      call ga_zero(g_a)
c
c     Check that it is indeed zero
c
      status = .true.     
      call ga_get(g_a, 1, n, 1, n, b, n)
      call ga_sync()
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. 0.0d0) then
               status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_zero .......................... OK'
        else
        write(6,'(A)') ' ga_zero .......................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Each node fills in disjoint sections of the array
c
      call ga_sync()
c
      status = .true.     
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
*               write(6,4) me, ilo, ihi, jlo, jhi
* 4             format(' node ',i2,' checking put ',4i4)
*               call ffflush(6)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      call util_dfill(n*n, 0.0d0, b, 1)
*     call ga_print(g_a,1)
      call ga_get(g_a, 1, n, 1, n, b, n)
*      write(6,*) ' after get'
*      call output(b, 1, n, 1, n, n, n, 1)
c
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. a(i,j)) then
               status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_put ........................... OK'
        else
        write(6,'(A)') ' ga_put ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Now check nloop random gets from each node
c
      call ga_sync()
c
      nwords = 0
c
      status = .true.
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      do loop = 1, nloop
         ilo = iran(loop)
         ihi = iran(loop)
         if (ihi.lt. ilo) then
            itmp = ihi
            ihi = ilo
            ilo = itmp
         endif
         jlo = iran(loop)
         jhi = iran(loop)
         if (jhi.lt. jlo) then
            itmp = jhi
            jhi = jlo
            jlo = itmp
         endif
c
         nwords = nwords + (ihi-ilo+1)*(jhi-jlo+1)
c     
         call util_dfill(n*n, 0.0d0, b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
         sum1 = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum1 = sum1 + b(i,j)
               if (b(i,j) .ne. a(i,j)) then
                  status = .false.
               endif
            enddo
         enddo
c
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_get ........................... OK'
        else
        write(6,'(A)') ' ga_get ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Each node accumulates into disjoint sections of the array
c
      call ga_sync()
c     
      crap = util_drand(12345)       ! Same seed for each process
      do j = 1, n
         do i = 1, n
c           b(i,j) = util_drand(0)
            b(i,j) = i+j 
         enddo
      enddo
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
c           x = util_drand(0)
            x = 10.
            ilo = i
            ihi = min(i+inc-1, n)
            if(ihi.eq.n-1)ihi=n
c           ihi = min(i+inc, n)
            jlo = j
            jhi = min(j+inc-1, n)
            if(jhi.eq.n-1)jhi=n
c           jhi = min(j+inc-1, n)
*               call ffflush(6)
            if (mod(ij,nproc) .eq. me) then
c               print *, me, 'checking accumulate ',ilo,ihi,jlo,jhi,x
* 11            format(' node ',i2,' checking accumulate ',4i4)
*               call ffflush(6)
               call ga_acc(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n, x)
            endif
            ij = ij + 1
c     
c           Each process applies all updates to its local copy
c
            do jj = jlo, jhi
               do ii = ilo, ihi
                  a(ii,jj) = a(ii,jj) + x * b(ii,jj)
               enddo
            enddo
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      status = .true.
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if(MISMATCH(b(i,j),a(i,j)))then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_acc (disjoint) ................ OK'
        else
        write(6,'(A)') ' ga_acc (disjoint) ................ Failed'
        endif
        call ffflush(6)
      endif
c
c     overlapping accumulate
      status = .true.
      call ga_sync()
      if (.not. ga_create(MT_DBL, n, n, 'b', 0, 0, g_b)) then
         status = .false.
      endif
c
      call ga_zero(g_b)
      call ga_acc(g_b, n/2, n/2, n/2, n/2, 1d0, 1, 1d0)
      call ga_sync()
      if (me.eq.0) then
         call ga_get(g_b, n/2, n/2, n/2, n/2, b(1,1), 1)
         x = abs(b(1,1) -1d0*nproc)
         if(x.gt. 1d-10)then
            status = .false.
         endif
      endif
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_acc (overlap) ................. OK'
        else
        write(6,'(A)') ' ga_acc (overlap) ................. Failed'
        endif
        call ffflush(6)
      endif
c
c     Check the ga_add function
c
      crap = util_drand(12345)       ! Everyone has same seed
      do j = 1, n
         do i = 1, n
            b(i,j) = util_drand(0)
            a(i,j) = 0.1d0*a(i,j) + 0.9d0*b(i,j)
         enddo
      enddo
      
      status = .true.
      if (me.eq.0) call ga_put(g_b, 1, n, 1, n, b, n)
      call ga_add(0.1d0, g_a, 0.9d0, g_b, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if(MISMATCH(b(i,j), a(i,j)))then
               status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_add ........................... OK'
        else
        write(6,'(A)') ' ga_add ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Check the ddot function
c
      crap = util_drand(12345)       ! Everyone has same seed
      sum1 = 0.0d0
      do j = 1, n
         do i = 1, n
            b(i,j) = util_drand(0)
            sum1 = sum1 + a(i,j)*b(i,j)
         enddo
      enddo
      if (me.eq.0) then 
         call ga_put(g_b, 1, n, 1, n, b, n)
         call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_sync()
      sum2 = ga_ddot(g_a,g_b)
      status = .true.
      if(MISMATCH(sum1, sum2))then
         status = .false.
      endif
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_ddot .......................... OK'
        else
        write(6,'(A)') ' ga_ddot .......................... Failed'
        endif
        call ffflush(6)
      endif
c
c     Check the ga_scale function
c
      call ga_scale(g_a, 0.123d0)
      call ga_get(g_a, 1, n, 1, n, b, n)
      status = .true.
      do j = 1, n
         do i = 1, n
            a(i,j) = a(i,j)*0.123d0
            if (MISMATCH(b(i,j), a(i,j)))then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scale ......................... OK'
        else
        write(6,'(A)') ' ga_scale ......................... Failed'
        endif
        call ffflush(6)
      endif
c
c     Check the ga_copy function
c
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      call ga_copy(g_a, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      status = .true.
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_copy .......................... OK'
        else
        write(6,'(A)') ' ga_copy .......................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_sync()
c
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      status = .true.
      do j = 1, 10 
       call ga_sync()
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
         do loop = 1,m 
           ilo = iran(n) 
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
         enddo 
         call ga_gather(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,v(loop),1)
           if(v(loop)  .ne. a(ilo,jlo))then
             status = .false.
           endif
         enddo
       endif
      enddo
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_gather ........................ OK'
        else
        write(6,'(A)') ' ga_gather ........................ Failed'
        endif
        call ffflush(6)
      endif
c
      status = .true.
      do j = 1,10 
       call ga_sync()
       if(me.eq.iran(ga_nnodes())-1) then
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
c          v(loop) = DSIN(a(ilo,jlo)+b(ilo,jlo)) 
           v(loop) = 1d0 *(ilo+jlo)
         enddo 
         call ga_scatter(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop) 
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
c          if(v(loop) .ne. w(loop))then 
           if(w(loop) .ne. 1d0 *(ilo+jlo) )then 
             status = .false.
           endif
         enddo 
       endif
       call ga_sync()
      enddo
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scatter ....................... OK'
        else
        write(6,'(A)') ' ga_scatter ....................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_sync()
c
c     scatter-acc available in GA ver. 3.0
#ifdef GA3
c
      crap = util_drand(1234)
      call ga_zero(g_a)
c
      do j = 1, n
         do i = 1, n
            b(i,j) =0.
         enddo
      enddo
c
      status = .true.
      x = .1d0
      ii =n
      do jj = 1,1 
	  call ga_sync()
        do loop = 1, ii 

c          generate unique i,j pairs
10         continue
              i = iran(n)
	      j=iran(n)
           if (found(i,j, iv, jv, loop-1) ) goto 10

           iv(loop) = i
           jv(loop) = j
           v(loop) = 1d0 *(i+j)
           b(i,j) = b(i,j) + nproc*x*v(loop) ! update local ref. array
        enddo
	call ga_scatter_acc(g_a,v,iv,jv, ii,x)
c
        call ga_sync()
c     
c       check the result
c
        call ga_get(g_a, 1, n, 1,n, a, n)

        do loop = 1,ii
          i = iv(loop)
          j = jv(loop)
          if(MISMATCH(a(i,j),b(i,j)))then
             status = .false.
*            if(me.eq.0)then
*              do ii=1,loop
*                   print *,'element',ii, iv(ii),jv(ii)
*              enddo
*            endif
             status = .false.
          endif
        enddo
        call ga_sync()
      enddo

      call ga_sync()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scatter_acc ................... OK'
        else
        write(6,'(A)') ' ga_scatter_acc ................... Failed'
        endif
        call ffflush(6)
      endif
#endif
c
c     Delete the global arrays
c
      status = ga_destroy(g_b)
      status = status .and. ga_destroy(g_a)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_destroy ....................... OK'
        else
        write(6,'(A)') ' ga_destroy ....................... Failed'
        endif
        call ffflush(6)
      endif
c     
      end

c-----------------------------------------------------------------
      subroutine check_complex()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 60)
      parameter (m = 2*n)
      double complex a(n,n), b(n,n), v(m),w(m)
      integer iv(m), jv(m)
      logical status
      integer g_a, g_b
      integer i,j, loop,nloop,maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, ij, inc, ii, jj
      parameter (maxloop = 100)
      integer maxproc
      parameter (maxproc = 128)
      double precision crap
      double precision nwords
      double complex   x, sum1, sum2, factor
      integer iran
      external iran
c
      nproc = ga_nnodes()
      me    = ga_nodeid()
      nloop = Min(maxloop,n)
c
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = cmplx(dble(i-1), dble((j-1)*n))
            b(i,j) = cmplx(-1d0,1d0)
         enddo
      enddo
c
c     Create a global array
c
c     print *,ga_nodeid(), ' creating array'
      call ffflush(6)
c     call setdbg(1)
      status = ga_create(MT_DCPL, n, n, 'a', 0, 0, g_a)
      status = status .and. ga_create(MT_DCPL, n, n, 'b', 0, 0, g_b)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_create ........................ OK'
        else
        write(6,'(A)') ' ga_create ........................ Failed'
        endif
        call ffflush(6)
      endif
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
      call ga_sync()
c
c     Zero the array
c
      call ga_zero(g_a)
c
c     Check that it is indeed zero
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      call ga_sync()
      status = .true.
      do i = 1, n
         do j = 1, n
            if(b(i,j).ne.(0d0,0d0)) then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_zero .......................... OK'
        else
        write(6,'(A)') ' ga_zero .......................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Each node fills in disjoint sections of the array
c
      call ga_sync()
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      call util_qfill(n*n, (0d0,0d0), b, 1)
      call ga_get(g_a, 1, n, 1, n, b, n)
c
      status = .true.
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. a(i,j)) then
               status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_put ........................... OK'
        else
        write(6,'(A)') ' ga_put ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Now check nloop random gets from each node
c
      call ga_sync()
c
      nwords = 0
c
      status = .true.
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      do loop = 1, nloop
         ilo = iran(loop)
         ihi = iran(loop)
         if (ihi.lt. ilo) then
            itmp = ihi
            ihi = ilo
            ilo = itmp
         endif
         jlo = iran(loop)
         jhi = iran(loop)
         if (jhi.lt. jlo) then
            itmp = jhi
            jhi = jlo
            jlo = itmp
         endif
c
         nwords = nwords + (ihi-ilo+1)*(jhi-jlo+1)
c
         call util_qfill(n*n, (0.0d0,0d0), b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
         do j = jlo, jhi
            do i = ilo, ihi
               if (b(i,j) .ne. a(i,j)) then
                 status = .false.
               endif
            enddo
         enddo
c
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_get ........................... OK'
        else
        write(6,'(A)') ' ga_get ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Each node accumulates into disjoint sections of the array
c
      call ga_sync()
c
      status = .true.
      crap = util_drand(12345)       ! Same seed for each process
      do j = 1, n
         do i = 1, n
            b(i,j) = cmplx(util_drand(0),util_drand(1))
         enddo
      enddo
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
c           x = cmplx(util_drand(0),0.333d0)
c           x = cmplx(0.333d0,0)
*           x = cmplx(0d0,0d0)
            x = 0 
            ilo = i
            ihi = min(i+inc-1, n)
            if(ihi.eq.n-1)ihi=n
            jlo = j
            jhi = min(j+inc-1, n)
            if(jhi.eq.n-1)jhi=n
            if (mod(ij,nproc) .eq. me) then
               call ga_acc(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n, x)
            endif
            ij = ij + 1
c
c           Each process applies all updates to its local copy
c
            do jj = jlo, jhi
               do ii = ilo, ihi
                  a(ii,jj) = a(ii,jj) + x * b(ii,jj)
               enddo
            enddo
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (MISMATCH(b(i,j), a(i,j)))then
               status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_acc (disjoint) ................ OK'
        else
        write(6,'(A)') ' ga_acc (disjoint) ................ Failed'
        endif
        call ffflush(6)
      endif
c
c     overlapping accumulate
c
      call ga_zero(g_b)
      call ga_acc(g_b, n/2, n/2, n/2, n/2, (1d0,-1d0), 1, (1d0,0d0))
      call ga_sync()
      status = .true.
      if (me.eq.0) then
         call ga_get(g_b, n/2, n/2, n/2, n/2, x, 1)
c        error = abs(x -(1d0,-1d0)*nproc)
         if (MISMATCH(x, ((1d0,-1d0)*nproc)))then
c        if(error.gt. (1d-8))then
           status = .false.
         endif
      endif
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_acc (overlap) ................. OK'
        else
        write(6,'(A)') ' ga_acc (overlap) ................. Failed'
        endif
        call ffflush(6)
      endif
c
c     Check the ga_copy function
c
      status = .true.
      call ga_sync()
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      call ga_copy(g_a, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_copy .......................... OK'
        else
        write(6,'(A)') ' ga_copy .......................... Failed'
        endif
        call ffflush(6)
      endif
c
c
c     Check the ga_scale function
c
      factor = (1d0,-1d0)
      call ga_scale(g_a, factor)
      call ga_get(g_a, 1, n, 1, n, b, n)
      status = .true.
      do j = 1, n
         do i = 1, n
            a(i,j) = a(i,j)*factor
            if (MISMATCH(b(i,j), a(i,j)))then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scale ......................... OK'
        else
        write(6,'(A)') ' ga_scale ......................... Failed'
        endif
        call ffflush(6)
      endif
c
c     Check scatter&gather
c
      call ga_sync()
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
c
      status = .true.
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      do j = 1, 10
       call ga_sync()
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
         enddo
         call ga_gather(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,v(loop),1)
           if(v(loop)  .ne. a(ilo,jlo))then
             status = .false.
           endif
         enddo
       endif
      enddo
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_gather ........................ OK'
        else
        write(6,'(A)') ' ga_gather ........................ Failed'
        endif
        call ffflush(6)
      endif
c
      status = .true.
      do j = 1,10
       call ga_sync()
       if(me.eq.iran(ga_nnodes())-1) then
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
           v(loop) = (1d0,-1d0) *(ilo+jlo)
         enddo
         call ga_scatter(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
           if(w(loop)  .ne. (1d0,-1d0) *(ilo+jlo) )then
             status = .false.
           endif
         enddo
       endif
       call ga_sync()
      enddo
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scatter ....................... OK'
        else
        write(6,'(A)') ' ga_scatter ....................... Failed'
        endif
        call ffflush(6)
      endif
c
c     Check ga_add
c
      call ga_get(g_a, 1, n, 1, n, a, n)
      crap = util_drand(12345)       ! Everyone has same seed
      do j = 1, n
         do i = 1, n
            b(i,j) = cmplx(util_drand(0), util_drand(1))
            a(i,j) = (0.1d0,-.1d0)*a(i,j) + (.9d0,-.9d0)*b(i,j)
         enddo
      enddo
      status = .true.
      if (me.eq.0) call ga_put(g_b, 1, n, 1, n, b, n)
      call ga_add((0.1d0,-.1d0), g_a, (0.9d0,-.9d0), g_b, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (MISMATCH(b(i,j), a(i,j)))then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_add ........................... OK'
        else
        write(6,'(A)') ' ga_add ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Check the zdot function
c
      crap = util_drand(12345)       ! Everyone has same seed
      sum1 = (0.0d0,0.d0)
      do j = 1, n
         do i = 1, n
            b(i,j) = cmplx(util_drand(0), util_drand(1))
            sum1 = sum1 + a(i,j)*b(i,j)
         enddo
      enddo
      if (me.eq.0) then
         call ga_put(g_b, 1, n, 1, n, b, n)
         call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_sync()
      sum2 =  ga_zdot(g_a,g_b)
      status = .true.
      if (MISMATCH(sum1, sum2))then
        status = .false.
      endif
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_zdot .......................... OK'
        else
        write(6,'(A)') ' ga_zdot .......................... Failed'
        endif
        call ffflush(6)
      endif
c
c     Delete the global arrays
c
      status = ga_destroy(g_b)
      status = status .and. ga_destroy(g_a)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_destroy ....................... OK'
        else
        write(6,'(A)') ' ga_destroy ....................... Failed'
        endif
        call ffflush(6)
      endif
c
      end
c-----------------------------------------------------------------




      subroutine check_int()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n
      parameter (n = 128)
      integer a(n,n), b(n,n)
      logical status
      integer g_a
      integer i, j, loop, nloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, ij, inc, dimi,dimj,iproc, ii, jj
      double precision nwords
      parameter (nloop = 100)
      integer maxproc
      parameter (maxproc = 128)
      integer map(5,maxproc), found, np,k
      double precision crap, sum1
      integer buf
      integer iran
      external iran
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*1000
         enddo
      enddo
c     
c     Create a global array
c
      status = ga_create(MT_INT, n, n, 'a', 0, 0, g_a)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_create ........................ OK'
        else
        write(6,'(A)') ' ga_create ........................ Failed'
        endif
        call ffflush(6)
      endif
c     
c     Zero the array
c
      call ga_zero(g_a)
c
c     Check that it is indeed zero
c
      status = .true.
      call ga_get(g_a, 1, n, 1, n, b, n)
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. 0) then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_zero .......................... OK'
        else
        write(6,'(A)') ' ga_zero .......................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Each node fills in disjoint sections of the array
c
      call ga_sync()
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
c              write(6,4) me, ilo, ihi, jlo, jhi
c4             format(' node ',i2,' checking put ',4i4)
c              call ffflush(6)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      status = .true.
      if(me.eq.0)then
      call ga_get(g_a, 1, n, 1, n, b, n)
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. a(i,j)) then
              status = .false.
            endif
         enddo
      enddo
      endif
      call ga_sync()
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_put ........................... OK'
        else
        write(6,'(A)') ' ga_put ........................... Failed'
        endif
        call ffflush(6)
      endif
c
c     Now check nloop random gets from each node
c
      call ga_sync()
c
      nwords = 0
c
      status = .true.
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      do loop = 1, nloop
         ilo = iran(loop)
         ihi = iran(loop)
         if (ihi.lt. ilo) then
            itmp = ihi
            ihi = ilo
            ilo = itmp
         endif
         jlo = iran(loop)
         jhi = iran(loop)
         if (jhi.lt. jlo) then
            itmp = jhi
            jhi = jlo
            jlo = itmp
         endif
c     
         nwords = nwords + (ihi-ilo+1)*(jhi-jlo+1)
c
         call util_ifill(n*n, 0, b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
c
         sum1 = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum1 = sum1 + b(i,j)
               if (b(i,j) .ne. a(i,j)) then
                 status = .false.
               endif
            enddo
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_get ........................... OK'
        else
        write(6,'(A)') ' ga_get ........................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_sync()
c
      status = .true.
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      inc =5 
c     every processor will be operating on somebody elses data
c
      iproc = ga_nnodes()-me-1
c
      call ga_distribution(g_a,iproc,ilo,ihi,jlo,jhi)
c
      dimi = ihi-ilo
      dimj = jhi-jlo
c     write(6,*) me,'..',ilo,ihi,jlo,jhi,'.',dimi,dimj
c     call ffflush(6)
      call ga_sync()
      if(ilo .gt.0 .and. jhi .gt. 0)then
       do loop = 1,nloop 
         ii= IABS(iran(dimi))
         jj= IABS(iran(dimj))
         i=ilo + Mod(ii,dimi)
         j=jlo + Mod(jj,dimj)
c
c        write(6,*) me,'..',ilo,ihi,jlo,jhi,'.',dimi,dimj,'..',i,j
c        call ffflush(6)
         buf = ga_read_inc(g_a,i,j,inc)
         if(a(i,j).ne. buf)then   
           status = .false.
         endif
         call ga_get(g_a, i,i,j,j, buf,1) 
         a(i,j) = a(i,j)+inc 
         if(a(i,j).ne.  buf)then   
           status = .false.
         endif
       enddo
      endif
      call ga_sync()
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_read_inc ...................... OK'
        else
        write(6,'(A)') ' ga_read_inc ...................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_zero(g_a)
c
c***  use ga_read_inc and elements g_a(1:2,1) to implement a lock
c***  compute g_a(:,n) = sum (1(:)) for P processors
c     
      status = ga_create_mutexes(1)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_create_mutexes ................ OK'
        else
        write(6,'(A)') ' ga_create_mutexes ................ Failed'
        endif
        call ffflush(6)
      endif

      if ((n.lt.2).and.(me.eq.0)) then
        write(6,'(A)') ' ga_fence ........................ N/A'
        write(6,'(A)') ' ga_lock ......................... N/A'
        call ffflush(6)
      endif

      call ga_lock(0)
c     call my_lock(g_a)

c          get original values  g_a(:,n)
           call ga_get(g_a, 1,n, n,n, b,n)
c          add my contribution
           do i =1,n
              b(i,1)= b(i,1)+1
           enddo
c
c          need to use fence to assure that coms complete before leaving
c          Critical Section
c
           call ga_init_fence()
                call ga_put(g_a, 1,n, n,n, b,n)
           call ga_fence()
      call ga_unlock(0)
c     call my_unlock(g_a)
c
333   status = ga_destroy_mutexes()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_destroy_mutexes ............... OK'
        else
        write(6,'(A)') ' ga_destroy_mutexes ............... Failed'
        endif
        call ffflush(6)
      endif
c
      status = .true.
      call ga_sync()
      if (me.eq.0) then
         call ga_get(g_a, 1,n, n,n, b,n)
         do i =1,n
              if(b(i,1).ne. nproc)then
                status = .false.
              endif
         enddo
      endif 
c
      status = ga_locate_region(g_a, 1, n, 1,n, map,np)
      found = 0
      do j=1,n
         do i=1,n
            b(i,j)=-1
         enddo
      enddo
      if(me.eq.0)call ga_put(g_a,1,n,1,n,b,n)
      call ga_sync()
      do k = 1, np
         if(map(5,k).eq.me)then
                if(found.eq.1) then
                    write(6,*)'double entry in map for proc ',me
                    call ffflush(6)
                endif
                do j= map(3,k), map(4,k)
                        do i= map(1,k), map(2,k)
                        b(i,j)=1*me
                        enddo
                enddo
                call ga_put(g_a, map(1,k),map(2,k),map(3,k),map(4,k),
     &                 b(map(1,k),map(3,k)),n)
                found = 1
         endif
      enddo
      call ga_sync()
c
      do k = 1, np
         if(map(5,k).eq.me)then
                call ga_get(g_a, map(1,k),map(2,k),map(3,k),map(4,k),
     &                 a(map(1,k),map(3,k)),n)
                do j= map(3,k), map(4,k)
                      do i= map(1,k), map(2,k)
                        if(b(i,j).ne.a(i,j)) then
                           write(6,*)
     &                     'proc ',me, 'overlap with ',a(i,j)
                           call ffflush(6)
                        endif
                      enddo
                enddo
         endif
      enddo
      call ga_sync()
c
      if(me.eq.0)then
        call ga_get(g_a,1,n,1,n,a,n)
        do j=1,n
         do i=1,n
            if(a(i,j).eq.-1)then
              status = .false.
            endif
         enddo
        enddo
      endif
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_locate_region ................. OK'
        else
        write(6,'(A)') ' ga_locate_region ................. Failed'
        endif
        call ffflush(6)
      endif
c     
c     Delete the global array
c
      status = ga_destroy(g_a)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_destroy ....................... OK'
        else
        write(6,'(A)') ' ga_destroy ....................... Failed'
        endif
        call ffflush(6)
      endif
c     
      end

c---------------------------------------------------------------------

      subroutine check_flt()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer n, m
      parameter (n =10) 
      parameter (m=2*n)
      real a(n,n), b(n,n), v(m), w(m)
      integer iv(m), jv(m)
      logical status
      integer g_a, g_b
      integer i, j, loop, nloop, maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, ij, inc, ii, jj
      double precision nwords
      parameter (maxloop = 100)
      integer maxproc
      parameter (maxproc = 128)
      double precision crap
      real x, sum1, sum2
      integer iran
      external iran

      nproc = ga_nnodes()
      me    = ga_nodeid()
      nloop = Min(maxloop,n)
c     
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) = -1.
         enddo
      enddo
c     
c     Create a global array
c
      status = ga_create(MT_REAL, n, n, 'a', 0, 0, g_a)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_create ........................ OK'
        else
        write(6,'(A)') ' ga_create ........................ Failed'
        endif
        call ffflush(6)
      endif
c
c     check if handle is valid. be quiet unless error
c
      status = .true.
      if(.not.ga_valid_handle(g_a)) status = .false.
c
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
      call ga_sync()      
c     
c     Zero the array
c
      call ga_zero(g_a)
c
c     Check that it is indeed zero
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      call ga_sync()
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. 0.0) then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_zero .......................... OK'
        else
        write(6,'(A)') ' ga_zero .......................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Each node fills in disjoint sections of the array
c
      status = .true.
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
c               write(6,4) me, ilo, ihi, jlo, jhi
c 4             format(' node ',i2,' checking put ',4i4)
c               call ffflush(6)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. a(i,j)) then
              status = .false.
            endif
         enddo
      enddo
      call ga_sync()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_put ........................... OK'
        else
        write(6,'(A)') ' ga_put ........................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_sync()
c
      nwords = 0
c
      status = .true.
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      do loop = 1, nloop
         ilo = iran(loop)
         ihi = iran(loop)
         if (ihi.lt. ilo) then
            itmp = ihi
            ihi = ilo
            ilo = itmp
         endif
         jlo = iran(loop)
         jhi = iran(loop)
         if (jhi.lt. jlo) then
            itmp = jhi
            jhi = jlo
            jlo = itmp
         endif
c     
         nwords = nwords + (ihi-ilo+1)*(jhi-jlo+1)
c
         call util_rfill(n*n, 0.0, b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
c
         sum1 = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum1 = sum1 + b(i,j)
               if (b(i,j) .ne. a(i,j)) then
                  status = .false.
               endif
            enddo
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_get ........................... OK'
        else
        write(6,'(A)') ' ga_get ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
c     Each node accumulates into disjoint sections of the array
c
      call ga_sync()
c
      crap = util_drand(12345)       ! Same seed for each process
      do j = 1, n
         do i = 1, n
c           b(i,j) = real(util_drand(0))
            b(i,j) = i+j
         enddo
      enddo
c
      status = .true.
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
c           x = real(util_drand(0))
            x = 10.
            ilo = i
            ihi = min(i+inc-1, n)
            if(ihi.eq.n-1)ihi=n
c           ihi = min(i+inc, n)
            jlo = j
            jhi = min(j+inc-1, n)
            if(jhi.eq.n-1)jhi=n
c           jhi = min(j+inc-1, n)
*               call ffflush(6)
            if (mod(ij,nproc) .eq. me) then
c               print *, me, 'checking accumulate ',ilo,ihi,jlo,jhi,x
* 11            format(' node ',i2,' checking accumulate ',4i4)
*               call ffflush(6)
               call ga_acc(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n, x)
            endif
            ij = ij + 1
c
c           Each process applies all updates to its local copy
c
            do jj = jlo, jhi
               do ii = ilo, ihi
                  a(ii,jj) = a(ii,jj) + x * b(ii,jj)
               enddo
            enddo
         enddo
      enddo
      call ga_sync()                       
c
c     All nodes check all of a
      call ga_get(g_a, 1, n, 1, n, b, n)
c
      do j = 1, n
         do i = 1, n
            if(MISMATCHF(b(i,j),a(i,j)))then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_acc (disjoint) ................ OK'
        else
        write(6,'(A)') ' ga_acc (disjoint) ................ Failed'
        endif
        call ffflush(6)
      endif
c
c     overlapping accumulate
      call ga_sync()
      status = .true.
      if (.not. ga_create(MT_REAL, n, n, 'b', 0, 0, g_b)) then
         status = .false.
      endif
c
      call ga_zero(g_b)
      call ga_acc(g_b, n/2, n/2, n/2, n/2, 1.0, 1, 1.0)
      call ga_sync()
      if (me.eq.0) then
         call ga_get(g_b, n/2, n/2, n/2, n/2, b(1,1), 1)
         x = abs(b(1,1) -1*nproc)
         if(x.gt. 1e-10)then
           status = .false.
         endif
      endif
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_acc (overlap) ................. OK'
        else
        write(6,'(A)') ' ga_acc (overlap) ................. Failed'
        endif
        call ffflush(6)
      endif
c
c     Check the ga_add function
c
      crap = util_drand(12345)       ! Everyone has same seed
      status = .true.
      do j = 1, n
         do i = 1, n
            b(i,j) = real(util_drand(0)*real(i)) + 1
            a(i,j) = 0.1*a(i,j) + 0.9*b(i,j)
         enddo
      enddo
      if (me.eq.0) call ga_put(g_b, 1, n, 1, n, b, n)
      call ga_add(0.1, g_a, 0.9, g_b, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if(MISMATCHF(b(i,j), a(i,j)))then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_add ........................... OK'
        else
        write(6,'(A)') ' ga_add ........................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
      status = .true.
      crap = util_drand(12345)       ! Everyone has same seed
      sum1 = 0.0
      do j = 1, n
         do i = 1, n
            b(i,j) = util_drand(0)
            sum1 = sum1 + a(i,j)*b(i,j)
         enddo
      enddo
      if (me.eq.0) then
         call ga_put(g_b, 1, n, 1, n, b, n)
         call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_sync()
      sum2 = ga_sdot(g_a,g_b)
      if(MISMATCHF(sum1, sum2))then
        status = .false.
      endif
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_sdot .......................... OK'
        else
        write(6,'(A)') ' ga_sdot .......................... Failed'
        endif
        call ffflush(6)
      endif
c
      status = .true.
      call ga_scale(g_a, 0.123)
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            a(i,j) = a(i,j)*0.123
            if (MISMATCHF(b(i,j), a(i,j)))then
              status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scale ......................... OK'
        else
        write(6,'(A)') ' ga_scale ......................... Failed'
        endif
        call ffflush(6)
      endif
c                                
      status = .true.
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      call ga_copy(g_a, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               status = .false.
            endif
         enddo
      enddo
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_copy .......................... OK'
        else
        write(6,'(A)') ' ga_copy .......................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_sync()
      status = .true.
      crap = util_drand(ga_nodeid()*51 + 1) ! Different seed for each process
      do j = 1, 10
       call ga_sync()
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
         enddo
         call ga_gather(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,v(loop),1)
           if(v(loop)  .ne. a(ilo,jlo))then
             status = .false.
           endif                                     
         enddo
       endif
      enddo
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_gather ........................ OK'
        else
        write(6,'(A)') ' ga_gather ........................ Failed'
        endif
        call ffflush(6)
      endif
c
      status = .true.
      do j = 1,10
       call ga_sync()
       if(me.eq.iran(ga_nnodes())-1) then
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
c          v(loop) = DSIN(a(ilo,jlo)+b(ilo,jlo))
           v(loop) = 1.0 *(ilo+jlo)
         enddo
         call ga_scatter(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
c          if(v(loop)  .ne. w(loop))then
           if(w(loop)  .ne. 1.0 *(ilo+jlo) )then
             status = .false.
           endif
         enddo
       endif
       call ga_sync()
      enddo
c
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scatter ....................... OK'
        else
        write(6,'(A)') ' ga_scatter ....................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_sync()
c
c     scatter-acc available in GA ver. 3.0
#ifdef GA3                                         
      status = .true.
      crap = util_drand(1234)
      call ga_zero(g_a)
c
      do j = 1, n
         do i = 1, n
            b(i,j) =0.
         enddo
      enddo
c
      x = .1d0
      ii =n
      do jj = 1,1
          call ga_sync()
        do loop = 1, ii
 
c          generate unique i,j pairs
10         continue
              i = iran(n)
              j=iran(n)
           if (found(i,j, iv, jv, loop-1) ) goto 10

           iv(loop) = i
           jv(loop) = j
           v(loop) = 1.0 *(i+j)
           b(i,j) = b(i,j) + nproc*x*v(loop) ! update local ref. array
        enddo
 
        call ga_scatter_acc(g_a,v,iv,jv, ii,x)
        call ga_sync()
c
c       check the result
c
        call ga_get(g_a, 1, n, 1,n, a, n)
        do loop = 1,ii
          i = iv(loop)
          j = jv(loop)
          if(MISMATCH(a(i,j),b(i,j)))then
            status = .false.
          endif
        enddo                                 
        call ga_sync()
      enddo
      call ga_sync()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_scatter_acc ................... OK'
        else
        write(6,'(A)') ' ga_scatter_acc ................... Failed'
        endif
        call ffflush(6)
      endif
#endif                      
c
c     Delete the global array
c
      status = ga_destroy(g_a)
      status = status .and. ga_destroy(g_b)
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_destroy ....................... OK'
        else
        write(6,'(A)') ' ga_destroy ....................... Failed'
        endif
        call ffflush(6)
      endif
c     
      end
c_____________________________________________________________

      subroutine check_wrappers
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      double precision sum
      integer isum, ibuf(2)
      integer me, nproc
      logical status
      real fsum
c
      nproc = ga_nnodes()
      me = ga_nodeid()
c
      status = .true.
      call ga_sync()
      ibuf(1) = 1
      ibuf(2) =  me
      call ga_igop(10000, ibuf, 2, '+')
      if(ibuf(1).ne.nproc)then
        status = .false.
      endif
      if(ibuf(2).ne.((nproc-1)*nproc/2))then
        status = .false.
      endif
      call ga_sync()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_igop .......................... OK'
        else
        write(6,'(A)') ' ga_igop .......................... Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
c
      status = .true.
      sum = 1d0 * me
      call ga_dgop(10000, sum, 1, '+')
      if(Int(sum).ne.((nproc-1)*nproc/2))then
        status = .false. 
      endif
      call ga_sync()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_dgop .......................... OK'
        else
        write(6,'(A)') ' ga_dgop .......................... Failed'
        endif
        call ffflush(6)
      endif
c
      call ga_sync()
      status = .true.
      fsum = 1.0 * me
      call ga_sgop(10000, fsum, 1, '+')
      if(Int(sum).ne.((nproc-1)*nproc/2))then
        status = .false.
      endif
      call ga_sync()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_sgop .......................... OK'
        else
        write(6,'(A)') ' ga_sgop .......................... Failed'
        endif
        call ffflush(6)
      endif
c                                   
      call ga_sync()
      status = .true.
      if(me.eq.nproc-1)then
        ibuf(1) = me
        ibuf(2) = nproc
      endif
      call ga_brdcst(1000,ibuf,util_mitob(2),nproc-1)
      if(ibuf(1).ne.nproc-1) status = .false.
      if(ibuf(2).ne.nproc) status = .false.
      call ga_sync()
      if (me.eq.0) then
        if (status) then
        write(6,'(A)') ' ga_brdcst ........................ OK'
        else
        write(6,'(A)') ' ga_brdcst ........................ Failed'
        endif
        call ffflush(6)
      endif
      call ga_sync()
      end


      subroutine check_mem
      implicit none
      integer mem_size
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,nmax,left,need, me,procs,g_a, g_b 
      integer stack, heap, global
      logical status, overify, ohardfail
c
      call input_mem_size(stack, heap, global, overify, ohardfail)
      mem_size = ma_sizeof(mt_dbl,global,mt_byte)
      write(*,*) 'mem_size = ',mem_size
      me = ga_nodeid()
      procs = ga_nnodes()
      nmax = int(dsqrt(dble(mem_size/util_mitob(1))))
      left =  mem_size/procs - ga_inquire_memory()
      n = nmax/2
      need = util_mdtob(n*n)/procs
c
      if(me.eq.0)then
        write(6,*)' '
        if(ga_uses_ma())then
           write(6,*)' CHECKING GA MEMORY RESTRICTIONS (MA used)'
        else
           write(6,*)' CHECKING GA MEMORY RESTRICTIONS (MA not used)'
        endif
        write(6,*)' '
        write(6,*)' '
        call print_mem_info(n,left, need, mem_size/procs)
      endif
c
      status = ga_create(MT_DBL, n, n, 'a', 0, 0, g_a)
c
      if(me.eq.0) then
        if(status) then
              write(6,*) '  success'
        else
              write(6,*) '  failure'
        endif
        call ffflush(6)
      endif
c
      n = nmax
      left =  mem_size/procs - ga_inquire_memory()
      need = util_mdtob(n*n)/procs
      if(me.eq.0)then
        call print_mem_info(n,left, need, mem_size/procs)
      endif
c
      status = ga_create(MT_DBL, n, n, 'b', 0, 0, g_b)
c
      if(me.eq.0) then
        if(status) then
              write(6,*) '  success'
        else
              write(6,*) '  failure'
        endif
        write(6,*)' '
        write(6,*)' '
        call ffflush(6)
      endif
      status = ga_destroy(g_a)
      end

     

      subroutine print_mem_info(n,left, need, total)
      implicit none
      integer n,left, need, total
c
      write(6,*)' '
      if(left - need .ge. 0) then
        write(6,1)n,n 
1       format('> Creating array ',i4,' by ',i4,' -- should succeed')
      else
        write(6,2)n,n 
2       format('> Creating array ',i4,' by ',i4,' -- SHOULD FAIL')
      endif
      write(6,3) need, left, total
3     format('  (need ',i7,' and  ',i7,' out of ',i7,' bytes are left)')
      write(6,*)' '
      call ffflush(6)
c
      end


      
      subroutine my_lock(g_b)
      implicit none
#include "global.fh"
      integer g_b, val, flag, i
      logical first_time
      double precision  dummy 
      common /lock/ val 
      common /dum/ dummy
      data first_time /.true./
c
c     this awkward initialization is to avoid a weird problem 
C     with block data on SUN
      if(first_time)then
        first_time = .false.
        dummy = .0
      endif
c
      val = ga_read_inc(g_b,1,1, 1)
10    call ga_get(g_b, 2,2,1,1, flag, 1)
           if(flag.eq.val) return 
c
c          to reduce memory stress,  wait a while before retrying
           do i = 1, 100
              dummy = dummy + .1 
           enddo 
      goto 10
      end


      subroutine my_unlock(g_b)
      implicit none
#include "global.fh"
      integer g_b, val
      common /lock/ val
c
      call ga_put(g_b, 2,2,1,1, val+1, 1) 
      end


      logical function found(i,j, iv, jv, n)
      integer n
      integer i,j, iv(n), jv(n)
      integer loop

      found = .false.
      do loop = 1, n
               if(i .eq. iv(loop) .and. j .eq.jv(loop))then
                 found = .true.
                 goto 99
               endif
            enddo
99    continue
      return
      end 
      
            
      subroutine proc_remap()
      implicit none
#include "global.fh"
      integer proc(100),nproc,i
      nproc = ga_nnodes()
      if(nproc.gt.100)
     $   call ga_error("remap requires<=100 processes",nproc)
      do i = 1, nproc
         proc(i) = nproc-i
      enddo
c     call ga_register_proclist(proc,nproc)
      end

      
      subroutine util_rfill(n,val,a,ia)
      implicit none
      real  a(*), val
      integer n, ia, i
c
c     initialise real array to scalar value
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end

      
      subroutine util_dfill(n,val,a,ia)
      implicit none
      double precision  a(*), val
      integer n, ia, i
c
c     initialise double precision array to scalar value
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end

      subroutine util_ifill(n,val,a,ia)
      implicit none
      integer n, ia, i, a(*),val
c
c     initialise integer array to scalar value
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end

      subroutine util_qfill(n,val,a,ia)
      implicit none
      double  complex  a(*), val
      integer n, ia, i
c
c     initialise double complex array to scalar value
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end


      integer function iran(i)
      implicit none
      double precision util_drand
      external util_drand
      integer i
      iran = int(util_drand(0)*dble(i))+1
      return
      end
