c $Id$
      subroutine util_test_eig
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c     integer heap, stack
c     
c
c***  Intitialize a message passing library
c
#ifdef MPI
c     integer ierr
c     call mpi_init(ierr)
#else
c     call pbeginf
#endif
c
c***  initialize PEIGS
#ifdef PAR_DIAG
c      call mxinit()   ! PEIGS needs mxinit
#endif
c
c     Intitialize the GA package
c
c     call ga_initialize()
c
c     Initialize the MA package
c
c     heap = 190000
c     stack= 190000
c     if (.not. ma_init(MT_DBL, heap, stack))
c    $    call ga_error("ma init failed",heap+stack)
c
      call testit5()
c     call ga_terminate()
c
#ifdef MPI
c     call mpi_finalize(ierr)
#else
c     call pend()
#endif
      end


c-----------------

      subroutine testit5()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c     
      integer n
      parameter (n = 100)
      double precision a(n,n), b(n,n), c(n,n), evals(n)
      double precision eva(n), evb(n)
      integer g_a,g_b,g_c,g_d
      integer  i, j, index(n),ind(n)
      integer nproc, me
      double precision dsin, sum
      logical status
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c***  a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = 1d0 * (i+j)  
            b(i,j) = DSIN(1d0* (i+j))   
	    if(i.eq.j) then
               b(i,j) = 2d0 *n
               a(i,j) = i
            endif
	    if(i.le.j)then
               c(i,j) = a(i,j)
            else
               c(i,j) = 0d0
            endif
         enddo
      enddo
c
c***  Create global arrays
c
      if (.not. ga_create(MT_DBL, n, n, 'A', 1, 1, g_a))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'B', 1, 1, g_b))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'C', 1, 1, g_c))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'D', 1, 1, g_d))
     $     call ga_error(' ga_create failed ',2)
c     
c***  Fill in arrays A & B
c
      if (me .eq. 0) then
c        write(6,21) 
c21      format(/' filling in ... ')
c        call ffflush(6)
         do j = 1, n
	    call ga_put(g_a, 1,n, j,j, a(1,j),n)
	    call ga_put(g_b, 1,n, j,j, b(1,j),n)
	    call ga_put(g_c, 1,n, j,j, c(1,j),n)
         enddo
c	print *,'A'
        do j = 1, n
       	  call GA_GET(g_a, 1,n, j,j, eva,1)
c         write(6,'(10e8.2)')(eva(i),i=1,n)
        enddo
      endif
c
c***  check symmetrization
c
c     if (me .eq. 0) then
c       print *,' '
c	print *,'>checking ga_symmetrize '
c       print *,' '
c       call ffflush(6)
c     endif
      call ga_symmetrize(g_c)
c
      call GA_GET(g_c,  1,n, 1,n,c,n)
      do j = ga_nodeid()+1, n, ga_nnodes()
         do i = j+1, n
            if(c(i,j).ne.c(j,i))then
                 print *, me, ' symmetrize ',i,j,c(i,j),c(j,i)
                 call ffflush(6)
                 call ga_error('exiting',-1)
            endif
         enddo
      enddo
      call ga_sync()
      if (me .eq. 0) then
        write(6,'(A)') ' ga_symmetrize .................... OK'
c       print *,' '
c       print *,' ga_symmetrize is OK'
c       print *,' '
        call ffflush(6)
      endif

c
c***  check symmetrization 
c
c     if (me .eq. 0) then
c	print *,' '
c       print *,'>checking ga_transpose '
c	print *,' '
c       call ffflush(6)
c     endif
c
      call ga_transpose(g_c,g_d)
*     call ga_print(g_c)
*     call ga_print(g_d)
      call GA_GET(g_d,  1,n, 1,n,a,n)
      do j = ga_nodeid()+1, n, ga_nnodes()
         call GA_GET(g_d,  1,n, j,j, a,n)
         do i = 1, n
            if(a(i,1).ne.c(j,i))then
                 print *, me, ' transpose ',i,j,a(i,1),c(j,i) 
                 call ffflush(6)
                 call ga_error('exiting',-1)
            endif
         enddo
      enddo
      call ga_sync()
      if (me .eq. 0) then
        write(6,'(A)') ' ga_transpose ..................... OK'
c	print *,' '
c	print *,' ga_transpose is OK'
c	print *,' '
        call ffflush(6)
      endif
c
c
c***  solve the eigenproblem
c     if (me .eq. 0)then
c	print *,' '
c       write(6,*) '>checking the generalized eigensolver ... '
c	print *,' '
c       call ffflush(6)
c     endif
      call ga_sync()
#ifndef PAR_DIAG
      call ga_diag_seq(g_a,g_b,g_c,evals)
#else
      call ga_diag(g_a,g_b,g_c,evals)
#endif
c     if (me .eq. 0) then
c	print *,' '
c	print *,' checking multiplication'
c	print *,' '
c       call ffflush(6)
c     endif
c
      call ga_sync()
      call ga_dgemm('t','n',n,n,n, 1d0, g_c, g_a, 0d0, g_d)
      call ga_dgemm('n','n',n,n,n, 1d0, g_d, g_c, 0d0, g_a)
c
      call ga_sync()
      if (me .eq. 0) then
         do j = 1, n
	    call GA_GET(g_a, j,j, j,j, eva(j),1)
         enddo
      endif

      call ga_sync()
      call ga_dgemm('t','n',n,n,n, 1d0, g_c, g_b, 0d0, g_d)
      call ga_sync()
      call ga_dgemm('n','n',n,n,n, 1d0, g_d, g_c, 0d0, g_b)
c
      call ga_sync()
      if (me .eq. 0) then
         do j = 1, n
	    call GA_GET(g_b, j,j, j,j, evb(j),1)
         enddo
c        write(6,*)'  j   lambda      eva      evb      eva/evb' 
c        write(6,*)'  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' 
         call ffflush(6)
         do j = 1, n
            if(ABS(evals(j) - eva(j)/evb(j)).gt.1d-5)
     &         write(6,'(i4,1h_,6(e10.3,1h,))')
     @            j,evals(j), eva(j), evb(j),eva(j)/evb(j)
         enddo
c        write(6,*)'       OK         OK       OK        OK' 
        write(6,'(A)') ' ga_dgemm ......................... OK'
        write(6,'(A)') ' ga_diag .......................... OK'
         call ffflush(6)
      endif
      if (me .eq. 0) then
c	print *,' '
c	print *,'  eigensolver & multiplication are OK'
c	print *,' '
c	print *,' '
        call ffflush(6)
      endif
c
c..................................................................
c
c***  solve the std eigenproblem
c     if (me .eq. 0)then
c       print *,' '
c       write(6,*) '>checking the standard eigensolver ... '
c       print *,' '
c       call ffflush(6)
c     endif
      do j =1,n
         index(j)=j
         ind(j)=j
      enddo

      call ga_sync()
#ifdef PAR_DIAG
      call ga_diag_std(g_a,g_c,evals)
#else
      call ga_diag_std_seq(g_a,g_c,evals)
#endif
c
      call ga_sync()

      call ga_zero(g_b)
      
      call ga_dgemm('n','n',n,n,n, 1d0, g_a, g_c, 0d0, g_d) ! d := a*c
c
c
      if (me .eq. 0) call ga_scatter(g_b, evals, index, ind, n)

      call ga_sync()
      call ga_dgemm('n','n',n,n,n, 1d0, g_c, g_b, 0d0, g_a) ! a := c*b
 
      call ga_sync()
      call ga_add(1d0, g_d, -1d0, g_a, g_c)
      sum = ga_ddot(g_c,g_c)
      if (me .eq. 0) then
c       print *,' '
        if(dsqrt(sum)/n.lt.1d-11)then
        write(6,'(A)') ' ga_diag_std ...................... OK'
c          print *, ' std eigensolver is OK '
        else
           print *, ' test failed: norm = ', dsqrt(sum)/n
        endif
c       print *,' '
        call ffflush(6)
      endif
c     status =  MA_summarize_allocated_blocks()
      status =  ga_destroy(g_d)
      status =  ga_destroy(g_c)
      status =  ga_destroy(g_b)
      status =  ga_destroy(g_a)
      end

