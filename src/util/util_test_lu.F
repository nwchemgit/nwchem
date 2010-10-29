c $Id$
      subroutine util_test_lu
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap, stack
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
c     Intitialize the GA package
c
c     call ga_initialize()
c     if(ga_nodeid().eq.0)print *,ga_nnodes(),' nodes'
c
c     Initialize the MA package
c
      heap = 190000
      stack= 190000
c     if (.not. ma_init(MT_DBL, heap, stack))
c    $    call ga_error("ma init failed",heap+stack) 
c
c
      call testit6()
c     call ga_terminate()
c
#ifdef MPI
c     call mpi_finalize(ierr)
#else
c     call pend()
#endif
      end


c-----------------

      subroutine testit6()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c     
      integer n
      parameter (n = 100)
      double precision a(n,n), b(n,n), c(n,n)
      integer g_a,g_b,g_c,g_d, g_e, g_f, g_g
      integer  i, j
      integer nproc, me
      double precision dsin, sum
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c     a() is a local copy of what the global array should start as
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
      if (.not. ga_create(MT_DBL, n, n, 'a', 1, 1, g_a))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'b', 1, 1, g_b))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'c', 1, 1, g_c))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'd', 1, 1, g_d))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, 1, 'e', 1, 1, g_e))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, 1, 'f', 1, 1, g_f))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, 1, 'g', 1, 1, g_g))
     $     call ga_error(' ga_create failed ',2)
c     
c     
c***  Fill in arrays A & B
      if (me .eq. 0) then
c        print *,  ' filling in A and B  '
         call ffflush(6)
	 call ga_put(g_e, 1,n, 1,1, b(1,1),n)
	 call ga_put(g_f, 1,n, 1,1, b(1,1),n)
      endif
      do j = 1+me, n, nproc 
	    call ga_put(g_a, 1,n, j,j, a(1,j),n)
	    call ga_put(g_b, 1,n, j,j, b(1,j),n)
	    call ga_put(g_c, 1,n, j,j, b(1,j),n)
      enddo
c
c     call ga_copy(g_b,g_c)
c
c     if (me .eq. 0) then
c       print *,' '
c       print *, '>Test of the LU-based solver with nxn rhs '
c       print *,' '
c       call ffflush(6)
c     endif
#ifndef SCALAPACK
      call ga_lu_solve_seq('n', g_a, g_b)
#else
      call ga_lu_solve('n', g_a, g_b)
#endif
      call ga_dgemm('n','n',n,n,n, 1d0, g_a, g_b, 0d0, g_d) ! d := a*b
      call ga_add(1d0, g_d, -1d0, g_c, g_c) 
      sum = ga_ddot(g_c,g_c)
      if (me .eq. 0) then
c	print *,' '
c       print *, ' norm = ', dsqrt(sum)/n
        if(dsqrt(sum)/n.lt.1d-10) 
     1  write(6,'(A)') ' ga_lu_solve (multiple RHS) ....... OK'
c	print *,' '
        call ffflush(6)
      endif
c
c     if (me .eq. 0) then
c       print *,' '
c       print *,'>Test of the LU-based solver with a single vector rhs'
c       print *,' '
c       call ffflush(6)
c     endif
c
#ifndef SCALAPACK
      call ga_lu_solve_seq('n', g_a, g_e)
#else
      call ga_lu_solve('n', g_a, g_e)
#endif
c
      call ga_dgemm('n','n',n,1,n, 1d0, g_a, g_e, 0d0, g_g) ! g := a*e
      call ga_add(1d0, g_g, -1d0, g_f, g_f) 
      sum = ga_ddot(g_f,g_f)
      if (me .eq. 0) then
c       print *,' '
c       print *, ' norm = ', dsqrt(sum)/n
        if(dsqrt(sum)/n.lt.1d-10)
     1  write(6,'(A)') ' ga_lu_solve (single RHS) ......... OK'
c       print *,' '
        call ffflush(6)
      endif
      if (.not. ga_destroy(g_a))
     $     call ga_error(' ga_destroy failed ',2)
      if (.not. ga_destroy(g_b))
     $     call ga_error(' ga_destroy failed ',2)
      if (.not. ga_destroy(g_c))
     $     call ga_error(' ga_destroy failed ',2)
      if (.not. ga_destroy(g_d))
     $     call ga_error(' ga_destroy failed ',2)
      if (.not. ga_destroy(g_e))
     $     call ga_error(' ga_destroy failed ',2)
      if (.not. ga_destroy(g_f))
     $     call ga_error(' ga_destroy failed ',2)
      if (.not. ga_destroy(g_g))
     $     call ga_error(' ga_destroy failed ',2)
      end

