c $Id$
C>
C> \defgroup kain Krylov subspace Accelerated Inexact Newton method (KAIN)
C>
C> \ingroup kain
C> @{
C>
C> \file ga_it2.F
C> \brief Implementation of Krylov-subspace nonlinear equation solver
C>
C> This files provides an implementation of the KAIN method to solve
C> linear and nonlinear systems of equations [1]. The solvers 
C> `ga_kain` and `ga_lkain` are general purpose routines. They require
C> a routine to calculate a matrix-vector product and a routine
C> to apply a preconditioner as arguments.
C>
C> [1] R.J. Harrison,
C>     <i>"Krylov Subspace Accelerated Inexact Newton Method for Linear
C>     and Nonlinear Systems",</i>
C>     Journal of Computational Chemistry (2004) <b>25</b>, pp 328-334,
C>     DOI: <a href="https://doi.org/10.1002/jcc.10108">10.1002/jcc.10108</a>
C>
C> @}
c
C> \ingroup kain
C> @{
C>
C> \brief A trivial matrix-vector product routine
C>
C> This routine is simply there to test the code for system solvers.
C> The matrix-vector product it calculates is simply a product involving
C> an explicitly stored matrix. This matrix is passed through a common
C> block.
C>
      subroutine test_product(acc,g_x, g_Ax)
      implicit none
#include "errquit.fh"
      double precision acc !< [Input] Accuracy (not used here)
      integer g_x  !< [Input] Global array with vector
      integer g_Ax !< [Output] Global array with matrix-vector product
c
      integer g_a, g_b
      common /testme/ g_A, g_B
      integer n, nvec, typex, typeax
      double complex one, zero
c
      one  = cmplx(1.0d0,0.0d0)
      zero = cmplx(0.0d0,0.0d0)
c
      call ga_inquire(g_Ax, typeax, n, nvec)
      call ga_inquire(g_x, typex, n, nvec)
      if (typex.ne.typeax)
     $  call errquit("test_product: g_x not same type as g_Ax",
     $               typex,UERR)
      call ga_zero(g_Ax)
      call ga_dgemm('n', 'n', n, nvec, n, one, g_A, g_x, zero, g_Ax)
c
      end
C>
C> \brief A trivial non-linear matrix-vector product routine
C>
C> This routine is simply there to test the code for system solvers.
C> The matrix-vector product it calculates is simply a product involving
C> an explicitly stored matrix and it adds another vector to achieve
C> something non-linear. The additional matrix and vector are passed
C> through a common block.
C>
      subroutine test_nlproduct(acc,g_x, g_Ax)
      implicit none
#include "errquit.fh"
      double precision acc
      integer g_x, g_Ax
c
      integer g_a, g_b
      common /testme/ g_A, g_B
      integer n, nvec, typex, typeax, typeb
      double complex mone, one, zero
c
      mone = cmplx(-1.0d0,0.0d0)
      one  = cmplx(1.0d0,0.0d0)
      zero = cmplx(0.0d0,0.0d0)
c
      call ga_inquire(g_b, typeb, n, nvec)
      call ga_inquire(g_Ax, typeax, n, nvec)
      call ga_inquire(g_x, typex, n, nvec)
      if (typex.ne.typeax)
     $  call errquit("test_nlproduct: g_x not same type as g_Ax",
     $               typex,UERR)
      if (typex.ne.typeb)
     $  call errquit("test_nlproduct: g_x not same type as g_b",
     $               typex,UERR)
      call ga_zero(g_Ax)
      call ga_dgemm('n', 'n', n, nvec, n, one, g_A, g_x, zero, g_Ax)
      call ga_add(one, g_AX, mone, g_B, g_AX)
c
      end
C>
C> \brief A trivial preconditioner
C>
C> This routine is simply there to test the code for system solvers.
C> It simply scales the elements of a vector by
C> \f{eqnarray*}{
C>    x_i = \frac{x_i}{i - \epsilon}
C> \f}
C> where \f$ \epsilon \f$ is the shift.
C>
      subroutine test_precond(g_x,shift)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
      integer g_x
      double precision shift
      integer n, nvec, type
      integer i, ivec
      double precision x
      double complex zx
c
      call ga_inquire(g_x, type, n, nvec)
      if (ga_nodeid() .eq. 0) then
         if (type.eq.MT_DBL) then
            do ivec = 1, nvec
               do i = 1, n
                  call ga_get(g_x, i, i, ivec, ivec, x, 1)
                  x = x / (dble(i) - shift)
                  call ga_put(g_x, i, i, ivec, ivec, x, 1)
               end do
            end do
         else if (type.eq.MT_DCPL) then
            do ivec = 1, nvec
               do i = 1, n
                  call ga_get(g_x, i, i, ivec, ivec, zx, 1)
                  zx = zx / (dble(i) - shift)
                  call ga_put(g_x, i, i, ivec, ivec, zx, 1)
               end do
            end do
         else
            call errquit('test_precond: illegal type',type,UERR)
         endif
      end if
      call ga_sync()
c
      end
C>
C> \brief test driver for ga_lkain
C>
C> This routine drives a standard test of the ga_lkain 
C> and the ga_kain solver.
C> Because the test is fixed the subroutine takes no inputs
C> and it prints a summary on standard output. 
C>
      subroutine ga_lkain_test()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer n, nvec, maxsub, maxiter
c     parameter (n=300, nvec=1)
      parameter (n=10, nvec=1)
      double complex angle
      integer g_x
      double precision util_random
      external test_precond, test_product, test_nlproduct, util_random
c
      integer g_tx, g_ta, g_tb
      integer g_a, g_b
      common /testme/ g_a, g_b
c
      integer i
      logical  ga_copy_dz
      external ga_copy_dz
****      integer info
****      double precision a(n,n), b(n), w(n)
c
      maxsub = 4*nvec
      maxiter = 100
c
      if (.not. ga_create(MT_DBL, n, nvec, 'testx', 0, 0, g_x))
     $     call errquit('test kain', 1, GA_ERR)
      if (.not. ga_create(MT_DBL, n, nvec, 'testb', 0, 0, g_b))
     $     call errquit('test kain', 2, GA_ERR)
      if (.not. ga_create(MT_DBL, n, n, 'testA', 0, 0, g_A))
     $     call errquit('test kain', 3, GA_ERR)
c
      call ga_ran_fill(g_A, 1, n, 1, n)
      call ga_ran_fill(g_b, 1, n, 1, nvec)
      if (ga_nodeid() .eq. 0) then
         do i = 1, n
            call ga_put(g_a, i, i, i, i, 0.5*dble(i), 1)
         end do
      end if
      call ga_sync()
c
      call ga_copy(g_b, g_x)
      call test_precond(g_x,0.0d0)
      call util_flush(6)
      call ga_print(g_b)
      call ga_print(g_x)
      call util_flush(6)
c
****      call ga_get(g_a, 1, n, 1, n, a, n)
****      call ga_get(g_b, 1, n, 1, nvec, b, n)
****      call dgesv(n, nvec, a, n, w, b, n, info)
****      write(6,*) ' info ', info
****      call ga_put(g_x, 1, n, 1, nvec, b, n)
c
C     This should have something other than zero
      call ga_lkain(0,g_x, g_b,test_product,test_precond,1d-6,maxsub,
     $     maxiter,.true.,.true.)
      call util_flush(6)
      call ga_print(g_x)
      call util_flush(6)
      call ga_summarize(0)
      call ma_summarize_allocated_blocks()
c
      write(6,*)
      write(6,*)
      write(6,*)
      write(6,*) ' DOING NON LINEAR '
      write(6,*)
      write(6,*)
      write(6,*)
      call ga_copy(g_b, g_x)
      call test_precond(g_x,0.0d0)
      maxsub = 10
c
      call ga_kain(g_x, 
     $     test_nlproduct, test_precond, 
     $     1d-6, 
     $     10.0d0, 10.0d0,
     $     maxsub, maxiter, 
     $     .true.)
c
      call util_flush(6)
      call ga_print(g_x)
      call util_flush(6)
      call ga_summarize(0)
      call ma_summarize_allocated_blocks()
c
c     Same but now COMPLEX
c
      g_tx = g_x
      g_tb = g_b
      g_tA = g_A
      maxsub = 4*nvec
      maxiter = 100
      if (.not. ga_create(MT_DCPL, n, nvec, 'testx', 0, 0, g_x))
     $     call errquit('test kain', 1, GA_ERR)
      if (.not. ga_create(MT_DCPL, n, nvec, 'testb', 0, 0, g_b))
     $     call errquit('test kain', 2, GA_ERR)
      if (.not. ga_create(MT_DCPL, n, n, 'testA', 0, 0, g_A))
     $     call errquit('test kain', 3, GA_ERR)
c
c     call ga_ran_fill(g_A, 1, n, 1, n)
c     call ga_ran_fill(g_b, 1, n, 1, nvec)
c     if (ga_nodeid() .eq. 0) then
c        do i = 1, n
c           call ga_put(g_a, i, i, i, i, 0.5*cmplx(dble(i)), 1)
c        end do
c     end if
      if (.not. ga_copy_dz(g_tb,g_b)) 
     +  call errquit('test kain: copy g_b failed',0,GA_ERR)
      if (.not. ga_copy_dz(g_tA,g_A))
     +  call errquit('test kain: copy g_A failed',0,GA_ERR)
      call ga_sync()
      angle = cmplx(1.0d0,util_random(0))
      angle = angle/abs(angle)
      call ga_scale(g_b,angle)
      angle = cmplx(1.0d0,util_random(0))
      angle = angle/abs(angle)
      call ga_scale(g_A,angle)
c
      if (.not. ga_destroy(g_tx))
     $     call errquit('test kain', 1, GA_ERR)
      if (.not. ga_destroy(g_tb))
     $     call errquit('test kain', 2, GA_ERR)
      if (.not. ga_destroy(g_tA))
     $     call errquit('test kain', 3, GA_ERR)
c
      call ga_copy(g_b, g_x)
      call test_precond(g_x,0.0d0)
      call util_flush(6)
      call ga_print(g_b)
      call ga_print(g_x)
      call util_flush(6)
c
****      call ga_get(g_a, 1, n, 1, n, a, n)
****      call ga_get(g_b, 1, n, 1, nvec, b, n)
****      call dgesv(n, nvec, a, n, w, b, n, info)
****      write(6,*) ' info ', info
****      call ga_put(g_x, 1, n, 1, nvec, b, n)
c
C     This should have something other than zero
      call ga_lkain(0,g_x, g_b,test_product,test_precond,1d-6,maxsub,
     $     maxiter,.true.,.true.)
      call util_flush(6)
      call ga_print(g_x)
      call util_flush(6)
      call ga_summarize(0)
      call ma_summarize_allocated_blocks()
c
      write(6,*)
      write(6,*)
      write(6,*)
      write(6,*) ' DOING NON LINEAR '
      write(6,*)
      write(6,*)
      write(6,*)
      call ga_copy(g_b, g_x)
      call test_precond(g_x,0.0d0)
      maxsub = 10
c
      call ga_kain(g_x, 
     $     test_nlproduct, test_precond, 
     $     1d-6, 
     $     10.0d0, 10.0d0,
     $     maxsub, maxiter, 
     $     .true.)
c
      call util_flush(6)
      call ga_print(g_x)
      call util_flush(6)
      call ga_summarize(0)
      call ma_summarize_allocated_blocks()
      call errquit('done',0, MEM_ERR)
c
      if (.not. ga_destroy(g_x))
     $     call errquit('test kain', 1, GA_ERR)
      if (.not. ga_destroy(g_b))
     $     call errquit('test kain', 2, GA_ERR)
      if (.not. ga_destroy(g_A))
     $     call errquit('test kain', 3, GA_ERR)
c
      end
C>
C> \brief The Linear System Solver
C>
C> The routine solves a linear system of equations using a Krylov
C> subspace method:
C> \f{eqnarray*}{
C>   F(x) &=& b \\\\
C>   F(x) &=& Ax
C> \f}
C> The Right-Hand Side (RHS) matrix and the solution
C> matrix store the data pertaining to a single linear system in 
C> columns. Using multiple columns the solver can solve multiple
C> equations simultaneously. Internally one Krylov subspace is used
C> for all systems of equations, and subspace vectors contribute to
C> the solutions of all equations.
C>
C> An special argument `odiff` is provided that allows one to use a
C> difference strategy for the required matrix-vector products. I.e.
C> rather than computing 
C> \f{eqnarray*}{
C>     y &=& Ax_i
C> \f}
C> each iteration one uses a formulation that requires
C> \f{eqnarray*}{
C>     y' &=& A(x_i-x_{i-1})
C> \f}
C> As \f$ x \f$ converges \f$ x_i - x_{i-1} \f$ approaches \f$ 0 \f$.
C> This may be exploited in the matrix-vector product to achieve 
C> considerable savings.
C>
C> Furthermore the routine takes two subroutines as arguments. One 
C> provides the matrix-vector product and has the interface
C> ~~~~
C> subroutine product(acc,g_x,g_Ax)
C> ~~~~
C> which evaluates
C> \f{eqnarray*}{
C>   g_{Ax} &=& A g_x
C> \f}
C> where
C>
C> - acc: the accuracy required for the matrix-vector elements
C>
C> - g_x: is a global array containing the input vectors in columns
C> 
C> - g_Ax: is a global array that will contain the matrix-vector product
C>   results in columns
C>
C> Note that `ga_lkain` may pass global arrays with a different number
C> of vectors than the number of linear systems being solved. Hence
C> the product routine should check the actual dimension of g_x and
C> g_Ax.
C>
C> The second subroutine provides a preconditioner and has the 
C> interface
C> ~~~~
C> subroutine precond(g_x,shift)
C> ~~~~
C> where 
C> 
C> - g_x: is a global array holding a set of vectors that will have
C>   the preconditioner applied to them 
C>
C> - shift: is a shift to ensure the preconditioner is non-singular
C>
C> Note that the number of vectors may also here be different from the
C> number of linear systems being solved. Hence the dimension of g_x
C> should be tested. Finally, the shift is not used in this linear
C> system solver. The argument is provided to ensure that the same
C> preconditioner can be used in eigenvalue solvers where the shift is
C> used.
C> 
      subroutine ga_lkain(rtdb,g_x, g_b, product, precond, 
     $     tol, mmaxsub, maxiter, odiff, oprint)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer g_x          !< [Input/Output] Initial guess/solution
      integer g_b          !< [Input] Right-hand side vectors
      external product     !< [Input] product routine
      external precond     !< [Input] preconditioner routine
      double precision tol !< [Input] convergence threshold
      integer mmaxsub      !< [Input] maximum subspace dimension
      integer maxiter      !< [Input] maximum no. of iterations
      logical odiff        !< [Input] use differences in product
      logical oprint       !< [Input] print flag
      integer rtdb         !< [Input] the RTDB handle
c
c     Solves the linear equations A(X)-B=0 for multiple vectors.
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . optional shift (not used here but used by the diagonalizer)
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c     at least while the product() routine is being executed.
c
      integer dimi, dimj
      integer iter, n, nvec, nsub, isub, typex, typeb, maxsub
      integer g_y, g_Ay, g_Ax, g_r, g_a, g_bb, g_c, g_xold, g_Axold
      double precision rmax, acc, ga_svd_tol
      double complex zero, one, mone
      logical converged
      logical odebug
      character*255 filestub,filesoln
c
      logical  file_write_ga, file_read_ga
      external file_write_ga, file_read_ga
c
      logical solver_restart
      external solver_restart
c 
      logical do_restart
c
      zero = cmplx(0.0d0,0.0d0)
      one  = cmplx(1.0d0,0.0d0)
      mone = cmplx(-1.0d0,0.0d0)
c
      odebug = util_print('debug lsolve', print_never) .and. 
     $     ga_nodeid().eq.0
      if (.not.rtdb_get(rtdb,'cphf:acc',mt_dbl,1,acc)) acc=1d-4*tol
      call ga_inquire(g_b, typeb, dimi, dimj)
      call ga_inquire(g_x, typex, n, nvec)
      if (typeb.ne.typex) 
     $  call errquit("ga_lkain: solution and right-hand-side vectors not
     $ of same type",0,UERR)
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec
c
      if (.not.rtdb_get(rtdb,'cphf:ga_svd_tol',mt_dbl,1,
     &                  ga_svd_tol)) then
c       See comment just before the ga_svd_solve_seq call to
c       understand these choices.
        if ((100*maxsub).lt.n) then
          ga_svd_tol = 1d-7
        else
          ga_svd_tol = 1d-14
        endif 
      endif
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) n, nvec, maxsub, maxiter, tol, util_wallsec()
 1       format(//,'Iterative solution of linear equations',/,
     $        '  No. of variables', i9,/,
     $        '  No. of equations', i9,/,
     $        '  Maximum subspace', i9,/,
     $        '        Iterations', i9,/,
     $        '       Convergence', 1p,d9.1,/,
     $        '        Start time', 0p,f9.1,/)
         call util_flush(6)
      end if
c
      if (.not. ga_create(typex, n, maxsub, 'lkain: Y', 
     $     0, 0, g_y))
     $     call errquit('lkain: failed allocating subspace', maxsub,
     &       GA_ERR)
      if (.not. ga_create(typex, n, maxsub, 'lkain: Ay', 
     $     0, 0, g_Ay))
     $     call errquit('lkain: failed allocating subspace2', maxsub,
     &       GA_ERR)
      if (.not. ga_create(typex, n, nvec, 'lkain: Ax',
     $     0, 0, g_Ax))
     $     call errquit('lkain: failed allocating subspace3', nvec,
     &       GA_ERR)
      if (.not. ga_create(typex, n, nvec, 'lkain: r',
     $     0, 0, g_r))
     $     call errquit('lkain: failed allocating subspace4', nvec,
     &       GA_ERR)
      if (odiff) then
         if (.not. ga_create(typex, n, nvec, 'lkain: xold',
     $        0, 0, g_xold))
     $        call errquit('lkain: failed allocating subspace5', nvec,
     &       GA_ERR)
         if (.not. ga_create(typex, n, nvec, 'lkain: Axold',
     $        0, 0, g_Axold))
     $        call errquit('lkain: failed allocating subspace6', nvec,
     &       GA_ERR)
         call ga_zero(g_xold)
         call ga_zero(g_Axold)
      end if
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_Ax)
      call ga_zero(g_r)
      call ga_sync()
c
c     Solution file
c
      if (.not. rtdb_cget(rtdb, 'solver:filestub', 1, filestub))
     &       filestub = 'lkain_soln'
      if (.not. rtdb_cget(rtdb, 'solver:filesoln', 1, filesoln))
     &       filesoln = 'lkain_soln'
#if 0
      call util_file_name(filestub,.false.,.false.,filesoln)
#else
      call cphf_fname(filestub,filesoln)
#endif
c     if (ga_nodeid().eq.0) write(luout,*) "ga_lkain filestub:",filestub
c     if (ga_nodeid().eq.0) write(luout,*) "ga_lkain filesoln:",filesoln
c
c     Check if this is a restart
c
      if (solver_restart(rtdb)) then
        do_restart = .true.
        if(.not.file_read_ga(filesoln,g_x)) do_restart = .false.
        if (do_restart) then
          if (ga_nodeid().eq.0) 
     &     write(luout,*) "Restarting solution from: ", filesoln
        else
          if (ga_nodeid().eq.0)
     &     write(luout,*) "Error in restart solution: ", filesoln
          goto 50
        end if  ! do_restart
      end if  ! solver_restart
c
 50   continue
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,2)
         call util_flush(6)
 2       format(/
     $        '   iter   nsub   residual    time',/,
     $        '   ----  ------  --------  ---------')
      end if
      nsub = 0
      converged = .false.
      do iter = 1, maxiter
         if (odiff) then
            call ga_add(one, g_x, mone, g_xold,  g_x)
         end if
         call product(acc,g_x, g_Ax)
         if (odiff) then
            call ga_add(one, g_Ax, one, g_Axold, g_Ax)
            call ga_add(one, g_x,  one, g_xold,  g_x)
            call ga_copy(g_x, g_xold)
            call ga_copy(g_Ax, g_Axold)
         end if
         call ga_zero(g_r)
         call ga_sync()
         call ga_add(one, g_b, mone, g_Ax, g_r) ! The residual
         call ga_sync()
         call ga_maxelt(g_r, rmax)
c
         if (.not. rtdb_put(rtdb, 'lkain:rmax', mt_dbl, 1, rmax))
     $     call errquit('ga_lkain: rmax put failed', 1, RTDB_ERR)
c
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,3) iter, nsub+nvec, rmax, util_wallsec()
            call util_flush(6)
 3          format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1)
         end if
c 
c        Check convergence
c
         if (rmax .lt. tol) then
            converged = .true.
            goto 100
         end if
         call precond(g_Ax,0.0d0)
         call precond(g_r,0.0d0)
         call ga_sync()
c
c        Copy the vectors to the subspace work area
c
         call ga_copy_patch('n', 
     $        g_Ax, 1, n, 1, nvec, 
     $        g_Ay, 1, n, nsub+1, nsub+nvec)
         call ga_copy_patch('n', 
     $        g_x, 1, n, 1, nvec, 
     $        g_y, 1, n, nsub+1, nsub+nvec)
         nsub = nsub + nvec
c
c        Form and solve the subspace equations using SVD in order
c        to manage near linear dependence in the subspace.
c     
         if (.not. ga_create(typex, nsub, nsub, 'lkain: A', 0, 0, g_a))
     $        call errquit('lkain: allocating g_a?', nsub, GA_ERR)
         if (.not. ga_create(typex, nsub, nvec, 'lkain: B', 0, 0,g_bb))
     $        call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
         if (.not. ga_create(typex, nsub, nvec, 'lkain: C', 0, 0, g_c))
     $        call errquit('lkain: allocating g_c?', nsub, GA_ERR)
         call ga_zero(g_a)
         call ga_zero(g_bb)
         call ga_zero(g_c)
         call ga_sync()
         call ga_dgemm('c','n',nsub,nsub,n,one,g_y,g_Ay,zero,g_a)
         call ga_sync()
         call ga_dgemm('c','n',nsub,nvec,n,one,g_y,g_r,zero,g_bb)
         call ga_sync()
         if (odebug) then
           call util_flush(6)
           call ga_print(g_a)
           call ga_print(g_bb)
           call util_flush(6)
         endif
c
c     The threshold used here should reflect the accuracy in the
c     products.  If very accurate products are used, then there is big
c     advantage for small cases (maxsub close to n) in using a very
c     small threshold in the SVD solve (e.g., 1e-14), but for more
c     realistic examples (maxsub << n) there is only a little
c     advantage and in the precence of real noise in the products
c     screening with a realistic threshold is important.
c
         call ga_svd_solve_seq(g_a,g_bb,g_c,ga_svd_tol)
         if (odebug) then
           call util_flush(6)
           call ga_print(g_c)
           call util_flush(6)
         endif
c
c     Form and add the correction, in parts, onto the solution
c
         call ga_sync()
         call ga_dgemm('n','n',n,nvec,nsub,mone,g_Ay,g_c,one,g_r)
         if (odebug) then
            call util_flush(6)
            write(6,*) ' The update in the complement '
            call util_flush(6)
            call ga_print(g_r)
         end if
         call ga_sync()
         call ga_add(one, g_r, one, g_x, g_x)
         call ga_sync()
         call ga_dgemm('n','n',n,nvec,nsub,one,g_y,g_c,zero,g_r)
         if (odebug) then
            call util_flush(6)
            write(6,*) ' The update in the subspace '
            call util_flush(6)
            call ga_print(g_r)
         end if
         call ga_sync()
         call ga_add(one, g_r, one, g_x, g_x)
         call ga_sync()
c
c        Save intermediate solution
c
         if(.not.file_write_ga(filesoln,g_x)) call errquit
     $     ('ga_lkain:could not write solution',1, DISK_ERR)
c
         if (.not. ga_destroy(g_a)) call errquit('lkain: a',0, GA_ERR)
         if (.not. ga_destroy(g_bb))call errquit('lkain: b',0, GA_ERR)
         if (.not. ga_destroy(g_c)) call errquit('lkain: c',0, GA_ERR)
c
c     Reduce the subspace as necessary
c
         if (nsub .eq. maxsub) then
            do isub = nvec+1, maxsub, nvec
               call ga_copy_patch('n', 
     $              g_Ay, 1, n, isub, isub+nvec-1, 
     $              g_Ax, 1, n, 1, nvec)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, nvec,
     $              g_Ay, 1, n, isub-nvec, isub-1)
c
               call ga_copy_patch('n', 
     $              g_y, 1, n, isub, isub+nvec-1, 
     $              g_Ax, 1, n, 1, nvec)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, nvec,
     $              g_y, 1, n, isub-nvec, isub-1)
            end do
            nsub = nsub - nvec
            call ga_sync()
         end if
c
      end do
 100  continue
c
      if (odiff) then
         if (.not. ga_destroy(g_xold)) call errquit('lkain: destroy',1,
     &       GA_ERR)
         if (.not. ga_destroy(g_Axold)) call errquit('lkain: destroy',2,
     &       GA_ERR)
      end if
      if (.not. ga_destroy(g_Ax)) call errquit('lkain: destroy',20,
     &       GA_ERR)
      if (.not. ga_destroy(g_Ay)) call errquit('lkain: destroy',3,
     &       GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit('lkain: destroy',4,
     &       GA_ERR)
      if (.not. ga_destroy(g_r)) call errquit('lkain: destroy',5,
     &       GA_ERR)
c
      if (.not. converged) call errquit('lkain: not converged',0,
     &       CALC_ERR)
c
      end
C>
C> \brief The Non-Linear System Solver
C>
C> The routine solves a non-linear system of equations using a Krylov
C> subspace method.
C> \f{eqnarray*}{
C>   F(x) &=& 0 
C> \f}
C> For linear systems the vector function in the above equation would
C> be \f$ F(x) = Ax - b \f$ but this routine also allows for non-linear
C> vector functions.
C> The Right-Hand Side (RHS) matrix and the solution
C> matrix store the data pertaining to a single linear system in 
C> columns. Using multiple columns the solver can solve multiple
C> equations simultaneously. Internally one Krylov subspace is used
C> for all systems of equations, and subspace vectors contribute to
C> the solutions of all equations.
C>
C> Furthermore the routine takes two subroutines as arguments. One 
C> provides the vector residual and has the interface
C> ~~~~
C> subroutine residual(acc,g_x,g_Ax)
C> ~~~~
C> where
C>
C> - acc: the accuracy required for the matrix-vector elements
C>
C> - g_x: is a global array containing the input vectors in columns
C> 
C> - g_Ax: is a global array that will contain the residual-vector
C>   results in columns
C>
C> Note that `ga_lkain` may pass global arrays with a different number
C> of vectors than the number of linear systems being solved. Hence
C> the product routine should check the actual dimension of g_x and
C> g_Ax.
C>
C> The second subroutine provides a preconditioner and has the 
C> interface
C> ~~~~
C> subroutine precond(g_x,shift)
C> ~~~~
C> where 
C> 
C> - g_x: is a global array holding a set of vectors that will have
C>   the preconditioner applied to them 
C>
C> - shift: is a shift to ensure the preconditioner is non-singular
C>
C> Note that the number of vectors may also here be different from the
C> number of linear systems being solved. Hence the dimension of g_x
C> should be tested. Finally, the shift is not used in this linear
C> system solver. The argument is provided to ensure that the same
C> preconditioner can be used in eigenvalue solvers where the shift is
C> used.
C> 
      subroutine ga_kain(
     $     g_x, 
     $     residual, precond, 
     $     tol, 
     $     trustmin, trustmax, 
     $     maxsub, maxiter, 
     $     oprint)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer g_x               !< [Input/Output] Initial guess/solution
      external residual         !< [Input] residual routine
      external precond          !< [Input] preconditioner routine
      double precision tol      !< [Input] convergence threshold
      double precision trustmin !< [Input] range to constrain trust radius
      double precision trustmax !< [Input] range to constrain trust radius
      integer maxsub            !< [Input] maximum subspace dimension
      integer maxiter           !< [Input] maximum no. of iterations
      logical oprint            !< [Input] print flag
c
c     Solves the non-linear equations f(X)=0 for multiple vectors.
c
c     call residual(acc, g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the residual
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the residual vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_kain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . optional shift (not used here but used by the diagonalizer)
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
      integer maxmaxsub
      parameter (maxmaxsub = 20)
      integer iter, n, nvec, nsub, isub, jsub, typex
      integer g_y, g_Ay, g_Ax, g_delta, g_a, g_b, g_c
      integer mdtob_nsub1
      double precision rmax, acc, trust
      double precision a(maxmaxsub,maxmaxsub), b(maxmaxsub), 
     $                 c(maxmaxsub), csum
      double complex za(maxmaxsub,maxmaxsub), zb(maxmaxsub), 
     $               zc(maxmaxsub), zcsum
      double complex zero, one, mone
      logical converged
      logical odebug
c
      zero = cmplx(0.0d0,0.0d0)
      one  = cmplx(1.0d0,0.0d0)
      mone = cmplx(-1.0d0,0.0d0)
c
      trust = trustmin
cold      acc = 0.01d0*tol
      acc = 0.0001d0*tol
      if (maxsub .gt. maxmaxsub) maxsub = maxmaxsub
      odebug = util_print('debug lsolve', print_never) .and. 
     $     ga_nodeid().eq.0
c
      call ga_inquire(g_x, typex, n, nvec)
      if (nvec .ne. 1) call errquit('kain: nvec?', nvec, GA_ERR)
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) n, maxsub, tol, trustmin, trustmax, util_wallsec()
 1       format(//,'Iterative solution of non-linear equations',/,
     $        '  No. of variables', i9,/,
     $        '  Maximum subspace', i9,/,
     $        '       Convergence', 1p,d9.1,/,
     $        '     Trust min/max', 0p,2f6.2,
     $        '        Start time', 0p,f9.1,/)
         call util_flush(6)
      end if
c
      if (.not. ga_create(typex, n, maxsub, 'kain: Y', 
     $     0, 0, g_y))
     $     call errquit('kain: failed allocating subspace', maxsub,
     &       GA_ERR)
      if (.not. ga_create(typex, n, maxsub, 'kain: Ay', 
     $     0, 0, g_Ay))
     $     call errquit('kain: failed allocating subspace2', maxsub,
     &       GA_ERR)
      if (.not. ga_create(typex, n, 1, 'kain: Ax',
     $     0, 0, g_Ax))
     $     call errquit('kain: failed allocating subspace3', 1, GA_ERR)
      call ga_zero(g_y)
      call ga_zero(g_Ay)
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,2)
 2       format(/
     $        '   iter   nsub   residual    time',/,
     $        '   ----  ------  --------  ---------')
      end if
      nsub = 0
      converged = .false.
      do iter = 1, maxiter
         call ga_zero(g_Ax)
         call residual(acc, g_x, g_Ax)
         call ga_maxelt(g_Ax, rmax)
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,3) iter, nsub+1, rmax, util_wallsec()
 3          format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1)
         end if
         if (rmax .lt. tol) then
            converged = .true.
            goto 100
         end if
c
c     Copy the vectors to the subspace work area
c
         call precond(g_Ax,0.0d0)
         call ga_copy_patch('n', 
     $        g_Ax, 1, n, 1, 1, 
     $        g_Ay, 1, n, nsub+1, nsub+1)
         call ga_copy_patch('n', 
     $        g_x, 1, n, 1, 1, 
     $        g_y, 1, n, nsub+1, nsub+1)
         nsub = nsub + 1
c
c     Not converged ... make the update
c
         g_delta = g_Ax         ! A reminder that these two are aliased
         call ga_scale(g_delta, mone)
c
         if (iter .gt. 1) then
c
c     Form the reduced space matrix and RHS
c
           if (typex.eq.MT_DBL) then
             call ga_local_mdot(n,nsub,nsub,a,maxmaxsub,g_y,g_Ay)
             do isub = 1, nsub-1
                b(isub) = -(a(isub,nsub) - a(nsub,nsub))
             end do
             do isub = 1, nsub-1
                do jsub = 1, nsub-1
                   a(isub,jsub) = a(isub,jsub)
     $                  - a(nsub,jsub) - a(isub,nsub) + a(nsub,nsub)
                end do
             end do
           else if (typex.eq.MT_DCPL) then
             call ga_local_zmdot(n,nsub,nsub,za,maxmaxsub,g_y,g_Ay)
             do isub = 1, nsub-1
                zb(isub) = -(za(isub,nsub) - za(nsub,nsub))
             end do
             do isub = 1, nsub-1
                do jsub = 1, nsub-1
                   za(isub,jsub) = za(isub,jsub)
     $                  - za(nsub,jsub) - za(isub,nsub) + za(nsub,nsub)
                end do
             end do
           else
             call errquit("ga_kain: illegal type",typex,UERR)
           endif
c
c     Solve the subspace equations (lazily using existing GA routine)
c     
            if (.not. ga_create(typex,nsub-1,nsub-1,'kain: A',
     $           nsub-1,nsub-1,g_a))
     $           call errquit('kain: allocating g_a?', nsub, GA_ERR)
            if (.not. ga_create(typex,nsub-1,1,'kain: B',nsub-1,1,g_b))
     $           call errquit('kain: allocating g_bb?', nsub, GA_ERR)
            if (.not. ga_create(typex,nsub-1,1,'kain: C',nsub-1,1,g_c))
     $           call errquit('kain: allocating g_c?', nsub, GA_ERR)
            if (ga_nodeid() .eq. 0) then
              if (typex.eq.MT_DBL) then
                call ga_put(g_a, 1, nsub-1, 1, nsub-1, a, maxmaxsub)
                call ga_put(g_b, 1, nsub-1, 1, 1, b, 1)
              else if (typex.eq.MT_DCPL) then
                call ga_put(g_a, 1, nsub-1, 1, nsub-1, za, maxmaxsub)
                call ga_put(g_b, 1, nsub-1, 1, 1, zb, 1)
              else
                call errquit("ga_kain: illegal type",typex,UERR)
              endif
            end if
            call ga_sync
c
            call ga_svd_solve_seq(g_a,g_b,g_c,1d-14)
c
            if (odebug) then
              call util_flush(6)
              call ga_print(g_c)
              call util_flush(6)
           endif
           mdtob_nsub1=MA_sizeof(MT_DBL,nsub-1,MT_BYTE)
            if (typex.eq.MT_DBL) then
              if (ga_nodeid() .eq. 0)
     $           call ga_get(g_c, 1, nsub-1, 1, 1, c, 1)
              call ga_brdcst(1, c, mdtob_nsub1, 0)
              write(6,*) ' KAIN SUBSPACE COEFFS'
              call output(c, 1, nsub-1, 1, 1, nsub-1, 1, 1)
            else if (typex.eq.MT_DCPL) then
              if (ga_nodeid() .eq. 0)
     $           call ga_get(g_c, 1, nsub-1, 1, 1, zc, 1)
              call ga_brdcst(1, zc, 2*mdtob_nsub1, 0)
              write(6,*) ' KAIN SUBSPACE COEFFS'
              call zoutput(zc, 1, nsub-1, 1, 1, nsub-1, 1, 1)
            endif
            call ga_sync
            if (.not. ga_destroy(g_a)) call errquit('kain: a',0, GA_ERR)
            if (.not. ga_destroy(g_b)) call errquit('kain: b',0, GA_ERR)
            if (.not. ga_destroy(g_c)) call errquit('kain: c',0, GA_ERR)
c     
c     Form the correction
c     
            if (typex.eq.MT_DBL) then
              csum = 0.0d0
              do isub = 1, nsub-1
                 csum = csum + c(isub)
                 call ga_add_patch( c(isub),  g_y, 1, n, isub, isub, 
     $                one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
                 call ga_add_patch(-c(isub), g_Ay, 1, n, isub, isub, 
     $                one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
              end do
              call ga_add_patch(-csum,  g_y, 1, n, nsub, nsub, 
     $             one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
              call ga_add_patch( csum,  g_Ay, 1, n, nsub, nsub, 
     $             one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
            else if (typex.eq.MT_DCPL) then
              zcsum = 0.0d0
              do isub = 1, nsub-1
                 zcsum = zcsum + zc(isub)
                 call ga_add_patch( zc(isub),  g_y, 1, n, isub, isub, 
     $                one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
                 call ga_add_patch(-zc(isub), g_Ay, 1, n, isub, isub, 
     $                one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
              end do
              call ga_add_patch(-zcsum,  g_y, 1, n, nsub, nsub, 
     $             one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
              call ga_add_patch( zcsum,  g_Ay, 1, n, nsub, nsub, 
     $             one, g_delta, 1, n, 1, 1, g_delta, 1, n, 1, 1)
            endif
         endif
c
c     Step restriction
c
         call ga_maxelt(g_delta, rmax)
         if (rmax .gt. trust) then
            if (oprint) write(6,*) ' RESTRICTION ', rmax, trust
            call ga_scale(g_delta, trust/rmax)
         end if
c
         call ga_add(one, g_delta, one, g_x, g_x)
c     
c     Reduce the subspace as necessary (note g_delta=g_Ax destroyed)
c     
         if (nsub .eq. maxsub) then
            do isub = 2, maxsub
               call ga_copy_patch('n', 
     $              g_Ay, 1, n, isub, isub, 
     $              g_Ax, 1, n, 1, 1)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, 1,
     $              g_Ay, 1, n, isub-1, isub-1)
c     
               call ga_copy_patch('n', 
     $              g_y, 1, n, isub, isub, 
     $              g_Ax, 1, n, 1, 1)
               call ga_copy_patch('n', 
     $              g_Ax, 1, n, 1, 1,
     $              g_y, 1, n, isub-1, isub-1)
            end do
            nsub = nsub - 1
         end if
c
      end do
 100  continue
c
      if (.not. ga_destroy(g_Ax)) call errquit('kain: destroy',20,
     &       GA_ERR)
      if (.not. ga_destroy(g_Ay)) call errquit('kain: destroy',3,
     &       GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit('kain: destroy',4, GA_ERR)
c
      if (.not. converged) call errquit('kain: not converged',0,
     &       CALC_ERR)
c
      end
C>
C> \brief A direct linear system solver
C>
C> Solve for X from the linear equations
C> \f{eqnarray*}{
C>    A*X &=& B
C> \f}
C> or more explicitly
C> \f{eqnarray*}{
C>    A(m,n)*X(n,nvec) = B(m,nvec)
C> \f}
C> Where \f$ A \f$ is a general real matrix (not necessarily square, or
C> symmetric, or full rank) and \f$ X \f$ and \f$ B \f$ are matrices with one or more
C> columns representing the solutions and right hand sides.  Singular
C> values of \f$ A \f$ less than \f$ tol\f$  are neglected. \f$ X \f$ is returned.
C>
C> If the SVD of \f$ A \f$ is \f$ U*values*VT \f$, then the solution
C> is of the form
C> \f{eqnarray*}{
C>    V*(1/values)*UT*B
C> \f}
C> where the reciprocal of `values` less than `tol` are neglected.
C     
      subroutine ga_svd_solve_seq(g_a, g_b, g_x, tol)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer g_a !< [Input] the matrix stored explicitly
      integer g_b !< [Input] the right-hand-sides
      integer g_x !< [Input/Output] the guess/solution
      double precision tol !< [Input] the tolerance
c
c     Solve for X from the linear equations
c
c     A*X = B
c
c     A(m,n)*X(n,nvec) = B(m,nvec)
c
c     Where A is a general real matrix (not necessarily square, or
c     symmetric, or full rank) and X and B are matrices with one or more
c     columns representing the solutions and right hand sides.  Singular
c     values of A less than tol are neglected.  X is returned.
c
c     If the SVD of A is U*values*VT, then the solution
c     is of the form
c
c     V*(1/values)*UT*B
c
c     where the reciprocal of values less than tol are neglected.
c     
      integer m,n,nn,type,nvec,nsing,l_val, k_val,g_u,g_vt,i,g_tmp
      double complex one, zero
      logical oprint
c
      oprint = util_print('debug svdsolve', print_high) .and.
     $     ga_nodeid().eq.0
c
      one  = cmplx(1.0d0,0.0d0)
      zero = cmplx(0.0d0,0.0d0)
c
      call ga_inquire(g_a, type, m, n)
      call ga_inquire(g_b, type, nn, nvec)
      if (nn .ne. n) call errquit('gasvdsol: b does not conform',nn,
     &       GA_ERR)
      nsing = min(m,n)
      if (.not. ma_push_get(MT_DBL, nsing, 'gasvdsol', l_val, k_val))
     $     call errquit('gasvdsol: val',nsing, MA_ERR)
      if (.not. ga_create(type,m,nsing,'gasvd',0,0,g_u))
     $     call errquit('gasvdsol: u',m*nsing, GA_ERR)
      if (.not. ga_create(type,nsing,n,'gasvd',0,0,g_vt))
     $     call errquit('gasvdsol: u',nsing*n, GA_ERR)
      if (.not. ga_create(type,nsing,nvec,'gasvd',0,0,g_tmp))
     $     call errquit('gasvdsol: tmp',nsing*nvec, GA_ERR)
      call ga_zero(g_tmp)
c
      call ga_svd_seq(g_a, g_u,g_vt,dbl_mb(k_val))
c
      do i = 0, nsing-1
        if (dbl_mb(k_val+i) .lt. tol) then
          if (ga_nodeid() .eq. 0 .and. oprint) then
            write(6,*) ' neglecting ', i+1, dbl_mb(k_val+i)
          endif
          dbl_mb(k_val+i) = 0.0d0
        else
          dbl_mb(k_val+i) = 1.0d0/dbl_mb(k_val+i)
        end if
      end do
c
      call ga_dgemm('c','n',nsing,nvec,m,one,g_u,g_b,zero,g_tmp)
      call ga_scale_lh(g_tmp,dbl_mb(k_val))
      call ga_zero(g_x)
      call ga_dgemm('c','n',n,nvec,nsing,one,g_vt,g_tmp,zero,g_x)
c
      if (.not. ga_destroy(g_tmp)) call errquit('gasvdsol: des',1,
     &       GA_ERR)
      if (.not. ga_destroy(g_u)) call errquit('gasvdsol: des',2,
     &       GA_ERR)
      if (.not. ga_destroy(g_vt)) call errquit('gasvdsol: des',3,
     &       GA_ERR)
      if (.not. ma_pop_stack(l_val)) call errquit('gasvdsol: pop',4,
     &       GA_ERR)
c
      end
C>
C> \brief Perform SVD on rectangular matrix
C>
      subroutine ga_svd_seq(g_a, g_u, g_vt, values)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_a, g_u, g_vt
      double precision values(*)
c     
c     Perform SVD on rectangular matrix
c
c     nsing = min(n,m)
c     g_a(m,n)      --- input matrix
c     g_u(m,nsing)  --- left singular vectors (output)
c     g_vt(nsing,n) --- right singular vectors transposed (output)
c     values(nsing) --- singular values (output)
c
c     A = U*values*VT
c
c     A possible parallel algorithm is to diagonalize ATA to get
c     V and AAT to get U --- both have values**2 as eigenvalues.
c
      integer n, m, type, l_a, k_a, l_u, k_u, l_vt, k_vt, 
     $     l_work, k_work, lwork, info, nsing
      integer l_rwork, k_rwork
c     
      call ga_inquire(g_a, type, m, n)
      nsing = min(m,n)
      if (ga_nodeid() .eq. 0) then
         lwork = 10*max(m,n)
         if (.not. ma_push_get(type, m*n, 'gasvd', l_a, k_a))
     $        call errquit('gasvd: a',m*n, MA_ERR)
         if (.not. ma_push_get(type, m*nsing, 'gasvd', l_u, k_u))
     $        call errquit('gasvd: u',m*nsing, MA_ERR)
         if (.not. ma_push_get(type, nsing*n, 'gasvd', l_vt, k_vt))
     $        call errquit('gasvd: vt',nsing*n, MA_ERR)
         if (.not. ma_push_get(type, lwork, 'gasvd', l_work, k_work))
     $        call errquit('gasvd: work',lwork, MA_ERR)
         if (type.eq.MT_DCPL) then
            if (.not. ma_push_get(MT_DBL, lwork, 'gasvd',
     $                            l_rwork, k_rwork))
     $           call errquit('gasvd: work',lwork, MA_ERR)
         endif
c
         if (type.eq.MT_DBL) then
c
           call ga_get(g_a, 1, m, 1, n, dbl_mb(k_a), m)
c
           call dgesvd('s','s',m,n,dbl_mb(k_a),m,values,
     $          dbl_mb(k_u),m,dbl_mb(k_vt),nsing,
     $          dbl_mb(k_work),lwork,info)
           if (info .ne. 0) call errquit('gasvd:d: failed',info,MEM_ERR)
c
           call ga_put(g_u,  1, n,     1, nsing, dbl_mb(k_u),  n)
           call ga_put(g_vt, 1, nsing, 1, m,     dbl_mb(k_vt), n)
c
         else if (type.eq.MT_DCPL) then
c
           call ga_get(g_a, 1, m, 1, n, dcpl_mb(k_a), m)
c
           call zgesvd('s','s',m,n,dcpl_mb(k_a),m,values,
     $          dcpl_mb(k_u),m,dcpl_mb(k_vt),nsing,
     $          dcpl_mb(k_work),lwork,dbl_mb(k_rwork),info)
           if (info .ne. 0) call errquit('gasvd:z: failed',info,MEM_ERR)
c
           call ga_put(g_u,  1, n,     1, nsing, dcpl_mb(k_u),  n)
           call ga_put(g_vt, 1, nsing, 1, m,     dcpl_mb(k_vt), n)
c
         else
c
           call errquit('gasvd: illegal data type',type,UERR)
c
         endif
c
         if (.not. ma_chop_stack(l_a)) call errquit('gasvd ma',0,
     &       MA_ERR)
      end if
      call ga_sync()
      call ga_brdcst(1,values,n*8,0)
      call ga_sync()
c     
      end
c
      logical function solver_restart(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer rtdb
c
      integer restr
c
c     Check for the restart flag
c
      solver_restart = .false.
      if (.not.rtdb_get(rtdb,'solver:restart',mt_int,1,restr))
     &     restr= 0
      if (restr.gt.0) solver_restart = .true.
c
      return
      end  
c
c     some parameters for the cphf solver (ga_lkain())
c
      subroutine solver_setup(rtdb,restr)
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "util.fh"
#include "inp.fh"
c
      integer rtdb
      integer restr
c
      character*256 cphf_sol
c
c     set up the files for the cphf solution
c
      call cphf_fname('cphf_sol',cphf_sol)
      if (.not. rtdb_cput(rtdb, 'solver:filestub', 1, 'cphf_sol'))
     *  call errquit('solver_setup: file stub',555,RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'solver:filesoln', 1,
     C     cphf_sol(1:inp_strlen(cphf_sol))))
     *  call errquit('solver_setup: intermediate solution file',555,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'solver:restart', MT_int, 1, restr))
     *  call errquit('solver_setup: solver restart flag',555, RTDB_ERR)
c
      return
      end
c
C> @}
