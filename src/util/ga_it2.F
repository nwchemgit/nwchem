      subroutine product(g_x, g_Ax)
      implicit none
      integer g_x, g_Ax
c
      integer g_a
      common /testme/ g_A
      integer n, nvec, type
c
      call ga_inquire(g_x, type, n, nvec)
      call ga_zero(g_Ax)
      call ga_dgemm('n', 'n', n, nvec, n, 1.0d0, g_A, g_x, 0.0d0, g_Ax)
c
      end
      subroutine precond(g_x)
      implicit none
#include "global.fh"
      integer g_x
      integer n, nvec, type
      integer i, ivec
      double precision x
c
      call ga_inquire(g_x, type, n, nvec)
      if (ga_nodeid() .eq. 0) then
         do ivec = 1, nvec
            do i = 1, n
               call ga_get(g_x, i, i, ivec, ivec, x, 1)
               x = x / dble(i)
               call ga_put(g_x, i, i, ivec, ivec, x, 1)
            end do
         end do
      end if
      call ga_sync()
c
      end
      subroutine ga_lkain_test()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer n, nvec, maxsub, maxiter
      parameter (n=5000, nvec=100)
      integer g_x, g_b
      external precond, product
c
      integer g_a
      common /testme/ g_a
c
      integer i
****      integer info
****      double precision a(n,n), b(n), w(n)
c
      maxsub = 4*nvec
      maxiter = 100
c
      if (.not. ga_create(MT_DBL, n, nvec, 'testx', 0, 0, g_x))
     $     call errquit('test kain', 1)
      if (.not. ga_create(MT_DBL, n, nvec, 'testx', 0, 0, g_b))
     $     call errquit('test kain', 2)
      if (.not. ga_create(MT_DBL, n, n, 'testx', 0, 0, g_A))
     $     call errquit('test kain', 3)
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
      call precond(g_x)
c
****      call ga_get(g_a, 1, n, 1, n, a, n)
****      call ga_get(g_b, 1, n, 1, nvec, b, n)
****      call dgesv(n, nvec, a, n, w, b, n, info)
****      write(6,*) ' info ', info
****      call ga_put(g_x, 1, n, 1, nvec, b, n)
c
      call ga_lkain(g_x, g_b, product, precond, 1d-6, maxsub, 
     $     maxiter,.true.,.true.)
      call ga_summarize(0)
      call ma_summarize_allocated_blocks()
      call errquit('done',0)
c
      end
      subroutine ga_lkain(
     $     g_x, g_b,
     $     product, precond, 
     $     tol, mmaxsub, maxiter, odiff, oprint)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer g_x               ! [input/output] Initial guess/solution
      integer g_b               ! [input] Right-hand side vectors
      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      integer mmaxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical odiff             ! [input] use differences in product
      logical oprint            ! [input] print flag
c
c     Solves the linear equations A(X)=0 for multiple vectors.
c
c     call product(g_x, g_Ax)
c     . where g_x contains the vectors and g_Ax should be filled
c     . with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x)
c     . apply preconditioning directly to the vectors in g_x
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c
      integer iter, n, nvec, nsub, isub, type, maxsub
      integer g_y, g_Ay, g_Ax, g_r, g_a, g_bb, g_c, g_xold, g_Axold
      double precision rmax
      logical converged
      logical odebug
c
      odebug = util_print('debug lsolve', print_never) .and. 
     $     ga_nodeid().eq.0
c
      call ga_inquire(g_x, type, n, nvec)
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) n, nvec, maxsub, tol, util_wallsec()
 1       format(//,'Iterative solution of linear equations',/,
     $        '  No. of variables', i9,/,
     $        '  No. of equations', i9,/,
     $        '  Maximum subspace', i9,/,
     $        '       Convergence', 1p,d9.1,/,
     $        '        Start time', 0p,f9.1,/)
         call util_flush(6)
      end if
c
      if (.not. ga_create(MT_DBL, n, maxsub, 'lkain: Y', 
     $     0, 0, g_y))
     $     call errquit('lkain: failed allocating subspace', maxsub)
      if (.not. ga_create(MT_DBL, n, maxsub, 'lkain: Ay', 
     $     0, 0, g_Ay))
     $     call errquit('lkain: failed allocating subspace2', maxsub)
      if (.not. ga_create(MT_DBL, n, nvec, 'lkain: Ax',
     $     0, 0, g_Ax))
     $     call errquit('lkain: failed allocating subspace3', nvec)
      if (.not. ga_create(MT_DBL, n, nvec, 'lkain: r',
     $     0, 0, g_r))
     $     call errquit('lkain: failed allocating subspace4', nvec)
      if (odiff) then
         if (.not. ga_create(MT_DBL, n, nvec, 'lkain: xold',
     $        0, 0, g_xold))
     $        call errquit('lkain: failed allocating subspace5', nvec)
         if (.not. ga_create(MT_DBL, n, nvec, 'lkain: xold',
     $        0, 0, g_Axold))
     $        call errquit('lkain: failed allocating subspace6', nvec)
         call ga_zero(g_xold)
         call ga_zero(g_Axold)
      end if
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_Ax)
      call ga_zero(g_r)
c
      nsub = 0
      converged = .false.
      do iter = 1, maxiter
         if (odiff) then
            call ga_add(1.0d0, g_x, -1.0d0, g_xold,  g_x)
         end if
         call product(g_x, g_Ax)
         if (odiff) then
            call ga_add(1.0d0, g_Ax, 1.0d0, g_Axold, g_Ax)
            call ga_add(1.0d0, g_x,  1.0d0, g_xold,  g_x)
            call ga_copy(g_x, g_xold)
            call ga_copy(g_Ax, g_Axold)
         end if
         call ga_zero(g_r)
         call ga_add(1.0d0, g_b, -1.0d0, g_Ax, g_r) ! The residual
         call ga_maxelt(g_r, rmax)
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,2) iter, nsub+nvec, rmax, util_wallsec()
 2          format(/
     $           '   iter   nsub   residual    time',/,
     $           '   ----  ------  --------  ---------',/,
     $           ' ', i5, i7, 3x,1p,d9.2,0p,f10.1/)
         end if
         if (rmax .lt. tol) then
            converged = .true.
            goto 100
         end if
         call precond(g_Ax)
         call precond(g_r)
c
c     Copy the vectors to the subspace work area
c
         call ga_copy_patch('n', 
     $        g_Ax, 1, n, 1, nvec, 
     $        g_Ay, 1, n, nsub+1, nsub+nvec)
         call ga_copy_patch('n', 
     $        g_x, 1, n, 1, nvec, 
     $        g_y, 1, n, nsub+1, nsub+nvec)
         nsub = nsub + nvec
c
c     Form and solve the subspace equations using SVD in order
c     to manage near linear dependence in the subspace.
c     
         if (.not. ga_create(MT_DBL, nsub, nsub, 'lkain: A', 0, 0, g_a))
     $        call errquit('lkain: allocating g_a?', nsub)
         if (.not. ga_create(MT_DBL, nsub, nvec, 'lkain: B', 0, 0,g_bb))
     $        call errquit('lkain: allocating g_bb?', nsub)
         if (.not. ga_create(MT_DBL, nsub, nvec, 'lkain: C', 0, 0, g_c))
     $        call errquit('lkain: allocating g_c?', nsub)
         call ga_zero(g_a)
         call ga_zero(g_bb)
         call ga_zero(g_c)
         call ga_dgemm('t','n',nsub,nsub,n,1.0d0,g_y,g_Ay,0.0d0,g_a)
         call ga_dgemm('t','n',nsub,nvec,n,1.0d0,g_y,g_r,0.0d0,g_bb)
         if (odebug) call ga_print(g_a)
         if (odebug) call ga_print(g_c)
c
c     The threshold used here should reflect the accuracy in the
c     products.  If very accurate products are used, then there is big
c     advantage for small cases (maxsub close to n) in using a very
c     small threshold in the SVD solve (e.g., 1e-14), but for more
c     realistic examples (maxsub << n) there is only a little
c     advantage and in the precence of real noise in the products
c     screening with a realistic threshold is important.
c
         call ga_svd_solve_seq(g_a,g_bb,g_c,tol*1d-2)
         if (odebug) call ga_print(g_c)
c
c     Form and add the correction, in parts, onto the solution
c
         call ga_dgemm('n','n',n,nvec,nsub,-1.0d0,g_Ay,g_c,1.0d0,g_r)
         if (odebug) then
            write(6,*) ' The update in the complement '
            call ga_print(g_r)
         end if
         call ga_add(1.0d0, g_r, 1.0d0, g_x, g_x)
         call ga_dgemm('n','n',n,nvec,nsub,1.0d0,g_y,g_c,0.0d0,g_r)
         if (odebug) then
            write(6,*) ' The update in the subspace '
            call ga_print(g_r)
         end if
         call ga_add(1.0d0, g_r, 1.0d0, g_x, g_x)
c
         if (.not. ga_destroy(g_a)) call errquit('lkain: a',0)
         if (.not. ga_destroy(g_bb))call errquit('lkain: b',0)
         if (.not. ga_destroy(g_c)) call errquit('lkain: c',0)
c
c     Reduce the subspace as necessary
c
         if (nsub .eq. maxsub) then
            do isub = nsub+1, maxsub, nvec
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
         end if
c
      end do
 100  continue
c
      if (odiff) then
         if (.not. ga_destroy(g_xold)) call errquit('lkain: destroy',1)
         if (.not. ga_destroy(g_Axold)) call errquit('lkain: destroy',2)
      end if
      if (.not. ga_destroy(g_Ax)) call errquit('lkain: destroy',20)
      if (.not. ga_destroy(g_Ay)) call errquit('lkain: destroy',3)
      if (.not. ga_destroy(g_y)) call errquit('lkain: destroy',4)
      if (.not. ga_destroy(g_r)) call errquit('lkain: destroy',5)
c
      if (.not. converged) call errquit('lkain: not converged',0)
c
      end
      subroutine ga_lu_solve_seq(op, g_a, g_b)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      character*(*) op
      integer g_a, g_b
c     
c     Solve A*X=B using LU replacing B with X
c
      integer n, nn, nvec, type, l_a, k_a, k_b, l_b, k_piv, l_piv, info
c     
      if (op .ne. 'n') call errquit('luseq op',0)
      if (ga_nodeid() .eq. 0) then
         call ga_inquire(g_a, type, n, nn)
         if (n .ne. nn) call errquit('luseq: matrix not square',0)
         call ga_inquire(g_b, type, nn, nvec)
         if (n .ne. nn) call errquit('luseq: matrices do not conform',0)
c     
         if (.not. ma_push_get(MT_DBL, n*n, 'luseq', l_a, k_a))
     $        call errquit('luseq: a',n*n)
         if (.not. ma_push_get(MT_DBL, n*nvec, 'luseq', l_b, k_b))
     $        call errquit('luseq: b',n*nvec)
         if (.not. ma_push_get(mt_int, n, 'luseq', l_piv, k_piv))
     $        call errquit('luseq: piv', n)
         call ga_get(g_a, 1, n, 1, n, dbl_mb(k_a), n)
         call ga_get(g_b, 1, n, 1, nvec, dbl_mb(k_b), n)
         call dgesv(n, nvec, dbl_mb(k_a), n, int_mb(k_piv), 
     $        dbl_mb(k_b), n, info)
         if (info .ne. 0) call errquit('luseq: singular', info)
         call ga_put(g_b, 1, n, 1, nvec, dbl_mb(k_b), n)
         if (.not. ma_chop_stack(l_a)) call errquit('luseq ma',0)
      end if
      call ga_sync()
c     
      end
      subroutine ga_svd_solve_seq(g_a, g_b, g_x, tol)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_a, g_b, g_x
      double precision tol
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
c
      call ga_inquire(g_a, type, m, n)
      call ga_inquire(g_b, type, nn, nvec)
      if (nn .ne. n) call errquit('gasvdsol: b does not conform',nn)
      nsing = min(m,n)
      if (.not. ma_push_get(MT_DBL, nsing, 'gasvdsol', l_val, k_val))
     $     call errquit('gasvdsol: val',nsing)
      if (.not. ga_create(MT_DBL,m,nsing,'gasvd',0,0,g_u))
     $     call errquit('gasvdsol: u',m*nsing)
      if (.not. ga_create(MT_DBL,nsing,n,'gasvd',0,0,g_vt))
     $     call errquit('gasvdsol: u',nsing*n)
      if (.not. ga_create(MT_DBL,nsing,nvec,'gasvd',0,0,g_tmp))
     $     call errquit('gasvdsol: tmp',nsing*nvec)
      call ga_zero(g_tmp)
c
      call ga_svd_seq(g_a, g_u,g_vt,dbl_mb(k_val))
c
      do i = 0, nsing-1
         if (dbl_mb(k_val+i) .lt. tol) then
            write(6,*) ' neglecting ', i+1, dbl_mb(k_val+i)
            dbl_mb(k_val+i) = 0.0d0
         else
            dbl_mb(k_val+i) = 1.0d0/dbl_mb(k_val+i)
         end if
      end do
c
      call ga_dgemm('t','n',nsing,nvec,m,1.0d0,g_u,g_b,0.0d0,g_tmp)
      call ga_scale_lh(g_tmp,dbl_mb(k_val))
      call ga_zero(g_x)
      call ga_dgemm('t','n',n,nvec,nsing,1.0d0,g_vt,g_tmp,0.0d0,g_x)
c
      if (.not. ga_destroy(g_tmp)) call errquit('gasvdsol: des',1)
      if (.not. ga_destroy(g_u)) call errquit('gasvdsol: des',2)
      if (.not. ga_destroy(g_vt)) call errquit('gasvdsol: des',3)
      if (.not. ma_pop_stack(l_val)) call errquit('gasvdsol: pop',4)
c
      end
      subroutine ga_svd_seq(g_a, g_u, g_vt, values)
      implicit none
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
c     
      call ga_inquire(g_a, type, m, n)
      nsing = min(m,n)
      if (ga_nodeid() .eq. 0) then
         lwork = 10*max(m,n)
         if (.not. ma_push_get(MT_DBL, m*n, 'gasvd', l_a, k_a))
     $        call errquit('gasvd: a',m*n)
         if (.not. ma_push_get(MT_DBL, m*nsing, 'gasvd', l_u, k_u))
     $        call errquit('gasvd: u',m*nsing)
         if (.not. ma_push_get(MT_DBL, nsing*n, 'gasvd', l_vt, k_vt))
     $        call errquit('gasvd: vt',nsing*n)
         if (.not. ma_push_get(MT_DBL, lwork, 'gasvd', l_work, k_work))
     $        call errquit('gasvd: work',lwork)
c
         call ga_get(g_a, 1, m, 1, n, dbl_mb(k_a), m)
c
         call dgesvd('s','s',m,n,dbl_mb(k_a),m,values,
     $        dbl_mb(k_u),m,dbl_mb(k_vt),nsing,
     $        dbl_mb(k_work),lwork,info)
         if (info .ne. 0) call errquit('gasvd: failed', info)
c
         call ga_put(g_u,  1, n,     1, nsing, dbl_mb(k_u),  n)
         call ga_put(g_vt, 1, nsing, 1, m,     dbl_mb(k_vt), n)
c
         if (.not. ma_chop_stack(l_a)) call errquit('gasvd ma',0)
      end if
      call ga_sync()
      call ga_brdcst(1,values,n*8,0)
      call ga_sync()
c     
      end
