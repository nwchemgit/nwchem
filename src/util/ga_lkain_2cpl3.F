      subroutine ga_lkain_2cpl3(rtdb, g_x, g_b, g_x_im, g_b_im,
     &   product, precond, 
     $   tol, mmaxsub, maxiter, odiff, oprint, omega, limag,
     &   lifetime, gamwidth, ncomp)

c     $Id: ga_lkain_2cpl3.F,v 1.4 2007-10-31 18:58:18 bert Exp $

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
c
      integer rtdb              ! [input] database handle
      integer g_x(2)            ! [input/output] Initial guess/solution
      integer g_x_im(2)         ! not used
      integer g_b(2)            ! [input] Right-hand side vectors 
      integer g_b_im(2)         ! not used
      double precision omega    ! [input] coupling parameter
      logical limag             ! [input] imaginary perturbation?
      logical lifetime          ! [input] consider damping or not?
      double precision gamwidth ! [input] damping parameter
      integer ncomp             ! [input] no. of components to treat
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
c ... jochen:
c     This is a modified version of ga_lkain from file ga_it2.F
c     This version allows to solve a coupled set of equations, i.e.
c     there are two right-hand vectors and two initial guesses and two
c     solutions which are coupled. The coupling is mediated by a
c     parameter omega in the call to the preconditioner 
c     (elsewhere, omega is simply called "frequency")
c
c     note: when called from cphf_solve3, odiff = .false. on input
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . coupling parameter omega
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
      integer iter, n, n2, nvec, nsub, isub, type, maxsub, ipm,
     &   ntmp1, ntmp2

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      integer g_y, g_Ay, g_Ax(2), g_r(2), g_r2, g_a, g_bb,
     &   g_c, g_xold(2), g_Axold(2), g_Ax_im(2)
      double precision rmax, rmax1, rmax2, acc
      logical converged
      logical odebug, debug, converge_precond
c
c     =================================================================
c
      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

c     check input key if we should check for convergence
c     after the preconditioner has been applied to the residual
      if (.not. rtdb_get(rtdb, 'aoresponse:precond',    mt_log, 1,
     &                            converge_precond))
     &  converge_precond = .false.
      
      if (debug) write (6,*) 'ga_lkain_2cpl3 omega =',omega
      if (debug) write (6,*) 'ga_lkain_2cpl3 limag =',limag
      if (debug) write (6,*) 'ga_lkain_2cpl3 lifetime,gamwidth',
     &   lifetime,gamwidth
      if (debug) write (6,*) 'ga_lkain_2cpl3 converge_precond',
     &   converge_precond
c
      if (lifetime) call errquit('ga_lkain_2cpl3 called with damping',
     &   0,UNKNOWN_ERR)
c     
      odebug = util_print('debug lsolve', print_never) .and. 
     $   ga_nodeid().eq.0
      if (.not. rtdb_get(rtdb, 'cphf:acc',    mt_dbl, 1,
     &                            acc)) acc = 0.01d0*tol
c     
      call ga_inquire(g_x(1), type, n, nvec)

      if (ncomp.gt.1) then
        call ga_inquire(g_x(2), type, ntmp1, ntmp2)      
c       ... jochen: do a sanity check on the array dimensions
        if (ntmp1.ne.n .or. ntmp2.ne.nvec) call errquit
     &     ('ga_lkain_2cpl:inconsistent dimensions of g_x components',
     &     nvec,CALC_ERR)
      endif
      
c     later we combine the two components to vecors of double
c     length if we have two components, otherwise not:
      n2 = n
      if (ncomp.gt.1) n2 = n+n                  

      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec
c     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,1) n2, nvec, maxsub, tol, util_wallsec()
    1   format(//,'Iterative solution of linear equations',/,
     $     '  No. of variables', i9,/,
     $     '  No. of equations', i9,/,
     $     '  Maximum subspace', i9,/,
     $     '       Convergence', 1p,d9.1,/,
     $     '        Start time', 0p,f9.1,/)
        call util_flush(6)
      end if
c     
      do ipm = 1,ncomp
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: Ax',
     $     0, 0, g_Ax(ipm)))
     $     call errquit('lkain: failed allocating Ax', nvec,
     &     GA_ERR)
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: r',
     $     0, 0, g_r(ipm)))
     $     call errquit('lkain_2cpl: failed allocating r', nvec,
     &     GA_ERR)
        if (odiff) then
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_xold(ipm)))
     $       call errquit('lkain: failed allocating xold', nvec,
     &       GA_ERR)
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_Axold(ipm)))
     $       call errquit('lkain: failed allocating Axold', nvec,
     &       GA_ERR)
          call ga_zero(g_xold(ipm))
          call ga_zero(g_Axold(ipm))
        end if                  ! odiff
        call ga_zero(g_Ax(ipm))
        call ga_zero(g_r(ipm))
c       
      enddo                     ! ipm = 1,ncomp
      
c     allocate g_y, g_Ay, and g_r2 with dimension n2 to hold
c     the number of components 
      if (.not. ga_create(MT_DBL, n2, maxsub, 'lkain_2cpl: Y', 
     $   0, 0, g_y))
     $   call errquit('lkain: failed allocating subspace', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n2, maxsub, 'lkain_2cpl: Ay', 
     $   0, 0, g_Ay))
     $   call errquit('lkain: failed allocating subspace2', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n2, nvec, 'lkain_2cpl: r2',
     $   0, 0, g_r2))
     $   call errquit('lkain_2cpl: failed allocating r2', nvec,
     &   GA_ERR)
      
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_r2)
      call ga_sync()
c     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,2)
        call util_flush(6)
    2   format(/
     $     '   iter   nsub   residual    time ',/,
     $     '   ----  ------  --------  --------- ')
      end if
c     
      nsub = 0
      converged = .false.
c     
c     ---------------------
c     start interation loop
c     ---------------------
c     
      do iter = 1, maxiter
        
c       
c ... jochen: here in the iteration loops we keep track
c       of two components of the solution vector, ipm = 1 and 2
c       (ipm stands for + (plus) and - (minus) components)
c       
        if (odiff) then
          do ipm = 1,ncomp   
            call ga_add(1.0d0, g_x(ipm), -1.0d0,
     &         g_xold(ipm),  g_x(ipm))
            call ga_sync()
          enddo
        endif
c       
c ... jochen: call product routine with initial or intermediate
c       solution vector: g_x and g_Ax MUST have two components here
        
        if (debug) write (6,*) 'calling product from ga_lkain_2cpl'
        call product(acc, g_x, g_Ax, g_x_im, g_Ax_im, omega, limag,
     &     lifetime, gamwidth, ncomp)
        call ga_sync()
        if (debug) write (6,*) 'returning product from ga_lkain_2cpl'

c       g_r is zeroed below so we should make sure to do the same
c       with g_r2 here
        call ga_zero(g_r2)
        
        do ipm = 1,ncomp
          
          if (odiff) then
            call ga_add(1.0d0, g_Ax(ipm), 1.0d0,
     &         g_Axold(ipm), g_Ax(ipm))
            call ga_add(1.0d0, g_x(ipm),  1.0d0,
     &         g_xold(ipm),  g_x(ipm))
            call ga_sync()
            call ga_copy(g_x(ipm), g_xold(ipm))
            call ga_copy(g_Ax(ipm), g_Axold(ipm))
          end if
          call ga_zero(g_r(ipm))
          call ga_sync()
c
c         g_Ax = g_b if the system is solved. During the first cycle,
c         g_Ax is calculated from the initial guess
          call ga_add(1.0d0, g_b(ipm),
     &       -1.0d0, g_Ax(ipm), g_r(ipm)) ! The residual
          call ga_sync()
c         
        enddo                   ! ipm = 1,ncomp
        call ga_sync()

c       convergence checking:
c       find the largest element of the residual either 
c       before or after the call to the preconditioner

        if (converge_precond) then
          call precond(g_r(1),  -omega)
          if (ncomp.gt.1) then
            call precond(g_r(2),   omega)
          endif
        endif

        call ga_maxelt(g_r(1), rmax1)
        if (ncomp.gt.1) then
          call ga_maxelt(g_r(2), rmax2)
        else
          rmax2 = 0d0
        endif
        rmax = max(rmax1, rmax2)  
      
        if (oprint .and. ga_nodeid().eq.0) then
          write(6,3) iter, nsub+nvec, rmax, util_wallsec()
          call util_flush(6)
    3     format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
        end if

c       stop iterations if residual is smaller than criterion
        if (rmax .lt. tol) then
          converged = .true.
          goto 100
        end if
c       
c ... jochen: changed 0d0 to omega in the calls below
c       (there were only two in ga_lkain, one with g_aX, one with g_r)
c       for array g_r the preconditioner call is only necessary
c       in case converge_precond is .false.

        call precond(g_Ax(1), -omega)
        if (.not.converge_precond) call precond(g_r(1), -omega)
        if (ncomp.gt.1) then
          call precond(g_Ax(2),  omega)
          if (.not.converge_precond) call precond(g_r(2), omega)
        endif       
        call ga_sync()
        
c       Copy the vectors to the subspace work area
c       
        call ga_copy_patch('n', 
     $     g_Ax(1), 1, n, 1, nvec, 
     $     g_Ay, 1, n, nsub+1, nsub+nvec)
        call ga_sync()
        if (ncomp.gt.1) then
          call ga_copy_patch('n', 
     $       g_Ax(2), 1, n, 1, nvec, 
     $       g_Ay, n+1, n2, nsub+1, nsub+nvec)
        endif
        call ga_copy_patch('n', 
     $     g_x(1), 1, n, 1, nvec, 
     $     g_y, 1, n, nsub+1, nsub+nvec)
        call ga_sync()
        if (ncomp.gt.1) then
          call ga_copy_patch('n', 
     $       g_x(2), 1, n, 1, nvec, 
     $       g_y, n+1, n2, nsub+1, nsub+nvec)
        endif

c       g_r2 is needed below for multiplication
        call ga_copy_patch('n', 
     $     g_r(1), 1, n, 1, nvec, 
     $     g_r2, 1, n, 1, nvec)
        if (ncomp.gt.1) then
          call ga_copy_patch('n', 
     $       g_r(2), 1, n, 1, nvec, 
     $       g_r2, n+1, n2, 1, nvec)
        endif
        
        nsub = nsub + nvec
c       
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
c       
        if (.not. ga_create(MT_DBL, nsub, nsub,
     &     'lkain_2cpl: A', 0, 0, g_a))
     $     call errquit('lkain: allocating g_a?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl: B', 0, 0,g_bb))
     $     call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl: C', 0, 0, g_c))
     $     call errquit('lkain: allocating g_c?', nsub, GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_bb)
        call ga_zero(g_c)
        call ga_sync()
        call ga_dgemm('t','n',nsub,nsub,n2,1.0d0,
     &     g_y,g_Ay,0.0d0,g_a)
        call ga_dgemm('t','n',nsub,nvec,n2,1.0d0,
     &     g_y,g_r2,0.0d0,g_bb)
        call ga_sync()
        if (odebug) call ga_print(g_a)
        if (odebug) call ga_print(g_c)
c       
c       The threshold used here should reflect the accuracy in the
c       products.  If very accurate products are used,
c       then there is big
c       advantage for small cases (maxsub close to n) in using a very
c       small threshold in the SVD solve (e.g., 1e-14), but for more
c       realistic examples (maxsub << n) there is only a little
c       advantage and in the precence of real noise in the products
c       screening with a realistic threshold is important.
c       
        call ga_svd_solve_seq(g_a,g_bb,g_c,1d-14)
        if (odebug) call ga_print(g_c)
c       
c       Form and add the correction, in parts, onto the solution
c       
        call ga_sync()
        call ga_dgemm('n','n',n2,nvec,nsub,-1.0d0,
     &     g_Ay,g_c,1.0d0,g_r2)
        if (odebug) then
          write(6,*) ' The update in the complement '
          call ga_print(g_r2)
        end if
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x
c
        call ga_sync()
        call ga_copy_patch('n', 
     $     g_r2, 1, n, 1, nvec, 
     $     g_r(1), 1, n, 1,nvec)
        if (ncomp.gt.1) then
          call ga_copy_patch('n', 
     $       g_r2, n+1, n2, 1, nvec, 
     $       g_r(2), 1, n, 1,nvec)
        endif
        call ga_sync()
c       
        do ipm = 1,ncomp
          call ga_add(1.0d0, g_r(ipm), 1.0d0, g_x(ipm), g_x(ipm))
        enddo
        
        call ga_sync()
        call ga_dgemm('n','n',n2,nvec,nsub,1.0d0,
     &     g_y,g_c,0.0d0,g_r2)
        call ga_sync()
        if (odebug) then
          write(6,*) ' The update in the subspace '
          call ga_print(g_r2)
        end if
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x
c
        call ga_copy_patch('n', 
     $     g_r2, 1, n, 1, nvec, 
     $     g_r(1), 1, n, 1,nvec)
        if (ncomp.gt.1) then
          call ga_copy_patch('n', 
     $       g_r2, n+1, n2, 1, nvec, 
     $       g_r(2), 1, n, 1,nvec)
        endif
        call ga_sync()
        do ipm = 1,ncomp
          call ga_add(1.0d0, g_r(ipm), 1.0d0, g_x(ipm), g_x(ipm))
        enddo
        call ga_sync()
c        
        if (.not. ga_destroy(g_a)) call errquit
     &     ('lkain_2cpl: a',0, GA_ERR)
        if (.not. ga_destroy(g_bb)) call errquit
     &     ('lkain_2cpl: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('lkain_2cpl: c',0, GA_ERR)
c       
c       Reduce the subspace as necessary
c       
        if (nsub .eq. maxsub) then
          do isub = nvec+1, maxsub, nvec
c           component 1: 
            call ga_copy_patch('n', 
     $         g_Ay,    1, n, isub, isub+nvec-1, 
     $         g_Ax(1), 1, n, 1, nvec)
            call ga_sync()
            call ga_copy_patch('n', 
     $         g_Ax(1), 1, n, 1, nvec,
     $         g_Ay,    1, n, isub-nvec, isub-1)
            call ga_sync()
c           component 2: 
            if (ncomp.gt.1) then
              call ga_copy_patch('n', 
     $           g_Ay,    n+1, n2, isub, isub+nvec-1, 
     $           g_Ax(2), 1,   n,  1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax(2), 1,   n,  1, nvec,
     $           g_Ay,    n+1, n2, isub-nvec, isub-1)            
              call ga_sync()
            endif
c           
c           component 1:
            call ga_copy_patch('n', 
     $         g_y,     1, n, isub, isub+nvec-1, 
     $         g_Ax(1), 1, n, 1, nvec)
            call ga_sync()
            call ga_copy_patch('n', 
     $         g_Ax(1), 1, n, 1, nvec,
     $         g_y,     1, n, isub-nvec, isub-1)            
            call ga_sync()
c           component 2:
            if (ncomp.gt.1) then
              call ga_copy_patch('n', 
     $           g_y,     n+1, n2, isub, isub+nvec-1, 
     $           g_Ax(2), 1,   n,  1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax(2), 1,   n,  1, nvec,
     $           g_y,     n+1, n2, isub-nvec, isub-1)
              call ga_sync()
            endif
c           
          end do                ! isub = nvec+1, maxsub, nvec
          nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub) 
c       
      end do                    ! iter = 1,maxiter
  100 continue                  ! jump here if converged
c     
c     
c     deallocate workspace:
c     
      do ipm = 1,ncomp
        if (odiff) then
          if (.not. ga_destroy(g_xold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',1, GA_ERR)
          if (.not. ga_destroy(g_Axold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',2,GA_ERR)
        end if
        if (.not. ga_destroy(g_Ax(ipm))) call errquit
     &     ('lkain_2cpl: destroy',20, GA_ERR)
        if (.not. ga_destroy(g_r(ipm))) call errquit
     &     ('lkain_2cpl: destroy',5, GA_ERR)
c       
      enddo                     ! ipm = 1,2
      
      if (.not. ga_destroy(g_Ay)) call errquit
     &   ('lkain_2cpl: destroy Ay',3, GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit
     &   ('lkain_2cpl: destroy r',4, GA_ERR)
      if (.not. ga_destroy(g_r2)) call errquit
     &   ('lkain_2cpl: destroy r2',6, GA_ERR)
      
      
c ... jochen: disable this error exit during debuging phase
c     but print a warning instead

c     error exit if this hasn't converged :-(
c     if (.not. converged) call errquit('lkain_2cpl: not converged',0,
c     &    CALC_ERR)

      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
        endif
      endif
c     
      end

c     ******************************************************************
c     ******************************************************************

      subroutine ga_lkain_2cpl3_damp(rtdb, g_x, g_b, g_x_im, g_b_im,
     &   product, precond, 
     $   tol, mmaxsub, maxiter, odiff, oprint, omega, limag,
     &   lifetime, gamwidth, ncomp)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
c
      integer rtdb              ! [input] database handle
      integer g_x(2)            ! [input/output] Initial guess/solution Re
      integer g_x_im(2)         ! [input/output] Initial guess/solution Im
      integer g_b(2)            ! [input] Right-hand side vectors Re
      integer g_b_im(2)         ! [input] Right-hand side vectors Im
      double precision omega    ! [input] coupling parameter
      logical limag             ! [input] imaginary perturbation?
      logical lifetime          ! [input] consider damping or not?
      double precision gamwidth ! [input] damping parameter
      integer ncomp             ! [input] no. of components to treat
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
c ... jochen:
c     This is a modified version of ga_lkain from file ga_it2.F
c     This version allows to solve a coupled set of equations, i.e.
c     there are two right-hand vectors and two initial guesses and two
c     solutions which are coupled. The coupling is mediated by a
c     parameter omega in the call to the preconditioner 
c     (elsewhere, omega is simply called "frequency")
c
c ... jochen: the above comment is from ga_lkain_2cpl3. This here is
c     a modified version of that routine and takes care of a real and an
c     imaginary part for each frequency component. I.e. now arrays
c     have four components ... 
c
c     note: when called from cphf_solve3, odiff = .false. on input
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . coupling parameter omega
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
      integer iter, n, n2, nvec, nsub, isub, type, maxsub, ipm,
     &   ntmp1, ntmp2, n3, n4

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      integer g_y, g_Ay, g_Ax(2), g_r(2), g_r2, g_a, g_bb,
     &   g_c, g_xold(2), g_Axold(2)
      integer g_r_im(2), g_Ax_im(2)
      double precision rmax, rmax1, rmax2, acc
      logical converged
      logical odebug, debug, converge_precond
c
c     =================================================================
c
      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

c     check input key if we should check for convergence
c     after the preconditioner has been applied to the residual
      if (.not. rtdb_get(rtdb, 'aoresponse:precond',    mt_log, 1,
     &                            converge_precond))
     &  converge_precond = .false.
      
      if (debug) write (6,*) 'ga_lkain_2cpl_damp omega =',omega
      if (debug) write (6,*) 'ga_lkain_2cpl_damp limag =',limag
      if (debug) write (6,*) 'ga_lkain_2cpl_damp lifetime =',lifetime
      if (debug) write (6,*) 'ga_lkain_2cpl_damp gamwidth =',gamwidth
      if (debug) write (6,*) 'ga_lkain_2cpl_damp ncomp =', ncomp
      if (debug) write (6,*) 'ga_lkain_2cpl3 converge_precond',
     &   converge_precond
c
c     exit if this is the wrong routine to call (lifetime switch
c     must be set)
      if (.not.lifetime) call errquit
     &   ('ga_lkain_2cpl_damp but lifetime=.F.',0,UNKNOWN_ERR)

c     make sure odiff is false (never tested for odiff = .true.)
      if (odiff) call errquit
     &   ('ga_lkain_2cpl_damp odiff=.T.',0,UNKNOWN_ERR)
c     
      odebug = util_print('debug lsolve', print_never) .and. 
     $   ga_nodeid().eq.0
c     
      if (.not. rtdb_get(rtdb, 'cphf:acc',    mt_dbl, 1,
     &                            acc)) acc = 0.01d0*tol

      call ga_inquire(g_x(1), type, n, nvec)

      if (ncomp.gt.1) then
        call ga_inquire(g_x(2), type, ntmp1, ntmp2)
        
c       ... jochen: do a sanity check on the array dimensions
        if (ntmp1.ne.n .or. ntmp2.ne.nvec) call errquit
     &     ('ga_lkain_2cpl:inconsistent dimensions of g_x components',
     &     nvec,CALC_ERR)
      endif
      
c     later we combine the two components to vecors of double
c     length and combine again Re and Im, i.e. 
c     the dimension is up to 4*n

      n2 = n
      if (ncomp.gt.1) n2 = 2 * n   
      n3 = n
      if (lifetime .or. ncomp.gt.1) n3 = 2 * n
      if (lifetime .and. ncomp.gt.1) n3 = 3 * n
      n4 = n
      if (ncomp.gt.1 .or. lifetime) n4 = 2* n
      if (lifetime .and. ncomp.gt.1) n4 = 4 * n
      if (debug) write (6,*) 'n1n2n3n4',n,n2,n3,n4

      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec
c     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,1) n4, nvec, maxsub, tol, util_wallsec()
    1   format(//,'Iterative solution of linear equations',/,
     $     '  No. of variables', i9,/,
     $     '  No. of equations', i9,/,
     $     '  Maximum subspace', i9,/,
     $     '       Convergence', 1p,d9.1,/,
     $     '        Start time', 0p,f9.1,/)
        call util_flush(6)
      end if
c     
      do ipm = 1,ncomp
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: Ax',
     $     0, 0, g_Ax(ipm)))
     $     call errquit('lkain: failed allocating Ax', nvec,
     &     GA_ERR)
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: r',
     $     0, 0, g_r(ipm)))
     $     call errquit('lkain_2cpl: failed allocating r', nvec,
     &     GA_ERR)
        if (lifetime) then
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: Ax_im',
     $       0, 0, g_Ax_im(ipm)))
     $       call errquit('lkain: failed allocating Ax_im', nvec,
     &       GA_ERR)
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: r_im',
     $       0, 0, g_r_im(ipm)))
     $       call errquit('lkain_2cpl: failed allocating r_im', nvec,
     &       GA_ERR)
        endif                   ! lifetime
        
        if (odiff) then
c         jochen: this part and all subsequent "odiff" parts were
c         never adapted for the imaginary components
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_xold(ipm)))
     $       call errquit('lkain: failed allocating xold', nvec,
     &       GA_ERR)
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_Axold(ipm)))
     $       call errquit('lkain: failed allocating Axold', nvec,
     &       GA_ERR)
          call ga_zero(g_xold(ipm))
          call ga_zero(g_Axold(ipm))
        end if                  ! odiff

        call ga_zero(g_Ax(ipm))
        call ga_zero(g_r(ipm))
        if (lifetime) then
          call ga_zero(g_Ax_im(ipm))
          call ga_zero(g_r_im(ipm))
        endif
c       
      enddo                     ! ipm = 1,ncomp
      
c     allocate g_y, g_Ay, and g_r2 with dimension n4 to hold
c     all necessary components simultaneously
      if (.not. ga_create(MT_DBL, n4, maxsub, 'lkain_2cpl: Y', 
     $   0, 0, g_y))
     $   call errquit('lkain: failed allocating subspace', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n4, maxsub, 'lkain_2cpl: Ay', 
     $   0, 0, g_Ay))
     $   call errquit('lkain: failed allocating subspace2', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n4, nvec, 'lkain_2cpl: r2',
     $   0, 0, g_r2))
     $   call errquit('lkain_2cpl: failed allocating r2', nvec,
     &   GA_ERR)
      
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_r2)
      call ga_sync()
c     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,2)
        call util_flush(6)
    2   format(/
     $     '   iter   nsub   residual    time ',/,
     $     '   ----  ------  --------  --------- ')
      end if
c     
      nsub = 0
      converged = .false.
c     
c     ---------------------
c     start interation loop
c     ---------------------
c     
      do iter = 1, maxiter
        
c       
c ... jochen: here in the iteration loops we keep track
c       of two components of the solution vector, ipm = 1 and 2
c       (ipm stands for + (plus) and - (minus) components)
c       
        if (odiff) then
          do ipm = 1,ncomp   
            call ga_add(1.0d0, g_x(ipm), -1.0d0,
     &         g_xold(ipm),  g_x(ipm))
            call ga_sync()
          enddo
        endif
c       
c ... jochen: call product routine with initial or intermediate
c       solution vector: g_x and g_Ax MUST have dimension two here
c       even if only one of them is used
        
        if (debug) write (6,*)
     &     'calling product from ga_lkain_2cpl_damp'
        call product(acc, g_x, g_Ax,  g_x_im, g_Ax_im, omega, limag,
     &     lifetime, gamwidth, ncomp)
        if (debug) write (6,*)
     &     'returning product from ga_lkain_2cpl_damp'

c       g_r is zeroed below so we should make sure to do the same
c       with g_r2 here
        call ga_zero(g_r2)
        call ga_sync()
        
        do ipm = 1,ncomp
          
          if (odiff) then
c           jochen: odiff stuff presently ignored
            call ga_add(1.0d0, g_Ax(ipm), 1.0d0,
     &         g_Axold(ipm), g_Ax(ipm))
            call ga_add(1.0d0, g_x(ipm),  1.0d0,
     &         g_xold(ipm),  g_x(ipm))
            call ga_sync()
            call ga_copy(g_x(ipm), g_xold(ipm))
            call ga_copy(g_Ax(ipm), g_Axold(ipm))
          end if                ! odiff

          call ga_zero(g_r(ipm))
          if (lifetime) call ga_zero(g_r_im(ipm))

c         g_r will be the quantity -Ax + b, i.e. if the equation system
c         Ax = b is solved then this vector will be zero
c
c         During the first cycle,
c         g_Ax is calculated from the initial guess for which the
c         preconditioner has already been applied (to be more clear:
c         we have divided the perturbation matrix elements by orbital
c         energy denominators, including the frequency term, 
c         and assigned real and imaginary parts accordingly)

          call ga_sync()
          call ga_add(1.0d0, g_b(ipm),
     &       -1.0d0, g_Ax(ipm), g_r(ipm)) ! The residual, Real part

          if (lifetime) then
            call ga_add(1.0d0, g_b_im(ipm),
     &         -1.0d0, g_Ax_im(ipm), g_r_im(ipm)) ! The residual, Im part
          endif
c         
        enddo                   ! ipm = 1,ncomp
        call ga_sync()

c       convergence checking:
c       find the largest element of the residual either 
c       before or after the call to the preconditioner

        if (converge_precond) then
          call precond(g_r(1),  g_r_im(1),  -omega, gamwidth)
          if (ncomp.gt.1) then
            call precond(g_r(2),  g_r_im(2),   omega, gamwidth)
          endif          
        endif
        call ga_sync()
c
        call ga_maxelt(g_r(1), rmax1)
        if (ncomp.gt.1) then
          call ga_maxelt(g_r(2), rmax2)
        else
          rmax2 = 0d0
        endif
        rmax = max(rmax1, rmax2)
        if (lifetime) then
          call ga_maxelt(g_r_im(1), rmax1)
          if (ncomp.gt.1) call ga_maxelt(g_r_im(2), rmax2)
          rmax = max(rmax, rmax1)
          rmax = max(rmax, rmax2)
        endif
        
        if (oprint .and. ga_nodeid().eq.0) then
          write(6,3) iter, nsub+nvec, rmax, util_wallsec()
          call util_flush(6)
    3     format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
        end if
        
c       stop iterations if residual is smaller than criterion
        if (rmax .lt. tol) then
          converged = .true.
          goto 100
        end if
c       
c       Call the preconditioner with the residual as well as with
c       the quantity Ax
c ... jochen: changed 0d0 to omega in the calls below compared to
c       the original routine ga_lkain
c       (there were only two in ga_lkain, one with g_aX, one with g_r)
c       for array g_r the preconditioner call is only necessary
c       in case converge_precond is .false.

        call precond(g_Ax(1), g_Ax_im(1), -omega, gamwidth)
        if (.not.converge_precond) 
     &     call precond(g_r(1),  g_r_im(1),  -omega, gamwidth)
        if (ncomp.gt.1) then
          call precond(g_Ax(2), g_Ax_im(2),  omega, gamwidth)
          if (.not.converge_precond) 
     &       call precond(g_r(2),  g_r_im(2),   omega, gamwidth)
        endif
        call ga_sync()
        if (debug) write (6,*) 'lkain3_damp: back from precond'
                
c       Copy the vectors to the subspace work area
c       
        call ga_copy_patch('n', 
     $     g_Ax(1), 1, n, 1, nvec, 
     $     g_Ay, 1, n, nsub+1, nsub+nvec)
        if (ncomp.gt.1) call ga_copy_patch('n', 
     $     g_Ax(2), 1, n, 1, nvec, 
     $     g_Ay, n+1, n2, nsub+1, nsub+nvec)
        if (lifetime) then
          call ga_copy_patch('n', 
     $       g_Ax_im(1), 1, n, 1, nvec, 
     $       g_Ay, n2+1, n3, nsub+1, nsub+nvec)
          if (ncomp.gt.1) call ga_copy_patch('n', 
     $       g_Ax_im(2), 1, n, 1, nvec, 
     $       g_Ay, n3+1, n4, nsub+1, nsub+nvec)
        endif                   ! lifetime
        call ga_copy_patch('n', 
     $     g_x(1), 1, n, 1, nvec, 
     $     g_y, 1, n, nsub+1, nsub+nvec)
        if (ncomp.gt.1)  call ga_copy_patch('n', 
     $     g_x(2), 1, n, 1, nvec, 
     $     g_y, n+1, n2, nsub+1, nsub+nvec)
        if (lifetime) then
          call ga_copy_patch('n', 
     $       g_x_im(1), 1, n, 1, nvec, 
     $       g_y, n2+1, n3, nsub+1, nsub+nvec)
          if (ncomp.gt.1) call ga_copy_patch('n', 
     $       g_x_im(2), 1, n, 1, nvec, 
     $       g_y, n3+1, n4, nsub+1, nsub+nvec)
        endif                   ! lifetime
        call ga_sync()
        if (debug) write (6,*) 'lkain3_damp: vec to sub complete'

c       g_r2 is needed below for multiplication
        call ga_copy_patch('n', 
     $     g_r(1), 1, n, 1, nvec, 
     $     g_r2, 1, n, 1, nvec)
        if (ncomp.gt.1) call ga_copy_patch('n', 
     $     g_r(2), 1, n, 1, nvec, 
     $     g_r2, n+1, n2, 1, nvec)
        if (lifetime) then
          call ga_copy_patch('n', 
     $       g_r_im(1), 1, n, 1, nvec, 
     $       g_r2, n2+1, n3, 1, nvec)
          if (ncomp.gt.1) call ga_copy_patch('n', 
     $       g_r_im(2), 1, n, 1, nvec, 
     $       g_r2, n3+1, n4, 1, nvec)
        endif                   ! lifetime
        call ga_sync()
        if (debug) write (6,*) 'lkain3_damp: r2 complete'
        
        nsub = nsub + nvec
c       
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
c       
        if (.not. ga_create(MT_DBL, nsub, nsub,
     &     'lkain_2cpl3_damp: A', 0, 0, g_a))
     $     call errquit('lkain: allocating g_a?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl3_damp: B', 0, 0,g_bb))
     $     call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl3_damp: C', 0, 0, g_c))
     $     call errquit('lkain: allocating g_c?', nsub, GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_bb)
        call ga_zero(g_c)
        call ga_sync()
        call ga_dgemm('t','n',nsub,nsub,n4,1.0d0,
     &     g_y,g_Ay,0.0d0,g_a)
        call ga_sync()
        call ga_dgemm('t','n',nsub,nvec,n4,1.0d0,
     &     g_y,g_r2,0.0d0,g_bb)
        call ga_sync()
        if (odebug) call ga_print(g_a)
        if (odebug) call ga_print(g_c)
c       
c       The threshold used here should reflect the accuracy in the
c       products.  If very accurate products are used,
c       then there is big
c       advantage for small cases (maxsub close to n) in using a very
c       small threshold in the SVD solve (e.g., 1e-14), but for more
c       realistic examples (maxsub << n) there is only a little
c       advantage and in the precence of real noise in the products
c       screening with a realistic threshold is important.
c       
        call ga_svd_solve_seq(g_a,g_bb,g_c,1d-14)
        if (odebug) call ga_print(g_c)
c       
c       Form and add the correction, in parts, onto the solution
c       
        call ga_sync()
        call ga_dgemm('n','n',n4,nvec,nsub,-1.0d0,
     &     g_Ay,g_c,1.0d0,g_r2)
        call ga_sync()
        if (odebug) then
          write(6,*) ' The update in the complement '
          call ga_print(g_r2)
        end if
        if (debug) write (6,*) 'lkain3_damp: after dgemm'
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x
        call ga_copy_patch('n', 
     $     g_r2, 1, n, 1, nvec, 
     $     g_r(1), 1, n, 1,nvec)
        if (ncomp.gt.1) call ga_copy_patch('n', 
     $     g_r2, n+1, n2, 1, nvec, 
     $     g_r(2), 1, n, 1,nvec)
        if (lifetime) then
          call ga_copy_patch('n', 
     $       g_r2, n2+1, n3, 1, nvec, 
     $       g_r_im(1), 1, n, 1,nvec)
          if (ncomp.gt.1) call ga_copy_patch('n', 
     $       g_r2, n3+1, n4, 1, nvec, 
     $       g_r_im(2), 1, n, 1,nvec)
        endif                   ! lifetime
        call ga_sync()
        if (debug) write (6,*) 'lkain3_damp: r2 copied'
c               
        do ipm = 1,ncomp
          call ga_add(1.0d0, g_r(ipm), 1.0d0, g_x(ipm), g_x(ipm))
          if (lifetime)
     &       call ga_add(1.0d0, g_r_im(ipm),
     &       1.0d0, g_x_im(ipm), g_x_im(ipm))
        enddo
        call ga_sync()
c       
        call ga_dgemm('n','n',n4,nvec,nsub,1.0d0,
     &     g_y,g_c,0.0d0,g_r2)
        call ga_sync()
        if (odebug) then
          write(6,*) ' The update in the subspace '
          call ga_print(g_r2)
        end if
        if (debug) write (6,*) 'lkain3_damp: subspace updated'
        
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x
        call ga_copy_patch('n', 
     $     g_r2, 1, n, 1, nvec, 
     $     g_r(1), 1, n, 1,nvec)
        if (ncomp.gt.1) call ga_copy_patch('n', 
     $     g_r2, n+1, n2, 1, nvec, 
     $     g_r(2), 1, n, 1,nvec)
        if (lifetime) then
          call ga_copy_patch('n', 
     $       g_r2, n2+1, n3, 1, nvec, 
     $       g_r_im(1), 1, n, 1,nvec)
          if (ncomp.gt.1) call ga_copy_patch('n', 
     $       g_r2, n3+1, n4, 1, nvec, 
     $       g_r_im(2), 1, n, 1,nvec)
        endif                   ! lifetime
        call ga_sync()
        do ipm = 1,ncomp
          call ga_add(1.0d0, g_r(ipm), 1.0d0, g_x(ipm), g_x(ipm))
          if (lifetime)
     &       call ga_add(1.0d0, g_r_im(ipm),
     &       1.0d0, g_x_im(ipm), g_x_im(ipm))
        enddo
        call ga_sync()
        if (debug) write (6,*) 'lkain3_damp: gr to gx complete'
c        
        if (.not. ga_destroy(g_a)) call errquit
     &     ('lkain_2cpl: a',0, GA_ERR)
        if (.not. ga_destroy(g_bb)) call errquit
     &     ('lkain_2cpl: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('lkain_2cpl: c',0, GA_ERR)
c       
c       Reduce the subspace as necessary
c       
        if (nsub .eq. maxsub) then
          do isub = nvec+1, maxsub, nvec
c           Real part:
c           component 1: 
            call ga_copy_patch('n', 
     $         g_Ay,    1, n, isub, isub+nvec-1, 
     $         g_Ax(1), 1, n, 1, nvec)
            call ga_sync()
            call ga_copy_patch('n', 
     $         g_Ax(1), 1, n, 1, nvec,
     $         g_Ay,    1, n, isub-nvec, isub-1)
c           component 2: 
            if (ncomp.gt.1) then
              call ga_copy_patch('n', 
     $           g_Ay,    n+1, n2, isub, isub+nvec-1, 
     $           g_Ax(2), 1,   n,  1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax(2), 1,   n,  1, nvec,
     $           g_Ay,    n+1, n2, isub-nvec, isub-1)  
            endif
            if (lifetime) then
c             Imaginary part:
c             component 1: 
              call ga_copy_patch('n', 
     $           g_Ay,    n2+1, n3, isub, isub+nvec-1, 
     $           g_Ax_im(1), 1, n, 1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax_im(1), 1, n, 1, nvec,
     $           g_Ay,    n2+1, n3, isub-nvec, isub-1)
c             component 2: 
              if (ncomp.gt.1) then
                call ga_copy_patch('n', 
     $             g_Ay,    n3+1, n4, isub, isub+nvec-1, 
     $             g_Ax_im(2), 1,   n,  1, nvec)
                call ga_sync()
                call ga_copy_patch('n', 
     $             g_Ax_im(2), 1,   n,  1, nvec,
     $             g_Ay,    n3+1, n4, isub-nvec, isub-1)         
              endif
            endif               ! lifetime
c             
c           Real part:
c           component 1:
            call ga_sync()
            call ga_copy_patch('n', 
     $         g_y,     1, n, isub, isub+nvec-1, 
     $         g_Ax(1), 1, n, 1, nvec)
            call ga_sync()
            call ga_copy_patch('n', 
     $         g_Ax(1), 1, n, 1, nvec,
     $         g_y,     1, n, isub-nvec, isub-1)
c           component 2:
            if (ncomp.gt.1)  then
              call ga_copy_patch('n', 
     $           g_y,     n+1, n2, isub, isub+nvec-1, 
     $           g_Ax(2), 1,   n,  1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax(2), 1,   n,  1, nvec,
     $           g_y,     n+1, n2, isub-nvec, isub-1)
            endif
c
            if (lifetime) then
c             Imaginary Part:
c             component 1:
              call ga_copy_patch('n', 
     $           g_y,     n2+1, n3, isub, isub+nvec-1, 
     $           g_Ax_im(1), 1, n, 1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax_im(1), 1, n, 1, nvec,
     $           g_y,     n2+1, n3, isub-nvec, isub-1)
c             component 2:
              if (ncomp.gt.1) then
                call ga_copy_patch('n', 
     $             g_y,     n3+1, n4, isub, isub+nvec-1, 
     $             g_Ax_im(2), 1,   n,  1, nvec)
                call ga_sync()
                call ga_copy_patch('n', 
     $             g_Ax_im(2), 1,   n,  1, nvec,
     $             g_y,     n3+1, n4, isub-nvec, isub-1)
              endif
            endif               ! lifetime
c           
          end do                ! isub = nvec+1, maxsub, nvec
          nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub) 
        if (debug) write (6,*) 'lkain3_damp: sub reduction comp.'
c       
      enddo                     ! iter = 1,maxiter
  100 continue                  ! jump here if converged
c     
c     
c     deallocate workspace:
c     
      do ipm = 1,ncomp
        if (odiff) then
          if (.not. ga_destroy(g_xold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',1, GA_ERR)
          if (.not. ga_destroy(g_Axold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',2,GA_ERR)
        end if
        if (.not. ga_destroy(g_Ax(ipm))) call errquit
     &     ('lkain_2cpl: destroy',20, GA_ERR)
        if (.not. ga_destroy(g_r(ipm))) call errquit
     &     ('lkain_2cpl: destroy',5, GA_ERR)
        if (lifetime) then
          if (.not. ga_destroy(g_Ax_im(ipm))) call errquit
     &       ('lkain_2cpl: destroy',201, GA_ERR)
          if (.not. ga_destroy(g_r_im(ipm))) call errquit
     &       ('lkain_2cpl: destroy',51, GA_ERR)
        endif
c       
      enddo                     ! ipm = 1,2
      
      if (.not. ga_destroy(g_Ay)) call errquit
     &   ('lkain_2cpl: destroy Ay',3, GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit
     &   ('lkain_2cpl: destroy r',4, GA_ERR)
      if (.not. ga_destroy(g_r2)) call errquit
     &   ('lkain_2cpl: destroy r2',6, GA_ERR)

      if (debug) then
        write (6,*) 'ga_lkain_2cpl3_damp: solution vectors:'
        call ga_print(g_x(1))
        if (ncomp.gt.1) call ga_print(g_x(2))
        if (lifetime) then
          call ga_print(g_x_im(1))
          if (ncomp.gt.1) call ga_print(g_x_im(2))
        endif
      endif
      
      
c ... jochen: disable this error exit during debuging phase
c     but print a warning instead

c     error exit if this hasn't converged :-(
c     if (.not. converged) call errquit('lkain_2cpl: not converged',0,
c     &    CALC_ERR)

      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
        endif
      endif
c     
      end
